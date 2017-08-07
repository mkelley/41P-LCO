# Licensed under a MIT style license - see LICENSE
"""science - Comet science with LCO data."""

import os
from collections import OrderedDict
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle

class Science:
    """Comet science with LCO data.

    Parameters
    ----------
    logger : logging.Logger, optional
      Log to this `Logger` instance, otherwise log to the python
      default.
    rlevel : int, optional
      LCO reduction level to consider, or `None` for the most current.
    log_to_file : bool, optional
      Set to `False` to disable logging to a file, e.g., for debugging
      in interactive mode.

    """

    def __init__(self, rlevel=None, log_to_file=True):
        import logging
        from .core import setup_logger, open_log_file, config

        self.logger = setup_logger('tgk.science')
        if log_to_file:
            log_file = os.sep.join((config['science path'], 'tgk-science.log'))
            open_log_file(log_file, 'tgk.science')

        self.config = config

        assert isinstance(rlevel, (int, type(None))), 'rlevel must be integer or `None`.'
        self.rlevel = rlevel

        if rlevel is None:
            self.logger.info('Processing the most recent reduction level.')
        else:
            self.logger.info('Processing reduction level {} only.'.format(
                rlevel))

        self.observation_log = ObservationLog()
        self.geometry_table = GeometryTable()
        self.processing_history = ProcessingHistory()

    def find_data(self):
        """Find comet data."""
        from glob import glob

        if self.rlevel is None:
            rlevel_dir = 'e[0-9][0-9]'
        else:
            rlevel_dir = 'e{:02d}'.format(self.rlevel)

        pat = os.sep.join((self.config['download path'], rlevel_dir,
                           '2017*', '*fz'))
        all_files = sorted(glob(pat))
        msg = '{} files found'.format(len(all_files))

        # frame name, rlevel, full path
        all_files = sorted([(os.path.basename(f)[:-12], f[-11:-8], f)
                            for f in all_files])

        if self.rlevel is None:
            unique, duplicates = self._find_duplicates(all_files)

            # thanks to sort order and rlevel format, x[-1] will
            # always be the most current version:
            all_files = sorted(unique + [x[-1] for x in duplicates])

            msg += ', {} remain after removing duplicates'.format(len(all_files))

        self.logger.debug(msg + '.')
        self.files = all_files

    def _find_duplicates(self, files):
        """Find unique and duplicated frames in `files`.

        `files` must be sorted.

        """

        assert all([files[i][0] <= files[i + 1][0]
                    for i in range(len(files) - 1)]), 'files must be sorted'
        
        i = 0
        unique = []
        duplicates = []
        while i < len(files):
            # files are sorted by frame name; how many are the same?
            found = []
            basename = files[i][0]
            while i < len(files) and basename in files[i][0]:
                found.append(files[i])
                i += 1

            if len(found) == 1:
                unique.append(found[0])
            else:
                duplicates.append(found)

        return unique, duplicates
        
    def find_new_data(self):
        """Determine which data have not been processed or have an updated rlevel."""
        from .core import timestamp
        new_data = []
        for f in self.files:
            try:
                row = self.processing_history.get_frame(f[0])
                last_rlevel = row[1]
                if last_rlevel < f[1]:
                    new_data.append(f)
            except IndexError:
                new_data.append(f)
                continue

        self.logger.info(timestamp()[:-7] + ' {} frames are new or updated.'.format(len(new_data)))
        return new_data

    def continuous_process(self):
        """Continuously run the science pipeline."""
        
        import time
        from astropy.time import Time
        from .core import timestamp
        
        delay = 5 * u.min
        last_science = Time('2000-01-01')
        self.logger.info('{} Entering continuous science mode, checking local archive every {}.'.format(timestamp(), delay))

        try:
            while True:
                now = Time.now()
                if (now - last_science) > delay:
                    self.process()
                    last_science = Time.now()
                else:
                    dt = int((now - last_science).sec)
                    sleep = int(delay.to(u.s).value) - dt + 2
                    self.logger.debug(
                        'Last science: {} s ago.  Sleep {} s.'.format(dt, sleep))
                    time.sleep(sleep)
                    self.logger.debug('Awake!')
        except KeyboardInterrupt:
            self.logger.info('Caught interrupt signal.  Shutdown.')

    @classmethod
    def get_frame_data(cls, frame, filename):
        """Data as Image and Observation classes."""
        from astropy.io import fits

        with fits.open(filename, mode='readonly', lazy_load=False) as hdu:
            im = Image(hdu)
            obs = Observation(im.header.copy())

        return im, obs

    @classmethod
    def get_geometry(cls, frame, obs):
        """Geometry metadata as Geometry class."""
        # only call JPL/HORIZONS if needed
        try:
            g = cls.geometry_table
        except AttributeError:
            g = GeometryTable()

        try:
            data = g.get_frame(frame)
        except IndexError:
            data = None
            
        return Geometry(obs, data=data)

    @classmethod
    def get_minion_history(cls, frame):
        """Return any known minion processing history."""
        try:
            processing_history = cls.processing_history
        except AttributeError:
            processing_history = ProcessingHistory()

        try:
            minion_history = processing_history.get_frame(frame)[2]
            if minion_history is not np.ma.masked:
                minion_history = minion_history.split(';')
            else:
                minion_history = []
        except IndexError:
            minion_history = []

        return minion_history
            
    def process(self, reprocess=[]):
        """Run the science pipeline and post-science hook.

        Processing history and observation logs are updated.

        Parameters
        ----------
        reprocess : list, optional
          Default is to run only if there are new files.  Set to 'all'
          to rerun on all frames and tables, or 'tables' to rerun
          tables, even if no new frames are found.

        """

        import subprocess
        from . import minions
        from .core import timestamp

        self.find_data()
        
        # determine which files to process
        if len(reprocess) == 0:
            files = self.find_new_data()
        elif 'all' in reprocess:
            files = self.files
        elif any([r in minions.frame_minion_names for r in reprocess]):
            files = self.files
        else:
            files = []

        if len(files) == 0:
            if len(reprocess) == 0:
                # we're done
                return
    
        # First, minions that operate on individual frames
        n_remaining = len(files)
        for frame, rlevel, filename in files:
            n_remaining -= 1
            self.logger.info('  {} [{} of {} remaining]'.format(frame, n_remaining, len(files)))

            try:
                im, obs = self.get_frame_data(frame, filename)
            except IOError as e:
                self.logger.error('{}: {}'.format(e, filename))
                
            if frame not in self.observation_log.tab['frame']:
                self.observation_log.update(obs.log_row())
                
            geom = self.get_geometry(frame, obs)

            if frame not in self.geometry_table.tab['frame']:
                self.geometry_table.update(geom.geometry_row())

            hist = minions.frame(self.config, im, obs, geom,
                                 reprocess=reprocess)
            minion_history = self.get_minion_history(frame)
            minion_history.extend(hist)

            row = (frame, rlevel, ';'.join(np.unique(minion_history)))
            self.processing_history.update(row)

        # Next, minions that operate on tables
        self.logger.debug('Running table minions.')
        minions.table(self.config, reprocess=reprocess)

        # Finally, the post-science hook
        if len(self.config['post-science hook']) > 0:
            self.logger.info(timestamp()[:-7] + ' Running post-science hook...')
            try:
                r = subprocess.check_output(self.config['post-science hook'])
                self.logger.info(r.decode())
            except Exception as e:
                err = '[Post-science hook] {}: {}'.format(type(e).__name__, e)
                self.logger.error(err)
            self.logger.info(timestamp()[:-7] + '...done.')

########################################################################
class Observation:
    """Telescope and image meta data for each LCO frame.

    Parameters
    ----------
    header : astropy.io.fits.Header
      The FITS header for the frame in question.

    Attributes
    ----------
    exptime : Quantity
      Exposure time.
    filter : string
      Filter used.
    rlevel : int
      Reduction level.

    time : Time
      Mid-point of observation.
    start : Time
      Start time of observation.
    stop : Time
      Stop time of observation.
    wcs : astropy.wcs.WCS
      The world coordinate system.

    gain : Quantity
      Detector gain.
    nonlinearity : Quantity
      Non-linearity limit.
    binning : list of int
      CCD pixel bin factor.
    pixel_scale : Quantity
      Nominal pixel scale, binned pixels.
    readnoise : Quantity
      Detector readnoise.
    trimsec : tuple of slices
      Section of useful data.

    site : tuple of strings
      The observatory site ID and name.
    enclosure : tuple of strings
      The enclosure ID and name.
    telescope : tuple of strings
      The telescope ID and name.
    instrument : string
      The instrument ID.

    longitude : Angle
      Observatory longitude (East).
    latitude : Angle
      Observatory latitude.
    elevation : Quantity
      Observatory elevation.

    lunar_alt : Angle
      Lunar altitude.
    lunar_elong : Angle
      Lunar elongation.
    solar_alt : Angle
      Solar altitude.
    solar_elong : Angle
      Solar elongation.

    """

    def __init__(self, header):        
        self.header = header

    def log_row(self):
        """Format data for adding to an observation log."""
        return (self.frame_name, self.time.iso[:10], self.time.iso[11:19],
                self.pixel_scale.value, '{0[0]}x{0[1]}'.format(self.binning),
                self.exptime.value, self.airmass, self.filter)

    ############################################################
    @property
    def frame_name(self):
        return self.header['ORIGNAME'][:-9]
    
    @property
    def exptime(self):
        return self.header['EXPTIME'] * u.s

    @property
    def filter(self):
        return self.header['FILTER']

    @property
    def rlevel(self):
        return self.header['RLEVEL']

    @property
    def airmass(self):
        return self.header['AIRMASS']

    ############################################################
    @property
    def time(self):
        from astropy.time import Time
        return Time(self.header['DATE-OBS'], scale='utc')

    @property
    def start(self):
        from astropy.time import Time
        return Time(self.header['DATE-OBS'], scale='utc')

    @property
    def stop(self):
        from astropy.time import Time
        t = self.header['DATE-OBS'][:11] + self.header['UTSTOP']
        return Time(t, scale='utc')

    @property
    def wcs(self):
        from astropy.wcs import WCS
        return WCS(self.header)

    ############################################################
    @property
    def site(self):
        return (self.header['SITEID'], self.header['SITE'])

    @property
    def enclosure(self):
        return (self.header['ENCID'], self.header['ENCLOSUR'])

    @property
    def telescope(self):
        return (self.header['TELID'], self.header['TELESCOP'])

    @property
    def instrument(self):
        return self.header['INSTRUME']
    
    ############################################################
    @property
    def gain(self):
        return self.header['GAIN'] * u.electron / u.adu

    @property
    def nonlinearity(self):
        return self.header['MAXLIN'] * u.adu

    @property
    def binning(self):
        return [int(x) for x in self.header['CCDSUM'].split()]
    
    @property
    def pixel_scale(self):
        return self.header['PIXSCALE'] * self.binning[0] * u.arcsec
    
    @property
    def readnoise(self):
        return self.header['RDNOISE'] * u.electron

    @property
    def trimsec(self):
        if self.site[0] == 'coj':
            # TRIMSEC is bad for coj, define our own
            return np.s_[31:4069, 1:4094]
        else:
            x, y = self.header['TRIMSEC'][1:-1].split(',')
            x1, x2 = [int(z) for z in x.split(':')]
            y1, y2 = [int(z) for z in y.split(':')]
            return np.s_[(y1 - 1):y2, (x1 - 1):x2]
    
    ############################################################
    @property
    def latitude(self):
        Angle(self.header['LATITUDE'] * u.deg)

    @property
    def longitude(self):
        Angle(self.header['LONGITUD'] * u.deg)

    @property
    def elevation(self):
        return self.header['HEIGHT'] * u.m

    ############################################################
    @property
    def lunar_alt(self):
        return Angle(self.header['MOONALT'] * u.deg)

    @property
    def lunar_elong(self):
        return Angle(self.header['MOONDIST'] * u.deg)

    @property
    def solar_alt(self):
        return Angle(self.header['SUNALT'] * u.deg)

    @property
    def solar_elong(self):
        return Angle(self.header['SUNDIST'] * u.deg)

class Geometry:
    """Comet and comet-observer geometrical circumstances.

    Parameters
    ----------
    obs : Observation
      The observation meta data.
    data : astropy.table.Row
      Get parameters from this geometry table row instead of
      JPL/HORIZONS.

    Attributes
    ----------
    radec_predict : astropy.coordinates.SkyCoord
      Predectied J2000 coordinates.
    mu : Quantity
      Proper motion.
    mu_ra : Qunatity
      Right ascention rate of change times cos(declination).
    mu_dec : Quantity
      Declination rate of change.

    rh : Quantity
      Heliocentric distance.
    delta : Quantity
      Observer-comet distance.
    phase : Angle
      Phase angle.

    sun_position_angle : Angle
      Projected comet->Sun vector.
    velocity_position_angle : Angle
      Projected comet velocity vector.

    """
    def __init__(self, obs, data=None):
        from astropy.table import Row
        from astropy.coordinates import SkyCoord

        self.obs = obs
        
        if data is None:
            self.get_ephemeris(obs)
        else:
            assert isinstance(data, Row)
            self._geom = data
            self._horizons_query = None

    def get_ephemeris(self, obs):
        """Get the comet ephemeris from JPL/HORIZONS."""
        import logging
        from astropy.coordinates import SkyCoord
        import callhorizons
        from . import lco

        logger = logging.getLogger('tgk.science')
        logger.debug('      Get geometry from HORIZONS.')
        q = callhorizons.query('41P')
        q.set_discreteepochs([obs.time.jd])
        obs_code = lco.mpc_codes[(obs.site[0], obs.enclosure[0], obs.telescope[0])]
        n = q.get_ephemerides(obs_code)
        if n != 1:
            raise EphemerisError('Bad return from JPL/HORIZONS, check URL: {}'.format(q.url))

        self._horizons_query = q

    def geometry_row(self):
        """Format data for adding to a geometry table."""
        return (self.obs.frame_name, self.time.iso[:10], self.time.iso[11:19],
                self.radec_predict.ra.deg, self.radec_predict.dec.deg,
                self.mu_ra.value, self.mu_dec.value, self.rh.value,
                self.delta.value, self.phase.value, self.sun_position_angle.value,
                self.velocity_position_angle.value)

    ############################################################
    @property
    def time(self):
        return self.obs.time

    ############################################################
    @property
    def radec_predict(self):
        from astropy.coordinates import SkyCoord
        if self._horizons_query is None:
            return SkyCoord(self._geom['ra predict'] * u.deg,
                            self._geom['dec predict'] * u.deg)
        else:
            return SkyCoord(self._horizons_query['RA'][0] * u.deg,
                            self._horizons_query['DEC'][0] * u.deg)

    @property
    def mu(self):
        return np.sqrt(self.mu_ra**2 + self.mu_dec**2)

    @property
    def mu_ra(self):
        if self._horizons_query is None:
            return self._geom['mu ra'] * u.arcsec / u.s
        else:
            return self._horizons_query['RA_rate'][0] * u.arcsec / u.s

    @property
    def mu_dec(self):
        if self._horizons_query is None:
            return self._geom['mu dec'] * u.arcsec / u.s
        else:
            return self._horizons_query['DEC_rate'][0] * u.arcsec / u.s

    ############################################################
    @property
    def rh(self):
        if self._horizons_query is None:
            return self._geom['rh'] * u.au
        else:
            return self._horizons_query['r'][0] * u.au

    @property
    def delta(self):
        if self._horizons_query is None:
            return self._geom['delta'] * u.au
        else:
            return self._horizons_query['delta'][0] * u.au

    @property
    def phase(self):
        if self._horizons_query is None:
            return self._geom['phase'] * u.deg
        else:
            return Angle(self._horizons_query['alpha'][0] * u.deg)

    ############################################################
    @property
    def sun_position_angle(self):
        if self._horizons_query is None:
            return self._geom['sun PA'] * u.deg
        else:
            a = (self._horizons_query['sunTargetPA'][0] + 180) % 360.0
            return Angle(a * u.deg)

    @property
    def velocity_position_angle(self):
        if self._horizons_query is None:
            return self._geom['velocity PA'] * u.deg
        else:
            a = (self._horizons_query['velocityPA'][0] + 180) % 360.0
            return Angle(a * u.deg)

class Image:
    """The image data and LCO photometry table.

    Parameters
    ----------
    hdu : astropy.io.fits.HDUList
      The FITS data for the frame in question.

    Attributes
    ----------
    header : astropy.io.fits.Header
      The science header.
    data : ndarray
      The science data in adu/s.
    bpm : ndarray
      The bad pixel mask.
    cat : astropy.table.Table
      The LCO pipeline photometry catalog.

    """

    def __init__(self, hdu):
        from astropy.table import Table
        from astropy.io import fits
        self._hdu = fits.HDUList()
        self._hdu.append(hdu[0].copy())
        self._hdu.append(hdu['sci'].copy())
        self._bpm = hdu['bpm'].data.astype(bool)
        self._cat = Table(hdu['cat'].data.copy())

    @property
    def header(self):
        return self._hdu['sci'].header

    @property
    def data(self):
        return self._hdu['sci'].data

    @property
    def bpm(self):
        return self._bpm

    @property
    def cat(self):
        return self._cat

class ScienceTableError(Exception):
    pass

class ScienceTable:
    """Science table convenience class.

    Parameters
    ----------
    filename : string
      The file to read/write table data to/from.  The science path is
      automatically prepended.
    verbose : bool, optional
      Log verbosity.

    Requires: _table_title, _table_columns, _table_dtypes,
      _table_meta, _table_formats, _table_sort

    """
    def __init__(self, filename, verbose=True):
        from .core import config
        self.verbose = verbose
        self.filename = os.sep.join([config['science path'], filename])
        self.read()

    def read(self):
        """Read in the table."""
        
        import logging
        from astropy.io import ascii
        from astropy.table import Table, vstack

        self.tab = Table(names=self._table_columns,
                         dtype=self._table_dtypes)
        self.tab.meta = self._table_meta

        assert isinstance(self._table_formats, (dict, tuple, list))
        if isinstance(self._table_formats, (tuple, list)):
            for i, col in enumerate(self.tab.colnames):
                self.tab[col].format = self._table_formats[i]
        else:
            for col, cformat in self._table_formats.items():
                self.tab[col].format = cformat

        logger = logging.getLogger('tgk.science')
        if os.path.exists(self.filename):
            try:
                tab = ascii.read(self.filename, format='ecsv')
            except ValueError as e:
                raise ScienceTableError("Error reading {}: {}".format(self.filename, e))
            
            tab.meta = OrderedDict()  # no need to duplicate meta data
            self.tab = vstack((self.tab, tab))

            if self.verbose:
                logger.debug('Read {} table from {}.'.format(
                    self._table_title, self.filename))
        else:
            if self.verbose:
                logger.debug('Initialized empty {} table.'.format(
                    self._table_title))

    def write(self, delimiter=','):
        """Write table to file."""
        if self._table_sort is not None:
            self.tab.sort(self._table_sort)
        self.tab.write(self.filename, delimiter=delimiter,
                       overwrite=True, format='ascii.ecsv')

    def _get_unique_row(self, column, value):
        """Get the row that matches value in column, or else `None`."""
        if value in self.tab[column]:
            i = np.flatnonzero(self.tab[column] == value)[0]
            return self.tab[i]
        else:
            raise IndexError('{} not in {} table.'.format(
                value, self._table_title))
        
    def _update_unique_column(self, unique_column, row):
        """Update table with new data.

        If `unique_column` is already present in the table, it is
        replaced.

        """

        assert np.unique(self.tab[unique_column]).size == len(self.tab), '{} entries are not unique in {} table.'.format(unique_column, self._table_title)

        i = self.tab.index_column(unique_column)
        j = self.tab[unique_column] == row[i]
        if any(j):
            self.tab.remove_row(np.flatnonzero(j)[0])

        self.tab.add_row(row)
        self.write()

    def update(self, row):
        """Add row to table."""
        self.tab.add_row(row)
        self.write()

class ProcessingHistory(ScienceTable):
    """Processing history table."""

    _table_title = 'processing history'
    _table_columns = ('frame', 'rlevel', 'minions')
    _table_dtypes = ('U64', 'U3', 'U256')
    _table_meta = OrderedDict()
    _table_formats = {}
    _table_sort = ['frame']

    def __init__(self):
        ScienceTable.__init__(self, 'processing-history.csv')

    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)

class ObservationLog(ScienceTable):
    """Observation log table."""

    _table_title = 'observation log'
    _table_columns = (
        'frame', 'date', 'time', 'pixel scale', 'binning', 'exptime',
        'airmass', 'filter'
    )
    _table_dtypes = ('U64', 'U10', 'U8', float, 'U3', float, float, 'U2')
    _table_meta = OrderedDict()
    _table_meta['date/time'] = 'UTC.'
    _table_meta['binning'] = 'Original CCD pixels bin factor.'
    _table_meta['pixel scale'] = 'Binned pixel scale, arcsec.'
    _table_meta['exptime'] = 'Frame exposure time, s.'
    _table_meta['filter'] = 'LCO filter name.'
    _table_formats = {
        'exptime': '{:.1f}',
        'airmass': '{:.3f}',
    }
    _table_sort = ['date', 'time']

    def __init__(self):
        ScienceTable.__init__(self, 'observation-log.csv')

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)

class GeometryTable(ScienceTable):
    """Geometry table."""
    # changes to this table should be reflected in the Geometry class
    _table_title = 'geometry'
    _table_columns = ('frame', 'date', 'time', 'ra predict', 'dec predict',
                      'mu ra', 'mu dec', 'rh', 'delta', 'phase', 'sun PA',
                      'velocity PA')
    _table_dtypes = ['U64', 'U10', 'U8'] + [float] * 9
    _table_meta = OrderedDict()
    _table_meta['date'] = 'UTC'
    _table_meta['time'] = 'midpoint, UTC'
    _table_meta['ra/dec predict'] = 'predicted for telescope by JPL/HORIZONS, deg'
    _table_meta['mu ra/dec'] = 'proper motion from JPL HORIZONS, arcsec/s'
    _table_meta['rh'] = 'au'
    _table_meta['delta'] = 'au'
    _table_meta['phase'] = 'Sun-comet-observer, deg'
    _table_meta['sun/velocity PA'] = 'Projected vector position angle, deg E of N'
    _table_formats = {
        'ra predict': '{:.6f}',
        'dec predict': '{:.6f}',
        'mu ra': '{:.6f}',
        'mu dec': '{:.6f}',
        'rh': '{:.4f}',
        'delta': '{:.4f}',
        'phase': '{:.2f}',
        'sun PA': '{:.2f}',
        'velocity PA': '{:.2f}',
    }
    _table_sort = ['date', 'time']

    def __init__(self):
        ScienceTable.__init__(self, 'geometry.csv')

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)

    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)
