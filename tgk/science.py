# Licensed under a MIT style license - see LICENSE
"""science - Comet science with LCO data."""

import os
from collections import OrderedDict
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from . import core
from .core import setup_logger, open_log_file, config

logger = setup_logger('tgk.science')
log_file = os.sep.join((config['science path'], 'tgk-science.log'))
open_log_file(log_file, 'tgk.science')

class Science:
    """Comet science with LCO data.

    Parameters
    ----------
    logger : logging.Logger, optional
      Log to this `Logger` instance, otherwise log to the python
      default.
    rlevel : int, optional
      LCO reduction level to consider, or `None` for the most current.

    """

    def __init__(self, rlevel=None):
        import logging
        from .core import config

        self.logger = logging.getLogger('tgk.science')
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
        self.read_processing_history()
        self.find_data()

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
        self.logger.info('{} files match.'.format(
            len(all_files)))

        # frame name, rlevel, full path
        all_files = sorted([(os.path.basename(f)[:-12], f[-11:-8], f)
                            for f in all_files])

        if self.rlevel is None:
            unique, duplicates = self._find_duplicates(all_files)
            self.logger.info('{} frames have multiple reduction levels.'.format(len(duplicates)))

            # thanks to sort order and rlevel format, x[-1] will
            # always be the most current version:
            all_files = sorted(unique + [x[-1] for x in duplicates])

            self.logger.info('{} files will be considered.'.format(
                len(all_files)))

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
        new_data = []
        for f in self.files:
            if f[0] in self.processing_history:
                last_rlevel = self.processing_history[f[0]][0]
                if last_rlevel < f[1]:
                    new_data.append(f)
            else:
                new_data.append(f)

        self.logger.info('{} frames are new or updated.'.format(len(new_data)))
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
                    self.logger.info(timestamp() + ' Checking local archive.')
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
    
    def process(self, all_files=False):
        """Run the science pipeline.

        Processing history and observation logs are updated.

        Parameters
        ----------
        all_files : bool, optional
          Default is to run on new files.  Set `all_files` to `True`
          to run on all files.

        """
        
        from astropy.io import fits
        from . import minions

        if all_files:
            files = self.files
        else:
            files = self.find_new_data()

        if len(files) == 0:
            self.logger.info('No files to process.')
            return

        # First, minions that operate on individual frames
        n_remaining = len(files)
        for frame, rlevel, filename in files:
            n_remaining -= 1
            self.logger.info('  {} [{} remaining]'.format(frame, n_remaining))

            minion_history = []
            with fits.open(filename) as hdu:
                im = Image(hdu)
                obs = Observation(im.header)

                # only call JPL/HORIZONS if needed
                if frame in self.geometry_table.tab['frame']:
                    i = np.flatnonzero(
                        self.geometry_table.tab['frame'] == frame)[0]
                    data = self.geometry_table.tab[i]
                else:
                    data = None
                geom = Geometry(obs, data=data)

                minion_history.extend(minions.frame(self.config, im, obs, geom))
                
                if frame not in self.observation_log.tab['frame']:
                    self.observation_log.update(obs.log_row())

                if frame not in self.geometry_table.tab['frame']:
                    self.geometry_table.update(geom.geometry_row())

            self.processing_history[frame] = (rlevel, minion_history)
            self.save_processing_history()

        # Next, minions that operate on derived data
        pass


    def read_processing_history(self):
        """Read in processing history file."""
        self.processing_history = {}
        
        fn = os.sep.join([self.config['science path'], 'processed-data.txt'])
        if os.path.exists(fn):
            with open(fn, 'r') as inf:
                for line in inf:
                    frame, rlevel, minions = line.strip().split(';')
                    minions = minions.split(',')
                    self.processing_history[frame] = (rlevel, minions)
            self.logger.info('Read processing history from {} .'.format(fn))
        else:
            self.logger.info('Processing history {} does not exist.'.format(fn))

    def save_processing_history(self):
        """Save the processing history."""
        fn = os.sep.join([self.config['science path'], 'processed-data.txt'])
        with open(fn, 'w') as outf:
            for frame, (rlevel, minions) in self.processing_history.items():
                outf.write(';'.join([frame, rlevel, ','.join(minions)]))
                outf.write('\n')

        self.logger.debug('Wrote processing history to {} .'.format(fn))

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
    trim : tuple of slices
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
                '{0[0]}x{0[1]}'.format(self.binning),
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
        x, y = self.header['TRIMSEC'][1:-1].split(',')
        x1, x2 = x.split(':')
        y1, y2 = y.split(':')
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
        logger.info('      Get geometry from HORIZONS.')
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
        self._hdu = hdu
        self._bpm = hdu['bpm'].data.astype(bool)
        self._cat = Table(hdu['cat'].data)

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
            tab = ascii.read(self.filename, format='ecsv')
            tab.meta = OrderedDict()  # no need to duplicate meta data
            self.tab = vstack((self.tab, tab))

            if self.verbose:
                logger.info('Read {} table from {}.'.format(
                    self._table_title, self.filename))
        else:
            if self.verbose:
                logger.info('Initialized empty {} table.'.format(
                    self._table_title))

    def write(self):
        """Write table to file."""
        if self._table_sort is not None:
            self.tab.sort(self._table_sort)
        self.tab.write(self.filename, overwrite=True, format='ascii.ecsv')

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

class ObservationLog(ScienceTable):
    """Observation log table."""

    _table_title = 'observation log'
    _table_columns = (
        'frame', 'date', 'time', 'binning', 'exptime', 'airmass', 'filter'
    )
    _table_dtypes = ('U64', 'U10', 'U8', 'U3', float, float, 'U2')
    _table_meta = OrderedDict()
    _table_meta['date/time'] = 'UTC.'
    _table_meta['binning'] = 'Original CCD pixels bin factor.'
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

class CometPhotometry(ScienceTable):
    """All comet photometry.

    Parameters
    ----------
    filename : string
      The comet photometry table to read and update.

    """
    _table_title = 'comet photometry'
    _table_columns = [
        'frame', 'filter', 'x', 'y', 'bg', 'bgsig', 'bgarea',
        'f1', 'ferr1', 'f2', 'ferr2', 'f3', 'ferr3',
        'm1', 'merr1', 'm2', 'merr2', 'm3', 'merr3',
    ]
    _table_dtypes = ['U64', 'U2'] + [float] * 29
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name.'
    _table_meta['x/y'] = 'Aperture center, 0-based index, pixels.'
    _table_meta['bg'] = 'Background estimate, ADU/s/pixel.'
    _table_meta['bgsig'] = 'Background standard deviation per pixel.'
    _table_meta['bgarea'] = 'Area used for background estimate.'
    _table_meta['fi, ferri'] = 'Background subtracted flux and error estimates for 1, 2, and 3" radius apertures, ADU/s.'
    _table_meta['mi, merri'] = 'Calibrated magnitudes for each aperture, AB mag.'
    _table_formats = {
        'x': '{:.2f}',
        'y': '{:.2f}',
        'bg': '{:.2f}',
        'bgsig': '{:.2f}',
        'f1': '{:.5g}',
        'f2': '{:.5g}',
        'f3': '{:.5g}',
        'ferr1': '{:.5g}',
        'ferr2': '{:.5g}',
        'ferr3': '{:.5g}',
        'm1': '{:.3f}',
        'm2': '{:.3f}',
        'm3': '{:.3f}',
        'merr1': '{:.3f}',
        'merr2': '{:.3f}',
        'merr3': '{:.3f}',
    }
    _table_sort = 'frame'

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)
