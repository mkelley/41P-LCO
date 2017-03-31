# Licensed under a MIT style license - see LICENSE
"""science - Comet science with LCO data."""

import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from .core import TGKMaster

class Science(TGKMaster):
    """Comet science with LCO data.

    Parameters
    ----------
    config_file : string
      The name of a configuration file to read.
    logger : logging.Logger, optional
      Log to this `Logger` instance, otherwise log to the python
      default.
    rlevel : int, optional
      LCO reduction level to consider, or `None` for the most current.

    """

    def __init__(self, config_file, logger=None, rlevel=None):
        TGKMaster.__init__(self, config_file, logger=logger,
                           log_file=('science path', 'tgk-science.log'))

        assert isinstance(rlevel, (int, type(None))), 'rlevel must be integer or `None`.'
        
        self.rlevel = rlevel
        if rlevel is None:
            self.logger.info('Processing the most recent reduction level.')
        else:
            self.logger.info('Processing reduction level {} only.'.format(
                rlevel))

        self.read_observing_log()
        self.read_geometry()
        self.read_processing_history()
        self.find_data()

    def find_data(self):
        """Find comet data."""
        import os
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
            if f[0] not in self.processing_history:
                new_data.append(f)

        self.logger.info('{} frames are new or updated.'.format(len(new_data)))
        return new_data

    def process(self, all_files=False):
        """Run the science pipeline.

        Processing history and observing logs are updated.

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
        for frame, rlevel, filename in files:
            self.logger.info('  ' + frame)

            minion_history = []
            with fits.open(filename) as hdu:
                im = Image(hdu)
                obs = Observation(im.header)

                # only call JPL/HORIZONS if needed
                if frame in self.geometry['frame']:
                    i = np.flatnonzero(self.geometry['frame'] == frame)[0]
                    data = self.geometry[i]
                else:
                    data = None
                geom = Geometry(obs, data=data)

                minion_history.append(minions.frame(self.config, im, obs, geom))
                
                if frame not in self.observing_log['frame']:
                    self.observing_log.add_row(obs.log_row())

                if frame not in self.geometry['frame']:
                    self.geometry.add_row(geom.geometry_row())

            self.processing_history[frame] = (rlevel, minion_history)

        # Next, minions that operate on derived data
        pass

        self.save_observing_log()
        self.save_geometry()
        self.save_processing_history()

    def read_geometry(self):
        """Read the geometry info."""
        import os
        from astropy.io import ascii
        from astropy.table import Table, vstack

        # changes to this table should be reflected in the Geometry class
        self.geometry = Table(
            names=('frame', 'date', 'time', 'ra', 'dec', 'mu ra', 'mu dec',
                   'rh', 'delta', 'phase', 'sun PA', 'velocity PA'),
            dtype=(['U64', 'U10', 'U8'] + [float] * 9))
        self.geometry.meta['comments'] = ['Units: UTC, deg, arcsec/s, au, deg']
        for col in ('ra', 'dec', 'mu ra', 'mu dec'):
            self.geometry[col].format = '{:.6f}'
        for col in ('rh', 'delta'):
            self.geometry[col].format = '{:.4f}'
        for col in ('phase', 'sun PA', 'velocity PA'):
            self.geometry[col].format = '{:.2f}'

        fn = os.sep.join([self.config['science path'], 'geometry.csv'])
        if os.path.exists(fn):
            self.geometry = vstack(
                (self.geometry, ascii.read(fn, format='ecsv')))
            self.logger.info('Read geomtery from {} .'.format(fn))
        else:
            self.logger.info('Geometry file does not exist.'.format(fn))

    def read_observing_log(self):
        """Read the observing log."""
        import os
        from astropy.io import ascii
        from astropy.table import Table, vstack

        # changes to this table should be reflected in the Observation
        # class
        self.observing_log = Table(
            names=('frame', 'date', 'time', 'exptime', 'airmass', 'filter'),
            dtype=('U64', 'U10', 'U8', float, float, 'U2'))
        self.observing_log.meta['comments'] = ['Units: UTC, s, au']
        self.observing_log['exptime'].format = '{:.1f}'
        self.observing_log['airmass'].format = '{:.3f}'

        fn = os.sep.join([self.config['science path'], 'observing-log.csv'])
        if os.path.exists(fn):
            self.observing_log = vstack(
                (self.observing_log, ascii.read(fn, format='ecsv')))
            self.logger.info('Read observing log from {} .'.format(fn))
        else:
            self.logger.info('Observing log {} does not exist.'.format(fn))

    def read_processing_history(self):
        """Read in processing history file."""
        import os
        
        self.processing_history = {}
        
        fn = os.sep.join([self.config['science path'], 'processed-data.txt'])
        if os.path.exists(fn):
            with open(fn, 'r') as inf:
                for line in inf:
                    frame, rlevel, minions = line.split(';')
                    minions = minions.split(',')
                    self.processing_history[frame] = (rlevel, minions)
            self.logger.info('Read processing history from {} .'.format(fn))
        else:
            self.logger.info('Processing history {} does not exist.'.format(fn))

    def save_geometry(self):
        """Save the geometry data."""
        import os

        self.geometry.sort(['date', 'time'])

        fn = os.sep.join([self.config['science path'], 'geometry.csv'])
        self.geometry.write(fn, overwrite=True, delimiter=',',
                                 format='ascii.ecsv')
        self.logger.info('Wrote geometry to {} .'.format(fn))

    def save_observing_log(self):
        """Save the observing log."""
        import os

        self.observing_log.sort(['date', 'time'])

        fn = os.sep.join([self.config['science path'], 'observing-log.csv'])
        self.observing_log.write(fn, overwrite=True, delimiter=',',
                                 format='ascii.ecsv')
        self.logger.info('Wrote observing log to {} .'.format(fn))

    def save_processing_history(self):
        """Save the processing history."""
        import os
        
        fn = os.sep.join([self.config['science path'], 'processed-data.txt'])
        with open(fn, 'w') as outf:
            for frame, (rlevel, minions) in self.processing_history.items():
                outf.write(';'.join([frame, rlevel, ','.join(minions)]) + '\n')

        self.logger.info('Wrote processing history to {} .'.format(fn))

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
    pixel_scale : Quantity
      Nominal pixel scale.
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
    def pixel_scale(self):
        return self.header['PIXSCALE'] * u.arcsec
    
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
        from astropy.coordinates import SkyCoord
        import callhorizons
        from . import lco

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
            return SkyCoord(self._geom['ra'] * u.deg, self._geom['dec'] * u.deg)
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
    phot : astropy.table.Table
      The LCO pipeline photometry table.

    """

    def __init__(self, hdu):
        from astropy.table import Table
        self._hdu = hdu
        self._bpm = hdu['bpm'].data.astype(bool)
        self._phot = Table(hdu['cat'].data)

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
    def phot(self):
        return self._phot
