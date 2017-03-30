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
        self.read_processing_history()
        self.find_data()

        new_data = self.find_new_data()
        if len(new_data) > 0:
            self.process(new_data)
            self.save_observing_log()

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

    def process(self, files):
        from astropy.io import fits
        
        for frame, rlevel, filename in files:
            self.logger.info(frame)
            with fits.open(filename) as hdu:
                obs = Observation(hdu['sci'].header)
                self.observing_log.add_row(obs.log_row())

    def read_observing_log(self):
        """Read the observing log."""
        import os
        from astropy.io import ascii
        from astropy.table import Table, vstack

        # changes to this table should be reflected in the Observation
        # class
        self.observing_log = Table(
            names=('frame', 'date', 'time', 'exptime', 'airmass',
                   'filter', 'rh', 'delta'),
            dtype=('U64', 'U10', 'U8', float, float, 'U2', float, float))
        self.observing_log.meta['comments'] = ['Units: UTC, s, au']
        self.observing_log['exptime'].format = '{:.1f}'
        self.observing_log['airmass'].format = '{:.3f}'
        for col in ('rh', 'delta'):
            self.observing_log[col].format = '{:.4f}'

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

    def save_observing_log(self):
        """Save the observing log."""
        import os

        self.observing_log.sort(['time', 'date'])

        fn = os.sep.join([self.config['science path'], 'observing-log.csv'])
        self.observing_log.write(fn, overwrite=True, delimiter=',',
                                 format='ascii.ecsv')
        self.logger.info('Wrote observing log to {} .'.format(fn))

class Observation:
    """Meta data for each LCO frame.

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

    radec_predict : astropy.coordinates.SkyCoord
      Predectied J2000 coordinates (JPL/HORIZONS).
    mu : Quantity
      Proper motion  (JPL/HORIZONS).
    mu_ra : Qunatity
      Right ascention rate of change times cos(declination) (JPL/HORIZONS).
    mu_dec : Quantity
      Declination rate of change (JPL/HORIZONS).

    rh : Quantity
      Heliocentric distance (JPL/HORIZONS).
    delta : Quantity
      Observer-comet distance (JPL/HORIZONS).
    phase : Angle
      Phase angle (JPL/HORIZONS).

    sun_position_angle : Angle
      Projected comet->Sun vector.
    velocity_position_angle : Angle
      Projected comet velocity vector.

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
        self.get_ephemeris()

    def get_ephemeris(self):
        """Get the comet ephemeris from JPL/HORIZONS."""
        import callhorizons
        from . import lco

        q = callhorizons.query('41P')
        q.set_discreteepochs([self.time.jd])
        obs_code = lco.mpc_codes[(self.site[0], self.enclosure[0], self.telescope[0])]
        n = q.get_ephemerides(obs_code)
        if n != 1:
            raise EphemerisError('Bad return from JPL/HORIZONS, check URL: {}'.format(q.url))
        self._horizons_query = q

    def log_row(self):
        return (self.frame_name, self.time.iso[:10], self.time.iso[11:19],
                self.exptime.value, self.airmass, self.filter, self.rh.value,
                self.delta.value)

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
    def radec_predict(self):
        from astropy.coordinates import SkyCoord
        ra = self._horizons_query['RA'][0] * u.deg
        dec = self._horizons_query['DEC'][0] * u.deg
        return SkyCoord(ra, dec)

    @property
    def mu(self):
        return np.sqrt(self.mu_ra**2 + self.mu_dec**2)

    @property
    def mu_ra(self):
        return self._horizons_query['RA_rate'][0] * u.arcsec / u.s

    @property
    def mu_dec(self):
        return self._horizons_query['DEC_rate'][0] * u.arcsec / u.s

    ############################################################
    @property
    def rh(self):
        return self._horizons_query['r'][0] * u.au

    @property
    def delta(self):
        return self._horizons_query['delta'][0] * u.au

    @property
    def phase(self):
        return Angle(self._horizons_query['alpha'][0] * u.deg)

    ############################################################
    @property
    def sun_position_angle(self):
        return Angle((self._horizons_query['sunTargetPA'][0] + 180) * u.deg)

    @property
    def velocity_position_angle(self):
        return Angle((self._horizons_query['velocityPA'][0] + 180) * u.deg)

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
