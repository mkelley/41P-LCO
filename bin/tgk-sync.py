# Licensed under a MIT style license - see LICENSE
"""tgk-sync/science - Master control programs for the 41P LCO Outburst project.

Thanks to Nestor Espinoza's lcogtDD for an example on syncing with the
LCO archive.

"""

import logging

class ConfigFileError(Exception):
    pass

class AuthorizationError(Exception):
    pass

class ArchiveFileAlreadyExists(Exception):
    pass

########################################################################
class Logger(logging.Logger):
    def __init__(self, debug=False):
        import os
        
        logging.Logger.__init__(self, 'TGK Sync')
        self.loglevel = logging.DEBUG if debug else logging.INFO
        self.setLevel(self.loglevel)

        formatter = logging.Formatter('%(levelname)s: %(message)s')
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(self.loglevel)
        console.setFormatter(formatter)
        self.addHandler(console)

    def open_file(self, log_file):
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        logfile = logging.FileHandler(log_file)
        logfile.setLevel(self.loglevel)
        logfile.setFormatter(formatter)
        self.addHandler(logfile)
        
        self.info('#' * 73)
        self.timestamp()
        self.info('Logging to {}'.format(log_file))

    def log_table(self, tab):
        import io
        with io.StringIO() as s:
            tab.write(s, format='ascii.fixed_width_two_line')
            s.seek(0)
            self.info("\n" + s.read(-1))
        
    def shutdown(self):
        self.info('End logging.')
        self.timestamp()
        logging.shutdown()

    def timestamp(self):
        from datetime import datetime
        self.info(datetime.now().isoformat())

########################################################################
class TGKMaster:
    """Master class for TGK programs.

    Parameters
    ----------
    config_file : string
      The name of a configuration file to read.
    logger : logging.Logger, optional
      Log to this `Logger` instance, otherwise log to the python
      default.
    log_file : string or tuple, optional
      Name of a file to log at, or the config parameter and a file
      name to join together.  Requires `logger`.

    """

    _defaults = {
        'username': 'your@email',
        'password': 'your_password',
        'proposal': 'LCO2016B-109',
        'download path': '/full/path/to/your/download/directory',
        'science path': '/full/path/to/your/science/directory'
    }

    def __init__(self, config_file, logger=None, log_file=None):
        self.config_file = config_file

        # read configuration file
        self._read_config()

        # begin logging
        if logger is None:
            self.logger = logging
        else:
            self.logger = logger
            if log_file is not None:
                assert isinstance(log_file, (tuple, str))
                if isinstance(log_file, tuple):
                    fn = os.sep.join((self.config[log_file[0]], log_file[1]))
                else:
                    fn = log_file[1]
                    
                self.logger.open_file(fn)

        self.logger.info('Loaded configuration file: {}'.format(
            self.config_file))

    @classmethod
    def show_config(cls, config_file):
        import os
        import sys
        import json

        if config_file is None:
            print(json.dumps(cls._defaults, indent=2))
        elif os.path.exists(config_file):
            with open(config_file) as inf:
                config = json.load(inf)
            print(json.dumps(config, indent=2))
        else:
            raise FileNotFoundError('Config file does not exist: {}'.format(config_file))

    def _read_config(self):
        """Read configuration file."""
        import os
        import json

        if os.path.exists(self.config_file):
            with open(self.config_file) as inf:
                try:
                    self.config = json.load(inf)
                except JSONDecodeError as e:
                    raise ConfigFileError(
                        'Error reading config file: {}\n{}'.format(
                            config_file, e))
        else:
            raise FileNotFoundError("""Configuration file not found: {}
Use --show-config for an example.""".format(self.config_file))

        # Verify all keys are present
        missing = [k for k in self._defaults.keys()
                   if k not in self.config.keys()]
        if len(missing) > 0:
            raise ConfigFileError('Missing {} from config file.  See\n --show-config for examples.'.format(missing))

        # verify download path
        if not os.path.isdir(self.config['download path']):
            raise OSError('Download path does not exist: {}'.format(self.config['download path']))

########################################################################
class TGKSync(TGKMaster):
    """Sync observations with LCO.

    Parameters
    ----------
    config_file : string
      The name of a configuration file to read.
    logger : logging.Logger, optional
      Log to this `Logger` instance, otherwise log to the python
      default.

    """

    def __init__(self, config_file, logger=None):
        import astropy.units as u
        from astropy.time import Time

        self.request_delay = 1 * u.s
        self.last_request = Time.now() - self.request_delay
        self.last_download = None

        TGKMaster.__init__(self, config_file, logger=logger,
                           log_file=('download path', 'tgk-sync.log'))

        # get http authorization token from LCO
        self._get_auth_token()

    def _download_frame(self, meta):
        """Download a frame described by metadata from an LCO frame payload.

        Target location:
          download_path/e(rlevel)/UTC_date/filename

        If the file already exists, the download is skipped.

        Raises
        ------
        ArchiveFileAlreadyExists

        """
        import os
        import requests
        
        # first, verify target path and create subdirectories if needed
        d = self.config['download path']
        for tail in ('e{}'.format(meta['RLEVEL']),
                     meta['DATE_OBS'][:10].replace('-', '')):
            d = os.sep.join([d, tail])
            if os.path.exists(d):
                assert os.path.isdir(d), \
                    '{} exists, but is not a directory'.format(d)
            else:
                os.mkdir(d)

        # archive file name format:
        # (site)(tel)-(instr)-YYYYMMDD-(frame)-(type)(red.level).fits
        filename = os.sep.join([d, meta['filename']])

        if os.path.exists(filename):
            self.logger.debug(
                '{} already exists, skipping download.'.format(filename))
            raise ArchiveFileAlreadyExists(filename)
        else:
            self.logger.info('Downloading to {}.'.format(filename))
            with open(filename, 'wb') as outf:
                outf.write(requests.get(meta['url']).content)
        
    def _get_auth_token(self):
        import requests

        data = {
            'username': self.config['username'],
            'password': self.config['password']
        }
        r = requests.post('https://archive-api.lco.global/api-token-auth/',
                          data=data).json()
        token = r.get('token')
        if token is None:
            raise AuthorizationError('Authorization token not returned.')

        self.auth = {'Authorization': 'Token ' + token}
        self.logger.info('Obtained authoriation token.')

    def continuous_sync(self, rlevels=[91], download=True):
        """Continuously check for new data.

        Parameters
        ----------
        rlevels : list of int, optional
          Reduction levels to check.
        download : bool, optional
          Download flag.

        """
        
        from astropy.time import Time
        import astropy.units as u
        import time

        window = 1 * u.day  # search window
        cadence = 2 * u.hr  # search cadence
        last_sync = Time('2000-01-01')
        self.logger.timestamp()
        self.logger.info('Entering continuous sync mode, checking LCO archive every {} with a {}-wide search window.'.format(cadence, window))

        try:
            while True:
                now = Time.now()
                if (now - last_sync) > cadence:
                    self.logger.timestamp()
                    self.logger.info('Sync with LCO archive.')
                    self.sync(now - window, rlevels=rlevels, download=download)
                    last_sync = Time.now()
                else:
                    self.logger.debug('Last sync: {:.0f} s ago.  Sleep 60 s.'.format((now - last_sync).sec))
                    time.sleep(60)
                    self.logger.debug('Awake!')
        except KeyboardInterrupt:
            self.logger.info('Caught interrupt signal.  Shutdown.')

    def request(self, url, query={}):
        """Send HTTP request and return the output.

        Limits the overall number of requests to slow down any runaway
        code and prevent from exceeding the request limit.

        Parameters
        ----------
        url : string
          The URL.
        param : dict, optional
          The HTTP get parameters.

        """
        import time
        from astropy.time import Time
        import requests

        while (Time.now() - self.last_request) < self.request_delay:
            time.sleep(1)

        self.logger.info('Request: {}, {}'.format(url, query))
        response = requests.get(url, params=query, headers=self.auth)
        self.logger.debug(response.url)

        data = response.json()
        return data

    def _summarize_payload(self, payload):
        """Summarize payload as a table."""
        from astropy.table import Table

        tab = Table(names=('filename', 'date_obs', 'filter', 'exptime'),
                    dtype=('U64', 'U32', 'U16', float))
        for meta in payload:
            tab.add_row((meta['filename'], meta['DATE_OBS'], meta['FILTER'],
                         float(meta['EXPTIME'])))
        return tab

    def sync(self, start, end=None, rlevels=[91], download=True):
        """Request frames list from LCO and download, if needed.

        Only whole days are checked.

        Parameters
        ----------
        start : Time
          Check for frames since `start`.
        end : Time, optional
          Check for frames no later than `end`.
        rlevels : list of int, optional
          Which reduction levels to check.
        download : bool, optional
          Flag to download data.

        """
        import os
        from astropy.time import Time

        for rlevel in rlevels:
            query = {
                'PROPID': self.config['proposal'],
                'limit': 50,
                'RLEVEL': rlevel,
                'start': start.iso[:10],
            }
            if end is not None:
                query['end'] = end.iso[:10]

            self.logger.timestamp()
            data = self.request('https://archive-api.lco.global/frames/',
                                query=query)
            self.logger.info('Found {} frames with reduction level {}.'.format(
                data['count'], rlevel))

            dl_count = 0
            skip_count = 0
            while True:  # loop over all payload sets
                payload = data['results']

                if data['count'] > 0:
                    tab = self._summarize_payload(payload)
                    if download:
                        downloaded = []
                        for i, meta in enumerate(payload):
                            try:
                                self._download_frame(meta)
                                dl_count += 1
                                downloaded.append(i)
                            except ArchiveFileAlreadyExists:
                                skip_count += 1
                                pass
                        if len(downloaded) > 0:
                            self.logger.log_table(tab[downloaded])
                    else:
                        self.logger.log_table(tab)

                if data['next'] is not None:
                    # get next payload set
                    data = self.request(data['next'])
                else:
                    break  # end while loop

            self.logger.timestamp()
            self.logger.info('Downloaded {} files, {} skipped.'.format(
                dl_count, skip_count))

def list_of(type):
    def to_list(s):
        return [type(x) for x in s.split(',')]
    return to_list
        
########################################################################
if __name__ == '__main__':
    import os
    import sys
    import argparse
    import astropy.units as u
    from astropy.time import Time

    default_config = os.sep.join([os.path.expanduser('~'), '.config',
                                  '41p-lco', 'config.json'])

    parser = argparse.ArgumentParser(description='Sync with LCO archive.')
    parser.add_argument('--no-download', dest='download', action='store_false', help='Check the archive, do not download any data.')
    parser.add_argument('--start', type=Time, default=Time.now() - 1 * u.day, help='Search for files taken on or after this date (UTC). (default: yesterday)')
    parser.add_argument('--end', type=Time, help='Search for files before this datem (UTC). (default: None)')
    parser.add_argument('--rlevels', type=list_of(int), default=[91], help='Check for frames with these reduction levels. (default: 91)')
    parser.add_argument('--continuous', action='store_true', help='Continously check LCO for new data.')
    parser.add_argument('--config', default=default_config, help='Use this configuration file. (default: {})'.format(default_config))
    parser.add_argument('--show-config', action='store_true', help='Read and print the configuration file.')
    parser.add_argument('--show-defaults', action='store_true', help='Print the default configuration file.')
    parser.add_argument('-v', action='store_true', help='Increase verbosity.')

    args = parser.parse_args()

    if args.show_config:
        TGKSync.show_config(args.config)
        sys.exit()

    if args.show_defaults:
        TGKSync.show_config(None)
        sys.exit()

    logger = Logger(debug=args.v)

    try:
        sync = TGKSync(args.config, logger=logger)
        if args.continuous:
            sync.continuous_sync(rlevels=args.rlevels, download=args.download)
        else:
            sync.sync(args.start, end=args.end, rlevels=args.rlevels,
                      download=args.download)
        logger.shutdown()
    except Exception as e:
        err = '{}: {}'.format(type(e).__name__, e)
        logger.error(err)
        logger.shutdown()

        if args.v:
            raise(e)
        else:
            print(err)
            sys.exit()

