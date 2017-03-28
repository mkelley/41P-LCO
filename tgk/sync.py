# Licensed under a MIT style license - see LICENSE
"""sync - Sync LCO data with local archive.

Thanks to Nestor Espinoza's lcogtDD for an example on syncing with LCO.

"""

from .core import TGKMaster
from .core import ArchiveFileAlreadyExists, AuthorizationError

class Sync(TGKMaster):
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
