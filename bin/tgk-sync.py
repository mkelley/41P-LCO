
class ConfigFileError(Exception):
    pass

class TGKSync:
    """Sync observations with LCO.

    """

    _defaults = {
        'username': 'your@email',
        'password': 'your_password',
        'proposals': 'LCO2016B-109',
        'download path': '/full/path/to/your/download/directory'
    }
    

    def __init__(self, config_file):
        import json
        
        self.config_file = config_file
        self.last_download = None

        try:
            self._read_config()
        except json.JSONDecodeError as e:
            raise ConfigFileError('Error reading config file: {}\n{}'.format(config_file, e))
        
        self._setup_logger()

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
        import os
        import json

        if os.path.exists(self.config_file):
            with open(self.config_file) as inf:
                self.config = json.load(inf)
        else:
            raise FileNotFoundError("""Configuration file not found: {}
Use --show-config for an example.""".format(self.config_file))

        # verify download path
        if not os.path.isdir(self.config['download path']):
            raise OSError('Download path does not exist: {}'.format(self.config['download path']))

    def _setup_logger(self):
        import os
        import logging
        from datetime import datetime
        
        fn = os.sep.join([self.config['download path'], 'tgk-sync.log'])
        
        self.logger = logging.Logger('TGK Sync')
        self.logger.setLevel(logging.DEBUG)
        if len(self.logger.handlers) == 0:
            formatter = logging.Formatter('%(levelname)s: %(message)s')
            console = logging.StreamHandler(sys.stdout)
            console.setLevel(logging.DEBUG)
            console.setFormatter(formatter)
            self.logger.addHandler(console)

            logfile = logging.FileHandler(fn)
            logfile.setLevel(logging.INFO)
            logfile.setFormatter(formatter)
            self.logger.addHandler(logfile)

        self.logger.info('#' * 73)
        self.logger.info(datetime.now().isoformat())
        self.logger.info('Logging to console and {}'.format(fn))
        self.logger.info('Loaded configuration file: {}'.format(
            self.config_file))
        
    def run(self):
        from datetime import datetime
        from astropy.time import Time
        import astropy.units as u

        now = Time(datetime.now())
        
    def sync(self, **kwargs):
        import requests

if __name__ == '__main__':
    import os
    import sys
    import argparse

    default_config = os.sep.join([os.path.expanduser('~'), '.config',
                                  '41p-lco', 'sync.cfg'])
    
    parser = argparse.ArgumentParser(description='Master control program for monitoring LCO images of 41P.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--config', default=default_config, help='Use this configuration file.')
    parser.add_argument('--show-config', action='store_true', help='Read and print the configuration file.')
    parser.add_argument('--show-defaults', action='store_true', help='Print the default configuration file.')

    args = parser.parse_args()

    if args.show_config:
        TGKSync.show_config(args.config)
        sys.exit()

    if args.show_defaults:
        TGKSync.show_config(None)
        sys.exit()

    try:
        sync = TGKSync(args.config)
        sync.run()
    except Exception as e:
        print(e)
        exit()

