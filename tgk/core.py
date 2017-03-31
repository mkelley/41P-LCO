# Licensed under a MIT style license - see LICENSE
import logging

########################################################################
# Exceptions
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
        import sys
        
        logging.Logger.__init__(self, 'TGK')
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
        'science path': '/full/path/to/your/science/directory',
        'calibrate match radius': 1,
    }

    def __init__(self, config_file, logger=None, log_file=None):
        import os
        
        self.config_file = config_file

        # read configuration file
        self._read_config()

        # verify target directories
        for k in ['download path', 'science path']:
            assert os.path.isdir(self.config[k]), (
                '{} is not a directory or does not exist.'.format(
                    self.config[k]))

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
                except json.JSONDecodeError as e:
                    raise ConfigFileError(
                        'Error reading config file: {}\n{}'.format(
                            self.config_file, e))
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
# for command line parsing
def list_of(type):
    def to_list(s):
        return [type(x) for x in s.split(',')]
    return to_list
