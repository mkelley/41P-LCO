# Licensed under a MIT style license - see LICENSE
import os
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
def setup_logger(logging_name, debug=False):
    import os
    import sys
    
    logger = logging.getLogger(logging_name)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    formatter = logging.Formatter('%(levelname)s: %(message)s')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    logger.addHandler(console)

def open_log_file(log_file, logger_name):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    logfile = logging.FileHandler(log_file)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

    logger.info('#' * 73)
    logger.info(timestamp())
    logger.info('Logging to {}'.format(log_file))

def log_table(tab, logger_name):
    import io
    logger = logging.getLogger(logger_name)
    with io.StringIO() as s:
        tab.write(s, format='ascii.fixed_width_two_line')
        s.seek(0)
        logger.info("\n" + s.read(-1))

def shutdown_logging(logger_name):
    logger = logging.getLogger(logger_name)
    logger.info('End logging: ' + timestamp())
    logging.shutdown()

def timestamp():
    from datetime import datetime
    return datetime.now().isoformat()

########################################################################
config_defaults = {
    'username': 'your@email',
    'password': 'your_password',
    'proposal': 'LCO2016B-109',
    'download path': '/full/path/to/your/download/directory',
    'science path': '/full/path/to/your/science/directory',
    'calibrate match radius': 1.5,
    'mpl backend': 'PDF',
    'post-science hook': ""
}

def show_config(defaults=False):
    import os
    import sys
    import json

    if defaults:
        print(json.dumps(config_defaults, indent=2))
    elif os.path.exists(config_file):
        with open(config_file) as inf:
            config = json.load(inf)
        print(json.dumps(config, indent=2))
    else:
        raise FileNotFoundError('Config file does not exist: {}'.format(config_file))

def configure():
    """Read configuration file."""
    import os
    import json

    if os.path.exists(config_file):
        with open(config_file) as inf:
            try:
                config = json.load(inf)
            except json.JSONDecodeError as e:
                raise ConfigFileError(
                    'Error reading config file: {}\n{}'.format(
                        config_file, e))
    else:
        raise FileNotFoundError("""Configuration file not found: {}
Use --show-config for an example.""".format(config_file))

    # Verify all keys are present
    missing = [k for k in config_defaults.keys()
               if k not in config.keys()]
    if len(missing) > 0:
        raise ConfigFileError('Missing {} from config file.  See\n --show-config for examples.'.format(missing))

    # verify a couple directories
    for k in ['download path', 'science path']:
        assert os.path.isdir(config[k]), (
            '{} is not a directory or does not exist.'.format(
                config[k]))

    return config

########################################################################
def get_filename(frame, rlevel=None):
    """Get the full path to a single frame.

    Parameters
    ----------
    frame : string
      The frame name (i.e., everything up to -e*.fits.fz).
    rlevel : string, optional
      The reduction level desired, nor `None` to return the highest
      level.

    Returns
    -------
    filename : string

    """

    from glob import glob

    if rlevel is None:
        pat = os.sep.join([config['download path'], 'e*', '*', frame + '*'])
    else:
        pat = os.sep.join([config['download path'], 'e{}'.format(rlevel),
                           '*', frame + '*'])

    # for multiple matches, use sort order to determine highest rlevel
    return sorted(glob(pat))[-1]

########################################################################
# for command line parsing
def list_of(type):
    def to_list(s):
        return [type(x) for x in s.split(',')]
    return to_list

########################################################################

config_file = os.sep.join([os.path.expanduser('~'), '.config',
                           '41p-lco', 'config.json'])
config = configure()
