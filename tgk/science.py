# Licensed under a MIT style license - see LICENSE
"""science - Comet science with LCO data."""

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
        
        self.find_data()
        self.find_new_data()
        self.process(self.new_data)

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
        self.logger.info('{} files match requested reduction level.'.format(
            len(all_files)))

        if self.rlevel is None:
            files = []
            multi_versions = 0
            i = 0
            while i < len(all_files):
                # remove e??.fits.fz
                basename = os.path.basename(all_files[i]).rstrip('.fits.fz')[:-4]
                m = [s for s in all_files if basename in s]
                # Thanks to sort order and the two-digit pipeline
                # suffix, the following will always be the highest
                # level pipeline version:
                files.append(m[-1])
                multi_versions += len(m) > 1
                i += len(m)

            all_files = files
            self.logger.info('{} frames have multiple reduction levels.'.format(multi_versions))

        self.logger.info('{} files will be processed.'.format(len(all_files)))

    def find_new_data(self):
        self.new_data = []
        pass

    def process(self, files):
        pass
