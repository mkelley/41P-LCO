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

        self.read_processing_history()
        self.find_data()
        new_data = self.find_new_data()
        self.process(new_data)

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
                found.append(files[i][1])
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
        pass

    def read_processing_history(self):
        """Read in processing history file."""
        import os
        
        self.processing_history = {}
        
        f = os.sep.join([self.config['science path'], 'processed-data.txt'])
        if os.path.exists(f):
            with open(f, 'r') as inf:
                for line in inf:
                    frame, rlevel, minions = line.split(';')
                    minions = minions.split(',')
                    self.processing_history[frame] = (rlevel, minions)
