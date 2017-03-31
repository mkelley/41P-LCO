# Licensed under a MIT style license - see LICENSE
"""minions - Science minions execute pipeline tasks."""

class MinionError(Exception):
    pass

class Minion:
    def __init__(self, config):
        import os
        import logging

        self.config = config
        self.logger = logging.getLogger('tgk.science')

        d = self.minion_file('')
        if not os.path.exists(d):
            os.mkdir(d)

        assert os.path.isdir(d), 'Expected {} to be a directory.'.format(d)

    def run(self):
        pass

    @property
    def name(self):
        """The name of the minion, used for creating files."""
        pass

    def minion_file(self, basename):
        """Get the name of a minion file."""
        import os
        fn = os.sep.join([self.config['science path'], self.name, basename])
        return fn

    def append_to_table(self, row, basename=None, unique=None):
        """Append a row to a minion data table.

        Parameters
        ----------
        row : array-like
          The row to append.
        basename : string, optional
          The table file name.  The default is `self.name + '.csv'`.
        unique : string, optional
          This row field must be unique in the table.  If a match is
          found with `row`, the table entry will be replaced.

        """

        import os
        from collections import OrderedDict
        from astropy.io import ascii
        from astropy.table import Table, vstack

        if basename is None:
            basename = self.name + '.csv'
        
        try:
            tab = Table(names=self._table_names, dtype=self._table_dtype)
            for col, cformat in zip(self._table_names, self._table_format):
                tab[col].format = cformat
            tab.meta = self._table_meta
            assert isinstance(self._table_sort, (type(None), list, tuple))
        except NameError as e:
            raise MinionError('Cannot create minion table; missing definitions.')

        fn = self.minion_file(basename)
        if os.path.exists(fn):
            _tab = ascii.read(fn, format='ecsv')
            _tab.meta = OrderedDict()  # no need to duplicate meta data
            tab = vstack((tab, _tab))

        # enforce unique fields
        if unique is not None:
            i = tab.index_column[unique]
            j = tab[unique] == row[i]
            if any(j):
                tab.remove_row(np.flatnonzero(j)[0])

        tab.add_row(row)
        if self._table_sort is not None:
            tab.sort(self._table_sort)
        tab.write(fn, overwrite=True, delimiter=',', format='ascii.ecsv')

class FrameMinion(Minion):
    def __init__(self, config, im, obs, geom):
        self.config = config
        self.im = im
        self.obs = obs
        self.geom = geom
        
        Minion.__init__(self, config)

def frame(config, im, obs, geom):
    """Run minions on an individual frame.

    Parameters
    ----------
    config : dict
      Configuration parameters.
    im : Image
      Frame data.
    obs : Observation
      Frame meta data.
    geom : Geometry
      Comet geometric circumstances.

    Returns
    -------
    history : list of strings
      Names of the executed minions.

    """

    import logging
    from .calibrate import Calibrate, CalibrationFailure

    logger = logging.getLogger('tgk.science')
    history = []
    for minion in (Calibrate,):
        try:
            m = minion(config, im, obs, geom)
            m.run()
            history.append(m.name)
        except CalibrationFailure as e:
            err = '{}: {}'.format(type(e).__name__, e)
            logger.error(err)
            return history

    return history
