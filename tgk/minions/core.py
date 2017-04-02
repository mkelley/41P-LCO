# Licensed under a MIT style license - see LICENSE
"""minions.core."""

class MinionError(Exception):
    pass

class Minion:
    def __init__(self, config, minion_directory=False):
        import os
        import logging

        self.config = config
        self.logger = logging.getLogger('tgk.science')

        self.minion_directory = minion_directory
        if minion_directory:
            d = self.minion_file('')
            if not os.path.exists(d):
                os.mkdir(d)

            assert os.path.isdir(d), 'Expected {} to be a directory.'.format(d)

    def run(self):
        pass

    # The name of the minion, used for creating files and reprocessing
    name = ''

    def minion_file(self, basename):
        """Get the name of a minion file."""
        import os
        if self.minion_directory:
            fn = os.sep.join([self.config['science path'], self.name, basename])
        else:
            fn = os.sep.join([self.config['science path'], basename])
        return fn

class FrameMinion(Minion):
    def __init__(self, config, im, obs, geom, minion_directory=False):
        self.config = config
        self.im = im
        self.obs = obs
        self.geom = geom
        
        Minion.__init__(self, config, minion_directory=minion_directory)

class TableMinion(Minion):
    pass

def frame(config, im, obs, geom, reprocess=[]):
    """Run minions on an individual frames.

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
    reprocess : list, optional
      Only run these minions (implies their requirements have already
      been met).

    Returns
    -------
    history : list of strings
      Names of the executed minions.

    """

    import logging
    from . import frame_minions
    from ..core import timestamp

    if len(reprocess) == 0:
        reprocess = [m.name for m in frame_minons]

    logger = logging.getLogger('tgk.science')
    history = []
    for minion in frame_minions:
        try:
            if minion.name not in reprocess:
                continue
            m = minion(config, im, obs, geom)
            m.run()
            history.append(m.name)
        except MinionError as e:
            err = '   {} {}: {}'.format(timestamp()[:-7], type(e).__name__, e)
            logger.error(err)
            break

    return history

def table(config, reprocess=[]):
    """Table minions run on tables after all frames have been processed.

    Parameters
    ----------
    config : dict
      Configuration parameters.

    Returns
    -------
    history : list of strings
      Names of the executed minions.
    reprocess : list, optional
      Only run these minions (implies their requirements have already
      been met).

    """

    import logging
    from . import table_minions
    from ..core import timestamp

    if len(reprocess) == 0:
        reprocess = [m.name for m in table_minions]
    
    logger = logging.getLogger('tgk.science')
    history = []
    for minion in table_minions:
        try:
            if minion.name not in reprocess:
                continue
            m = minion(config)
            m.run()
            history.append(m.name)
        except MinionError as e:
            err = '   {} {}: {}'.format(timestamp()[:-7], type(e).__name__, e)
            logger.error(err)
            break

    return history
