# Licensed under a MIT style license - see LICENSE
"""minions - Science minions execute pipeline tasks."""

class Minion:
    def __init__(self, config):
        import os
        import logging

        self.config = config
        self.logger = logging.getLogger('TGK')

        d = self.minion_file('')
        if not os.path.exists(d):
            os.mkdir(d)

        assert os.path.isdir(d), 'Expected {} to be a directory.'.format(d)

    def run(self):
        pass

    @property
    def name(self):
        pass

    def minion_file(self, basename):
        import os
        fn = os.sep.join([self.config['science path'], self.name, basename])
        return fn

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
    
    from .calibrate import Calibrate

    history = []

    for minion in (Calibrate,):
        m = minion(config, im, obs, geom)
        m.run()
        history.append(m.name)

    return history
