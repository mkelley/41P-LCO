# Licensed under a MIT style license - see LICENSE
"""calcomet - Calibrate comet photometry."""

from . import FrameMinion, MinionError

class CalCometFailure(MinionError):
    pass

class CalComet(FrameMinion):
    """Calibrate comet photometry.

    Requires image calibration.

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

    """

    name = 'calcomet'

    def run(self):
        import logging
        import numpy as np
        import astropy.units as u
        from .cometphot import CometPhotometry
        from .calibrate import CalibrationTable

        logger = logging.getLogger('tgk.science')
        logger.info('    Calibrate comet photometry.')

        try:
            cal = CalibrationTable().get_frame(self.obs.frame_name)
            comet = CometPhotometry().get_frame(self.obs.frame_name)
        except IndexError as e:
            raise CalCometFailure(e)

        for rap in [1, 2, 3]:
            f = 'f{}'.format(rap)
            ferr = 'ferr{}'.format(rap)
            m = 'm{}'.format(rap)
            merr = 'merr{}'.format(rap)
            comet[m] = -2.5 * np.log10(comet[f]) + cal['scmean(dm)']
            comet[merr] = np.sqrt((1.0857 * comet[ferr] / comet[f])**2
                                  + cal['scstdev(dm)']**2)
            
        CometPhotometry().update(comet)
