# Licensed under a MIT style license - see LICENSE
"""plotcometloc - Show estimated comet position."""

from ..core import config
import matplotlib
matplotlib.use(config['mpl backend'])

from . import FrameMinion, MinionError
from .cometphot import CometPhotometry

class PlotCometLocFailure(MinionError):
    pass

class PlotCometLoc(FrameMinion):
    """Plot the estimated comet positions.

    Requires CometPhot.

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

    def __init__(self, *args):
        FrameMinion.__init__(self, *args, minion_directory=True)
    
    name = 'plotcometloc'

    def run(self):
        import os
        import logging
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.colors import SymLogNorm
        from astropy.visualization import ZScaleInterval, ImageNormalize
        from astropy.wcs.utils import skycoord_to_pixel

        logger = logging.getLogger('tgk.science')
        logger.info('    Plot comet location.')
        try:
            comet = CometPhotometry().get_frame(self.obs.frame_name)
        except IndexError as e:
            raise PlotCometLocFailure(e)

        yc, xc = int(comet['y']), int(comet['x'])

        fig = plt.figure(figsize=(8, 8))
        fig.clear()
        axes = [fig.add_subplot(gs)
                for gs in plt.GridSpec(2, 2, wspace=0, hspace=0,
                                       bottom=0, top=1, left=0, right=1)]

        opts = dict(origin='lower', cmap='viridis',
                    norm=ImageNormalize(self.im.data[self.obs.trimsec],
                                        interval=ZScaleInterval()))
        axes[0].imshow(self.im.data, **opts)
        axes[2].imshow(self.im.data, **opts)

        im = self.im.data - comet['bg']
        opts['norm'] = SymLogNorm(comet['bgsig'], vmin=-comet['bgsig'],
                                  vmax=im[yc, xc])
        axes[1].imshow(im, **opts)
        axes[3].imshow(im, **opts)

        xjpl, yjpl = skycoord_to_pixel(self.geom.radec_predict, self.obs.wcs)
        opts = dict(color='k', linewidths=1)
        for i in range(4):
            axes[i].scatter(self.im.cat['X'], self.im.cat['Y'], s=24,
                            marker='o', color='w', linewidth=0.5, label=None,
                            facecolor='none')
            axes[i].scatter(self.im.cat['X'], self.im.cat['Y'], s=20,
                            marker='o', color='k', linewidth=0.5, label=None,
                            facecolor='none')
            axes[i].scatter([xjpl], [yjpl], marker='+', label='JPL/HORIZONS',
                            **opts)
            axes[i].scatter([xc], [yc], marker='x', label='Centroid',
                            **opts)

        axes[0].legend(numpoints=1, prop=dict(size='medium'), loc='upper left')
        plt.setp(axes[2:], xlim=[xc - 100, xc + 100], ylim=[yc - 100, yc + 100])
        plt.setp(axes[:2], xlim=[0, self.im.data.shape[1] - 1],
                 ylim=[0, self.im.data.shape[0] - 1])
        plt.setp(axes, frame_on=False, xticks=[], yticks=[])
        fig.canvas.draw()

        date = self.obs.time.iso[:10].replace('-', '')
        fn = self.minion_file('{}/{}.png'.format(date, self.obs.frame_name))
        d = os.path.split(fn)[0]
        if not os.path.exists(d):
            os.mkdir(d)
            self.logger.info('Created directory {}.'.format(d))
            
        fig.savefig(fn, dpi=75)
        plt.close(fig)
