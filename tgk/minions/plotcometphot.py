# Licensed under a MIT style license - see LICENSE
"""plotcometphot - Plot comet photometry."""

from . import TableMinion, MinionError

class PlotCometPhotFailure(MinionError):
    pass

class PlotCometPhot(TableMinion):
    """Plot comet photometry.

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

    def __init__(self, config):
        TableMinion.__init__(self, config, minion_directory=True)

    @property
    def name(self):
        return 'plotcometphot'

    def run(self):
        import logging
        import numpy as np
        import matplotlib.pyplot as plt
        from .cometphot import CometPhotometry
        from ..science import GeometryTable

        logger = logging.getLogger('tgk.science')
        logger.info('    Plot comet photometry.')

        comet = CometPhotometry()
        geom = GeometryTable()

        fig = plt.figure()
        fig.clear()
        ax = fig.add_subplot(111)

        date = ['T'.join((geom.get_frame(frame)['date'],
                          geom.get_frame(frame)['time']))
                for frame in comet.tab['frame']]
        date = np.array(date)
        
        rh = [geom.get_frame(frame)['rh'] for frame in comet.tab['frame']]
        rh = np.array(rh)
        
        rh[date < '2017-04-12T00:45:07'] *= -1
        
        i = comet.tab['filter'] == 'rp'
        opts = dict(ls='none', marker='.', ecolor='0.5')
        ax.errorbar(rh[i], comet.tab['m2'][i], comet.tab['merr2'][i],
                    color='b', label='2" radius', **opts)
        ax.errorbar(rh[i], comet.tab['m3'][i], comet.tab['merr3'][i],
                    color='r', label='3" radius', **opts)

        plt.setp(ax, xlabel=r'$r_h$ (au)', ylabel=r"$r'$ (mag)",
                 ylim=ax.get_ylim()[::-1])
        plt.legend(numpoints=1, prop=dict(size='medium'))
        
        fig.canvas.draw()
        fig.savefig(self.minion_file('m-vs-rh.png'), dpi=200)
        fig.savefig(self.minion_file('m-vs-rh.pdf'), dpi=200)
        plt.close(fig)
