# Licensed under a MIT style license - see LICENSE
"""plotcometphot - Plot comet photometry."""

from ..core import config
import matplotlib
matplotlib.use(config['mpl backend'])

from . import TableMinion, MinionError

class PlotCometPhotFailure(MinionError):
    pass

class PlotCometPhot(TableMinion):
    """Plot comet photometry.

    Parameters
    ----------
    config : dict
      Configuration parameters.

    """

    def __init__(self, config):
        TableMinion.__init__(self, config, minion_directory=True)

    name = 'plotcometphot'

    def run(self):
        import logging
        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.time import Time
        from .cometphot import CometPhotometry
        from ..science import GeometryTable

        logger = logging.getLogger('tgk.science')
        logger.debug('  Plot comet photometry.')

        comet = CometPhotometry()
        geom = GeometryTable()


        try:
            date = ['T'.join((geom.get_frame(frame)['date'],
                              geom.get_frame(frame)['time']))
                    for frame in comet.tab['frame']]
        except IndexError as e:
            raise PlotCometPhotFailure(e)  # missing geometry
        
        date = np.array(date)
        tmtp = Time(date, scale='utc').jd - 2457856.252
        rh = [geom.get_frame(frame)['rh'] for frame in comet.tab['frame']]
        delta = [geom.get_frame(frame)['delta'] for frame in comet.tab['frame']]
        d_correction = 2.5 * np.log10(delta)
        opts = dict(ls='none', marker='.', ecolor='0.5', alpha=0.5)

        ######################################################################
        fig = plt.figure(figsize=(8, 8))
        fig.clear()
        axes = [fig.add_subplot(gs) for gs in
                plt.GridSpec(2, 1, wspace=0, hspace=0)]

        i = comet.tab['filter'] == 'rp'
        if i.sum() != 0:
            for ap, color in zip('26', 'rb'):
                t = tmtp[i]
                m = comet.tab['m' + ap][i]
                merr = comet.tab['merr' + ap][i]
                kwargs = dict(color=color, label=ap + '" radius', **opts)
                axes[0].errorbar(t, m, merr, **kwargs)
                axes[1].errorbar(t, m - d_correction[i], merr, **kwargs)

            plt.setp(axes[0], ylabel=r"$r'$ (mag)", ylim=axes[0].get_ylim()[::-1])
            plt.setp(axes[1], xlabel=r'$T-T_p$ (days)',
                     ylim=axes[1].get_ylim()[::-1],
                     ylabel=r"$r' - 2.5 \log{\Delta}$ (mag)")
            axes[0].legend(numpoints=1, prop=dict(size='medium'))

            fig.canvas.draw()
            fig.savefig(self.minion_file('rp-vs-rh.png'), dpi=200)
            fig.savefig(self.minion_file('rp-vs-rh.pdf'), dpi=200)

        ######################################################################
        fig.clear()
        axes = [fig.add_subplot(gs) for gs in
                plt.GridSpec(2, 1, wspace=0, hspace=0)]

        i = comet.tab['filter'] == 'gp'
        if i.sum() != 0:
            for ap, color in zip('26', 'rb'):
                t = tmtp[i]
                m = comet.tab['m' + ap][i]
                merr = comet.tab['merr' + ap][i]
                kwargs = dict(color=color, label=ap + '" radius', **opts)
                axes[0].errorbar(t, m, merr, **kwargs)
                axes[1].errorbar(t, m - d_correction[i], merr, **kwargs)

            plt.setp(axes[0], ylabel=r"$g'$ (mag)", ylim=axes[0].get_ylim()[::-1])
            plt.setp(axes[1], xlabel=r'$T-T_p$ (days)',
                     ylim=axes[1].get_ylim()[::-1],
                     ylabel=r"$g' - 2.5 \log{\Delta}$ (mag)")

            axes[0].legend(numpoints=1, prop=dict(size='medium'))

            fig.canvas.draw()
            fig.savefig(self.minion_file('gp-vs-rh.png'), dpi=200)
            fig.savefig(self.minion_file('gp-vs-rh.pdf'), dpi=200)
            plt.close(fig)
