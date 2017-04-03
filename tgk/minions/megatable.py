# Licensed under a MIT style license - see LICENSE
"""megatable - Join all frame results into one table."""

from . import TableMinion, MinionError

class MegaTableFailure(MinionError):
    pass

class MegaTable(TableMinion):
    """Join frame science tables.

    Parameters
    ----------
    config : dict
      Configuration parameters.

    """

    name = 'megatable'

    def run(self):
        import logging
        from astropy.table import join
        from .cometphot import CometPhotometry
        from .background import BackgroundTable
        from .calibrate import CalibrationTable
        from ..science import GeometryTable, ProcessingHistory, ObservationLog

        logger = logging.getLogger('tgk.science')
        logger.info('  Create mega-table.')

        try:
            obslog = ObservationLog()
            hist = ProcessingHistory()
            geom = GeometryTable()
            cal = CalibrationTable()
            background = BackgroundTable()
            phot = CometPhotometry()

            tab = join(obslog.tab, hist.tab)
            tab = join(tab, geom.tab)
            tab = join(tab, cal.tab)
            tab = join(tab, background.tab)
            tab = join(tab, phot.tab)

            tab.write(self.minion_file('science.csv'), delimiter=' ',
                      format='ascii.ecsv', overwrite=True)
        except Exception as e:
            raise MegaTableFailure('{}: {}'.format(type(e).__name__, e))

