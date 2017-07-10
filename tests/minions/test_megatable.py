import numpy as np
import tgk.core
tgk.core.configure()

def test_megatable():
    from astropy.table import join
    from tgk.minions.cometphot import CometPhotometry
    from tgk.minions.background import BackgroundTable
    from tgk.minions.calibrate import CalibrationTable
    from tgk.science import GeometryTable, ProcessingHistory, ObservationLog

    obslog = ObservationLog()
    hist = ProcessingHistory()
    geom = GeometryTable()
    cal = CalibrationTable()
    background = BackgroundTable()
    phot = CometPhotometry()

    tab = join(obslog.tab, hist.tab)
    assert len(tab) > 0

    tab = join(tab, geom.tab)
    assert len(tab) > 0

    tab = join(tab, cal.tab)
    assert len(tab) > 0

    tab = join(tab, background.tab)
    assert len(tab) > 0

    tab = join(tab, phot.tab)
    assert len(tab) > 0

