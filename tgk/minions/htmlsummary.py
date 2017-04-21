# Licensed under a MIT style license - see LICENSE
"""htmlsummary - Summarize results into an HTML format."""

from . import TableMinion, MinionError

class HTMLSummaryFailure(MinionError):
    pass

class HTMLSummary(TableMinion):
    """Summary.

    Parameters
    ----------
    config : dict
      Configuration parameters.

    """

    name = 'htmlsummary'

    def run(self):
        import logging
        from astropy.table import join
        from .cometphot import CometPhotometry
        from .background import BackgroundTable
        from .calibrate import CalibrationTable
        from ..science import GeometryTable, ProcessingHistory, ObservationLog

        logger = logging.getLogger('tgk.science')
        logger.info('  Create HTML summary.')

        try:
            obslog = ObservationLog()
            hist = ProcessingHistory()
            geom = GeometryTable()
            cal = CalibrationTable()
            background = BackgroundTable()
            phot = CometPhotometry()

            with open(self.minion_file('summary.html'), 'w') as outf:
                outf.write('''<html lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta charset="utf-8">
    <base href="/~msk/">
    <link href="css/default.css" rel="stylesheet">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
    <script type="text/javascript" src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha256-ZosEbRLbNQzLpnKIkEdrPv7lOy9C27hHQ+Xp8a4MxAQ=" crossorigin="anonymous"></script><script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script><meta name="description" content="The astronomical homepage of Michael S. P. Kelley.">
    <title>41P LCO Pipeline Science Summary</title>
  </head>
  <body>
''')

                outf.write('\n    <h2>Observing log</h2>\n')
                obslog.tab.write(outf, format='html')

                outf.write('\n    <h2>Geometry</h2>\n')
                geom.tab.write(outf, format='html')

                outf.write('\n    <h2>Processing history</h2>\n')
                tab = join(hist.tab, obslog.tab['date', 'time'],
                           join_type='outer')
                tab.sort(['date', 'time'])
                reorder = tab.colnames
                reorder.remove('frame')
                reorder.remove('date')
                reorder.remove('time')
                reorder = ['frame', 'date', 'time'] + reorder
                tab = tab[reorder]
                tab.write(outf, format='html')

                outf.write('\n    <h2>Calibration</h2>\n')
                cal.tab.add_column(obslog.tab['time'], index=1)
                cal.tab.add_column(obslog.tab['date'], index=1)
                cal.tab.sort(['date', 'time'])
                cal.tab.write(outf, format='html')

                outf.write('\n    <h2>Background</h2>\n')
                background.tab.add_column(obslog.tab['time'], index=1)
                background.tab.add_column(obslog.tab['date'], index=1)
                background.tab.sort(['date', 'time'])
                background.tab.write(outf, format='html')

                outf.write('\n    <h2>Photometry</h2>\n')
                phot.tab.add_column(obslog.tab['time'], index=1)
                phot.tab.add_column(obslog.tab['date'], index=1)
                phot.tab.sort(['date', 'time'])
                phot.tab.write(outf, format='html')

        except Exception as e:
            raise HTMLSummaryFailure('{}: {}'.format(type(e).__name__, e))
