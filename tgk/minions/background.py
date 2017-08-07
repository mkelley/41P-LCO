# Licensed under a MIT style license - see LICENSE
"""background - Estiate image background."""

from collections import OrderedDict
from . import FrameMinion, MinionError
from ..science import ScienceTable

class BackgroundFailure(MinionError):
    pass

class Background(FrameMinion):
    """Estimate comet background.

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

    name = 'background'

    def run(self):
        import os
        import logging
        import numpy as np
        from astropy.io import fits
        
        logging.getLogger('tgk.science').debug('    Measuring background.')

        fn = '{}/source_mask/{}.fits'.format(self.config['science path'],
                                             self.obs.frame_name)
        mask = fits.getdata(fn).astype(bool)
        sky = self.im.data[~mask]
        
        row = [self.obs.frame_name, self.obs.filter]
        row.extend(self._est(sky)[1:])
        
        BackgroundTable().update(row)

    def _est(self, sky, sigma_lower=3, sigma_upper=3, **kwargs):
        import numpy as np
        from astropy import stats

        clipped = stats.sigma_clip(sky, sigma_lower=sigma_lower,
                                   sigma_upper=sigma_upper)
        bg = np.ma.median(clipped)
        bgsig = clipped.std()
        bgarea = np.sum(~clipped.mask, dtype=int)
        return clipped, bg, bgsig, bgarea

class BackgroundTable(ScienceTable):
    """Background measurements.

    Parameters
    ----------
    filename : string
      The background table to read and update.

    """

    _table_title = 'background'
    _table_columns = [
        'frame', 'filter', 'bg', 'bgsig', 'bgarea'
    ]
    _table_dtypes = ['U64', 'U2', float, float, int]
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name.'
    _table_meta['bg'] = 'Sigma-clipped background estimate, ADU/pixel.'
    _table_meta['bgsig'] = 'Background standard deviation per pixel.'
    _table_meta['bgarea'] = 'Area used for background estimate, pixels.'
    _table_formats = [None, None, '{:.2f}', '{:.2f}', None]
    _table_sort = 'frame'

    def __init__(self, verbose=False):
        ScienceTable.__init__(self, 'background.csv', verbose=verbose)
    
    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)
