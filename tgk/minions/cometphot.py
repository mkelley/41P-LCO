# Licensed under a MIT style license - see LICENSE
"""cometphot - Measure comet photometry."""

from collections import OrderedDict
from . import FrameMinion, MinionError
from ..science import ScienceTable

class CometPhotFailure(MinionError):
    pass

class CometPhot(FrameMinion):
    """Comet photometry.

    Requires background estimate.

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

    name = 'cometphot'

    def run(self):
        import logging
        import warnings
        import numpy as np
        from numpy.ma.core import MaskedArrayFutureWarning
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astropy.wcs.utils import skycoord_to_pixel
        from .background import BackgroundTable
        from ..utils import apphot

        warnings.simplefilter('ignore', MaskedArrayFutureWarning)
            
        logger = logging.getLogger('tgk.science')
        logger.info('    Comet photometry.')

        # centroid around the JPL/HORIZONS prediction
        yxg = skycoord_to_pixel(self.geom.radec_predict, self.obs.wcs)[::-1]
        yxc = gcentroid(self.im.data, yxg, box=21, niter=5)
        sep = np.sqrt(np.sum((np.array(yxg) - np.array(yxc))**2))

        # Get background estimate
        try:
            bg = BackgroundTable().get_frame(self.obs.frame_name)
        except IndexError as e:
            raise CometPhotFailure(e)
        
        # photometry
        rap = np.array((2, 4, 6)) / self.obs.pixel_scale
        area, flux = apphot(self.im.data - bg['bg'], yxc, rap, subsample=1)
        flux /= self.obs.exptime.value
        bgvar = area * bg['bgsig']**2 * (1 + area / bg['bgarea'])
        ferr = np.sqrt(flux / self.obs.gain.value + bgvar)
        
        row = [self.obs.frame_name, self.obs.filter,
               sep, yxc[1], yxc[0], bg['bg'], bg['bgsig'], bg['bgarea']]
        row.extend(zip(flux, ferr))
        row.extend(np.zeros(6))  # magnitude columns

        CometPhotometry().update(row)

class CometPhotometry(ScienceTable):
    """All comet photometry.

    Parameters
    ----------
    filename : string
      The comet photometry table to read and update.

    """
    _table_title = 'comet photometry'
    _table_columns = [
        'frame', 'filter', 'match sep', 'x', 'y', 'bg', 'bgsig', 'bgarea',
        'f2', 'ferr2', 'f4', 'ferr4', 'f6', 'ferr6',
        'm2', 'merr2', 'm4', 'merr4', 'm6', 'merr6',
    ]
    _table_dtypes = ['U64', 'U2'] + [float] * 5 + [int] + [float] * 12
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name.'
    _table_meta['match sep'] = 'Distance between HORIZONS prediction and object centroid, pixels.'
    _table_meta['x/y'] = 'Aperture center, 0-based index, pixels.'
    _table_meta['bg'] = 'Background estimate, ADU/s/pixel.'
    _table_meta['bgsig'] = 'Background standard deviation per pixel.'
    _table_meta['bgarea'] = 'Area used for background estimate.'
    _table_meta['fi, ferri'] = 'Background subtracted flux and error estimates for 2, 4, and 6 radius apertures, ADU/s.'
    _table_meta['mi, merri'] = 'Calibrated magnitudes for each aperture, AB mag.'
    _table_formats = {
        'match sep': '{:.2f}',
        'x': '{:.2f}',
        'y': '{:.2f}',
        'bg': '{:.2f}',
        'bgsig': '{:.2f}',
        'f2': '{:.5g}',
        'f4': '{:.5g}',
        'f6': '{:.5g}',
        'ferr2': '{:.5g}',
        'ferr4': '{:.5g}',
        'ferr6': '{:.5g}',
        'm2': '{:.3f}',
        'm4': '{:.3f}',
        'm6': '{:.3f}',
        'merr2': '{:.3f}',
        'merr4': '{:.3f}',
        'merr6': '{:.3f}',
    }
    _table_sort = 'frame'

    def __init__(self, verbose=False):
        ScienceTable.__init__(self, 'cometphot.csv', verbose=verbose)

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)
    
    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)
