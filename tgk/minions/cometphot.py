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

    @property
    def name(self):
        return 'cometphot'

    def run(self):
        import logging
        import warnings
        import numpy as np
        from numpy.ma.core import MaskedArrayFutureWarning
        import astropy.units as u
        from astropy.coordinates import SkyCoord, match_coordinates_sky
        from .background import BackgroundTable

        warnings.simplefilter('ignore', MaskedArrayFutureWarning)
            
        logger = logging.getLogger('tgk.science')
        logger.info('    Comet photometry.')

        # find comet in LCO catalog
        lco = SkyCoord(ra=self.im.cat['RA'],
                       dec=self.im.cat['DEC'],
                       unit='deg')
        c = self.geom.radec_predict
        match, sep = match_coordinates_sky(c, lco)[:2]
        logger.info('      Matched comet from HORIZONS coordinates to LCO object {:.2f} away.'.format(sep[0].to(u.arcsec)))
        
        comet = self.im.cat[match]
        
        # Take comet photometry and use our own background estimate
        bg = BackgroundTable().get_frame(self.obs.frame_name)
        
        # first index: aperture, second: flux and ferr
        flux = np.empty((3, 2))
        for i, col in enumerate(['FLUXAPER2', 'FLUXAPER4', 'FLUXAPER6']):
            area = 3.1416 * (1 / self.obs.pixel_scale.value)**2
            flux[i, 0] = comet[col] + comet['BACKGROUND'] * area

            bgvar = area * bg['bgsig']**2 * (1 + area / bg['bgarea'])
            flux[i, 1] = np.sqrt(flux[i, 0] / self.obs.gain.value + bgvar)
            flux[i, 0] -= bg['bg'] * area

        flux /= self.obs.exptime.value
        
        row = [self.obs.frame_name, self.obs.filter,
               sep[0].to(u.arcsec).value, comet['X'], comet['Y'],
               bg['bg'], bg['bgsig'], bg['bgarea']]
        row.extend(flux.ravel())
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
        'f1', 'ferr1', 'f2', 'ferr2', 'f3', 'ferr3',
        'm1', 'merr1', 'm2', 'merr2', 'm3', 'merr3',
    ]
    _table_dtypes = ['U64', 'U2'] + [float] * 5 + [int] + [float] * 12
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name.'
    _table_meta['match sep'] = 'Distance between HORIZONS prediction and matched LCO object.'
    _table_meta['x/y'] = 'Aperture center, 0-based index, pixels.'
    _table_meta['bg'] = 'Background estimate, ADU/s/pixel.'
    _table_meta['bgsig'] = 'Background standard deviation per pixel.'
    _table_meta['bgarea'] = 'Area used for background estimate.'
    _table_meta['fi, ferri'] = 'Background subtracted flux and error estimates for 1, 2, and 3" radius apertures, ADU/s.'
    _table_meta['mi, merri'] = 'Calibrated magnitudes for each aperture, AB mag.'
    _table_formats = {
        'match sep': '{:.2f}',
        'x': '{:.2f}',
        'y': '{:.2f}',
        'bg': '{:.2f}',
        'bgsig': '{:.2f}',
        'f1': '{:.5g}',
        'f2': '{:.5g}',
        'f3': '{:.5g}',
        'ferr1': '{:.5g}',
        'ferr2': '{:.5g}',
        'ferr3': '{:.5g}',
        'm1': '{:.3f}',
        'm2': '{:.3f}',
        'm3': '{:.3f}',
        'merr1': '{:.3f}',
        'merr2': '{:.3f}',
        'merr3': '{:.3f}',
    }
    _table_sort = 'frame'

    def __init__(self, verbose=False):
        ScienceTable.__init__(self, 'cometphot.csv', verbose=verbose)

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)
    
    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)
