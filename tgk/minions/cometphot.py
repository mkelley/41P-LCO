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

        warnings.simplefilter('ignore', MaskedArrayFutureWarning)
            
        logger = logging.getLogger('tgk.science')
        logger.debug('    Comet photometry.')

        yxc, sep = self.centroid()

        # Get background estimate
        try:
            bg = BackgroundTable().get_frame(self.obs.frame_name)
        except IndexError as e:
            raise CometPhotFailure(e)
        
        # photometry
        rap = np.array((2, 4, 6)) / self.obs.pixel_scale
        area, flux, ferr = self.apphot(yxc, rap, bg)
        
        row = [self.obs.frame_name, self.obs.filter,
               sep, yxc[1], yxc[0],
               bg['bg'] / self.obs.exptime.value,
               bg['bgsig'] / self.obs.exptime.value,
               bg['bgarea']]
        
        row.extend([flux[0], ferr[0], flux[1], ferr[1], flux[2], ferr[2]])
        row.extend(np.zeros(6))  # magnitude columns

        CometPhotometry().update(row)

    def centroid(self):
        """Find the comet using the HORIZONS position as a guess.

        1) Smooth the image with a 1/ρ kernel.

        2) Find the peak pixel in a box near the ephemeris position.

        3) Centroid about this point.

        """

        import logging
        import numpy as np
        import astropy.units as u
        from astropy.convolution import convolve
        from astropy.wcs.utils import skycoord_to_pixel
        from ..utils import cutout, gcentroid, UnableToCenter

        logger = logging.getLogger('tgk.science')

        # pre-computed 1/ρ kernel
        K = np.array(
            [[ 0.14142136,  0.15617376,  0.17149859,  0.18569534,  0.19611614, 0.2       ,  0.19611614,  0.18569534,  0.17149859,  0.15617376, 0.14142136],
             [ 0.15617376,  0.1767767 ,  0.2       ,  0.2236068 ,  0.24253563, 0.25      ,  0.24253563,  0.2236068 ,  0.2       ,  0.1767767 , 0.15617376],
             [ 0.17149859,  0.2       ,  0.23570226,  0.2773501 ,  0.31622777, 0.33333333,  0.31622777,  0.2773501 ,  0.23570226,  0.2       , 0.17149859],
             [ 0.18569534,  0.2236068 ,  0.2773501 ,  0.35355339,  0.4472136 , 0.5       ,  0.4472136 ,  0.35355339,  0.2773501 ,  0.2236068 , 0.18569534],
             [ 0.19611614,  0.24253563,  0.31622777,  0.4472136 ,  0.70710678, 1.        ,  0.70710678,  0.4472136 ,  0.31622777,  0.24253563, 0.19611614],
             [ 0.2       ,  0.25      ,  0.33333333,  0.5       ,  1.        , 1.        ,  1.        ,  0.5       ,  0.33333333,  0.25      , 0.2       ],
             [ 0.19611614,  0.24253563,  0.31622777,  0.4472136 ,  0.70710678, 1.        ,  0.70710678,  0.4472136 ,  0.31622777,  0.24253563, 0.19611614],
             [ 0.18569534,  0.2236068 ,  0.2773501 ,  0.35355339,  0.4472136 , 0.5       ,  0.4472136 ,  0.35355339,  0.2773501 ,  0.2236068 , 0.18569534],
             [ 0.17149859,  0.2       ,  0.23570226,  0.2773501 ,  0.31622777, 0.33333333,  0.31622777,  0.2773501 ,  0.23570226,  0.2       , 0.17149859],
             [ 0.15617376,  0.1767767 ,  0.2       ,  0.2236068 ,  0.24253563, 0.25      ,  0.24253563,  0.2236068 ,  0.2       ,  0.1767767 , 0.15617376],
             [ 0.14142136,  0.15617376,  0.17149859,  0.18569534,  0.19611614, 0.2       ,  0.19611614,  0.18569534,  0.17149859,  0.15617376, 0.14142136]]
        )
        
        yxg = skycoord_to_pixel(self.geom.radec_predict, self.obs.wcs)[::-1]
        # smooth and find peak pixel
        if self.geom.delta < 0.25 * u.au:
            cut = cutout(np.array(yxg, int), 30, self.im.data.shape)
        else:
            cut = cutout(np.array(yxg, int), 15, self.im.data.shape)

        subim = self.im.data[cut]
        if subim.shape[0] == 0 or subim.shape[1] == 0:
            raise CometPhotFailure('Ephemeris position outside of image.')
            
        sim = convolve(subim, K, boundary='fill', fill_value=0)
        y, x = np.unravel_index(sim.argmax(), sim.shape)
        y = float(y + cut[0].start)
        x = float(x + cut[1].start)

        box = 13
        while box >= 3:
            try:
                yxc = gcentroid(self.im.data, (y, x), box=13, niter=3)
                break
            except UnableToCenter:
                box -= 2
        else:
            raise CometPhotFailure('Unable to center target.')
                
        sep = np.sqrt(np.sum((np.array(yxg) - np.array(yxc))**2))
        return yxc, sep

    def apphot(self, yxc, rap, bg):
        """Standard aperture photometry on the image data."""
        import numpy as np
        from ..utils import apphot
        
        area, flux = apphot(self.im.data - bg['bg'], yxc, rap, subsample=1)
        bgvar = area * bg['bgsig']**2 * (1 + area / bg['bgarea'])
        ferr = np.sqrt(flux / self.obs.gain.value + bgvar)
        flux = flux / self.obs.exptime.value
        ferr = ferr / self.obs.exptime.value
        return area, flux, ferr

class CometPhotometry(ScienceTable):
    """All comet photometry.

    Parameters
    ----------
    filename : string
      The comet photometry table to read and update.

    """
    _table_title = 'comet photometry'
    _table_columns = [
        'frame', 'filter', 'match sep', 'x', 'y',
        'comet bg', 'comet bgsig', 'bgarea',
        'f2', 'ferr2', 'f4', 'ferr4', 'f6', 'ferr6',
        'm2', 'merr2', 'm4', 'merr4', 'm6', 'merr6',
    ]
    _table_dtypes = ['U64', 'U2'] + [float] * 5 + [int] + [float] * 12
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name.'
    _table_meta['match sep'] = 'Distance between HORIZONS prediction and object centroid, pixels.'
    _table_meta['x/y'] = 'Aperture center, 0-based index, pixels.'
    _table_meta['comet bg'] = 'Background estimate, ADU/s/pixel.'
    _table_meta['comet bgsig'] = 'Background standard deviation per pixel.'
    _table_meta['comet bgarea'] = 'Area used for background estimate.'
    _table_meta['fi, ferri'] = 'Background subtracted flux and error estimates for 2, 4, and 6 radius apertures, ADU/s.'
    _table_meta['mi, merri'] = 'Calibrated magnitudes for each aperture, AB mag.'
    _table_formats = {
        'match sep': '{:.2f}',
        'x': '{:.2f}',
        'y': '{:.2f}',
        'comet bg': '{:.2f}',
        'comet bgsig': '{:.2f}',
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
