# Licensed under a MIT style license - see LICENSE
"""source_mask - Create image source masks."""

from collections import OrderedDict
from . import FrameMinion, MinionError

class SourceMaskFailure(MinionError):
    pass

class SourceMask(FrameMinion):
    """Mask sources in photometry table, mask comet, directly measure bg.

    Based on a 60s frame (elp1m008-fl05-20170402-0100):

      Counts  Mask radius
      ------  -----------
      1e7     36
      7e5     16
      1e5     11
      9e3      4

    >>> linefit(np.log10([1e7, 7e5, 1e5, 9e3]), np.log10([36, 16, 11, 4]),
                         None, [1, 0])

    Log-log fit: log10(radius) = 0.30 * log10(counts) - 0.56

    Use this fit +2 pixels, minimum of 5.


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

    name = 'source mask'

    def run(self):
        import os
        import logging
        from astropy.io import fits

        logger = logging.getLogger('tgk.science')
        
        fn = self.minion_file('{}/{}.fits'.format(
            'source_mask', self.obs.frame_name))
        d = os.path.split(fn)[0]
        if not os.path.exists(d):
            os.mkdir(d)
            logger.debug('Created directory {}.'.format(d))

        if os.path.exists(fn):
            logger.debug('    Using previous source mask.')
            return
        else:
            logger.debug('    Creating source mask.')

        mask = self._mask()

        hdu = fits.PrimaryHDU(mask)
        hdu.writeto(fn)

    def _mask(self):
        import numpy as np
        from scipy.ndimage import binary_dilation, median_filter
        from astropy.wcs import WCS
        from astropy.wcs.utils import skycoord_to_pixel

        mask = np.zeros_like(self.im.data, bool)
        for source in self.im.cat:
            m = np.zeros_like(mask, bool)
            # need sqrt(2) to get binary dilation to cover the full circle
            ap = int(10**(0.30 * np.log10(source['FLUX']) - 0.56) * np.sqrt(2))
            ap = max(3, ap) + 2
            m[int(np.floor(source['Y'])), int(np.floor(source['X']))] = True
            mask += binary_dilation(m, iterations=ap)

        ap = int(300 / self.obs.pixel_scale.value * np.sqrt(2))  # 10 arcmin

        xy = skycoord_to_pixel(self.geom.radec_predict, self.obs.wcs)
        m = np.zeros_like(mask, bool)
        m[int(np.floor(xy[1])), int(np.floor(xy[0]))] = True
        mask += binary_dilation(m, iterations=ap)
        mask[~np.isfinite(self.im.data)] = True

        m = np.ones_like(mask, bool)
        m[self.obs.trimsec] = False
        mask += m

        return mask.astype(int)

