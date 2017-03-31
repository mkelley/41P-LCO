# Licensed under a MIT style license - see LICENSE
"""calibrate - Derive zero-point magnitude with PanSTARRS."""

from . import FrameMinion

class Calibrate(FrameMinion):
    """Calibrate with PanSTARRS

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
        return 'calibrate'

    def run(self):
        import os
        import astropy.units as u
        from astropy.io import votable
        import logging
        import requests

        fn = self.minion_file('{}.xml'.format(self.obs.frame_name))
        if not os.path.exists(fn):
            self.logger.info('    Retrieving PS1 catalog for {} .'.format(self.obs.frame_name))
            r = max(self.im.data.shape) * self.obs.pixel_scale / 2
            c = self.geom.radec_predict
            params = dict(RA=c.ra.deg, DEC=c.dec.deg, SR=r.to(u.deg).value,
                          max_records=1000)
            q = requests.get('https://archive.stsci.edu/panstarrs/search.php',
                             params=params)

            with open(fn, 'w') as outf:
                outf.write(q.text)

        self.logger.debug('    Reading PS1 catalog from {} .'.format(fn))
        cat = votable.parse(fn)

