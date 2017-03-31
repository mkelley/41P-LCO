# Licensed under a MIT style license - see LICENSE
"""calibrate - Derive zero-point magnitude with PanSTARRS."""

from . import FrameMinion

class Calibrate(FrameMinion):
    """Calibrate with PanSTARRS

    VO Simple Cone Search parameters at:
      https://archive.stsci.edu/vo/general_params.html

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

    _table_names = ['frame', 'filter', 'N cat', 'min(sep)', 'max(sep)',
                    'scmean(sep)', 'scmedian(sep)', 'scstdev(sep)',
                    'N match', 'min(dm)', 'max(dm)', 'scmean(dm)',
                    'scmedian(dm)', 'scstdev(dm)']
    _table_dtype = ['U64', 'U2'] + ([int] + [float] * 5) * 2
    _table_format = ([None, None, None] + ['{:.2f}'] * 5
                     + [None] + ['{:.3f}'] * 5)
    _table_comments = ['Units: arcsec, magnitude']

    @property
    def name(self):
        return 'calibrate'

    def run(self):
        import os
        import logging
        import requests
        import numpy as np
        import astropy.units as u
        from astropy.coordinates import SkyCoord, match_coordinates_sky
        from astropy import stats
        from astropy.io import votable
        from .. import lco

        log = [self.obs.frame_name, self.obs.filter]  # calibration log
        fn = self.minion_file('{}.xml'.format(self.obs.frame_name))
        ps1filter = lco.filter2PS1[self.obs.filter]
        
        if not os.path.exists(fn):
            self.logger.info('    Retrieving PS1 catalog for {} .'.format(self.obs.frame_name))
            r = max(self.im.data.shape) * self.obs.pixel_scale * 2 / 3
            c = self.geom.radec_predict
            columns = ','.join([
                'objname',
                'objid',
                'ramean',
                'decmean',
                'rameanerr',
                'decmeanerr',
                'ndetections',
                'randomid',
                'projectionid',
                'skycellid',
                'objinfoflag',
                'qualityflag',
                'rastack',
                'decstack',
                'rastackerr',
                'decstackerr',
                'epochmean',
                'nstackdetections',
                'ng',
                'nr',
                'gqfperfect',
                'gmeanpsfmag',
                'gmeanpsfmagerr',
                'gmeankronmag',
                'gmeankronmagerr',
                'gmeanapmag',
                'gmeanapmagerr',
                'gflags',
                'rqfperfect',
                'rmeanpsfmag',
                'rmeanpsfmagerr',
                'rmeankronmag',
                'rmeankronmagerr',
                'rmeanapmag',
                'rmeanapmagerr',
                'rflags',
            ])

            params = dict(RA=c.ra.deg, DEC=c.dec.deg, SR=r.to(u.deg).value,
                          max_records=3000, ordercolumn1='ndetections',
                          descending1='on', selectedColumnsCsv=columns)
            
            q = requests.get('https://archive.stsci.edu/panstarrs/search.php',
                             params=params)

            with open(fn, 'w') as outf:
                outf.write(q.text)

        self.logger.info('    Reading PS1 catalog from {} .'.format(fn))
        cat = votable.parse_single_table(fn)

        # match catalogs
        ps1 = SkyCoord(ra=cat.array['ramean'],
                       dec=cat.array['decmean'],
                       unit='deg')
        lco = SkyCoord(ra=self.im.cat['RA'],
                       dec=self.im.cat['DEC'],
                       unit='deg')
        match, sep = match_coordinates_sky(lco, ps1)[:2]

        # match stats
        minmax = min(sep.to(u.arcsec).value), max(sep.to(u.arcsec).value)
        mms = stats.sigma_clipped_stats(sep.to(u.arcsec).value)
        self.logger.info('    Matched {0} objects to PS1.'.format(len(match)))
        self.logger.debug('''      - Best/worst match separation distance: {0[0]:.2f}/{0[1]:.2f}
      - Sigma-clipped separation: mean/median/stddev = {1[0]:.2f}/{1[1]:.2f}/{1[2]:.2f} arcsec.'''.format(minmax, mms))
        log.append(len(match))
        log.extend(minmax)
        log.extend(mms)

        # user's match radius cutoff
        i = sep < self.config['calibrate match radius'] * u.arcsec

        # Remove PS1 catalog missing data
        ps1magfield = '{}meanapmag'.format(ps1filter)
        i *= cat.array[ps1magfield][match] > -999
        
        # measure magnitude offset, use ADU/s
        m = -2.5 * np.log10(self.im.cat['FLUX'][i] / self.obs.exptime.value)
        merr = 1.0857 * self.im.cat['FLUXERR'][i] / self.im.cat['FLUX'][i]
        dm = cat.array[ps1magfield][match][i] - m
        dm_err = np.sqrt(merr**2 + cat.array[ps1magfield + 'err'][match][i]**2)

        # cal stats
        minmax = min(dm), max(dm)
        mms = stats.sigma_clipped_stats(dm)
        self.logger.info('    Calibrated using {0} objects.'.format(i.sum()))
        self.logger.debug('''      - delta-mag range: {minmax[0]:.3f} - {minmax[1]:.3f}
      - Sigma-clipped mean/median/stddev = {mms[0]:.3f}/{mms[1]:.3f}/{mms[2]:.3f} arcsec.'''.format(minmax=minmax, mms=mms))
        log.append(i.sum())
        log.extend(minmax)
        log.extend(mms)

        # save to calibration log
        self.append_to_table(log)
        stop
