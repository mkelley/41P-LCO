# Licensed under a MIT style license - see LICENSE
"""calibrate - Derive zero-point magnitude with PanSTARRS."""

from . import FrameMinion

class CalibrationFailure(Exception):
    pass

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

    _table_names = ['frame', 'filter', 'airmass', 'N cat', 'N match',
                    'mean(FWHM)', 'mean(background)', 'min(sep)', 'max(sep)',
                    'scmean(sep)', 'scmedian(sep)', 'scstdev(sep)', 'min(dm)',
                    'max(dm)', 'scmean(dm)', 'scmedian(dm)', 'scstdev(dm)']
    _table_dtype = (['U64', 'U2', float, int, int] + [float] * 12)
    _table_format = ([None, None, '{:.3f}', None, None] + ['{:.2f}'] * 7
                     + ['{:.3f}'] * 5)
    _table_sort = ['frame']
    _table_meta = {
        'filter': 'LCO filter name',
        'N cat': 'Number of objects in LCO photometry catalog.',
        'N match': 'Number of objects matched with PS1 catalog, after distance and non-zero flux tests.',
        'mean(FWHM)': 'Mean LCO FWHM of matched objects.',
        'mean(background)': 'Mean LCO photometry catalog background for matched objects, ADU/s.',
        'sep': 'Spherical distance between matched objects, arcsec.',
        'dm': 'Difference in magnitudes between matched objects.',
        'comments': 'scmean/scmedian/scstdev are sigma-clipped mean/median/standard deviation; LCO catalog magnitudes are based on flux in ADU/s.',
    }

    @property
    def name(self):
        return 'calibrate'

    def run(self):
        import os
        import logging
        import warnings
        import requests
        import numpy as np
        import astropy.units as u
        from astropy.coordinates import SkyCoord, match_coordinates_sky
        from astropy import stats
        from astropy.io import votable
        from .. import lco

        # calibration log
        log = [self.obs.frame_name, self.obs.filter, self.obs.airmass]
        
        fn = self.minion_file('{}.xml'.format(self.obs.frame_name))
        ps1filter = lco.filter2PS1[self.obs.filter]
        
        if not os.path.exists(fn):
            return
            self.logger.info('    Retrieving PS1 catalog.'.format())
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

        # suppress VOTable version warning
        warnings.simplefilter('ignore', category=votable.VOWarning)
        cat = votable.parse_single_table(fn)
        warnings.resetwarnings()

        # match catalogs
        ps1 = SkyCoord(ra=cat.array['ramean'],
                       dec=cat.array['decmean'],
                       unit='deg')
        lco = SkyCoord(ra=self.im.cat['RA'],
                       dec=self.im.cat['DEC'],
                       unit='deg')
        match, sep = match_coordinates_sky(lco, ps1)[:2]

        # user's match radius cutoff
        i = sep < self.config['calibrate match radius'] * u.arcsec

        # Remove PS1 catalog missing data
        ps1magfield = '{}meanapmag'.format(ps1filter)
        i *= cat.array[ps1magfield][match] > -999
        
        # Remove bogus photometry
        i *= (self.im.cat['FLUX'] * self.im.cat['FLUXERR']) > 0

        # any objects left?
        if i.sum() == 0:
            raise CalibrationFailure('{} has no good objects.'.format(
                self.obs.frame_name))
        
        # match stats
        minmax = min(sep[i].to(u.arcsec).value), max(sep[i].to(u.arcsec).value)
        mms = stats.sigma_clipped_stats(sep[i].to(u.arcsec).value)
        self.logger.info('    {} of {} objects matched to PS1 after distance and non-zero flux tests.'.format(i.sum(), len(lco)))
        self.logger.debug('''      - Best/worst match separation distance: {0[0]:.2f}/{0[1]:.2f}
      - Sigma-clipped separation: mean/median/stddev = {1[0]:.2f}/{1[1]:.2f}/{1[2]:.2f} arcsec.'''.format(minmax, mms))

        log.append(len(lco))
        log.append(i.sum())
        log.append(np.mean(self.im.cat['FWHM'][i]) * self.obs.pixel_scale)
        log.append(np.mean(self.im.cat['BACKGROUND'][i])
                   / self.obs.exptime.value)
        log.extend(minmax)
        log.extend(mms)
        
        # measure magnitude offset, use ADU/s
        m = -2.5 * np.log10(self.im.cat['FLUX'][i] / self.obs.exptime.value)
        merr = 1.0857 * self.im.cat['FLUXERR'][i] / self.im.cat['FLUX'][i]
        dm = cat.array[ps1magfield][match][i] - m
        dm_err = np.sqrt(merr**2 + cat.array[ps1magfield + 'err'][match][i]**2)

        if any(~np.isfinite(dm * dm_err)):
            raise CalibrationFailure('Not all values are finite: {}.'.format(
                self.obs.frame_name))
        
        # cal stats
        minmax = min(dm), max(dm)
        mms = stats.sigma_clipped_stats(dm)
        self.logger.debug('''      - delta-mag range: {minmax[0]:.3f} - {minmax[1]:.3f}
      - Sigma-clipped mean/median/stddev = {mms[0]:.3f}/{mms[1]:.3f}/{mms[2]:.3f} arcsec.'''.format(minmax=minmax, mms=mms))
        log.extend(minmax)
        log.extend(mms)

        # save to calibration log
        self.append_to_table(log)
