# Licensed under a MIT style license - see LICENSE
"""calibrate - Derive zero-point magnitude with PanSTARRS."""

from collections import OrderedDict
from . import FrameMinion, MinionError
from ..science import ScienceTable

class CalibrationFailure(MinionError):
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

    def __init__(self, *args):
        FrameMinion.__init__(self, *args, minion_directory=True)
    
    name = 'calibrate'

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
        from astropy.table import Table
        from .. import lco

        # calibration log
        log = [self.obs.frame_name, self.obs.filter, self.obs.airmass]
        
        fn = self.minion_file('{}.xml'.format(self.obs.frame_name))
        ps1filter = lco.filter2PS1[self.obs.filter]

        cols = self.im.cat.colnames
        if 'RA' not in cols or 'DEC' not in cols:
            raise CalibrationFailure('Missing RA and DEC in source table.')
        
        if not os.path.exists(fn):
            self.logger.info('    Retrieving PS1 catalog.'.format())
            r = max(self.im.data.shape) * self.obs.pixel_scale * 2 / 3
            r = min(r, 30 * u.arcmin)  # STScI has a 30 arcmin max radius
            c = self.geom.radec_predict
            columns = ','.join([    
                'objname', 'objid', 'ramean',
                'decmean', 'rameanerr', 'decmeanerr', 'ndetections',
                'randomid', 'projectionid', 'skycellid',
                'objinfoflag', 'qualityflag', 'rastack', 'decstack',
                'rastackerr', 'decstackerr', 'epochmean',
                'nstackdetections', 'ng', 'nr', 'gqfperfect',
                'gmeanpsfmag', 'gmeanpsfmagerr', 'gmeankronmag',
                'gmeankronmagerr', 'gmeanapmag', 'gmeanapmagerr',
                'gflags', 'rqfperfect', 'rmeanpsfmag',
                'rmeanpsfmagerr', 'rmeankronmag', 'rmeankronmagerr',
                'rmeanapmag', 'rmeanapmagerr', 'rflags'
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
        m = -2.5 * np.log10(self.im.cat['FLUX'] / self.obs.exptime.value)
        merr = 1.0857 * self.im.cat['FLUXERR'] / self.im.cat['FLUX']
        dm = cat.array[ps1magfield][match] - m
        dm_err = np.sqrt(merr**2 + cat.array[ps1magfield + 'err'][match]**2)

        if any(~np.isfinite(dm[i] * dm_err[i])):
            raise CalibrationFailure('Not all values are finite: {}.'.format(
                self.obs.frame_name))
        
        # cal stats
        minmax = min(dm[i]), max(dm[i])
        mms = stats.sigma_clipped_stats(dm[i])
        self.logger.debug('''      - delta-mag range: {minmax[0]:.3f} - {minmax[1]:.3f}
      - Sigma-clipped mean/median/stddev = {mms[0]:.3f}/{mms[1]:.3f}/{mms[2]:.3f} arcsec.'''.format(minmax=minmax, mms=mms))
        log.extend(minmax)
        log.extend(mms)

        # save to calibration log
        CalibrationTable().update(log)

        # save match catalog
        tab = Table(data=(cat.array['objname'][match],
                          cat.array['ramean'][match],
                          cat.array['decmean'][match],
                          self.im.cat['RA'] - cat.array['ramean'][match],
                          self.im.cat['DEC'] - cat.array['decmean'][match],
                          sep.to(u.arcsec).value,
                          cat.array['gflags'][match].astype(int),
                          cat.array['gmeanapmag'][match],
                          cat.array['gmeanapmagerr'][match],
                          cat.array['rflags'][match].astype(int),
                          cat.array['rmeanapmag'][match],
                          cat.array['rmeanapmagerr'][match],
                          self.im.cat['X'],
                          self.im.cat['Y'],
                          self.im.cat['BACKGROUND'] / self.obs.exptime.value,
                          self.im.cat['FWHM'] * self.obs.pixel_scale,
                          self.im.cat['FLUX'] / self.obs.exptime.value,
                          self.im.cat['FLUXERR'] / self.obs.exptime.value,
                          dm,
                          i),
                    names=('objname',
                           'RAmean',
                           'DECmean',
                           'RA_LCO - RA_PS1',
                           'DEC_LCO - DEC_PS1',
                           'sep',
                           'gflags',
                           'gmeanapmag',
                           'gmeanapmagerr',
                           'rflags',
                           'rmeanapmag',
                           'rmeanapmagerr',
                           'x',
                           'y',
                           'bg',
                           'FWHM',
                           'flux',
                           'fluxerr',
                           'dm',
                           'valid match'),
                    dtype=(['U64'] + [float] * 5
                           + [int, float, float] * 2
                           + [float] * 7
                           + [bool])
        )
        cformats = ([None] + ['{:.6f}'] * 4
                    + ['{:.2f}']
                    + [None, '{:.3f}', '{:.3f}'] * 2
                    + ['{:.2f}'] * 6 + ['{:.3f}', None])
        tab.meta = OrderedDict()
        tab.meta['objname'] = 'PS1 object name.'
        tab.meta['ramean/decmean'] = 'PS1 mean RA DEC (deg).'
        tab.meta['x_LCO - x_PS1'] = 'Offset between LCO object and PS1 match (deg).'
        tab.meta['sep'] = 'Radial offset magnitude (arcsec).'
        tab.meta['gflags'] = 'PS1 catalog g\' flags.'
        tab.meta['gmeanapmag/err'] = 'PS1 g\' aperture photometry magnitude and err (mag).'
        tab.meta['rflags'] = 'PS1 catalog r\' flags.'
        tab.meta['rmeanapmag/err'] = 'PS1 r\' aperture photometry magnitude and err (mag).'
        tab.meta['x/y'] = 'LCO source position on CCD (pixels).'
        tab.meta['bg'] = 'LCO source background (ADU/s/pixel).'
        tab.meta['FWHM'] = 'LCO source FWHM (arcsec).'
        tab.meta['flux'] = 'LCO source flux (ADU/s).'
        tab.meta['fluxerr'] = 'LCO source flux uncertainty (ADU/s).'
        tab.meta['dm'] = 'Difference between calibrated magnitude and instrumental magnitude (mag).'
        tab.meta['valid match'] = 'Flag that indicates if the match and photometry was presumed to be valid.'
        for col, cformat in zip(tab.colnames, cformats):
            tab[col].format = cformat
        fn = self.minion_file('{}-matches.csv'.format(self.obs.frame_name))
        tab.write(fn, format='ascii.ecsv', overwrite=True, delimiter=',')

class CalibrationTable(ScienceTable):
    _table_title = 'calibration'
    _table_columns = [
        'frame', 'filter', 'airmass', 'N cat', 'N match',
        'mean(FWHM)', 'mean(background)', 'min(sep)', 'max(sep)',
        'scmean(sep)', 'scmedian(sep)', 'scstdev(sep)', 'min(dm)',
        'max(dm)', 'scmean(dm)', 'scmedian(dm)', 'scstdev(dm)'
    ]
    _table_dtypes = ['U64', 'U2', float, int, int] + [float] * 12
    _table_formats = ([None, None, '{:.3f}', '{:4d}', '{:04d}']
                      + ['{:.2f}'] * 7 + ['{:.3f}'] * 5)
    _table_meta = OrderedDict()
    _table_meta['filter'] = 'LCO filter name'
    _table_meta['N cat'] = 'Number of objects in LCO photometry catalog.'
    _table_meta['N match'] = 'Number of objects matched with PS1 catalog, after distance and non-zero flux tests.'
    _table_meta['mean(FWHM)'] = 'Mean LCO FWHM of matched objects, arcsec.'
    _table_meta['mean(background)'] = 'Mean LCO photometry catalog background for matched objects, ADU/s.'
    _table_meta['sep'] = 'Spherical distance between matched objects, arcsec.'
    _table_meta['dm'] = 'Difference in magnitudes between matched objects.'
    _table_meta['comments'] = 'scmean/scmedian/scstdev are sigma-clipped mean/median/standard deviation; LCO catalog magnitudes are based on flux in ADU/s.'
    _table_sort = ['frame']

    def __init__(self, verbose=False):
        ScienceTable.__init__(self, 'calibrate.csv', verbose=verbose)

    def update(self, row):
        """Add row to table."""
        self._update_unique_column('frame', row)
    
    def get_frame(self, frame):
        return self._get_unique_row('frame', frame)
