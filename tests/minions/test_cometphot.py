import numpy as np
import tgk.core
tgk.core.configure()

class TestCometPhot:
    def test_centroid(self):
        import os
        from tgk import minions, core
        from tgk.science import Science

        frame = 'elp1m008-fl05-20170329-0170'
        fn = os.path.join(core.config['download path'],
                          'e91', '20170330',
                          'elp1m008-fl05-20170329-0170-e91.fits.fz')
        im, obs = Science.get_frame_data(frame, fn)
        geom = Science.get_geometry(frame, obs)

        m = minions.CometPhot(core.config, im, obs, geom)
        yxc, sep = m.centroid()
        
        yx0 = np.array([2222.7,  2103.0])
        d = np.sqrt(np.sum((yxc - yx0)**2))
        assert d < 0.5

    def test_photometry(self):
        import os
        from tgk import minions, core
        from tgk.science import Science

        frame = 'elp1m008-fl05-20170329-0170'
        fn = os.path.join(core.config['download path'],
                          'e91', '20170330',
                          'elp1m008-fl05-20170329-0170-e91.fits.fz')
        im, obs = Science.get_frame_data(frame, fn)
        geom = Science.get_geometry(frame, obs)

        im._hdu['sci'].data = 2 * np.ones_like(im.data)
        
        m = minions.CometPhot(core.config, im, obs, geom)
        bg = dict(bg=0, bgsig=0.5, bgarea=1000)
        area, flux, ferr = m.apphot([100, 100], 6, bg)
        
        assert np.isclose(area, 109)
        assert np.isclose(flux, 109 * 2 / obs.exptime.value)
        assert np.isclose(ferr, 15.75500714058867 / obs.exptime.value)
