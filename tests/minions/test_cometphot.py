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
