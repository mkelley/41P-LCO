#!/usr/bin/env python
from numpy.distutils.core import setup

if __name__ == "__main__":
    setup(name='tgk',
          version='0.2.0-dev',
          description='41P/Tuttle-Giacobini-Kresak Outburst Project Pipeline',
          author="Michael S. P. Kelley",
          author_email="msk@astro.umd.edu",
          url="https://github.com/mkelley/41P-LCO",
          packages=['tgk', 'tgk.minions'],
          scripts=['scripts/tgk-science', 'scripts/tgk-sync'],
          requires=['numpy', 'astropy', 'requests', 'astroquery'],
          license='MIT',
      )
