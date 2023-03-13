# -*- coding: utf-8 -*-
try:
    import setuptools
except ImportError:
    pass

import os
import glob
from distutils.core import setup

setup(name = "infrapy",
      license='LANL-MIT',
      version = '0.4.0.1',
      description = "A tool for processing infrasound and seismic array data.",
      keywords=['infrasound', 'geophysics', 'seismic', 'array'],
      author = "LANL Seismoacoustics Infrasound (LANL-SA) Team",
      author_email = 'pblom@lanl.gov',
      packages = ['infrapy',
                  'infrapy.association',
                  'infrapy.cli',
                  'infrapy.database',
                  'infrapy.database.taskbase',
                  'infrapy.detection',
                  'infrapy.location',
                  'infrapy.performance',
                  'infrapy.propagation',
                  'infrapy.utils',
                  'InfraView',
                  'InfraView.widgets',
                  'InfraView.graphics'],

      entry_points = {'console_scripts':['infrapy=infrapy.cli.__main__:main'],
                      'gui_scripts':['infraview = InfraView.__main__:main']},

      scripts=['scripts/detection_list.py',
               'scripts/plot1_rfk.py',
               'scripts/plot_fd.py',
               'scripts/plot_network.py',
               'scripts/plot_psd_array.py',
               'scripts/print_fk.py',
               'scripts/print_rfd.py',
               'scripts/read_pfd.py',
               'scripts/read_pfk.py',
               'scripts/update_refsta.py',
               'scripts/update_calib.py',
               'scripts/read_rassoc_numdet.py',
               'scripts/read_rassoc.py',
               'scripts/read_rfd_fast.py',
               'scripts/read_rlocBISL.py',
               'scripts/read_rloc.py',
               'scripts/update_calib.py',
               'scripts/update_chan.py',
               'scripts/update_refsta.py'],

      package_dir={'infrapy.propagation' : 'infrapy/propagation'},
      package_data={'infrapy.propagation' : ['compass.png',
                                             'ak135_1st_arrivals.dat'],
                    'infrapy.resources' :['default.config']},

      install_requires = ['numpy',
                          'scipy',
                          'obspy',
                          'pisces',
                          'scikit-learn',
                          'pathos',
                          'numba',
                          'pyqtgraph',
                          'pyproj',
                          'ipython',
                          'matplotlib',
                          'cx_Oracle']
     )
