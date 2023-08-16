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
      version = '0.4.0.6',
      description = "A tool for processing infrasound and seismic array data.",
      keywords=['infrasound', 'geophysics', 'seismic', 'array'],
      author = "LANL Seismoacoustics Infrasound (LANL-SA) Team",
      author_email = 'pblom@lanl.gov',
      packages = ['infrapy',
                  'infrapy.association',
                  'infrapy.characterization',
                  'infrapy.cli',
                  'infrapy.detection',
                  'infrapy.location',
                  'infrapy.performance',
                  'infrapy.propagation',
                  'infrapy.resources',
                  'infrapy.resources.travelTimeTables',
                  'infrapy.utils',
                  'InfraView',
                  'InfraView.widgets',
                  'InfraView.graphics'],

      entry_points = {'console_scripts':['infrapy=infrapy.cli.__main__:main'],
                      'gui_scripts':['infraview = InfraView.__main__:main']},

      
      package_dir={'infrapy.propagation' : 'infrapy/propagation'},
      package_data={'infrapy.propagation' : ['compass.png'],
                    'infrapy.resources' : ['default.config'],
                    'infrapy.resources.travelTimeTables' : ['ak135_1st_arrivals.dat']},
                    

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
