import os
import sys
import re

try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup


setup(
    name='galstreams',
    version='1.0.0',
    author='C. Mateu',
    author_email='cmateu@fisica.edu.uy',
    packages=['galstreams','galstreams/lib','galstreams/tracks'],
    package_data={'galstreams/lib':['*.dat',
                                   '*.log',
                                   'master*',
                                   'streams_lib_notes.ipynb',
                                   'globular_cluster_params.harris2010.tableI.csv'],
		  'galstreams/tracks':['*.ecsv'],
		  'galstreams/notebooks':['*ipy*']},
#    scripts=['bin/'],
    url='https://github.com/cmateu/galstreams',
    license='LICENSE',
    description='MW stream library toolkit',
    long_description=open('README.md').read(),
    install_requires=[
      "numpy",
      "scipy",
      "astropy",
    ],
)
