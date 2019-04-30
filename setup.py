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
    packages=['bovy_coords','galstreams','gcutils','galstreams/lib'],
    package_data={'galstreams/lib':['*.dat',
                                   '*.log',
                                   'streams_lib_notes.ipynb',
                                   'globular_cluster_params.harris2010.tableI.csv']},
#    scripts=['bin/'],
    url='https://github.com/cmateu/galstreams',
    license='LICENSE',
    description='MW stream library toolkit',
    long_description=open('README.md').read(),
    install_requires=[
      "numpy",
      "scipy"
    ],
)
