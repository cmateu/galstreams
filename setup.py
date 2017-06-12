import os
import sys
import re

try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup


setup(
    name='mwstreams',
    version='1.0.0',
    author='C. Mateu',
    author_email='cmateu.cida@gmail.com,cmateu@cida.gob.ve',
    packages=['bovy_coords','mwstreams','gcutils','mwstreams/lib'],
    package_data={'mwstreams/lib':['*.dat',
                                   '*.log',
                                   'streams_lib_notes.ipynb',
                                   'globular_cluster_params.harris2010.tableI.csv']},
    scripts=['bin/'],
    url='https://github.com/cmateu/mwstreams',
    license='LICENSE',
    description='MW stream library toolkit',
    long_description=open('README.md').read(),
    install_requires=[
      "numpy",
      "scipy"
    ],
)
