#!/usr/bin/env python

from distutils.core import setup

setup(name='swalign',
      version='0.3.1',
      description='Smith-Waterman local aligner',
      author='Marcus Breese',
      author_email='marcus@breese.com',
      url='http://github.com/mbreese/swalign/',
      packages=['swalign'],
      scripts=['bin/swalign']
     )
