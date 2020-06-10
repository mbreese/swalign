from setuptools import setup

setup(name='swalign',
      version='0.3.5',
      description='Smith-Waterman local aligner',
      author='Marcus Breese',
      author_email='marcus@breese.com',
      url='http://github.com/mbreese/swalign/',
      packages=['swalign'],
      scripts=['bin/swalign'],
      python_requires='>=3.1',

     )
