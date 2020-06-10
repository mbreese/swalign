from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='swalign',
      version='0.3.6',
      description='Smith-Waterman local aligner',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Marcus Breese',
      author_email='marcus@breese.com',
      url='http://github.com/mbreese/swalign/',
      packages=['swalign'],
      scripts=['bin/swalign'],
      python_requires='>=3.1',

     )
