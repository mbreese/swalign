from setuptools import setup, Extension
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

ext_modules = [
    Extension(
        r'swalign',
        [r'swalign/__init__.py']
    ),
]
    
setup(name='swalign',
      version='0.3.7',
      description='Smith-Waterman local aligner',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Marcus Breese',
      author_email='marcus@breese.com',
      packages=['swalign'],
      scripts=['bin/swalign'],
      python_requires='>=3.1',
      ext_modules=cythonize(ext_modules),
     )
