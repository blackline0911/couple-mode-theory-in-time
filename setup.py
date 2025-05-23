from setuptools import setup
from Cython.Build import cythonize
import numpy

"""
python setup.py build_ext --inplace
"""
setup(ext_modules=cythonize("cython_test.pyx"),
      include_dirs=numpy.get_include())

setup(ext_modules=cythonize("raise_cosine.pyx"),
      include_dirs=numpy.get_include())