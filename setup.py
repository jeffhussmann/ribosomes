import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('trim_cython', ['trim_cython.pyx'], include_dirs=[numpy.get_include()]),
              ]

setup(
  name = 'cython stuff',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
