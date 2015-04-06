from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('trim_cython', ['trim_cython.pyx']),
               Extension('composition_cython', ['composition_cython.pyx']),
               Extension('find_polyA_cython', ['find_polyA_cython.pyx']),
               Extension('pausing_cython', ['pausing_cython.pyx']),
              ]

setup(
  name = 'cython stuff',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
