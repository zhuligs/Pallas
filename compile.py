from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules = [
    Extension("zgg",  ["zgraph.py"],extra_compile_args=['-O3'],), 
#   ... all your modules that need be compiled ...
]
setup(
    name = 'zgg',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
