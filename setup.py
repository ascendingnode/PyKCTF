from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = "PyKCTF",
        ext_modules = cythonize(('PyKCTF.pyx'))
        )
