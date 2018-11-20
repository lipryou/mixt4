from distutils.core import setup, Extension

import numpy

setup(
    name="mixtures",
    version="1.0",
    ext_modules=[Extension("_mixtures", ["../src/mixt4.c",
                                         "src/mixt4module.c"])],
    include_dirs=[numpy.get_include()]
)
