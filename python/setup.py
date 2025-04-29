from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("tdavec.tdavec_core", ["tdavec/tdavec_core.pyx"],
              include_dirs=[numpy.get_include()])
]

setup(
    name="tdavec",
    version="0.1.43",
    packages=["tdavec"],
    ext_modules=cythonize(extensions),
    zip_safe=False,
    setup_requires=["Cython", "numpy"]
)
