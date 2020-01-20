from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools
import pybind11
import os


__version__ = "0.1"
__author__ = "ykatsu111"
__description__ = "A fitting algorithm for PSVD data with the EM algorithm."
__long_description__ = ""


if os.environ("ENABLE_CXX") == "OFF":

    setup(
        name="pyempsvd",
        version=__version__,
        author=__author__,
        description=__description__,
        long_description=__long_description__,
        package_dir={"": "python"},
        py_modules=["pyempsvd"],
        requires=["numpy", "scipy"]
    )

else:

    ext_modules = [
        Extension(
            "pyempsvd",
            ["cpp/EmpsvdCore.cpp", "cpp/PyEmpsvd.cpp"],
            include_dirs=[
                pybind11.get_include(),
                pybind11.get_include(True)
            ],
            language="c++",
            requires=["pybind11"],
            extra_compile_args=["-std=c++11"]
        )
    ]

    setup(
        name="pyempsvd",
        version=__version__,
        author=__author__,
        description=__description__,
        long_description=__long_description__,
        ext_modules=ext_modules,
        cmdclass={"build_ext": build_ext}
    )