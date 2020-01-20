from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools
import pybind11
import os


__version__ = "0.1"


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
    author="ykatsu111",
    description="A fitting algorithm for PSVD data with the EM algorithm.",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext}
)