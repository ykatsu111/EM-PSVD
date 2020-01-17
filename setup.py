from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools
import pybind11


___version__ = "0.1"


ext_modules = [
    Extension(
        "pyempsvd",
        ["cpp/EmpsvdCore.cpp", "cpp/PyEmpsvd.cpp"],
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(True)
        ],
        language="c++"
    )
]