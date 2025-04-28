"""
EM-PSVD
(C) Y. Katsuyama (2020)
See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html
"""


from setuptools import setup
import setuptools
import os


__version__ = "1.0"
__author__ = "Y. Katsuyama"
__description__ = "A fitting algorithm for PSVD data with the EM algorithm."
__long_description__ = """
EM-PSVD
(C) Y. Katsuyama (2020)
See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html
"""


if "ENABLE_CXX" in os.environ and os.environ["ENABLE_CXX"] == "OFF":

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

    from setuptools import Extension, find_packages
    from setuptools.command.build_ext import build_ext
    import pybind11

    ext_modules = [
        Extension(
            "pyempsvd",
            ["cpp/EmpsvdCore.cpp", "cpp/PyEmpsvd.cpp"],
            include_dirs=[
                pybind11.get_include(),
                pybind11.get_include(True)
            ],
            language="c++",
            requires=["numpy", "pybind11"],
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
        packages=find_packages(),
        cmdclass={"build_ext": build_ext}
    )
