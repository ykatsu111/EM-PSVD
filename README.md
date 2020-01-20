# Installation

```
pip3 install pybind11
CC=/usr/local/src/gcc-4.9/bin/gcc4.9 CXX=/usr/local/src/gcc-4.9/bin/g++4.9 python3 setup.py build_ext
```

or

```
pip3 install pybind11
CC=/usr/local/src/gcc-4.9/bin/gcc4.9 CXX=/usr/local/src/gcc-4.9/bin/g++4.9 pip3 install . --user
```

if gcc is not supported c++11,

```
ENABLE_CXX=OFF pip3 install . --user
```