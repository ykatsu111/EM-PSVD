#!/usr/local/bin/python3

"""
EM-PSVD
(C) Y. Katsuyama (2020)
See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html
"""

import numpy as np
from pyempsvd import EmpsvdCore

dat = np.loadtxt("../data/psvd.csv", delimiter=",", dtype=float, comments="#")
em = EmpsvdCore(2, dat[:, 0], dat[:, 1], tol=1e-2)
em.fit()
print("----Fitting result----")
print(em.theta)
