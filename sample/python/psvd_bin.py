#!/usr/local/bin/python3

"""
EM-PSVD
(C) Y. Katsuyama (2020)
See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html
"""

import numpy as np
from pyempsvd import EmpsvdCore

dat = np.loadtxt("../data/psvd_bin.csv", delimiter=",", dtype=float, comments="#")
mask = dat[:, 2] > 0
em = EmpsvdCore(
    2, dat[mask, 0], dat[mask, 1], dat[mask, 2]
)
em.fit()
print("----Fitting result----")
print(em.theta)
