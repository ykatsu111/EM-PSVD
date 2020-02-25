#!/usr/local/bin/python3

import numpy as np
from pyempsvd import EmpsvdCore

dat = np.loadtxt("../data/psvd.csv", delimiter=",")
em = EmpsvdCore(2, dat[:, 0], dat[:, 1], tol=1e-2)
em.fit()
print("----Fitting result----")
print(em.theta)
