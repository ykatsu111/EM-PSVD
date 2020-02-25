from typing import Tuple, List, Union, Callable, Any, AnyStr, Optional, Iterable, Dict
from numpy import ndarray


class EmpsvdCore:
    M = 6

    def __init__(
        self,
        k: int,
        x: ndarray,
        y: ndarray,
        z: ndarray =None,
        theta0: ndarray =None,
        tol: float =1e-2,
        max_iter: int =1000,
        fix_ab: bool =False,
        fix_alpha: bool =False
    ) -> None:
        import numpy as np

        if x.shape != y.shape or x.ndim != 1:
            raise ValueError(
                "x and y must be same size 1-dimensional numpy array."
            )
        if not isinstance(k, int):
            raise ValueError("k must be integer larger than 1.")

        self._x = x  # 1-dim
        self._y = y  # 1-dim
        if z is None:
            self._z = np.ones(x.size, dtype=x.dtype)
        else:
            self._z = z

        self._n = x.size
        self._k = k  # scalar
        self._m = self.M
        self._tol = tol
        self._max_iter = max_iter
        self._fix_ab = fix_ab
        self._fix_alpha = fix_alpha

        if theta0 is None:
            self.theta = self.make_theta0(k, x, y)
        else:
            self.theta = theta0

        self.gamma = self._get_gamma(self.theta)
        
        self.niter = 0

        return
    
    @property
    def x(self) -> ndarray:
        import numpy as np
        return np.array(self._x)

    @property
    def y(self) -> ndarray:
        import numpy as np
        return np.array(self._y)

    @property
    def z(self) -> ndarray:
        import numpy as np
        return np.array(self._z)

    @property
    def k(self) -> int:
        return self._k

    @property
    def n(self) -> int:
        return self._n

    @property
    def m(self) -> int:
        return self._m

    @property
    def max_iter(self) -> int:
        return self._max_iter

    @property
    def tol(self) -> float:
        return self._tol

    @property
    def fix_alpha(self) -> bool:
        return self._fix_alpha

    @property
    def fix_ab(self) -> bool:
        return self._fix_ab

    def fit(self) -> None:
        """fit and estimate parameters"""
        l1 = self.get_loglikelihood()
        while 1:
            if self.niter >= self._max_iter:
                raise RuntimeError("Failed to converge!")
            self.e_step()
            self.m_step()
            self.niter += 1
            l2 = self.get_loglikelihood()
            if abs(l1 - l2) < self._tol:
                return

            l1 = l2

    def e_step(self):
        """E-step"""
        self.gamma = self._get_gamma()
        return

    def m_step(self):
        """M-step"""
        for j in range(self.k):
            if not self._fix_ab:
                self.theta[j, 2] = self._get_new_bk(j)
                self.theta[j, 1] = self._get_new_ak(j, self.theta[j, 2])

            self.theta[j, 3] = self._get_new_sigma2k(j, self.theta[j, 1], self.theta[j, 2])

            if not self._fix_alpha:
                self.theta[j, 4] = self._get_new_alphak(j)

            self.theta[j, 5] = self._get_new_lambdak(j, self.theta[j, 4])
            self.theta[j, 0] = self._get_new_omegak(j)

        return

    def get_aic(self) -> float:
        return self.get_loglikelihood() - (self.k * self.m)

    def get_bic(self) -> float:
        import numpy as np
        return self.get_loglikelihood() - (0.5 * self.k * self.m * np.log(self.n))

    def get_loglikelihood(self, theta: ndarray =None) -> float:
        import numpy as np
        return np.sum(
            self._get_logsum_pxy(theta) * self._z
        )

    @staticmethod
    def make_theta0(k: int, x: ndarray, y: ndarray, z: ndarray =None):
        import numpy as np

        if z is None:
            z = np.ones(x.size, dtype=x.dtype)

        theta = np.zeros([k, EmpsvdCore.M], dtype=x.dtype)
        n = self.z.sum()
        v_mean = (y * z).sum() / n
        v_var = (((y - v_mean) ** 2) * z).sum() / n
        d_mean = (x * z).sum() / n
        d_var = (((x - d_mean) ** 2) * z).sum() / n

        a1 = y[z > 0.].max() / d_mean
        a2 = v_mean / 4.

        if k > 1:
            for ik in range(k):
                theta[ik] = (
                    1 / k, 
                    np.exp(np.log(a2) + ((np.log(a1) - np.log(a2)) * ik) / (k - 1)),
                    ik / (k - 1),
                    v_var, (d_mean ** 2) / d_var, d_mean / d_var
                ) 
        else:
            theta[0] = (
                1.0, 
                np.exp((np.log(a1) + np.log(a2)) / 2),
                0.5,
                v_var, (d_mean ** 2) / d_var, d_mean / d_var
            ) 
        return theta

    @staticmethod
    def calc_log_pxy(
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        theta: ndarray
    ) -> ndarray:
        import numpy as np
        from scipy import special

        k = theta.shape[0]
        if isinstance(x, np.ndarray):
            p = np.zeros([k] + list(x.shape), dtype=theta.dtype)
        elif isinstance(x, float):
            p = np.zeros([k])
        else:
            raise TypeError("args x and y must be float or ndarray.")

        for j in range(k):
            theta_k = theta[j]
            p[j] = np.log(
                theta_k[0] * (theta_k[5] ** theta_k[4]) * (x ** (theta_k[4] - 1)) / (
                    np.sqrt(2 * np.pi * theta_k[3]) * special.gamma(theta_k[4])
                )
            ) - (
                ((y - (theta_k[1] * x ** theta_k[2])) ** 2) /
                (2 * theta_k[3]) +
                (theta_k[5] * x)
            )
        return np.moveaxis(p, 0, -1)

    @classmethod
    def calc_pxy(
        cls,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        theta: Union[float, ndarray]
    ) -> ndarray:
        import numpy as np
        return np.exp(cls.calc_log_pxy(x, y, theta))

    def _get_log_pxy(self, theta=None) -> ndarray:
        if theta is None:
            theta = self.theta
        return self.calc_log_pxy(self._x, self._y, theta)

    def _get_pxy(self, theta=None) -> ndarray:
        if theta is None:
            theta = self.theta
        return self.calc_pxy(self._x, self._y, theta=theta)

    def _get_logsum_pxy(self, theta=None) -> ndarray:
        import numpy as np

        log_pxy = self._get_log_pxy(theta=theta)
        max_log_pxy = np.max(log_pxy, axis=1)
        for j in range(self.k):
            log_pxy[:, j] -= max_log_pxy

        return np.log(
            np.sum(np.exp(log_pxy), axis=1)
        ) + max_log_pxy

    def _get_sum_pxy(self, theta=None) -> ndarray:
        import numpy as np
        return np.exp(self._get_logsum_pxy(theta=theta))

    def _get_gamma(self, theta=None) -> ndarray:
        gma = self._get_pxy(theta=theta)
        sum_pxy = self._get_sum_pxy(theta=theta)
        for j in range(self.k):
            gma[:, j] /= sum_pxy
        return gma

    def _get_new_omegak(self, ik: int) -> float:
        import numpy as np
        return np.sum(self.gamma[:, ik]) / self._z.sum()

    def _get_new_ak(self, ik: int, new_bk: float) -> float:
        import numpy as np
        gammak = self.gamma[:, ik]
        q1 = self._z * gammak * self._y * self._x ** new_bk
        q2 = self._z * gammak * self._x ** (2.0 * new_bk)
        return np.sum(q1) / np.sum(q2)

    def _get_new_bk(self, ik: int) -> float:
        import numpy as np
        from scipy.optimize import newton

        gammak = self.gamma[:, ik]

        def bkdot(bk):
            ak = self._get_new_ak(ik, bk)
            q = self._z * gammak * np.log(self._x) * (
                (self._y * self._x ** bk) - (ak * self._x ** (2 * bk))
            )
            return np.sum(q)

        return newton(bkdot, self.theta[ik, 2])

    def _get_new_sigma2k(self, ik: int, new_ak: float, new_bk: float) -> float:
        import numpy as np
        gammak = self.gamma[:, ik]
        q = self._z * gammak * (self._y - (new_ak * self._x ** new_bk)) ** 2
        return np.sum(q) / np.sum(self._z * gammak)

    def _get_new_alphak(self, ik: int):
        return self._get_new_alphak_by_invdigamma(ik, self.theta[ik, 4])

    def _get_new_lambdak(self, ik: int, new_alphak: float) -> float:
        import numpy as np
        gammak = self.gamma[:, ik]
        q = self._z * gammak * self._x
        return new_alphak * np.sum(self._z * gammak) / np.sum(q)

    def _get_new_alphak_by_invdigamma(self, ik: int, alphak0: float, max_iter: int =500, tol: float =1e-2) -> float:
        import numpy as np
        from scipy.special import digamma

        gammak = self.gamma[:, ik]

        def alky(x, gammak):
            q1 = np.sum(self._z * gammak)
            q2 = np.sum(self._z * gammak * np.log(x))
            x_mean = np.sum(self._z * gammak * x) / q1
            return np.log(x_mean) - (q2 / q1)

        y = alky(self._x, gammak)
        qalk0 = alphak0
        # qalk0 = (1. + np.sqrt(1. + (4. * y / 3.))) / (4. * y)
        
        for i in range(max_iter):
            dig = digamma(qalk0)
            qalk1 = qalk0 * (np.log(qalk0) - dig) / y
            if abs(qalk1 - qalk0) < tol:
                return qalk1
            qalk0 = qalk1
                
        raise RuntimeError(f"Failed to converge new alk with {max_iter:d} iteration")
