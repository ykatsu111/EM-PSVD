# Introduction

This repositry provides some classes and modules to fit the mixed joint probability density function (PDF) to precipitation particle diameter velocity distribution (PSVD) data using the expectation maximazation (EM) algorithm and their sample source codes.  
The supported languages are now Python, C++, and fortran.  

In the fitting algorithm, we assume a mixed joint PDF of diameter ("D") and velocity ("V"):  
<a href="https://www.codecogs.com/eqnedit.php?latex=P(V,D|\mathbf{\theta})&space;=&space;\sum_{k=1}^{K}&space;\omega_k&space;Normal(V|\mu=a_k&space;D^{b_k},\sigma^2_k)&space;Gamma(D|\mu_k,&space;\lambda_k)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P(V,D|\mathbf{\theta})&space;=&space;\sum_{k=1}^{K}&space;\omega_k&space;Normal(V|\mu=a_k&space;D^{b_k},\sigma^2_k)&space;Gamma(D|\mu_k,&space;\lambda_k)" title="P(V,D|\mathbf{\theta}) = \sum_{k=1}^{K} \omega_k Normal(V|\mu=a_k D^{b_k},\sigma^2_k) Gamma(D|\mu_k, \lambda_k)" /></a>  
, where  
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}=(\omega_1,&space;...,&space;\omega_K,&space;a_1,&space;...,&space;a_K,&space;b_1,&space;...,&space;b_K,&space;\sigma^2_1,&space;...,&space;\sigma^2_K,&space;\mu_1,&space;...,&space;\mu_K,&space;\lambda_1,&space;...,&space;\lambda_K)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}=(\omega_1,&space;...,&space;\omega_K,&space;a_1,&space;...,&space;a_K,&space;b_1,&space;...,&space;b_K,&space;\sigma^2_1,&space;...,&space;\sigma^2_K,&space;\mu_1,&space;...,&space;\mu_K,&space;\lambda_1,&space;...,&space;\lambda_K)" title="\mathbf{\theta}=(\omega_1, ..., \omega_K, a_1, ..., a_K, b_1, ..., b_K, \sigma^2_1, ..., \sigma^2_K, \mu_1, ..., \mu_K, \lambda_1, ..., \lambda_K)" /></a>.  
Here, "ω" is a mixing fraction, "a" and "b" are parameters of the velocity-diameter relationship  
<a href="https://www.codecogs.com/eqnedit.php?latex=V=aD^b" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V=aD^b" title="V=aD^b" /></a>  
, "σ2" is a variance of velocity distribution assumed as the Normal distribution, "μ" and "λ" is respectively the shape and slope parameter of the diameter distribution assumed as the Gamma distribution, and "K" is the number of PDF elements. The fitting algorithm in this repositry can provide the optimal parameter set of "θ" giving the number of PDF elements.

Please see [Katsuyama and Inatsu (2020)]() for details of the algorithm.

Katsuyama, Y., and M. Inatsu, 2020: Fitting precipitation particle size-velocity data to mixed joint probability density function with an expectation maximization algorithm, J. Atmos. Oceanic Tech., in revision

# Licence and Aggreement

The source codes of this repositry is distributed under the XXX Licence.  
You must cite [Katsuyama and Inatsu (2020)]() in an appropriate way when you present and/or publish scientific results, products, and so on, obtained by using this fitting algorithm.

# For python users

## Requirements

The python module requires followings.

1. Python3
2. Numpy
3. Scipy
4. pybind11 (if your environment supports C++11.)
5. Eigen3 (if your environment supports C++11.)

## Installation

If your environment supports C++11, you can install the python module implemented with C++ by following commands after installing Eigen3.

```
pip3 install pybind11
wget XXX
cd XXX
pip3 install . --user
```

C/C++ compiler can be optionally specified with the enviroment variables of `CC` and `CXX` as follow.

```
CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++ pip3 install . --user 
```

If your environment does not support C++11, the python module implemented with pure python can be installed by following commands.

```
wget XXX
cd XXX
ENABLE_CXX=OFF pip3 install . --user
```

Or, you can use the module by simply putting [pyempsvd.py](python/pyempsvd.py) in a directory.

## Demonstration

A demonstration has been prepared in [sample/python](sample/python): [psvd.py](sample/python/psvd.py) and [psvd_bin.py](sample/python/psvd_bin.py) provides a demonstration for pure PSVD data and a demonstration for binned PSVD data, respectively. Here, we explain the usage of class following the sample demonstations.

### A case for pure PSVD data

Here, we demonstrate a case for fitting to pure PSVD data, enclosed at [sample/data/psvd.csv](sample/data/psvd.csv). This sample data contains diameter and velocity at first and second column, respectively. Therefore, the PSVD can be illustrated as follow.


```Python
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("./sample/data/psvd.csv", delimiter=",")
plt.scatter(
    data[:, 0], data[:, 1], c="k", s=5, marker="."
)
plt.xlim([0, 30])
plt.ylim([0, 5])
plt.xlabel("Diameter (mm)")
plt.ylabel("Velocity (m/s)")
plt.show()
```

<img src="fig/psvd.png" />

To fit this PSVD data, a EmpsvdCore object is firstly constructed.

```python
from pyempsvd import EmpsvdCore
em = EmpsvdCore(
    2, data[:, 0], data[:, 1]
)
```

In this case, the object has been constructed with the number of PDF elements `K=2`, specified at the first argument. The second and third agrument is respectively diameter and velocity of data. In this case, the initial parameters has been automatically generated using `EmpsvdCore.make_theta0`, a static method to generate an initial parameter sets based on the PSVD data given.  

If you would like to specify an initial parameter set explicitly, you can use a keyword argument `theta0`.

```python
theta0 = np.array([
    # omega, a, b, sigma2 (m^2 s^-2), mu+1, lambda (mm^-1) 
    [0.5, 0.3, 0.2, 0.04, 2.5, 0.1],  # 1st PDF element
    [0.5, 1.0, 0.8, 0.04, 2.5, 0.1]   # 2nd PDF element
])
em = EmpsvdCore(
    2, data[:, 0], data[:, 1], theta0=theta0
)
```

`theta0` must be an numpy ndarray with a shape of (K, 6) containing a parameter set for individual PDF element with an order of ω, a, b, σ2, μ+1, and λ with axis=1.  
As the other keyword arguments, the followings can be liste.

* tol: a tolerance parameter for the convergence judgement based on log-likelihood.
* max_iter: a maximum number of iteration in the fitting.
* fix_ab: fix parameters a and b at the initial state in the fitting.
* fix_alpha: fix a parameter μ at the initial state in the fitting.

Once the EmpsvdCore object has been constructed, the fitting can be easily done calling `fit` method.

```python
em.fit()
```

If the fitting has been successfully done, the optimal parameters set can be obtained by `em.theta` as a numpy ndarray.

```python
print(em.theta)
```

Otherwise, `RuntimeError` will be raise. The shape and element order are same as `theta0`.  
The fitting result can be illustrated as follow.

```python
plt.scatter(
    em.x, em.y, c="k", s=5, marker="."
)
x = np.linspace(0.01, 30, 100)
y = np.linspace(0.01, 5, 100)
xx, yy = np.meshgrid(x, y)
plt.contour(
    xx, yy, em.calc_pxy(xx, yy, em.theta).sum(axis=-1), 
    levels=[0.001, 0.01, 0.1], colors="r"
)
plt.xlim([0, 30])
plt.ylim([0, 5])
plt.xlabel("Diameter (mm)")
plt.ylabel("Velocity (m/s)")
plt.show()
```

<img src="fig/fitted_psvd.png" />

Here, `calc_pxy` method returns probability density of each PDF element corresponding to diemaeter (=`xx` in this case) and velocity (=`yy` in this case), so that we need to sum up the all PDF elements by `.sum(axis=-1)` in order to desplay the mixed PDF.

### A case for binned PSVD data

Here, we demonstrate a case for fitting to binned PSVD data, enclosed at [sample/data/psvd_bin.csv](sample/data/psvd_bin.csv). This sample data contains a number of particles (at third column) for a bin of diameter (at first column) and velocity (at second column), converted from [the pure PSVD data](sample/data/psvd.csv). The number of bins is 100 for diameter axis and 60 for velocity axis. Therefore, the PSVD can be illustrated as follow.

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
data = np.loadtxt("./sample/data/psvd_bin.csv", delimiter=",")
xx = data[:, 0].reshape(100, 60)
yy = data[:, 1].reshape(100, 60)
zz = data[:, 2].reshape(100, 60)
cmap = ListedColormap([[i/5] * 3 for i in range(5, 0, -1)])
bounds = [0, 1, 2, 3, 5, 10]
norm = BoundaryNorm(bounds,cmap.N)
plt.pcolor(
    xx, yy, zz, 
    cmap=cmap, norm=norm
)
plt.colorbar()
plt.xlim([0, 30])
plt.ylim([0, 5])
plt.xlabel("Diameter (mm)")
plt.ylabel("Velocity (m/s)")
plt.show()
```

<img src="fig/psvd_bin.png" />

For this binned PSVD data, the fitting can be done as following code.

```python
from pyempsvd import EmpsvdCore
mask = data[:, 2] > 0
em = EmpsvdCore(
    2, data[mask, 0], data[mask, 1], data[mask, 2]
)
em.fit()
```

Here, the EmpsvdCore object is constructed with the data where `z > 0` in order to avoid the zero division error. The fitting result is then illustrated as follow.

```python
xx = data[:, 0].reshape(100, 60)
yy = data[:, 1].reshape(100, 60)
zz = data[:, 2].reshape(100, 60)
cmap = ListedColormap([[i/5] * 3 for i in range(5, 0, -1)])
bounds = [0, 1, 2, 3, 5, 10]
norm = BoundaryNorm(bounds,cmap.N)
plt.pcolor(
    xx, yy, zz, 
    cmap=cmap, norm=norm
)
plt.colorbar()
x = np.linspace(0.01, 30, 100)
y = np.linspace(0.01, 5, 100)
xx, yy = np.meshgrid(x, y)
plt.contour(xx, yy, em.calc_pxy(xx, yy, em.theta).sum(axis=-1), levels=[0.001, 0.01, 0.1], colors="r")
plt.xlim([0, 30])
plt.ylim([0, 5])
plt.xlabel("Diameter (mm)")
plt.ylabel("Velocity (m/s)")
plt.show()
```

<img src="fig/fitted_psvd_bin.png" />

# For C++ users

## Requirements

The C++ class requires followings.

1. C++11
2. Eigen3

## Demonstration

A demonstration has been prepared in [sample/cpp](sample/cpp): [psvd.cpp](sample/cpp/psvd.cpp) and [psvd_bin.cpp](sample/cpp/psvd_bin.cpp) provides a demonstration for pure PSVD data and a demonstration for binned PSVD data, respectively.

```
cd sample/cpp
```

To compile the demonstrations, please modify `Makefile` following your environment and run `make pre` and `make` commands.

```
make pre
make
```

Then, you can see the executable files, `psvd.exe` and `psvd_bin.exe`. 

The usage is mostly same as the python. First, construct `EmpsvdCore` object after diameter (=`x`) and velocity (=`y`) data is obtained.

```cpp
constexpr int COLUMNS = 2;
constexpr int ROWS = 900;
io::CSVReader<COLUMNS> in("../data/psvd.csv");
Eigen::ArrayXd x(ROWS), y(ROWS);
double xi, yi;
Eigen::Index i = 0;	
while (in.read_row(xi, yi)) {
    x(i) = xi;
    y(i) = yi;
    i++;
}
Empsvd::EmpsvdCore em(2, x, y);
```

The fitting can be easily conduct by calling `em.fit()`.

```cpp
em.fit();
```

If the fitting has not been successfully done, `std::runtime_error` will be thrown. The fitting result can be contained in `em.theta` variable as well as the python.

If the PSVD data is binned data, the code is slightly modified.

```cpp
constexpr int COLUMNS = 3;
constexpr int ROWS = 100 * 60;
io::CSVReader<COLUMNS> in("../data/psvd_bin.csv");
Eigen::ArrayXd x_in(ROWS), y_in(ROWS), z_in(ROWS);
double xi, yi, zi;
Eigen::Index i = 0;

while (in.read_row(xi, yi, zi)) {
    x_in(i) = xi;
    y_in(i) = yi;
    z_in(i) = zi;
    i++;
}
Eigen::ArrayXd exist_data = (z_in.array() > 0.).select(					     
    Eigen::ArrayXd::Constant(z_in.size(), 1),
    Eigen::ArrayXd::Constant(z_in.size(), 0)
);
size_t n = exist_data.count();
Eigen::ArrayXd x(n), y(n), z(n);
Eigen::Index j = 0;
for (i = 0; i < x_in.size(); i++) {
        if (exist_data(i)) {
            x(j) = x_in(i);
            y(j) = y_in(i);
            z(j) = z_in(i);
            j++;
        }
}

Empsvd::EmpsvdCore em(2, x, y, z);
```

Please note that the elements where `z=0` should be removed in advance to construct `EmpsvdCore` object in order to avoid the zero division error.

# For fortran uesrs

## Requirements

The fortran module works under the normal fortran environment, fotran90/95.

## Demonstration

A demonstration has been prepared in [sample/fortran](sample/fortran): [psvd.f90](sample/fortran/psvd.f90) and [psvd_bin.f90](sample/fortran/psvd_bin.f90) provides a demonstration for pure PSVD data and a demonstration for binned PSVD data, respectively.

```
cd sample/fortran
```

To compile the demonstrations, please modify `Makefile` following your environment and run `make pre` and `make` commands.

```
make pre
make
```

Then, you can see the executable files, `psvd.exe` and `psvd_bin.exe`. 

The usage is mostly same as the python and C++. First, initialize the module by calling `init` function after diameter (=`x`) and velocity (=`y`) data is obtained.

```f90
use empsvd_core, only: init, fit, theta
implicit none
integer(8), parameter :: N=900, K=2
real(8) :: x(N), y(N), r
integer(8) :: i, info

open(10, file="../data/psvd.csv", status="old")
do i = 1, N
    read(10, *) x(i), y(i)
end do
close(10)

call init(K, x, y)
```

The fitting can be easily done by calling `fit` function.

```f90
call fit(info)
```

If the fitting has been successfully done, `info` will be `0`. The fitting result can be contained in `theta` variable as well as python and C++ demonstrations.

For the binned PSVD data, the code is slightly modified.

```f90
use empsvd_core, only: init, fit, theta, M, niter
implicit none
integer(8), parameter :: N=100 * 60, K=2, max_iter=1000
real(8) :: x_in(N), y_in(N), z_in(N), r
real(8), allocatable :: x(:), y(:), z(:)
integer(8) :: i, info, L

open(10, file="../data/psvd_bin.csv", status="old")
do i = 1, N
    read(10, *) x_in(i), y_in(i), z_in(i)
end do
close(10)

L = count( z_in > 0. )
allocate( x(L) )
allocate( y(L) )
allocate( z(L) )
x = pack( x_in, z_in > 0. )
y = pack( y_in, z_in > 0. )
z = pack( z_in, z_in > 0. )

call init(K, x, y, z)
```

Please note that the elements where `z=0` should be removed before calling `init` in order to avoid the zero division error.
