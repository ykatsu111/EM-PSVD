! EM-PSVD
! (C) Y. Katsuyama (2020)
! See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html

program main
  use empsvd_core, only: init, fit, theta, M, niter, get_loglikelihood, get_aic, get_bic
  implicit none
  integer(8), parameter :: N=900, K=2, max_iter=1000
  real(8) :: x(N), y(N), r
  integer(8) :: i, info

  ! read data and parameters
  open(10, file="../data/psvd.csv", status="old")
  read(10, *)
  do i = 1, N
     read(10, *) x(i), y(i)
  end do
  close(10)
  
  ! print initial condition
  call init(K, x, y, max_iter=max_iter, tol=1d-5)
  write(*, "(A)") "----initial condition----"
  write(*, "(A)") "                    omega              a              b         sigma2    alpha(mu+1)         lambda"
  do i = 1, K
     write(*, "(A,I0,A,6F15.5,A)") "theta(", i, ")=(", theta(i, :), ")"
  end do
  call get_loglikelihood(r)
  write(*, "(A,F15.5)") "log-likelihood: ", r
  call get_aic(r)
  write(*, "(A,F15.5)") "           aic: ", r
  call get_bic(r)
  write(*, "(A,F15.5)") "           bic: ", r

  ! fitting
  call fit(info)
  if ( info /= 0 ) then
     write(*, "(A)") "Failed to converge!"
     stop 99
  end if

  ! print result
  write(*, "(A)") "----Fitting result----"
  write(*, "(A)") "                    omega              a              b         sigma2    alpha(mu+1)         lambda"
  do i = 1, K
     write(*, "(A,I0,A,6F15.5,A)") "theta(", i, ")=(", theta(i, :), ")"
  end do
  write(*, "(A,I0)") "number of iteration: ", niter
  call get_loglikelihood(r)
  write(*, "(A,F15.5)") "log-likelihood: ", r
  call get_aic(r)
  write(*, "(A,F15.5)") "           aic: ", r
  call get_bic(r)
  write(*, "(A,F15.5)") "           bic: ", r
  
end program main
