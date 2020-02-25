program main
  use empsvd_core, only: init, fit, theta, M, niter, get_loglikelihood, get_aic, get_bic
  implicit none
  integer(8), parameter :: N=100 * 60, K=2, max_iter=1000
  real(8) :: x_in(N), y_in(N), z_in(N), r
  real(8), allocatable :: x(:), y(:), z(:)
  integer(8) :: i, info, L

  ! read data and parameters
  open(10, file="psvd_bin.csv", status="old")
  do i = 1, N
     read(10, *) x_in(i), y_in(i), z_in(i)
  end do
  close(10)

  ! mask out elements where z=0
  L = count( z_in > 0. )
  allocate( x(L) )
  allocate( y(L) )
  allocate( z(L) )
  x = pack( x_in, z_in > 0. )
  y = pack( y_in, z_in > 0. )
  z = pack( z_in, z_in > 0. )
  
  ! print initial condition
  call init(x, y, K, max_iter=max_iter, tol=1d-5, weight=z)
  write(*, "(A)") "----initial condition----"
  do i = 1, K
     write(*, "(A,I0,A,6F15.5,A))") "theta(", i, ")=(", theta(i, :), ")"
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
  do i = 1, K
     write(*, "(A,I0,A,6F15.5,A))") "theta(", i, ")=(", theta(i, :), ")"
  end do
  write(*, "(A,I0)") "number of iteration: ", niter
  call get_loglikelihood(r)
  write(*, "(A,F15.5)") "log-likelihood: ", r
  call get_aic(r)
  write(*, "(A,F15.5)") "           aic: ", r
  call get_bic(r)
  write(*, "(A,F15.5)") "           bic: ", r


  deallocate( x, y, z )

end program main
