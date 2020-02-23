program main
  use empsvd_bin_core, only: init, fit, theta, M
  implicit none
  integer(8), parameter :: N=100 * 60, K=2, max_iter=1000
  real(8) :: x(N), y(N), z(N)
  integer(8) :: i, info

  ! read data and parameters
  open(10, file="psvd_bin.csv", status="old")
  do i = 1, N
     read(10, *) x(i), y(i), z(i)
  end do
  close(10)
  
  ! EM Fitting
  call init(x, y, z, K, max_iter=max_iter, tol=1d-5)
  write(*, "(A)") "----initial condition----"
  do i = 1, K
     write(*, "(A,I0,A,6F15.5,A))") "theta(", i, ")=(", theta(i, :), ")"
  end do
  call fit(info)

  ! print result
  write(*, "(A)") "----Fitting result----"
  if ( info == 0 ) then
     do i = 1, K
        write(*, "(A,I0,A,6F15.5,A))") "theta(", i, ")=(", theta(i, :), ")"
     end do
  else
     write(*, "(A)") "Failed to converge!"
  end if
  
end program main
