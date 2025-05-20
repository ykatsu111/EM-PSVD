! EM-PSVD
! (C) Y. Katsuyama (2020)
! See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html

module empsvd_static
  implicit none

  public :: calc_log_pxy, calc_pxy
  public :: calc_new_ak
  public :: bkdot, calc_new_bk_by_newton
  public :: alky, calc_new_alk_by_invdigamma
  public :: calc_new_lk

  public :: stop_with_error, linspace, get_rand
  public :: make_theta0

  real(8), public, parameter :: pi = 3.1415926535897932384626433832795

  interface digamma
     module procedure digamma_0d
  end interface


contains


  subroutine stop_with_error(msg)
    implicit none
    character(*), intent(in) :: msg

    write(0, '(A)') trim(msg)
    stop 99

  end subroutine stop_with_error


  subroutine linspace(left, right, num, ret)
    implicit none
    real(8)   , intent(in)  :: left, right
    integer(8), intent(in)  :: num
    real(8)   , intent(out) :: ret(num)
    real(8)    :: step
    integer(8) :: i

    step = (right - left) / real(num - 1)
    ret(1) = left
    do i = 2, size(ret)
       ret(i) = ret(i - 1) + step
    end do

  end subroutine linspace


  subroutine get_rand(x)
    implicit none
    real(8), intent(out) :: x(:)
    integer(4), allocatable :: seed(:)
    integer(4) :: c, sz

    call random_seed( size=sz )
    allocate( seed(sz) )
    call random_seed( get=seed )
    call system_clock( count=c )
    seed = seed + c
    call random_seed( put=seed )
    call random_number( x )
    deallocate( seed )

  end subroutine get_rand


  subroutine make_theta0(N, K, x, y, z, theta0)
    implicit none
    integer(8), parameter :: M = 6
    integer(8), intent(in) :: N, K
    real(8)   , intent(in) :: x(N), y(N), z(N)
    real(8)   , intent(out) :: theta0(K, M)
    real(8) :: y_mean, y_var, x_mean, x_var, a_upp, a_low, nr, kr
    real(8) :: lins(K)

    nr = sum(z)
    kr = real(K)

    y_mean = sum(y * z) / nr
    y_var = sum( ((y - y_mean) ** 2) * z ) / nr
    x_mean = sum(x * z) / nr
    x_var = sum( ((x - x_mean) ** 2) * z ) / nr

    theta0(:, 1) = 1d0 / kr                 ! omega  : mixing fraction
    theta0(:, 4) = y_var                    ! sigma^2: variance of terminal velocity
    theta0(:, 5) = (x_mean ** 2) / x_var    ! alpha  : the shape parameter + 1
    theta0(:, 6) = x_mean / x_var           ! lambda : the slope parameter

    a_upp = maxval(y, mask=z>0.) / x_mean
    a_low = y_mean / 4d0

    if ( K > 1 ) then

       call linspace(0d0, 1d0, K, lins)
       theta0(:, 2) = exp( log(a_low) + (log(a_upp / a_low) * lins) ) ! a: the intercept parameter
       theta0(:, 3) = lins ! b: the slope parameter

    else

       theta0(:, 2) = exp( log(a_upp * a_low) / 2d0 )
       theta0(:, 3) = 0.5d0

    end if

  end subroutine make_theta0


  subroutine calc_log_pxy(x, y, theta, log_pxy)
    implicit none
    real(8), intent(in)  :: x(:), y(:), theta(:, :)
    real(8), intent(out) :: log_pxy(:, :)
    integer(8) :: shp(2)
    integer(8) :: j, N, K, M

    N = size(x)
    shp = shape(theta)
    K = shp(1)
    M = shp(2)
    if (M /= 6) stop 99  ! Invalid shape of theta
    if ( size(y) /= N ) stop 99  ! Invalid size of y

    do j = 1, K
       log_pxy(:, j) = log( theta(j, 1) * (theta(j, 6) ** theta(j, 5)) * (x ** (theta(j, 5) - 1d0)) /  &
            &               (sqrt( 2. * pi * theta(j, 4) ) * gamma( theta(j, 5) )) )                   &
            &          - ( ((y - (theta(j, 2) * x ** theta(j, 3))) ** 2) /                             &
            &              (2d0 * theta(j, 4)) +                                                       &
            &              (theta(j, 6) * x) )
    end do

  end subroutine calc_log_pxy


  subroutine calc_pxy(x, y, theta, pxy)
    implicit none
    real(8), intent(in)  :: x(:), y(:), theta(:, :)
    real(8), intent(out) :: pxy(:, :)

    call calc_log_pxy(x, y, theta, pxy)
    pxy = exp(pxy)
  end subroutine calc_pxy


  subroutine calc_new_ak(bk, x, y, z, gammak, new_ak)
    implicit none
    real(8), intent(in)  :: bk
    real(8), intent(in)  :: x(:), y(:), z(:), gammak(:)
    real(8), intent(out) :: new_ak
    real(8) :: q1, q2

    q1 = sum( z * gammak * y * x ** bk )
    q2 = sum( z * gammak * x ** (2d0 * bk) )
    new_ak = q1 / q2
  end subroutine calc_new_ak


  function bkdot(bk, x, y, z, gammak) result(r)
    implicit none
    real(8), intent(in) :: bk
    real(8), intent(in) :: x(:), y(:), z(:), gammak(:)
    real(8) :: r, ak

    call calc_new_ak(bk, x, y, z, gammak, ak)
    r = sum( z * gammak * log(x) * &
         &   ((y * x ** bk) - (ak * x ** (2d0 * bk))) )
  end function bkdot


  subroutine calc_new_bk_by_newton(bk0, x, y, z, gammak, new_bk, info, bk1)
    implicit none
    real(8), intent(in)  :: bk0
    real(8), intent(in)  :: x(:), y(:), z(:), gammak(:)
    real(8), intent(out) :: new_bk
    integer(8), intent(out) :: info
    real(8), intent(in), optional  :: bk1
    integer(8), parameter :: max_iter = 100
    real(8)   , parameter :: tol = 1.48d-08, offset = 0.001
    real(8) :: qbk0, qbk1, qbk2, fqbk0, fqbk1, fqbk2
    integer(8) :: i

    qbk0 = bk0
    if ( present(bk1) ) then
       qbk1 = bk1
    else
       qbk1 = bk0 + offset
    end if

    fqbk0 = bkdot(qbk0, x, y, z, gammak)
    fqbk1 = bkdot(qbk1, x, y, z, gammak)

    do i = 1, max_iter
       qbk2 = ((qbk0 * fqbk1) - (qbk1 * fqbk0)) / (fqbk1 - fqbk0)
       fqbk2 = bkdot(qbk2, x, y, z, gammak)

       if ( abs(fqbk2) < tol ) then
          new_bk = qbk2
          info = 0
          return
       end if

       qbk0 = qbk1
       fqbk0 = fqbk1
       qbk1 = qbk2
       fqbk1 = fqbk2
    end do

    new_bk = qbk2
    info = 2

  end subroutine calc_new_bk_by_newton


  function alky(x, z, gammak) result(y)
    implicit none
    real(8), intent(in) :: x(:), z(:), gammak(:)
    real(8) :: y, x_mean, q1, q2

    q1 = sum(z * gammak)
    q2 = sum(z * gammak * log(x))
    x_mean = sum(z * gammak * x) / q1
    y = log(x_mean) - (q2 / q1)
    
  end function alky


  subroutine calc_new_alk_by_invdigamma(qalk0, x, z, gammak, new_alk, info)
    implicit none
    real(8), intent(in)  :: qalk0
    real(8), intent(in)  :: x(:), z(:), gammak(:)
    real(8), intent(out) :: new_alk
    integer(8), intent(out) :: info
    integer(8), parameter :: max_iter = 500
    real(8)   , parameter :: tol = 1d-02
    real(8) :: dig, trig, alk0, alk1, y
    integer(8) :: i

    y = alky(x, z, gammak)
    alk0 = qalk0
    ! alk0 = (1d0 + sqrt(1d0 + (y * 4d0 / 3d0))) / (4d0 * y)
    ! alk0 = exp(y) + 0.5d0

    do i = 1, max_iter
       if (alk0 > 1d10) then
         info = 4
         new_alk = alk0
         return
       end if
       call digamma(alk0, dig)
       alk1 = alk0 * ((log(alk0) - dig) / y)
       ! call trigamma(alk0, trig, n=100)
       ! alk1 = alk0 + ((y - dig) / trig)
       if (abs(alk1 - alk0) < tol) then
          info = 0
          new_alk = alk1
          return
       end if
       alk0 = alk1
    end do

    info = 3
    new_alk = alk1

  end subroutine calc_new_alk_by_invdigamma


  subroutine calc_new_lk(alk, x, z, gammak, new_lk)
    implicit none
    real(8), intent(in)  :: alk
    real(8), intent(in)  :: x(:), z(:), gammak(:)
    real(8), intent(out) :: new_lk
    real(8)              :: q1, q2

    q1 = sum( z * gammak )
    q2 = sum( z * gammak * x )
    new_lk = alk * q1 / q2
    
  end subroutine calc_new_lk


  subroutine digamma_0d(a, dig) 
    implicit none
    real(8), intent(in)  :: a
    real(8), intent(out) :: dig
    real(8)    :: a1, a2
    integer(8) :: i
    real(8), parameter :: c(6) = (/ 0.64493313d0,  &
         &                         -0.20203181d0,  &
         &                          0.08209433d0,  &
         &                         -0.03591665d0,  &
         &                          0.01485925d0,  &
         &                         -0.00472050d0 /)
    real(8), parameter :: g = 0.57721566490153286061



    a2 = real(int(a))  ! a >= 1
    a1 = a - a2  ! 0 <= a < 1

    if (a1 > 0d0) then
       ! first compute digamma value for 0 < a < 1
       dig = (a1 / (a1 + 1d0)) - g + (0.5d0 * a1 ** 7)
       do i = 1, size(c)
          dig = dig + (c(i) * ((a1 ** i) - a1 ** 7))
       end do
       dig = dig - (1d0 / a1)
       ! then increment upto a >= 1
       if (int(a2) > 0) then
          do i = 1, int(a2)
             dig = dig + (1d0 / (a1 + real(i) - 1d0))
          end do
       end if

    else
       ! a1 = 0 means a2 is integer, then digamma value is special value
       dig = -g
       if (int(a2) > 1) then
          do i = 1, int(a2) - 1
             dig = dig + (1d0 / real(i))
          end do
       end if
    end if

  end subroutine digamma_0d


end module empsvd_static
