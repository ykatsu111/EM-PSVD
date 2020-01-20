module empsvd_core
  implicit none
  
  ! mandatory variables to be initialize by init subroutine
  real(8)   , protected, allocatable :: x(:), y(:)
  real(8)   , protected, allocatable :: theta(:, :)
  real(8)   , protected, allocatable :: gamm(:, :)
  integer(8), protected              :: N, K
  integer(8), public   , parameter   :: M = 6
  integer(8), protected              :: niter = 0

  real(8), private :: tol_ = 1d-2
  integer(8), private :: max_iter_ = 1000
  logical, private :: fix_ab_ = .false.
  logical, private :: fix_alpha_ = .true.

  ! mandatory subroutine
  public :: init, fit, e_step, m_step
  ! optional subroutine
  public :: get_loglikelihood, get_aic, get_bic

  private :: calc_new_pk, calc_new_ak, calc_new_bk, calc_new_sk, calc_new_lk
  private :: calc_log_pxy, calc_pxy, calc_logsum_pxy, calc_sum_pxy, calc_loglikelihood


contains


  subroutine init(x_in, y_in, k_in, theta0, max_iter, tol, fix_alpha, fix_ab)
    use empsvd_static, only: stop_with_error, make_theta0
    implicit none
    real(8)   , intent(in) :: x_in(:), y_in(:)
    integer(8), intent(in) :: k_in
    real(8)   , intent(in), optional :: theta0(:, :)
    integer(8), intent(in), optional :: max_iter
    real(8)   , intent(in), optional :: tol
    logical   , intent(in), optional :: fix_alpha, fix_ab
    integer(8) :: shp(2)

    if ( size(x_in) /= size(y_in) ) call stop_with_error("data length of x and y must be same.")
    if ( present(theta0) ) then
       shp = shape(theta0)
       if ( shp(1) /= k_in ) call stop_with_error("Invalid shape of theta0.")
       if ( shp(2) /= M ) call stop_with_error("Invalid shape of theta0.")
    end if
    
    if ( allocated(x) ) deallocate(x)
    if ( allocated(y) ) deallocate(y)
    if ( allocated(theta) ) deallocate(theta)
    if ( allocated(gamm) ) deallocate(gamm)

    N = size(x_in)
    K = k_in

    allocate( x(N) ) 
    allocate( y(N) )
    allocate( theta(K, M) )
    allocate( gamm(N, K) )

    if ( present(max_iter) ) max_iter_ = max_iter
    if ( present(tol) ) tol_ = tol
    if ( present(fix_alpha) ) fix_alpha_ = fix_alpha
    if ( present(fix_ab) ) fix_ab_ fix_ab

    x = x_in
    y = y_in

    if ( present(theta0) ) then
       theta = theta0
    else
       call make_theta0(theta)
    end if
       
    call calc_gamma(gamm)

    niter = 0

  end subroutine init


  subroutine fit(info)
    implicit none
    integer(8), intent(out) :: info
    real(8) :: l1, l2

    call get_loglikelihood(l1)
    do

       if ( niter > max_iter_ ) then
          info = 1
          return
       end if

       call e_step()
       call m_step(info)
       if (info > 0) return

       call get_loglikelihood(l2)

       if ( abs(l1 - l2) < tol_ ) then
          info = 0
          exit
       end if

       l1 = l2
       niter = niter + 1

    end do
    
  end subroutine fit


  subroutine get_bic(bic)
    implicit none
    real(8), intent(out) :: bic
    real(8) :: nn, kk, mm

    nn = real(N)
    kk = real(K)
    mm = real(M)
    
    call get_loglikelihood(bic)
    bic = bic - 0.5d0 * kk * mm* log(nn)

  end subroutine get_bic


  subroutine get_aic(aic)
    implicit none
    real(8), intent(out) :: aic
    real(8) :: kk, mm

    kk = real(K)
    mm = real(M)

    call get_loglikelihood(aic)
    aic = aic - kk * mm

  end subroutine get_aic


  subroutine get_loglikelihood(loglikelihood)
    implicit none
    real(8), intent(out) :: loglikelihood
    
    call calc_loglikelihood(loglikelihood)

  end subroutine get_loglikelihood


  subroutine calc_log_pxy(log_pxy, theta_in)
    use empsvd_static, only: calc_log_pxy_ => calc_log_pxy
    implicit none
    real(8), intent(out)           :: log_pxy(N, K)
    real(8), intent(in) , optional :: theta_in(K, M)

    if ( present(theta_in) ) then
       call calc_log_pxy_(x, y, theta_in, log_pxy)
    else
       call calc_log_pxy_(x, y, theta, log_pxy)
    end if
    
  end subroutine calc_log_pxy


  subroutine calc_pxy(pxy, theta_in)
    use empsvd_static, only: calc_pxy_ => calc_pxy
    implicit none
    real(8), intent(out) :: pxy(N, K)
    real(8), intent(in), optional :: theta_in(K, M)

    if ( present(theta_in) ) then
       call calc_pxy_(x, y, theta_in, pxy)
    else
       call calc_pxy_(x, y, theta, pxy)
    end if

  end subroutine calc_pxy


  subroutine calc_logsum_pxy(logsum_pxy, theta_in)
    implicit none
    real(8), intent(out) :: logsum_pxy(N)
    real(8), intent(in) , optional :: theta_in(K, M)
    real(8) :: log_pxy(N, K), dlog_pxy(N, K), max_log_pxy(N)
    integer(8) :: j

    if ( present(theta_in) ) then
       call calc_log_pxy(log_pxy, theta_in)
    else
       call calc_log_pxy(log_pxy)
    end if

    max_log_pxy = maxval( log_pxy, dim=2 )
    do j = 1, K
       dlog_pxy(:, j) = log_pxy(:, j) - max_log_pxy
    end do
    logsum_pxy = log( sum( exp(dlog_pxy), dim=2 ) ) + max_log_pxy

  end subroutine calc_logsum_pxy

  
  subroutine calc_sum_pxy(sum_pxy, theta_in)
    implicit none
    real(8), intent(out) :: sum_pxy(N)
    real(8), intent(in) , optional :: theta_in(K, M)
    
    if ( present(theta_in) ) then
       call calc_logsum_pxy(sum_pxy, theta_in)
    else
       call calc_logsum_pxy(sum_pxy)
    end if
    sum_pxy = exp(sum_pxy)

  end subroutine calc_sum_pxy


  subroutine calc_loglikelihood(loglikelihood, theta_in)
    implicit none
    real(8), intent(out) :: loglikelihood
    real(8), intent(in) , optional :: theta_in(K, M)
    real(8) :: logsum_pxy(N)

    if ( present(theta_in) ) then
       call calc_logsum_pxy(logsum_pxy, theta_in)
    else
       call calc_logsum_pxy(logsum_pxy)
    end if

    loglikelihood = sum( logsum_pxy )

  end subroutine calc_loglikelihood


  subroutine calc_gamma(g)
    implicit none
    real(8), intent(out) :: g(N, K)
    real(8) :: sum_pxy(N)
    integer(8) :: j

    call calc_sum_pxy(sum_pxy)
    call calc_pxy(g)
    
    do j = 1, K
       g(:, j) = g(:, j) / sum_pxy
    end do
    
  end subroutine calc_gamma


  subroutine e_step()
    implicit none

    call calc_gamma(gamm)

  end subroutine e_step


  subroutine m_step(info)
    implicit none
    integer(8), intent(out) :: info
    integer(8) :: j
    real(8)    :: tmp

    do j = 1, K

       if ( .not. fix_ab_ ) then
          call calc_new_bk(gamm(:, j), theta(j, 3), tmp, info)
          if ( info > 0 ) return
          theta(j, 3) = tmp
          call calc_new_ak(gamm(:, j), theta(j, 3), theta(j, 2))
       end if

       call calc_new_sk(gamm(:, j), theta(j, 2), theta(j, 3), theta(j, 4))

       if ( .not. fix_alpha_ ) then
          call calc_new_alk(gamm(:, j), theta(j, 5), tmp, info)
          if ( info > 0 ) return
          theta(j, 5) = tmp
       end if

       call calc_new_lk(gamm(:, j), theta(j, 5), theta(j, 6))
       call calc_new_pk(gamm(:, j), theta(j, 1))
       
    end do

    info = 0
    
  end subroutine m_step


  subroutine calc_new_pk(gammak, new_pk)
    implicit none
    real(8), intent(in)  :: gammak(N)
    real(8), intent(out) :: new_pk

    new_pk = sum( gammak ) / real(N)

  end subroutine calc_new_pk


  subroutine calc_new_ak(gammak, bk, new_ak)
    use empsvd_static, only: calc_new_ak_ => calc_new_ak
    implicit none
    real(8), intent(in)  :: gammak(N)
    real(8), intent(in)  :: bk
    real(8), intent(out) :: new_ak

    call calc_new_ak_(bk, x, y, gammak, new_ak)
    
  end subroutine calc_new_ak


  subroutine calc_new_bk(gammak, old_bk, new_bk, info)
    use empsvd_static, only: calc_new_bk_by_newton, get_rand
    implicit none
    real(8), intent(in)  :: gammak(N), old_bk
    real(8), intent(out) :: new_bk
    integer(8), intent(out) :: info
    integer(8), parameter :: retry = 30
    real(8)    :: bk1(retry)
    integer(8) :: i

    call calc_new_bk_by_newton(old_bk, x, y, gammak, new_bk, info)
    if ( info /= 0) then

       call get_rand(bk1)
       do i = 1, retry
!          write(*, "(A,f0.2)") "new_bk: retry from another bk1 ", old_bk + bk1(i)
          call calc_new_bk_by_newton(old_bk, x, y, gammak, new_bk, info, bk1(i))

          if ( info == 0 ) return
       end do
    end if

  end subroutine calc_new_bk


  subroutine calc_new_sk(gammak, ak, bk, new_sk)
    implicit none
    real(8), intent(in)  :: gammak(N)
    real(8), intent(in)  :: ak, bk
    real(8), intent(out) :: new_sk
    real(8) :: q1, q2

    q1 = sum( gammak * (y - (ak * x ** bk)) ** 2 )
    q2 = sum( gammak )
    new_sk = q1 / q2
    
  end subroutine calc_new_sk


  subroutine calc_new_alk(gammak, old_alk, new_alk, info)
    use empsvd_static, only: calc_new_alk_by_invdigamma
    implicit none
    real(8), intent(in)  :: gammak(N)
    real(8), intent(in)  :: old_alk
    real(8), intent(out) :: new_alk
    integer(8), intent(out) :: info

    call calc_new_alk_by_invdigamma(old_alk, x, gammak, new_alk, info)

  end subroutine calc_new_alk


  subroutine calc_new_lk(gammak, alk, new_lk)
    use empsvd_static, only: calc_new_lk_ => calc_new_lk
    implicit none
    real(8), intent(in)  :: gammak(N), alk
    real(8), intent(out) :: new_lk
    
    call calc_new_lk_(alk, x, gammak, new_lk)

  end subroutine calc_new_lk


end module empsvd_core
