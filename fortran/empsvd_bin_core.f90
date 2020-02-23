module empsvd_bin_core
  use empsvd_core, only: &
       ! override variables
       & x, y, theta, gamm, N, K, M, niter, &
       & max_iter_, tol_, fix_ab_, fix_alpha_
  implicit none

  ! mandatory variables to be niitialize by init subbroutine
  real(8), protected, allocatable :: z(:)

  ! mandatory subroutine
  public :: init, fit, e_step, m_step
  ! optional subroutine
  public :: get_loglikelihood, get_aic, get_bic

  private :: calc_gamma, calc_new_pk

contains

  
  subroutine init(x_in, y_in, z_in, k_in, theta0, max_iter, tol, fix_alpha, fix_ab)
    use empsvd_core, only: super_init => init 
    use empsvd_static, only: stop_with_error, make_theta0_bin
    real(8)   , intent(in) :: x_in(:), y_in(:), z_in(:)
    integer(8), intent(in) :: k_in
    real(8)   , intent(in), optional :: theta0(k_in, M)
    integer(8), intent(in), optional :: max_iter
    real(8)   , intent(in), optional :: tol
    logical   , intent(in), optional :: fix_alpha, fix_ab
    integer(8) :: shp(2), n_tmp
    real(8)    :: theta0_tmp(k_in, M)

    if ( size(x_in) /= size(y_in) ) call stop_with_error("data length of x and y must be same.")
    if ( size(x_in) /= size(z_in) ) call stop_with_error("data length of x, y, and z must be same.")

    if ( present(theta0) ) then
       theta0_tmp = theta0
    else
       n_tmp = size(x_in)
       call make_theta0_bin(n_tmp, k_in, x_in, y_in, z_in, theta0_tmp)
    end if

    call super_init(x_in, y_in, k_in, theta0_tmp, max_iter, tol, fix_alpha, fix_ab)
    if (allocated(z)) deallocate(z)
    allocate( z(N) )    

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
    use empsvd_core, only: calc_logsum_pxy
    implicit none
    real(8), intent(out) :: loglikelihood
    real(8) :: logsum_pxy(N)

    call calc_logsum_pxy(logsum_pxy)
    loglikelihood = sum(logsum_pxy * z)
  end subroutine get_loglikelihood


  subroutine calc_gamma(g)
    use empsvd_core, only: super_calc_gamma => calc_gamma
    implicit none
    real(8), intent(out) :: g(N, K)
    integer(8) :: j

    call super_calc_gamma(g)
    do j = 1, K
       where ( .not. z > 0.)
          g(:, j) = 0.
       end where
       g(:, j) = g(:, j) * z
    end do

  end subroutine calc_gamma


  subroutine e_step()
    implicit none

    call calc_gamma(gamm)

  end subroutine e_step


  subroutine m_step(info)
    use empsvd_core, only: &
         & calc_new_ak, calc_new_bk, calc_new_sk, &
         & calc_new_alk, calc_new_lk
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

    new_pk = sum( gammak ) / sum( z )
  end subroutine calc_new_pk

  
end module empsvd_bin_core
