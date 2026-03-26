!============================================================================
! RTP analytic P(x,t) from Eq. (10) — smooth curves like the paper (infinite line)
! Units: dimensionless (v = 1, gamma = 1). Only parameter is D.
! Integral evaluated with y = t*sin(theta), theta ∈ [-pi/2, pi/2] (Simpson).
!============================================================================
program rtp_eq10_analytic
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: nx = 1601     ! x grid points (odd helps plotting symmetry)
  integer, parameter :: ntheta = 2000 ! must be even for Simpson
  real(dp), parameter :: xmin = -2.5_dp, xmax = 2.5_dp

  real(dp), parameter :: Dlist(3) = (/ 0.03_dp, 0.10_dp, 0.20_dp /)  
  real(dp), parameter :: tlist(10) = (/ 0.2_dp,0.4_dp,0.6_dp,0.8_dp,1.0_dp, &
                                       1.2_dp,1.4_dp,1.6_dp,1.8_dp,2.0_dp /)

  real(dp) :: D, t
  integer :: iD, it, uout
  character(len=64) :: fname

  do iD = 1, 3
     D = Dlist(iD)
     do it = 1, 10
        t = tlist(it)
        write(fname,'("P_eq10_D",f0.2,"_t",f0.1,".txt")') D, t
        call sanitize_filename(fname)

        open(newunit=uout, file=trim(fname), status="replace", action="write")
        write(uout,*) "# x    P(x,t)    (analytic, D=",D,", t=",t,")"
        call write_pdf_line(uout, D, t, xmin, xmax, nx, ntheta)
        close(uout)
     end do
  end do
contains

  subroutine write_pdf_line(uout, D, t, xmin, xmax, nx, ntheta)
    implicit none
    integer, intent(in) :: uout, nx, ntheta
    real(dp), intent(in) :: D, t, xmin, xmax
    real(dp) :: dx, x, px
    integer :: ix
    dx = (xmax - xmin) / real(nx-1,dp)
    do ix = 0, nx-1
      x = xmin + dx*real(ix,dp)
      px = rtp_pdf_eq10(D, t, x, ntheta)
      write(uout,'(2(1X,ES20.12))') x, px
    end do
  end subroutine write_pdf_line

  function rtp_pdf_eq10(D, t, x, ntheta) result(px)
    ! P(x,t) = term_A + e^{-t}/2 * \int_{-pi/2}^{pi/2} G(x,t,theta) * [t cosθ I0(tc)+ t I1(tc)] dθ
    implicit none
    integer, intent(in) :: ntheta
    real(dp), intent(in) :: D, t, x
    real(dp) :: px, termA, theta, dth, s, c, y, G, sum
    integer :: j

    ! "ballistic shell" closed term:
    termA = dexp(-t - (x*x + t*t)/(4.0_dp*D*t)) * dcosh(x/(2.0_dp*D)) / dsqrt(4.0_dp*acos(-1.0_dp)*D*t)

    ! integral part (Simpson over theta ∈ [-π/2, π/2])
    dth = (acos(-1.0_dp)) / real(ntheta,dp)    ! (π) / ntheta ; interval width is π
    sum = 0.0_dp
    do j = 0, ntheta
       theta = -0.5_dp*acos(-1.0_dp) + dth*real(j,dp)
       s = dsin(theta); c = dcos(theta)          ! y = t*s ; sqrt(t^2-y^2)=t*c
       y = t*s
       ! Gaussian kernel
       G = dexp(-(x - y)*(x - y) / (4.0_dp*D*t)) / dsqrt(4.0_dp*acos(-1.0_dp)*D*t)
       ! integrand after substitution (no singularity):
       ! [t*c*I0(t*c) + t*I1(t*c)]
       ! weight for Simpson:
       if (j==0 .or. j==ntheta) then
          sum = sum + G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       else if (mod(j,2)==0) then
          sum = sum + 2.0_dp * G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       else
          sum = sum + 4.0_dp * G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       end if
    end do
    px = termA + 0.5_dp * dexp(-t) * (dth/3.0_dp) * sum
  end function rtp_pdf_eq10

  ! -------- Modified Bessel I0, I1 (accurate double-precision approximations) --------
  pure function i0_mod(x) result(y)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: y, ax, t
    ax = dabs(x)
    if (ax < 3.75_dp) then
      t = x/3.75_dp; t = t*t
      y = 1.0_dp + t*(3.5156229_dp + t*(3.0899424_dp + t*(1.2067492_dp + t*(0.2659732_dp + t*(0.0360768_dp + t*0.0045813_dp)))))
    else
      t = 3.75_dp/ax
      y = (dexp(ax)/dsqrt(ax)) * (0.39894228_dp + t*(0.01328592_dp + t*(0.00225319_dp + t*(-0.00157565_dp + t*(0.00916281_dp + t*(-0.02057706_dp + t*(0.02635537_dp + t*(-0.01647633_dp + t*0.00392377_dp))))))))
    end if
  end function i0_mod

  pure function i1_mod(x) result(y)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: y, ax, t
    ax = dabs(x)
    if (ax < 3.75_dp) then
      t = x/3.75_dp; t = t*t
      y = x*(0.5_dp + t*(0.87890594_dp + t*(0.51498869_dp + t*(0.15084934_dp + t*(0.02658733_dp + t*(0.00301532_dp + t*0.00032411_dp))))))
    else
      t = 3.75_dp/ax
      y = (dexp(ax)/dsqrt(ax)) * (0.39894228_dp + t*(-0.03988024_dp + t*(-0.00362018_dp + t*(0.00163801_dp + t*(-0.01031555_dp + t*(0.02282967_dp + t*(-0.02895312_dp + t*(0.01787654_dp - t*0.00420059_dp))))))))
      if (x < 0.0_dp) y = -y
    end if
  end function i1_mod

  ! replace '.' with 'p' for tidy filenames
  subroutine sanitize_filename(str)
    character(len=*), intent(inout) :: str
    integer :: i
    do i=1,len_trim(str)
      if (str(i:i)=='.') str(i:i)='p'
      if (str(i:i)==' ') str(i:i)='_'
    end do
  end subroutine sanitize_filename

end program rtp_eq10_analytic
