!============================================================================
! RTP Fig. 2 (analytic): P(0,t) and central curvature for D = 0.10, 0.175, 0.20
! Model (dimensionless): xdot = sigma(t) + sqrt(2D)*eta(t), sigma=±1 flips at rate 1
! Evaluates Eq. (10) via y = t*sin(theta), Simpson rule in theta.
!
! Outputs ('.' -> 'p' in names, as in your setup):
!   P0_vs_t_D0p10ptxt, P0_vs_t_D0p175ptxt, P0_vs_t_D0p20ptxt
!   curvature0_vs_t_D0p10ptxt, curvature0_vs_t_D0p175ptxt, curvature0_vs_t_D0p20ptxt
!   (columns noted below)
!============================================================================
program rtp_fig2_analytic
  implicit none
  ! ---------------- user knobs ----------------
  integer,  parameter :: ntheta = 2400   ! even number (Simpson); increase if you want
  real*8,   parameter :: tmin = 0.02d0   ! start > 0 to avoid singular t=0
  real*8,   parameter :: tmax = 2.00d0
  real*8,   parameter :: dt   = 0.02d0
  real*8,   parameter :: dxFD = 1.0d-3   ! finite-difference step for curvature at x=0
  ! --------------------------------------------

  integer :: nt, it
  real*8  :: Dlist(3)
  character(len=64) :: fP, fC

  ! Three D values used in the figure
  Dlist = (/ 0.10d0, 0.175d0, 0.20d0 /)

  nt = nint( (tmax - tmin)/dt ) + 1

  do it = 1, 3
     call write_one_D(Dlist(it))
  end do

contains

  subroutine write_one_D(D)
    implicit none
    real*8, intent(in) :: D
    integer :: it, uP, uC
    real*8 :: t, P0, Pm, Pp, d2P, d2P_scaled
    character(len=64) :: fP, fC

    write(fP,'("P0_vs_t_D",f0.3,".txt")') D
    write(fC,'("curvature0_vs_t_D",f0.3,".txt")') D
    call sanitize_filename(fP)   ! -> replace '.' by 'p'
    call sanitize_filename(fC)

    open(newunit=uP, file=trim(fP), status="replace", action="write")
    open(newunit=uC, file=trim(fC), status="replace", action="write")

    write(uP,*) "# t    P(0,t)         (D=",D,")"
    write(uC,*) "# t    d2P(0,t)       D^(5/2)*d2P(0,t)    (D=",D,")"

    t = tmin
    do it = 1, nint((tmax - tmin)/dt) + 1
       ! P at 0 and ±dx for finite-difference curvature
       P0 = rtp_pdf_eq10(D, t, 0.0d0, ntheta)
       Pp = rtp_pdf_eq10(D, t, +dxFD, ntheta)
       Pm = rtp_pdf_eq10(D, t, -dxFD, ntheta)
       d2P = (Pp - 2.0d0*P0 + Pm) / (dxFD*dxFD)
       d2P_scaled = d2P * D**(2.5d0)

       write(uP,'(2(1X,ES20.12))') t, P0
       write(uC,'(3(1X,ES20.12))') t, d2P, d2P_scaled

       t = t + dt
    end do

    close(uP); close(uC)
    print *,"Wrote:", trim(fP), "and", trim(fC)
  end subroutine write_one_D

  !==================== Analytic Eq. (10) evaluator ====================!
  function rtp_pdf_eq10(D, t, x, ntheta) result(px)
    implicit none
    real*8, intent(in) :: D, t, x
    integer, intent(in):: ntheta
    real*8 :: px, termA, theta, dth, s, c, y, G, sum, pi
    integer :: i

    pi = dacos(-1.0d0)

    ! Closed "shell" term (see paper’s Eq. 10; matches our earlier code)
    termA = dexp(-t - (x*x + t*t)/(4.0d0*D*t)) * dcosh(x/(2.0d0*D)) / dsqrt(4.0d0*pi*D*t)

    ! Integral over theta ∈ [-pi/2, pi/2], Simpson weights
    dth = pi / dble(ntheta)
    sum = 0.0d0
    do i = 0, ntheta
       theta = -0.5d0*pi + dth*dble(i)
       s = dsin(theta); c = dcos(theta)
       y = t*s
       G = dexp(- (x - y)*(x - y) / (4.0d0*D*t)) / dsqrt(4.0d0*pi*D*t)
       if (i==0 .or. i==ntheta) then
          sum = sum + G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       else if (mod(i,2)==0) then
          sum = sum + 2.0d0 * G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       else
          sum = sum + 4.0d0 * G * ( t*c*i0_mod(t*c) + t*i1_mod(t*c) )
       end if
    end do
    px = termA + 0.5d0 * dexp(-t) * (dth/3.0d0) * sum
  end function rtp_pdf_eq10

  !---------------- Bessel I0 / I1 (accurate double-precision) ----------------!
  pure function i0_mod(x) result(y)
    implicit none
    real*8, intent(in) :: x
    real*8 :: y, ax, t
    ax = dabs(x)
    if (ax < 3.75d0) then
      t = x/3.75d0; t = t*t
      y = 1.0d0 + t*(3.5156229d0 + t*(3.0899424d0 + t*(1.2067492d0 + t*(0.2659732d0 + t*(0.0360768d0 + t*0.0045813d0)))))
    else
      t = 3.75d0/ax
      y = (dexp(ax)/dsqrt(ax)) * (0.39894228d0 + t*(0.01328592d0 + t*(0.00225319d0 + t*(-0.00157565d0 + t*(0.00916281d0 + t*(-0.02057706d0 + t*(0.02635537d0 + t*(-0.01647633d0 + t*0.00392377d0))))))))
    end if
  end function i0_mod

  pure function i1_mod(x) result(y)
    implicit none
    real*8, intent(in) :: x
    real*8 :: y, ax, t
    ax = dabs(x)
    if (ax < 3.75d0) then
      t = x/3.75d0; t = t*t
      y = x*(0.5d0 + t*(0.87890594d0 + t*(0.51498869d0 + t*(0.15084934d0 + t*(0.02658733d0 + t*(0.00301532d0 + t*0.00032411d0))))))
    else
      t = 3.75d0/ax
      y = (dexp(ax)/dsqrt(ax)) * (0.39894228d0 + t*(-0.03988024d0 + t*(-0.00362018d0 + t*(0.00163801d0 + t*(-0.01031555d0 + t*(0.02282967d0 + t*(-0.02895312d0 + t*(0.01787654d0 - t*0.00420059d0))))))))
      if (x < 0.0d0) y = -y
    end if
  end function i1_mod

  !---------------- filename sanitizer: '.' -> 'p'; spaces -> '_' -------------!
  subroutine sanitize_filename(str)
    implicit none
    character(len=*), intent(inout) :: str
    integer :: i
    do i=1,len_trim(str)
      if (str(i:i)=='.') str(i:i)='p'
      if (str(i:i)==' ') str(i:i)='_'
    end do
  end subroutine sanitize_filename

end program rtp_fig2_analytic