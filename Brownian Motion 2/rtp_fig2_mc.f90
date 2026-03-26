!============================================================================
! RTP Fig. 2 via MONTE CARLO (infinite line, v=1, gamma=1; dimensionless)
! - Event-driven Poisson tumbles; exact Brownian increments
! - Histograms on x in [-2,2]; center bin sits at x=0 (nbins odd)
! - Outputs: P(0,t) and central curvature with scaled column D^(5/2)*curv
!
! Files written ('.'->'p', '.txt'->'ptxt'):
!   P0_vs_t_D0p10ptxt, P0_vs_t_D0p175ptxt, P0_vs_t_D0p20ptxt
!   curvature0_vs_t_D0p10ptxt, ... (columns: t, d2P, D^(5/2)d2P)
!
! Toggle options:
!   logical, parameter :: USE_QUADFIT = .true.   ! parabola fit on 5 bins
!   logical, parameter :: SMOOTH_TIME = .true.   ! 3-pt moving average over t
!
! Requires: mt.f90 providing  double-precision function grnd()
!============================================================================
program rtp_fig2_mc
  implicit none
  !---------------------- user knobs ----------------------
  integer,  parameter :: Ntraj  = 1000000      ! raise if you can (e.g. 2e6)
  real*8,   parameter :: xmin   = -2.0d0
  real*8,   parameter :: xmax   =  2.0d0
  integer,  parameter :: nbins  = 401          ! ODD → unique center bin
  real*8,   parameter :: tmin   = 0.02d0
  real*8,   parameter :: tmax   = 2.00d0
  real*8,   parameter :: dt_obs = 0.02d0
  logical,  parameter :: USE_QUADFIT = .true.
  logical,  parameter :: SMOOTH_TIME = .true.
  !-------------------------------------------------------

  integer :: ntimes, k, Didx, seed
  real*8  :: Dlist(3), D
  real*8, allocatable :: tgrid(:)
  real*8  :: binw
  integer :: idx0
  integer, allocatable :: hist(:,:)   ! (nbins, ntimes)

  character(len=64) :: fP, fC
  integer :: uP, uC

  ! RNG from mt.f90
  real*8, external :: grnd

  ! ----- build time grid -----
  ntimes = nint((tmax - tmin)/dt_obs) + 1
  allocate(tgrid(ntimes))
  do k = 1, ntimes
     tgrid(k) = tmin + dt_obs * dble(k-1)
  end do

  ! ----- histogram geometry -----
  binw = (xmax - xmin) / dble(nbins)
  idx0 = int((0.0d0 - xmin)/binw) + 1
  if (mod(nbins,2) == 0) then
     print *, "ERROR: nbins must be ODD so that center bin sits at x=0."
     stop
  end if
  if (idx0 < 3 .or. idx0 > nbins-2) then
     print *, "ERROR: center index too close to edge; adjust xrange/nbins."
     stop
  end if

  ! three D values
  Dlist = (/ 0.10d0, 0.175d0, 0.20d0 /)

  ! RNG seed
  seed = 20250821
  call sgrnd(seed)

  do Didx = 1, 3
     D = Dlist(Didx)
     print *,"[D=", D, "] MC trajectories:", Ntraj
     call simulate_one_D(D, tgrid, ntimes, binw, hist)
     call write_P0_and_curv(D, tgrid, ntimes, binw, idx0, hist)
     deallocate(hist)
  end do

contains

  subroutine simulate_one_D(D, tgrid, ntimes, binw, hist)
    implicit none
    real*8, intent(in)  :: D, tgrid(ntimes), binw
    integer, intent(in) :: ntimes
    integer, allocatable, intent(out) :: hist(:,:)

    integer :: n, k, idx
    real*8  :: x, tcur, sigma, next_flip, flip_wait, tobs, dt, g
    real*8, external :: grnd

    allocate(hist(nbins, ntimes)); hist = 0

    do n = 1, Ntraj
       ! initial state
       x    = 0.0d0
       tcur = 0.0d0
       sigma = merge(1.0d0, -1.0d0, grnd() >= 0.5d0)
       call exponential(flip_wait, 1.0d0)       ! rate = 1
       next_flip = tcur + flip_wait

       do k = 1, ntimes
          tobs = tgrid(k)

          ! march to observation time, handling flips exactly
          do
             if (next_flip < tobs) then
                dt = next_flip - tcur
                if (dt > 0.0d0) then
                   call gaussian(g)             ! N(0,1)
                   x = x + sigma*dt + dsqrt(2.0d0*D*dt)*g
                end if
                tcur = next_flip
                sigma = -sigma
                call exponential(flip_wait, 1.0d0)
                next_flip = tcur + flip_wait
             else
                dt = tobs - tcur
                if (dt > 0.0d0) then
                   call gaussian(g)
                   x = x + sigma*dt + dsqrt(2.0d0*D*dt)*g
                end if
                tcur = tobs
                exit
             end if
          end do

          ! bin
          if (x > xmin .and. x < xmax) then
             idx = int((x - xmin)/binw) + 1
             if (idx >= 1 .and. idx <= nbins) hist(idx,k) = hist(idx,k) + 1
          end if
       end do

       if (mod(n, 100000) == 0) print *,"  ...", n, "done"
    end do
  end subroutine simulate_one_D

  subroutine write_P0_and_curv(D, tgrid, ntimes, h, idx0, hist)
    implicit none
    real*8, intent(in)  :: D, tgrid(ntimes), h
    integer, intent(in) :: ntimes, idx0
    integer, intent(in) :: hist(nbins, ntimes)

    real*8, allocatable :: P0(:), curv(:), curv_scaled(:)
    real*8 :: Pm2,Pm1,P0c,Pp1,Pp2, a,b,c, x1,x2,x3,x4,x5, y1,y2,y3,y4,y5
    integer :: k, uP, uC
    character(len=64) :: fP, fC

    allocate(P0(ntimes), curv(ntimes), curv_scaled(ntimes))

    do k = 1, ntimes
       ! PDF heights from counts (density = counts / (Ntraj * binwidth))
       Pm2 = dble(hist(idx0-2,k)) / (dble(Ntraj)*h)
       Pm1 = dble(hist(idx0-1,k)) / (dble(Ntraj)*h)
       P0c = dble(hist(idx0  ,k)) / (dble(Ntraj)*h)
       Pp1 = dble(hist(idx0+1,k)) / (dble(Ntraj)*h)
       Pp2 = dble(hist(idx0+2,k)) / (dble(Ntraj)*h)
       P0(k) = P0c

       if (USE_QUADFIT) then
          ! Quadratic fit through 5 central bin centers: x = {-2h,-h,0,h,2h}
          x1=-2.0d0*h; x2=-1.0d0*h; x3=0.0d0; x4=1.0d0*h; x5=2.0d0*h
          y1=Pm2; y2=Pm1; y3=P0c; y4=Pp1; y5=Pp2
          call quadfit_5(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5, a,b,c)
          curv(k) = 2.0d0*a                      ! second derivative at 0
       else
          ! 5-point central difference (O(h^4))
          curv(k) = (-Pm2 + 16.0d0*Pm1 - 30.0d0*P0c + 16.0d0*Pp1 - Pp2) / (12.0d0*h*h)
       end if

       curv_scaled(k) = curv(k) * D**(2.5d0)
    end do

    if (SMOOTH_TIME) then
       call smooth3(P0, ntimes)
       call smooth3(curv, ntimes)
       call smooth3(curv_scaled, ntimes)
    end if

    write(fP,'("P0_vs_t_D",f0.3,".txt")') D
    write(fC,'("curvature0_vs_t_D",f0.3,".txt")') D
    call sanitize_filename(fP); call sanitize_filename(fC)

    open(newunit=uP, file=trim(fP), status="replace", action="write")
    write(uP,*) "# t    P(0,t)     (MC, D=",D,")"
    do k=1, ntimes
      write(uP,'(2(1X,ES20.12))') tgrid(k), P0(k)
    end do
    close(uP)

    open(newunit=uC, file=trim(fC), status="replace", action="write")
    write(uC,*) "# t    d2P(0,t)    D^(5/2)*d2P(0,t)   (MC, D=",D,")"
    do k=1, ntimes
      write(uC,'(3(1X,ES20.12))') tgrid(k), curv(k), curv_scaled(k)
    end do
    close(uC)

    print *,"Wrote:", trim(fP), "and", trim(fC)

    deallocate(P0, curv, curv_scaled)
  end subroutine write_P0_and_curv

  !-------------------- helpers --------------------!
  subroutine exponential(s, lambda)
    implicit none
    real*8, intent(out) :: s
    real*8, intent(in)  :: lambda
    real*8 :: u
    real*8, external :: grnd
    do
      u = grnd()
      if (u > 0.d0) exit
    end do
    s = -dlog(u) / lambda
  end subroutine exponential

  subroutine gaussian(z)
    real*8, intent(out) :: z
    real*8 :: x1, x2, w
    real*8, external :: grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
      x1 = 1.0d0 - 2.0d0*grnd()
      x2 = 1.0d0 - 2.0d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    z = dsqrt(-2.0d0 * dlog(w)/w) * x1
  end subroutine gaussian

  ! 5-point quadratic fit (least squares) on equispaced x: returns y ≈ a x^2 + b x + c
  subroutine quadfit_5(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5, a,b,c)
    implicit none
    real*8, intent(in)  :: x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    real*8, intent(out) :: a,b,c
    ! Because x are symmetric (-2h,-h,0,h,2h), closed-form normal equations simplify.
    ! We solve for a,b,c by linear least squares; with symmetry, b≈0, but we keep it general.
    real*8 :: h, Sx2, Sx4, Sy, Sxy, Sx2y
    ! Assume equal spacing: let h = x4 - x3
    h = x4 - x3
    ! Sums
    Sx2  = x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5
    Sx4  = x1**4 + x2**4 + x3**4 + x4**4 + x5**4
    Sy   = y1 + y2 + y3 + y4 + y5
    Sxy  = x1*y1 + x2*y2 + x3*y3 + x4*y4 + x5*y5
    Sx2y = (x1*x1)*y1 + (x2*x2)*y2 + (x3*x3)*y3 + (x4*x4)*y4 + (x5*x5)*y5
    ! Normal equations for [a b c] on basis [x^2, x, 1]
    ! M = [[Sx4,  0,   Sx2],
    !      [ 0 , Sx2,  0 ],
    !      [Sx2,  0,    5 ]]
    ! rhs = [Sx2y, Sxy, Sy]
    ! Solve exploiting zeros:
    if (Sx2 /= 0.0d0) then
       b = Sxy / Sx2
    else
       b = 0.0d0
    end if
    ! 2x2 for a,c:
    ! [Sx4  Sx2][a] = [Sx2y]
    ! [Sx2   5 ][c]   [Sy   ]
    call solve2x2(Sx4, Sx2, Sx2, 5.0d0, Sx2y, Sy, a, c)
  end subroutine quadfit_5

  subroutine solve2x2(a11,a12,a21,a22, b1,b2, x1,x2)
    implicit none
    real*8, intent(in)  :: a11,a12,a21,a22, b1,b2
    real*8, intent(out) :: x1,x2
    real*8 :: det
    det = a11*a22 - a12*a21
    if (det == 0.0d0) then
       x1 = 0.0d0; x2 = 0.0d0
    else
      x1 = ( b1*a22 - b2*a12) / det
      x2 = (-b1*a21 + b2*a11) / det
    end if
  end subroutine solve2x2

  subroutine smooth3(y, n)
    implicit none
    integer, intent(in) :: n
    real*8, intent(inout) :: y(n)
    real*8 :: tmp(n)
    integer :: i
    tmp(1) = y(1)
    do i=2, n-1
      tmp(i) = (y(i-1) + y(i) + y(i+1))/3.0d0
    end do
    tmp(n) = y(n)
    y = tmp
  end subroutine smooth3

  subroutine sanitize_filename(str)
    character(len=*), intent(inout) :: str
    integer :: i, L
    L = len_trim(str)
    do i=1, L
      if (str(i:i)=='.') str(i:i)='p'
      if (str(i:i)==' ') str(i:i)='_'
    end do
    if (L >= 4) then
      if (str(L-3:L)==".txt") str(L-3:L)="ptxt"
    end if
  end subroutine sanitize_filename

end program rtp_fig2_mc

! RNG (Mersenne Twister) must provide double-precision function grnd()
include 'mt.f90'
