!============================================================================
! RTP Fig. 3: Steady state in a finite interval (stick-until-flip walls)
! Dimensionless units: v=1, gamma=1; interval x in [-ell, ell], ell=1
! D ∈ {0.10, 0.50, 1.00}
!
! Key correctness points:
!   - Zero-flux BC via: stick at walls for outward drift, suppress outward
!     Brownian kicks only while at the wall.
!   - Time-average histogram over [Teq, Tfin] to get true steady-state.
!   - Endpoint-safe binning: include x = ±ell into edge bins.
!
! Output files ('.'->'p', '.txt'->'ptxt'):
!   fig3_D0p10ptxt, fig3_D0p50ptxt, fig3_D1p00ptxt
! Each line: x_center   P_MC(x)   P_exact(x)
!============================================================================
program rtp_fig3_reflecting
  implicit none
  ! ---------------- simulation knobs ----------------
  integer,  parameter :: Ntraj  = 400000     ! raise to 1e6 for publication smoothness
  real*8,   parameter :: ell    = 1.0d0
  real*8,   parameter :: dt     = 1.0d-3     ! reduce to 5e-4 if you want tighter wall behavior
  real*8,   parameter :: Teq    = 3.0d0      ! start time for steady averaging
  real*8,   parameter :: Tfin   = 5.0d0      ! end time for steady averaging
  real*8,   parameter :: dsamp  = 0.01d0     ! histogram every dsamp in [Teq, Tfin]
  integer,  parameter :: nbins  = 401        ! ODD → unique center bin at x=0
  ! --------------------------------------------------

  real*8  :: Dlist(3) = (/ 0.10d0, 0.50d0, 1.00d0 /)
  real*8  :: D, xmin, xmax, binw
  integer :: steps_eq, steps_fin, stride, nsamp_per_traj
  integer :: Didx
  integer(kind=8), allocatable :: hist(:)   ! large counts (steady averaging)
  character(len=64) :: fname
  integer :: uout

  ! RNG
  integer :: seed
  real*8, external :: grnd
  seed = 20250821
  call sgrnd(seed)

  ! Geometry / time discretization
  xmin = -ell; xmax =  ell
  if (mod(nbins,2) == 0) then
     print *,"ERROR: nbins must be ODD so that x=0 is a bin center."; stop
  end if
  binw = (xmax - xmin) / dble(nbins)

  steps_eq  = nint(Teq/dt)
  steps_fin = nint(Tfin/dt)
  stride    = max(1, nint(dsamp/dt))
  nsamp_per_traj = ((steps_fin - steps_eq) / stride) + 1

  do Didx = 1, 3
     D = Dlist(Didx)
     allocate(hist(nbins)); hist = 0_8

     call simulate_and_average(D, ell, dt, steps_eq, steps_fin, stride, nbins, xmin, binw, nsamp_per_traj, hist)

     ! ----- write MC (time-averaged) and analytic P_ss on identical x-grid -----
     call make_filename(D, fname)          ! -> fig3_D0p10ptxt etc.
     open(newunit=uout, file=trim(fname), status="replace", action="write")
     write(uout,*) "# x_center  P_MC(x)  P_exact(x)   (ell=1, D=",D,", Teq=",Teq,", Tfin=",Tfin,")"
     call dump_hist_plus_exact(uout, D, ell, nbins, xmin, binw, hist, Ntraj, nsamp_per_traj)
     close(uout)
     print *,"Wrote:", trim(fname)

     deallocate(hist)
  end do

contains

  subroutine simulate_and_average(D, ell, dt, steps_eq, steps_fin, stride, nbins, xmin, binw, nsamp_per_traj, hist)
    implicit none
    real*8, intent(in)        :: D, ell, dt, xmin, binw
    integer, intent(in)       :: steps_eq, steps_fin, stride, nbins, nsamp_per_traj
    integer(kind=8), intent(inout) :: hist(nbins)

    integer :: n, s, idx, samples
    real*8  :: x, sigma, pflip, w

    pflip = 1.0d0 - dexp(-dt)          ! flip prob (rate=1)

    do n = 1, Ntraj
       ! ---- equilibrate up to Teq ----
       x = 0.0d0
       if (grnd() < 0.5d0) then; sigma = -1.0d0; else; sigma = 1.0d0; end if

       do s = 1, steps_eq
          if (grnd() < pflip) sigma = -sigma

          ! ballistic: stick at outward wall, else move inside
          if (x >= ell - 1d-12 .and. sigma > 0.0d0) then
             x = ell
          else if (x <= -ell + 1d-12 .and. sigma < 0.0d0) then
             x = -ell
          else
             x = x + sigma*dt
             if (x > ell)  x = ell
             if (x < -ell) x = -ell
          end if

          ! diffusion: suppress outward kicks ONLY when at the wall
          call gaussian(w) ; w = dsqrt(2.0d0*D*dt) * w
          if (x >= ell - 1d-12 .and. w > 0.0d0)  w = 0.0d0
          if (x <= -ell + 1d-12 .and. w < 0.0d0) w = 0.0d0
          x = x + w
          if (x > ell)  x = ell
          if (x < -ell) x = -ell
       end do

       ! ---- steady averaging from Teq to Tfin ----
       samples = 0
       do s = steps_eq, steps_fin
          ! dynamics
          if (grnd() < pflip) sigma = -sigma

          if (x >= ell - 1d-12 .and. sigma > 0.0d0) then
             x = ell
          else if (x <= -ell + 1d-12 .and. sigma < 0.0d0) then
             x = -ell
          else
             x = x + sigma*dt
             if (x > ell)  x = ell
             if (x < -ell) x = -ell
          end if

          call gaussian(w) ; w = dsqrt(2.0d0*D*dt) * w
          if (x >= ell - 1d-12 .and. w > 0.0d0)  w = 0.0d0
          if (x <= -ell + 1d-12 .and. w < 0.0d0) w = 0.0d0
          x = x + w
          if (x > ell)  x = ell
          if (x < -ell) x = -ell

          ! sample on the stride
          if (mod(s - steps_eq, stride) == 0) then
             idx = int((x - xmin)/binw) + 1
             if (idx < 1)     idx = 1        ! endpoint-safe binning
             if (idx > nbins) idx = nbins
             hist(idx) = hist(idx) + 1_8
             samples = samples + 1
          end if
       end do
       ! samples should equal nsamp_per_traj for every trajectory
    end do
  end subroutine simulate_and_average

  subroutine dump_hist_plus_exact(uout, D, ell, nbins, xmin, binw, hist, Ntraj, nsamp_per_traj)
    implicit none
    integer, intent(in)              :: uout, nbins, Ntraj, nsamp_per_traj
    integer(kind=8), intent(in)      :: hist(nbins)
    real*8, intent(in)               :: D, ell, xmin, binw
    integer :: b
    real*8  :: x, Pmc, Pex, normchk

    normchk = 0.0d0
    do b = 1, nbins
       x   = xmin + (dble(b) - 0.5d0)*binw
       Pmc = dble(hist(b)) / ( dble(Ntraj)*dble(nsamp_per_traj)*binw )
       Pex = Pss_exact(D, ell, x)
       write(uout,'(3(1X,ES20.12))') x, Pmc, Pex
       normchk = normchk + Pmc*binw
    end do
    write(*,'("   MC normalization check (~1): ",F10.6)') normchk
  end subroutine dump_hist_plus_exact

  ! ---------------- analytic steady state (paper Eq. 16) ----------------
  function Pss_exact(D, ell, x) result(P)
    implicit none
    real*8, intent(in) :: D, ell, x
    real*8 :: P, mu, normA, denom
    mu    = dsqrt(2.0d0*D + 1.0d0) / D
    normA = 1.0d0 / ( dtanh(mu*ell)/dsqrt(2.0d0*D+1.0d0) + 2.0d0*ell )
    denom = 2.0d0*D*dcosh(mu*ell) + 1.0d0
    P = normA * dcosh(mu*x) / denom
  end function Pss_exact

  ! ---------------- RNG helpers ----------------
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

  subroutine make_filename(D, fname)
    implicit none
    real*8, intent(in)  :: D
    character(len=*), intent(out) :: fname
    integer :: Dint, Ddec
    ! Build a robust tag: D0p10, D0p50, D1p00
    Dint = int(D)                    ! integer part
    Ddec = nint(100.0d0*(D - dble(Dint)))
    write(fname,'("fig3_D",I1,"p",I2.2,".txt")') Dint, Ddec
    call sanitize_filename(fname)
  end subroutine make_filename

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

end program rtp_fig3_reflecting

! Mersenne Twister RNG: must provide double-precision function grnd()
include 'mt.f90'
