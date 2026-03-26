!===========================================================
! 2D Active Brownian Particle — Trajectory dumper
!  - Short-time path (Tshort << 1/DR)  → nearly ballistic
!  - Long-time path  (Tlong  >> 1/DR)  → diffusive-looking
! Files written:
!   traj_short.txt : n  x  y   (every step)
!   traj_long.txt  : n  x  y   (downsampled for size)

!===========================================================
program abp_traj_dump
  implicit none

  ! ---------- Physical parameters ----------
  real*8, parameter :: v0 = 1.0d0          ! self-propulsion speed
  real*8, parameter :: DR = 1.0d-2         ! rotational diffusivity (tau_R = 1/DR = 100)

  ! ---------- Integrator controls ----------
  real*8, parameter :: dt     = 2.0d-2     ! time step
  real*8, parameter :: Tshort = 5.0d0      ! 0.05 * tau_R  (t << 1/DR)
  real*8, parameter :: Tlong  = 1000.0d0   ! 10   * tau_R  (t >> 1/DR)
  integer, parameter :: stride_long = 10   ! write every 10th step for long run

  ! ---------- State & work ----------
  integer :: nS, nL, n, uS, uL
  real*8  :: x, y, phi, phi_new, phi_mid, g
  real*8  :: sq2DRdt, pi

  ! ---------- RNG (Mersenne Twister) ----------
  external :: sgrnd, grnd
  integer  :: seed

  ! ---------- constants ----------
  pi = dacos(-1.0d0)
  sq2DRdt = dsqrt(2.0d0*DR*dt)

  ! ---------- seed RNG ----------
  seed = 20250903
  call sgrnd(seed)

  ! =========================================================
  ! SHORT-TIME TRAJECTORY
  ! =========================================================
  nS  = nint(Tshort/dt)
  x   = 0.0d0
  y   = 0.0d0
  phi = 0.0d0           ! start along +x to show anisotropy

  open(newunit=uS, file="traj_short.txt", status="replace", action="write")
  write(uS,'(I10,1X,2(ES20.12,1X))') 0, x, y

  do n = 1, nS
     call gaussian(g)                             ! g ~ N(0,1)
     phi_new = phi + sq2DRdt * g                  ! Ito step for angle
     ! midpoint angle for smoother x,y increment:
     phi_mid = 0.5d0*(phi + phi_new)
     x = x + v0*dcos(phi_mid)*dt
     y = y + v0*dsin(phi_mid)*dt
     ! wrap to (-pi,pi] (hygiene only)
     if (phi_new >  pi)  phi_new = phi_new - 2d0*pi
     if (phi_new <= -pi) phi_new = phi_new + 2d0*pi
     phi = phi_new
     write(uS,'(I10,1X,2(ES20.12,1X))') n, x, y
  end do
  close(uS)

  ! =========================================================
  ! LONG-TIME TRAJECTORY (downsampled writes)
  ! =========================================================
  nL  = nint(Tlong/dt)
  x   = 0.0d0
  y   = 0.0d0
  phi = 0.0d0

  open(newunit=uL, file="traj_long.txt", status="replace", action="write")
  write(uL,'(I10,1X,2(ES20.12,1X))') 0, x, y

  do n = 1, nL
     call gaussian(g)
     phi_new = phi + sq2DRdt * g
     phi_mid = 0.5d0*(phi + phi_new)
     x = x + v0*dcos(phi_mid)*dt
     y = y + v0*dsin(phi_mid)*dt
     if (phi_new >  pi)  phi_new = phi_new - 2d0*pi
     if (phi_new <= -pi) phi_new = phi_new + 2d0*pi
     phi = phi_new
     if (mod(n, stride_long) == 0) then
        write(uL,'(I10,1X,2(ES20.12,1X))') n, x, y
     end if
  end do
  close(uL)

  print *, "Wrote: traj_short.txt , traj_long.txt"

contains

  ! --------- Standard polar Box–Muller using grnd() ----------
  subroutine gaussian(s)
    real*8 :: s, x1, x2, w
    real*8, external :: grnd
    w = 2d0
    do while (w > 1d0 .or. w == 0d0)
      x1 = 1d0 - 2d0*grnd()
      x2 = 1d0 - 2d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    s = dsqrt(-2d0*dlog(w)/w) * x1
  end subroutine gaussian

end program abp_traj_dump

! ---------- RNG backend (Mersenne Twister) ----------
include 'mt.f90'
