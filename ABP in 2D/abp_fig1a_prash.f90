!===========================================================================
! GOAL: Reproduce Fig. 1(a) for Active Brownian Particle in 2D
!       Show anisotropic MSD at short times (σx^2 ~ t^4, σy^2 ~ t^3)
!       and diffusive tail at late times with D_eff = v0^2/(2 DR).
!
! MODEL (overdamped, constant speed v0, diffusing orientation φ):
!   xdot = v0 * cos(phi)
!   ydot = v0 * sin(phi)
!   phidot = sqrt(2*DR) * eta(t),  <eta>=0, <eta(t)eta(t')>=δ(t-t')
!
! OUTPUTS:
!   1) msd_sim_abp_fig1.txt    : t, σx^2(sim), σy^2(sim), <x>   (for sanity)
!   2) msd_theory_abp_fig1.txt : t, σx^2(theory), σy^2(theory)  (exact, Eq. C4)
!
! NOTE:
!   Keep phi(0)=0 to preserve anisotropy (like the paper’s setup).
!   Use your mt.f90 RNG (sgrnd, grnd) + gaussian(s) for normals.
!===========================================================================

program abp_fig1a_prash
  implicit none

  !-------------------- Physical Parameters --------------------
  real*8, parameter :: v0   = 1.0d0
  real*8, parameter :: DR   = 1.0d-2   ! -> tau_R = 100

  !-------------------- Simulation Controls --------------------
  integer, parameter :: Ntraj = 20000       
  real*8,  parameter :: dt    = 1.0d-2
  real*8,  parameter :: tmax  = 1000.0d0
  integer            :: nsteps

  !-------------------- Accumulators / Work --------------------
  real*8, allocatable :: sumx(:), sumy(:), sumx2(:), sumy2(:), tgrid(:)
  real*8, allocatable :: thx(:), thy(:)
  real*8 :: x, y, phi, z, t, ex, ey, ex2, ey2, sx2, sy2, sq2DRdt
  integer :: i, istep, seed, u1, u2

  !-------------------- RNG Init (mt.f90) ----------------------
  seed = 24681357
  call sgrnd(seed)

  !-------------------- Setup Arrays ---------------------------
  nsteps = nint(tmax/dt)
  allocate(sumx(nsteps), sumy(nsteps), sumx2(nsteps), sumy2(nsteps), tgrid(nsteps))
  allocate(thx(nsteps), thy(nsteps))
  sumx=0d0; sumy=0d0; sumx2=0d0; sumy2=0d0
  sq2DRdt = dsqrt(2.0d0*DR*dt)

  !-------------------- Precompute Theory ----------------------
  do istep = 1, nsteps
     t = dt*dble(istep)
     tgrid(istep) = t
     call theory_sigma_xy(t, v0, DR, thx(istep), thy(istep))
  end do

  !=================================================================
  ! MAIN ENSEMBLE LOOP
  !=================================================================
  do i = 1, Ntraj
     ! ICs (important: phi(0)=0 to break x/y symmetry)
     x = 0d0; y = 0d0; phi = 0d0

     do istep = 1, nsteps
        ! ---- angle step: phi_{n+1} = phi_n + sqrt(2DR dt) * N(0,1)
        call gaussian(z)                  ! your subroutine below uses grnd()
        phi = phi + sq2DRdt * z

        ! optional hygiene: wrap to (-pi,pi]
        if (phi >  3.141592653589793d0) phi = phi - 6.283185307179586d0
        if (phi <=-3.141592653589793d0) phi = phi + 6.283185307179586d0

        ! ---- position step using current phi
        x = x + v0*dcos(phi)*dt
        y = y + v0*dsin(phi)*dt

        ! ---- on-the-fly accumulators
        sumx(istep)  = sumx(istep)  + x
        sumy(istep)  = sumy(istep)  + y
        sumx2(istep) = sumx2(istep) + x*x
        sumy2(istep) = sumy2(istep) + y*y
     end do
  end do

  !-------------------- Write Outputs ---------------------------
  open(newunit=u1, file="msd_sim_abp_fig1.txt",    status="replace")
  open(newunit=u2, file="msd_theory_abp_fig1.txt", status="replace")

  do istep = 1, nsteps
     t  = tgrid(istep)
     ex = sumx(istep)/dble(Ntraj)
     ey = sumy(istep)/dble(Ntraj)
     ex2 = sumx2(istep)/dble(Ntraj)
     ey2 = sumy2(istep)/dble(Ntraj)
     sx2 = ex2 - ex*ex
     sy2 = ey2 - ey*ey

     write(u1,'(4(1X,ES20.12))') t, sx2, sy2, ex
     write(u2,'(3(1X,ES20.12))') t, thx(istep), thy(istep)
  end do

  close(u1); close(u2)

  print *, "Done. Files: msd_sim_abp_fig1.txt, msd_theory_abp_fig1.txt"

contains

  !=================================================================
  ! Exact MSDs (Appendix C, Eq. C4)
  !=================================================================
  subroutine theory_sigma_xy(t, v0, DR, sx2, sy2)
    implicit none
    real*8, intent(in)  :: t, v0, DR
    real*8, intent(out) :: sx2, sy2
    real*8 :: e1, e2, e4
    e1 = dexp(-DR*t)
    e2 = dexp(-2.0d0*DR*t)
    e4 = dexp(-4.0d0*DR*t)
    sx2 = (v0*v0/DR)*t + (v0*v0/(12.0d0*DR*DR))*( e4 - 12.0d0*e2 + 32.0d0*e1 - 21.0d0 )
    
    sy2 = (v0*v0/DR)*t - (v0*v0/(12.0d0*DR*DR))*( e4 - 16.0d0*e1 + 15.0d0 )

  end subroutine theory_sigma_xy

end program abp_fig1a_prash

!---------------------------------------------------------------
! RNG: your Mersenne Twister (sgrnd, grnd) & gaussian below
!---------------------------------------------------------------
include 'mt.f90'

!========================================================
! Unit Gaussian using Polar Box–Muller (your style)
!========================================================
subroutine gaussian(s)
  implicit none
  real*8 :: x1, x2, w, s, grnd
  w = 2.0d0
  do while (w > 1.0d0 .or. w == 0.0d0)
     x1 = 1.0d0 - 2.0d0*grnd()
     x2 = 1.0d0 - 2.0d0*grnd()
     w  = x1*x1 + x2*x2
  end do
  s = dsqrt(-2.0d0 * dlog(w)/w) * x1
end subroutine gaussian
