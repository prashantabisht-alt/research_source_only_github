!===========================================================================
! GOAL: Fig. 3(a) — Survival Sx(t; x0) for ABP forward coordinate x.
! MODEL: xdot=v0 cos(phi), ydot=v0 sin(phi), phidot=sqrt(2 DR) * eta(t).
! ICs : x(0)=x0>0 (tiny), y(0)=0, phi(0)=0 (forward).
! SURVIVAL: Sx(t) = P[ x(s) > 0 for all s in [0,t] ].
!
! OUTPUTS (per DR):
!   Sx_DR<val>.txt           : t, Sx(t)
!   Sx_collapsed_DR<val>.txt : u=t*DR, Sx(t)     (for panel (b) collapse)
!
! NOTE:
!   - Uses linear interpolation to get the exact (in-step) first hit time.
!     This removes the step-bias and makes the t^{-1/2} tail match theory.
!   - Everything else (style/structure) left as you had it.
!===========================================================================

program abp_fig3a_Sx
  implicit none

  ! -------- physics --------
  real*8, parameter :: v0 = 1.0d0
  real*8, parameter :: x0 = 0.1d0    ! small positive offset (paper uses x0=0.1)

  ! -------- DR list (pick modest range for runtime; add more if you like) --------
  integer, parameter :: NDR = 3
  real*8,  parameter :: DRlist(NDR) = (/ 5.0d-3, 1.0d-2, 5.0d-2 /)  ! 0.005, 0.01, 0.05

  ! -------- numerics --------
  integer, parameter :: Ntraj = 10000       ! raise to 2e5+ for smoother tails
  real*8,  parameter :: dt    = 5.0d-3      ! small enough for all DR in list
  real*8,  parameter :: tmax  = 1000.0d0    ! bump to 10000 for long/clean tail
  integer :: nsteps

  ! -------- work arrays / vars --------
  integer :: iDR, i, istep, alive, k, cumd
  real*8  :: DR, tauR, t, phi, x, y, z, sq2DRdt
  real*8  :: x_prev, frac, Sx_val
  integer, allocatable :: deaths_per_step(:)
  real*8,  allocatable :: thit(:)            ! first-hit time per trajectory
  integer :: uid1, uid2
  character(len=64) :: fout_raw, fout_coll

  ! RNG
  integer :: seed
  seed = 4242421
  call sgrnd(seed)

  ! ---------------- setup ----------------
  nsteps = nint(tmax/dt)
  allocate(thit(Ntraj))              ! allocate once (reused per DR)

  do iDR = 1, NDR
    DR      = DRlist(iDR)
    tauR    = 1.0d0 / DR
    sq2DRdt = dsqrt(2.0d0*DR*dt)

    allocate(deaths_per_step(nsteps))
    deaths_per_step = 0

    ! ----- ensemble: compute precise first-hit times thit(i) -----
    do i = 1, Ntraj
      ! ICs
      x   = x0
      y   = 0.0d0
      phi = 0.0d0
      alive = 1
      thit(i) = tmax + dt     ! default: no hit within window

      if (mod(i, 1000) == 0) print *, "Completed", i, "of", Ntraj, "trajectories for DR=", DR

      do istep = 1, nsteps
        if (alive == 0) exit

        ! pre-step time and x
        t      = dble(istep-1)*dt
        x_prev = x

        ! angle step (Euler–Maruyama)
        call gaussian(z)
        phi = phi + sq2DRdt * z

        ! optional wrap for trig stability
        if (phi >  3.141592653589793d0) phi = phi - 6.283185307179586d0
        if (phi <=-3.141592653589793d0) phi = phi + 6.283185307179586d0

        ! position step
        x = x + v0*dcos(phi)*dt
        y = y + v0*dsin(phi)*dt

        ! detect first crossing inside this step: x_prev > 0, x <= 0
        if (x_prev > 0.0d0 .and. x <= 0.0d0) then
          frac     = x_prev / (x_prev - x)         ! fraction in (0,1]
          thit(i)  = t + frac*dt                   ! linear-in-step hit time
          alive    = 0
          exit
        end if
      end do
    end do

    ! ----- convert hit times -> survival on the integration grid -----
    ! histogram deaths to the first index >= t_hit
    do i = 1, Ntraj
      if (thit(i) <= tmax) then
        k = max(1, min(nsteps, ceiling(thit(i)/dt)))
        deaths_per_step(k) = deaths_per_step(k) + 1
      end if
    end do

    ! write outputs (raw t, and collapsed u=t*DR)
    write(fout_raw,  '("Sx_DR",F6.4,".txt")') DR
    write(fout_coll, '("Sx_collapsed_DR",F6.4,".txt")') DR
    open(newunit=uid1, file=adjustl(fout_raw),  status="replace")
    open(newunit=uid2, file=adjustl(fout_coll), status="replace")

    cumd = 0
    do istep = 1, nsteps
      t      = dble(istep)*dt
      cumd   = cumd + deaths_per_step(istep)
      Sx_val = (dble(Ntraj - cumd))/dble(Ntraj)
      write(uid1,'(2(1X,ES20.12))') t,    Sx_val
      write(uid2,'(2(1X,ES20.12))') t*DR, Sx_val   ! for panel (b)
    end do
    close(uid1); close(uid2)

    deallocate(deaths_per_step)
    print *, "Done DR=", DR, " -> ", trim(adjustl(fout_raw)), " & ", trim(adjustl(fout_coll))
  end do

  deallocate(thit)
  print *, "All DR cases finished."

contains
  !========================================================
  ! Unit Gaussian via Polar Box–Muller (uses mt.f90: grnd())
  !========================================================
  subroutine gaussian(s)
    implicit none
    real*8 :: s, x1, x2, w, grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
      x1 = 1.0d0 - 2.0d0*grnd()
      x2 = 1.0d0 - 2.0d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0*dlog(w)/w) * x1
  end subroutine gaussian
end program abp_fig3a_Sx

include 'mt.f90'
