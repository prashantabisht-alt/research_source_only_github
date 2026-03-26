!===========================================================================
! GOAL: Fig. 3(c,d) — Survival Sy(t; y0) for ABP transverse coordinate y.
! MODEL: xdot=v0 cos(phi), ydot=v0 sin(phi), phidot=sqrt(2 DR) * eta(t).
! ICs : y(0)=y0>0 (tiny), x(0)=0, phi(0)=0 (forward).
! SURVIVAL: Sy(t) = P[ y(s) > 0 for all s in [0,t] ].
!
! OUTPUTS (per DR):
!   Sy_DR<val>.txt                : t, Sy(t)                         (Fig. 3c)
!   Sy_collapsed_DR<val>.txt      : u=t*DR, (v0/(y0*DR))^(1/6)*Sy(t) (Fig. 3d)
!   Sy_RAP_asym_DR<val>.txt       : t, C_RAP*(y0*DR/v0)^(1/6)*(t*DR)^(-1/4) (guide)
!
! NOTE:
!   - Uses linear in-step interpolation for the first crossing y=0.
!   - DR list includes small values to expose a visible RAP (-1/4) regime.
!===========================================================================

program abp_fig3cd_Sy
  implicit none

  ! -------- physics --------
  real*8, parameter :: v0 = 1.0d0
  real*8, parameter :: y0 = 0.1d0   ! small positive offset (paper uses y0=0.1)

  ! -------- DR list (tune as you like) --------
  integer, parameter :: NDR = 3
  real*8,  parameter :: DRlist(NDR) = (/ 5.0d-4, 1.0d-3, 5.0d-3 /)  ! 0.0005, 0.001, 0.005

  ! -------- numerics --------
  integer, parameter :: Ntraj = 15000        ! raise to 1e5–5e5 for smooth tails
  real*8,  parameter :: dt    = 5.0d-3
  real*8,  parameter :: tmax  = 10000.0d0    ! long enough to see -1/4 then -1/2
  integer :: nsteps

  ! -------- work arrays / vars --------
  integer :: iDR, i, istep, alive, k, cumd
  real*8  :: DR, t, phi, x, y, z, sq2DRdt
  real*8  :: y_prev, frac, Sy_val, u, scale_y
  integer, allocatable :: deaths_per_step(:)
  real*8,  allocatable :: thit(:)            ! first-hit time per trajectory
  integer :: uid1, uid2, uid3
  character(len=64) :: f_raw, f_coll, f_rap

  ! RAP asymptotic constant C for Eq. (16) in paper (numerical value)
  ! Sy ~ C * (y0*DR/v0)^(1/6) * (t*DR)^(-1/4) for 1 << t << 1/DR
  real*8, parameter :: C_RAP = 0.8883354788623832d0

  ! RNG
  integer :: seed
  seed = 987654321
  call sgrnd(seed)

  ! ---------------- setup ----------------
  nsteps = nint(tmax/dt)
  allocate(thit(Ntraj))              ! reused per DR

  do iDR = 1, NDR
    DR      = DRlist(iDR)
    sq2DRdt = dsqrt(2.0d0*DR*dt)
    allocate(deaths_per_step(nsteps))
    deaths_per_step = 0

    ! ----- ensemble: compute precise first-hit times thit(i) -----
    do i = 1, Ntraj
      ! ICs
      x   = 0.0d0
      y   = y0
      phi = 0.0d0
      alive = 1
      thit(i) = tmax + dt     ! default: no hit within window

      if (mod(i, 2000) == 0) print *, "Completed", i, "of", Ntraj, "trajectories for DR=", DR

      do istep = 1, nsteps
        if (alive == 0) exit

        ! pre-step time and y
        t      = dble(istep-1)*dt
        y_prev = y

        ! angle step
        call gaussian(z)
        phi = phi + sq2DRdt * z
        if (phi >  3.141592653589793d0) phi = phi - 6.283185307179586d0
        if (phi <=-3.141592653589793d0) phi = phi + 6.283185307179586d0

        ! position step
        x = x + v0*dcos(phi)*dt
        y = y + v0*dsin(phi)*dt

        ! detect first crossing inside this step: y_prev > 0, y <= 0
        if (y_prev > 0.0d0 .and. y <= 0.0d0) then
          frac     = y_prev / (y_prev - y)         ! fraction in (0,1]
          thit(i)  = t + frac*dt                   ! linear-in-step hit time
          alive    = 0
          exit
        end if
      end do
    end do

    ! ----- convert hit times -> survival on the integration grid -----
    do i = 1, Ntraj
      if (thit(i) <= tmax) then
        k = max(1, min(nsteps, ceiling(thit(i)/dt)))
        deaths_per_step(k) = deaths_per_step(k) + 1
      end if
    end do

    ! ----- write outputs -----
    write(f_raw,  '("Sy_DR",F8.4,".txt")') DR
    write(f_coll, '("Sy_collapsed_DR",F8.4,".txt")') DR
    write(f_rap,  '("Sy_RAP_asym_DR",F8.4,".txt")') DR
    open(newunit=uid1, file=adjustl(f_raw),  status="replace")
    open(newunit=uid2, file=adjustl(f_coll), status="replace")
    open(newunit=uid3, file=adjustl(f_rap),  status="replace")

    cumd = 0
    do istep = 1, nsteps
      t      = dble(istep)*dt
      cumd   = cumd + deaths_per_step(istep)
      Sy_val = (dble(Ntraj - cumd))/dble(Ntraj)
      u      = t * DR
      scale_y = (v0/(y0*DR))**(1.0d0/6.0d0)

      ! raw: Sy(t)
      write(uid1,'(2(1X,ES20.12))') t, Sy_val

      ! collapsed: (v0/(y0*DR))^(1/6) * Sy, vs u=t*DR
      write(uid2,'(2(1X,ES20.12))') u, scale_y * Sy_val

      ! RAP asymptote (only meaningful for 1 << t << 1/DR; else write 0)
      if (t > 1.0d0 .and. t < 1.0d0/DR) then
        write(uid3,'(2(1X,ES20.12))') t, C_RAP * ((y0*DR)/v0)**(1.0d0/6.0d0) * u**(-0.25d0)
      else
        write(uid3,'(2(1X,ES20.12))') t, 0.0d0
      end if
    end do

    close(uid1); close(uid2); close(uid3)
    deallocate(deaths_per_step)

    print *, "Done DR=", DR, " -> ", trim(adjustl(f_raw)), ", ", trim(adjustl(f_coll)), ", ", trim(adjustl(f_rap))
  end do

  deallocate(thit)
  print *, "All DR cases finished (Fig. 3 c,d)."

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
end program abp_fig3cd_Sy

include 'mt.f90'
