!=====================================================================
! Inertial jerky particle + AOUP (active) -- Fig. 4(a)
!
! Mechanics:  λ x'''(t) + m x''(t) = γ_act * u(t)   (γ_fric=0, k=0)
! Define:     τI = λ / m
!
! Noise (AOUP): u_{n+1} = e^{-dt/τp} u_n + v0*sqrt(1 - e^{-2dt/τp})*N(0,1)
!
! Outputs (normalized to match paper axes):
!   1) msd_active_taupratio_1e-3.txt : τp = τI/1000  →  6→5→3
!   2) msd_active_taupratio_1e3.txt  : τp = 1000 τI  →  6→4→3
!
! Axis normalization for Fig. 4(a):
!   x: t/τI
!   y: MSD / (γ^2 v0^2 τI^6 / λ^2)
!
! Notes:
! - Uses trapezoid u_avg=(u_n+u_{n+1})/2 and exact in-step solution for a,v,x.
! - Choose dt as a small fraction of min(τI, τp) for speed+accuracy.
! 
!=====================================================================
program inertial_jerky_active_fig4a
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)

  ! ---------- physical parameters ----------
  real(dp), parameter :: lambda    = 1.0_dp     ! jerk coefficient λ
  real(dp), parameter :: gamma_act = 1.0_dp     ! RHS amplitude γ
  real(dp), parameter :: v0        = 1.0_dp     ! OU stationary std

  ! ---------- simulation controls ----------
  integer,  parameter :: Ntraj = 2000           ! ensemble size (raise later if needed)
  real(dp), parameter :: dt_frac = 2.0d-2       ! dt = dt_frac * min(τI, τp)

  integer :: seed
  seed = 24681357
  call sgrnd(seed)

  ! ----- Case A (red): τp = τI / 1000  → 6→5→3 -----
  call run_one_case("msd_active_taupratio_1e-3.txt", &
                    tauI=1.0_dp, taup=1.0e-3_dp, Tmax=120.0_dp)

  ! ----- Case B (blue): τp = 1000 τI → 6→4→3 -----
  ! Need Tmax ~ τp to see 4→3 crossover; here τI=1, τp=1000
  call run_one_case("msd_active_taupratio_1e3.txt", &
                    tauI=1.0_dp, taup=1000.0_dp, Tmax=5000.0_dp)

contains

  subroutine run_one_case(outfile, tauI, taup, Tmax)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: tauI, taup, Tmax
    real(dp) :: dt, tnow, normY
    integer  :: Nsteps, it, n, uout
    real(dp), allocatable :: msd(:)

    ! timestep chosen relative to min(τI, τp)
    dt = dt_frac * min(tauI, taup)
    Nsteps = nint(Tmax / dt)

    allocate(msd(0:Nsteps)); msd = 0.0_dp

    do n = 1, Ntraj
      call one_trajectory(msd, Nsteps, dt, tauI, taup)
      if (mod(n, 50) == 0) print *, "Completed", n, "of", Ntraj, "trajectories"
    end do

    msd = msd / real(Ntraj, dp)

    ! Fig. 4(a) normalization
    normY = (gamma_act**2 * v0**2 * tauI**6) / (lambda**2)

    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout,*) "# t/τI    MSD/(γ^2 v0^2 τI^6 / λ^2)"
    do it = 0, Nsteps
      tnow = it * dt
      write(uout,'(2(1X,ES20.12))') tnow/tauI, msd(it)/normY
    end do
    close(uout)

    deallocate(msd)
  end subroutine run_one_case

  subroutine one_trajectory(msd, Nsteps, dt, tauI, taup)
    implicit none
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Nsteps
    real(dp), intent(in)    :: dt, tauI, taup

    real(dp) :: x, v, a, u, u_new, u_avg
    real(dp) :: v_old, a_old, B, r, em
    real(dp) :: aou_e, aou_sig, g
    integer  :: it

    ! OU parameters
    aou_e  = exp(-dt/taup)
    aou_sig = v0 * sqrt(max(0.0_dp, 1.0_dp - aou_e*aou_e))

    ! initial conditions
    call gaussian(g); u = v0 * g   ! stationary OU
    x = 0.0_dp; v = 0.0_dp; a = 0.0_dp
    msd(0) = msd(0) + 0.0_dp

    r  = dt / tauI
    em = exp(-r)

    do it = 1, Nsteps
      ! OU step
      call gaussian(g)
      u_new = aou_e * u + aou_sig * g
      u_avg = 0.5_dp * (u + u_new)

      ! cache old state
      v_old = v
      a_old = a

      ! effective steady acceleration
      B = (gamma_act / lambda) * u_avg * tauI

      ! update acceleration
      a = a_old*em + B*(1.0_dp - em)

      ! update velocity
      v = v_old + a_old*tauI*(1.0_dp - em) + B*(dt - tauI*(1.0_dp - em))

      ! update position (use v_old)
      x = x + v_old*dt                                       &
            + a_old*(tauI*dt - tauI*tauI*(1.0_dp - em))      &
            + B*(0.5_dp*dt*dt - tauI*dt + tauI*tauI*(1.0_dp - em))

      u = u_new
      msd(it) = msd(it) + x*x
    end do
  end subroutine one_trajectory

  !---------------- RNG & Gaussian (uses mt.f90's grnd()) ----------------
  subroutine gaussian(z)
    implicit none
    real(dp), intent(out) :: z
    real(dp) :: x1, x2, w, grnd
    w = 2.0_dp
    do while (w > 1.0_dp .or. w == 0.0_dp)
      x1 = 1.0_dp - 2.0_dp*grnd()
      x2 = 1.0_dp - 2.0_dp*grnd()
      w  = x1*x1 + x2*x2
    end do
    z = sqrt(-2.0_dp * log(w) / w) * x1
  end subroutine gaussian

end program inertial_jerky_active_fig4a
include 'mt.f90'

