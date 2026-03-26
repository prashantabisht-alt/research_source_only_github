!=====================================================================
! Inertial jerky particle + AOUP (active) -- Fig. 4(a) slow-noise case
! Case B: τp = 1000 τI → scaling 6 → 4 → 3
!
! Output: msd_active_taupratio_1e3two.txt
!==================================================
program jerky_active_slow
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)

  ! Physical parameters
  real(dp), parameter :: lambda    = 1.0_dp
  real(dp), parameter :: gamma_act = 1.0_dp
  real(dp), parameter :: v0        = 1.0_dp
  real(dp), parameter :: tauI      = 1.0_dp
  real(dp), parameter :: taup      = 1000.0_dp

  ! Simulation controls (lighter but enough)
  integer,  parameter :: Ntraj = 10000
  real(dp), parameter :: dt_frac = 2.0d-2      ! dt = 0.02
  real(dp), parameter :: Tmax    = 50000.0_dp   ! ~2 τp, enough to see drift 4→3

  integer :: seed
  seed = 987654
  call sgrnd(seed)

  call run_case("msd_active_taupratio_1e3two.txt", tauI, taup, Tmax, Ntraj, dt_frac)

contains

  subroutine run_case(outfile, tauI, taup, Tmax, Ntraj_local, dt_frac_local)
    implicit none
    character(len=*), intent(in) :: outfile
    real(dp), intent(in) :: tauI, taup, Tmax, dt_frac_local
    integer,  intent(in) :: Ntraj_local
    real(dp) :: dt, tnow, normY
    integer  :: Nsteps, it, n, uout
    real(dp), allocatable :: msd(:)

    dt = dt_frac_local * min(tauI, taup)
    Nsteps = nint(Tmax / dt)

    allocate(msd(0:Nsteps)); msd = 0.0_dp

    do n = 1, Ntraj_local
      call one_traj(msd, Nsteps, dt, tauI, taup)
      if (mod(n, 100) == 0) print *, "Completed", n, "of", Ntraj_local, "trajectories"
    end do

    msd = msd / real(Ntraj_local, dp)

    normY = (gamma_act**2 * v0**2 * tauI**6) / (lambda**2)

    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout,*) "# t/τI    MSD/(γ^2 v0^2 τI^6 / λ^2)"
    do it = 0, Nsteps
      tnow = it * dt
      write(uout,'(2(1X,ES20.12))') tnow/tauI, msd(it)/normY
    end do
    close(uout)
  end subroutine run_case

  subroutine one_traj(msd, Nsteps, dt, tauI, taup)
    implicit none
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Nsteps
    real(dp), intent(in)    :: dt, tauI, taup

    real(dp) :: x, v, a, u, u_new, u_avg
    real(dp) :: v_old, a_old, B, r, em
    real(dp) :: aou_e, aou_sig, g
    integer  :: it

    aou_e  = exp(-dt/taup)
    aou_sig = v0 * sqrt(max(0.0_dp, 1.0_dp - aou_e*aou_e))

    call gaussian(g); u = v0 * g
    x=0.0_dp; v=0.0_dp; a=0.0_dp
    msd(0) = msd(0) + 0.0_dp

    r  = dt/tauI
    em = exp(-r)

    do it = 1, Nsteps
      call gaussian(g)
      u_new = aou_e*u + aou_sig*g
      u_avg = 0.5_dp*(u+u_new)

      v_old = v
      a_old = a
      B = (gamma_act/lambda) * u_avg * tauI

      a = a_old*em + B*(1.0_dp - em)
      v = v_old + a_old*tauI*(1.0_dp - em) + B*(dt - tauI*(1.0_dp - em))
      x = x + v_old*dt + a_old*(tauI*dt - tauI*tauI*(1.0_dp - em)) &
              + B*(0.5_dp*dt*dt - tauI*dt + tauI*tauI*(1.0_dp - em))

      u = u_new
      msd(it) = msd(it) + x*x
    end do
  end subroutine one_traj

  subroutine gaussian(z)
    real(dp), intent(out) :: z
    real(dp) :: x1, x2, w, grnd
    w = 2.0_dp
    do while (w>1.0_dp .or. w==0.0_dp)
      x1 = 1.0_dp-2.0_dp*grnd()
      x2 = 1.0_dp-2.0_dp*grnd()
      w  = x1*x1 + x2*x2
    end do
    z = sqrt(-2.0_dp*log(w)/w)*x1
  end subroutine gaussian

end program jerky_active_slow
include 'mt.f90'
