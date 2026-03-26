!=====================================================================
! Pure jerky particle + OU (persistent) drive -- Fig. 3(a), clean version
! Mechanics:  λ x'''(t) = γ_act * u(t)   (m=0, γ_fric=0, k=0)
! Noise:      OU/AOUP:  u_{n+1} = e^{-dt/τp} u_n + v0*sqrt(1-e^{-2dt/τp}) * N(0,1)
!
! Outputs (normalized to match paper axes):
!   1) msd_persistent_norm.txt :  s = t/τp ,  MSD / (γ^2 v0^2 τp^6 / λ^2)
!   2) theory_persistent_norm.txt : same, analytic Eq. (23) for overlay
!
! Small-t accuracy fixes:
!   - Stationary init for OU: u0 ~ N(0, v0^2)
!   - Per-step averaged jerk and exact kinematic update on [t, t+dt]
!=====================================================================
program aoup_fig3a
  implicit none

  !---------------- precision & params ----------------
  integer, parameter :: dp = selected_real_kind(15, 300)

  real(dp), parameter :: lambda    = 1.0_dp   ! jerk coefficient λ
  real(dp), parameter :: gamma_act = 1.0_dp   ! RHS amplitude γ (paper’s symbol)
  real(dp), parameter :: taup      = 1.0_dp   ! OU persistence time τp
  real(dp), parameter :: v0        = 1.0_dp   ! OU stationary std

  real(dp), parameter :: dt   = 0.0005_dp      ! time step (units of τp=1 here)
  real(dp), parameter :: Tmax = 20.0_dp       ! total time (~20 τp covers both regimes)
  integer,  parameter :: Nsteps = nint(Tmax/dt)

  integer,  parameter :: Ntraj  = 50000      ! ensemble size (raise for smoother curves)

  !---------------- arrays & vars ----------------
  real(dp), allocatable :: msd(:), tgrid(:)
  real(dp) :: normY
  integer  :: it, n, uout
  integer  :: seed

  !---------------- RNG init ----------------
  seed = 24681357
  call sgrnd(seed)

  allocate(msd(0:Nsteps)); msd = 0.0_dp
  allocate(tgrid(0:Nsteps))
  do it = 0, Nsteps
    tgrid(it) = it * dt
  end do
  ! y-axis normalization used in the paper
  normY = (gamma_act**2 * v0**2 * taup**6) / (lambda**2)

  !---------------- run ensemble ----------------
  call run_persistent_jerky(msd, Ntraj, Nsteps, dt, lambda, gamma_act, taup, v0)

  ! average
  do it = 0, Nsteps
    msd(it) = msd(it) / real(Ntraj, dp)
  end do

  !---------------- write normalized MSD ----------------
  open(newunit=uout, file="msd_persistent_norm.txt", status="replace", action="write")
  write(uout,*) "# s=t/taup   MSD/(γ^2 v0^2 τp^6 / λ^2)   (persistent OU, simulation)"
  do it = 0, Nsteps
    write(uout,'(2(1X,ES20.12))') tgrid(it)/taup, msd(it)/normY
  end do
  close(uout)

  !---------------- write analytic Eq.(23) overlay ----------------
  call write_theory_eq23("theory_persistent_norm.txt", Nsteps, dt, taup, v0, gamma_act, lambda, normY)

  print *, "Wrote:"
  print *, "  msd_persistent_norm.txt"
  print *, "  theory_persistent_norm.txt"
  print *, "Plot both on log-log axes and overlay t^6 and t^5 guides."

contains

  !-------------------------------------------------------------------
  ! Ensemble simulation: OU drive + pure jerky mechanics with
  ! stationary init and per-step averaged jerk + exact kinematics
  !-------------------------------------------------------------------
  subroutine run_persistent_jerky(msd, Ntraj, Nsteps, dt, lambda, gamma_act, taup, v0)
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Ntraj, Nsteps
    real(dp), intent(in)    :: dt, lambda, gamma_act, taup, v0

    real(dp) :: x, v, a, u, g
    real(dp) :: x0, v0loc, a0, u_old, u_new
    real(dp) :: aou_e, aou_sig, j_avg
    integer  :: n, it

    aou_e  = exp(-dt/taup)
    aou_sig = v0 * sqrt(max(0.0_dp, 1.0_dp - aou_e*aou_e))

    do n = 1, Ntraj
      ! stationary OU init for u; mechanics start at rest
      call gaussian(g)
      u = v0 * g
      x = 0.0_dp;  v = 0.0_dp;  a = 0.0_dp

      msd(0) = msd(0) + x*x

      do it = 1, Nsteps
        ! store pre-step values
        x0 = x; v0loc = v; a0 = a; u_old = u

        ! exact OU step
        call gaussian(g)
        u_new = aou_e * u_old + aou_sig * g

        ! averaged jerk over the step
        j_avg = (gamma_act / lambda) * 0.5_dp * (u_old + u_new)

        ! exact kinematics for constant jerk on [t, t+dt]
        a = a0 + dt * j_avg
        v = v0loc + dt * a0 + 0.5_dp * dt*dt * j_avg
        x = x0 + dt * v0loc + 0.5_dp * dt*dt * a0 + (1.0_dp/6.0_dp) * dt*dt*dt * j_avg

        ! commit OU state
        u = u_new

        msd(it) = msd(it) + x*x
      end do
    end do
  end subroutine run_persistent_jerky

  !-------------------------------------------------------------------
  ! Analytic OU MSD (Eq. 23), normalized like the data file
  !-------------------------------------------------------------------
  subroutine write_theory_eq23(fname, Nsteps, dt, taup, v0, gamma_act, lambda, normY)
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: Nsteps
    real(dp),         intent(in) :: dt, taup, v0, gamma_act, lambda, normY
    integer  :: it, uout
    real(dp) :: t, s, val

    open(newunit=uout, file=fname, status="replace", action="write")
    write(uout,*) "# s=t/taup   MSD_theory_norm (Eq.23)"
    do it = 0, Nsteps
      t = it * dt
      s = t / taup
      ! Eq.(23):  MSD = (γ^2 v0^2 τp^6 / λ^2) * [ s^5/10 - s^4/4 + s^3/3 - 2 + e^{-s}(s^2+2s+2) ]
      val = (gamma_act*gamma_act) * (v0*v0) * (taup**6) / (lambda*lambda) * &
            ( s**5/10.0_dp - s**4/4.0_dp + s**3/3.0_dp - 2.0_dp + exp(-s) * (s*s + 2.0_dp*s + 2.0_dp) )
      write(uout,'(2(1X,ES20.12))') s, val / normY
    end do
    close(uout)
  end subroutine write_theory_eq23

  !---------------- RNG & Gaussian ----------------
  
  subroutine gaussian(z)
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

end program aoup_fig3a
include 'mt.f90'
