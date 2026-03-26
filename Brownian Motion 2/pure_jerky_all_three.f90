!=====================================================================
! Pure jerky particle:  λ x'''(t) = γ_act * u(t)
!
! Three noise choices for u(t):
!   (1) ACTIVE WHITE:   <u(t)u(t')> = 2 D_act  δ(t-t')
!   (2) PERSISTENT OU:  u_{n+1} = e^{-dt/τp} u_n + v0*sqrt(1-e^{-2dt/τp})*N(0,1)
!   (3) PASSIVE WHITE:  <u(t)u(t')> = 2 D_pass δ(t-t'),  with D_pass = v0^2 τp
!
! Outputs (normalized like the paper’s Fig. 3a):
!   - msd_active_norm.txt
!   - msd_persistent_norm.txt
!   - msd_passive_norm.txt
!   - theory_persistent_norm.txt   (Eq. 23 overlay)
!
! Per-step kinematics (constant or averaged jerk on [t,t+dt]):
!   a_{n+1} = a_n + dt * j
!   v_{n+1} = v_n + dt*a_n + (1/2) dt^2 * j
!   x_{n+1} = x_n + dt*v_n + (1/2) dt^2 * a_n + (1/6) dt^3 * j
!=====================================================================
program pure_jerky_all_three
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! ---- physical/driver params ----
  real(dp), parameter :: lambda    = 1.0_dp
  real(dp), parameter :: gamma_act = 1.0_dp
  real(dp), parameter :: taup      = 1.0_dp      ! OU persistence time
  real(dp), parameter :: v0        = 1.0_dp      ! OU stationary std

  ! ---- time/discretization ----
  real(dp), parameter :: dt   = 5.0e-4_dp
  real(dp), parameter :: Tmax = 20.0_dp
  integer,  parameter :: Nsteps = nint(Tmax/dt)
  integer,  parameter :: Ntraj  = 50000

  ! ---- diffusion amplitudes for white-noise u(t) ----
  ! Passive (thermal) white via OU→white limit: D_pass = v0^2 * taup
  real(dp), parameter :: D_pass = v0*v0*taup
  ! Active white: free amplitude; choose equal to D_pass for apples-to-apples
  real(dp), parameter :: D_act  = D_pass

  ! ---- storage ----
  real(dp), allocatable :: msd_act(:), msd_pers(:), msd_pass(:), tgrid(:)
  real(dp) :: normY
  integer  :: it, uout

  ! ---- RNG ----
  integer :: seed
  seed = 13579
  call sgrnd(seed)

  allocate(msd_act(0:Nsteps));  msd_act  = 0.0_dp
  allocate(msd_pers(0:Nsteps)); msd_pers = 0.0_dp
  allocate(msd_pass(0:Nsteps)); msd_pass = 0.0_dp
  allocate(tgrid(0:Nsteps))
  do it = 0, Nsteps
    tgrid(it) = it*dt
  end do

  ! y-normalization used in the paper
  normY = (gamma_act**2 * v0**2 * taup**6) / (lambda**2)

  ! ---- run all three ensembles ----
  call run_active_white (msd_act , Ntraj, Nsteps, dt, lambda, gamma_act, D_act)
  call run_persistent_ou(msd_pers, Ntraj, Nsteps, dt, lambda, gamma_act, taup, v0)
  call run_passive_white(msd_pass, Ntraj, Nsteps, dt, lambda, gamma_act, D_pass)

  ! average
  do it = 0, Nsteps
    msd_act(it)  = msd_act(it)  / real(Ntraj,dp)
    msd_pers(it) = msd_pers(it) / real(Ntraj,dp)
    msd_pass(it) = msd_pass(it) / real(Ntraj,dp)
  end do

  ! ---- write normalized outputs ----
  open(newunit=uout, file="msd_active_norm.txt", status="replace", action="write")
  write(uout,*) "# s=t/taup   MSD_active / (γ^2 v0^2 τp^6 / λ^2)"
  do it=0,Nsteps; write(uout,'(2(1X,ES20.12))') tgrid(it)/taup, msd_act(it)/normY; end do
  close(uout)

  open(newunit=uout, file="msd_persistent_norm.txt", status="replace", action="write")
  write(uout,*) "# s=t/taup   MSD_persistent(OU) / (γ^2 v0^2 τp^6 / λ^2)"
  do it=0,Nsteps; write(uout,'(2(1X,ES20.12))') tgrid(it)/taup, msd_pers(it)/normY; end do
  close(uout)

  open(newunit=uout, file="msd_passive_norm.txt", status="replace", action="write")
  write(uout,*) "# s=t/taup   MSD_passive / (γ^2 v0^2 τp^6 / λ^2)"
  do it=0,Nsteps; write(uout,'(2(1X,ES20.12))') tgrid(it)/taup, msd_pass(it)/normY; end do
  close(uout)

  ! theory overlay for persistent (Eq. 23)
  call write_theory_eq23("theory_persistent_norm.txt", Nsteps, dt, taup, v0, gamma_act, lambda, normY)

  print *, "Wrote:"
  print *, "  msd_active_norm.txt"
  print *, "  msd_persistent_norm.txt"
  print *, "  msd_passive_norm.txt"
  print *, "  theory_persistent_norm.txt"
  print *, "Plot on log–log; overlay t^6 (1/36 x^6) and t^5 (0.1 x^5) guides."

contains

  !---------------- ACTIVE WHITE: u ~ N(0, 2D_act/dt), i.i.d. each step ---------------
  subroutine run_active_white(msd, Ntraj, Nsteps, dt, lambda, gamma_act, D)
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Ntraj, Nsteps
    real(dp), intent(in)    :: dt, lambda, gamma_act, D

    real(dp) :: x, v, a, a0, u_mid, j, g, sig
    integer  :: n, it

    sig = sqrt(2.0_dp*D/dt)   ! runtime init; must not be in declaration

    do n=1,Ntraj
      x=0.0_dp; v=0.0_dp; a=0.0_dp
      msd(0)=msd(0)+x*x
      do it=1,Nsteps
        call gaussian(g)
        u_mid = sig * g
        j = (gamma_act/lambda) * u_mid

        a0 = a
        a  = a0 + dt*j
        v  = v  + dt*a0 + 0.5_dp*dt*dt*j
        x  = x  + dt*v  + 0.5_dp*dt*dt*a0 + (1.0_dp/6.0_dp)*dt**3*j

        msd(it) = msd(it) + x*x
      end do
    end do
  end subroutine run_active_white

  !---------------- PERSISTENT (OU/AOUP): exact OU + averaged jerk -------------------
  subroutine run_persistent_ou(msd, Ntraj, Nsteps, dt, lambda, gamma_act, taup, v0)
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Ntraj, Nsteps
    real(dp), intent(in)    :: dt, lambda, gamma_act, taup, v0

    real(dp) :: x, v, a, a0, u_old, u_new, j_avg, g
    real(dp) :: aou_e, aou_sig
    integer  :: n, it

    aou_e  = exp(-dt/taup)
    aou_sig = v0 * sqrt(max(0.0_dp, 1.0_dp - aou_e*aou_e))

    do n=1,Ntraj
      call gaussian(g); u_old = v0*g    ! stationary init for OU
      x=0.0_dp; v=0.0_dp; a=0.0_dp
      msd(0)=msd(0)+x*x
      do it=1,Nsteps
        call gaussian(g)
        u_new = aou_e*u_old + aou_sig*g
        j_avg = (gamma_act/lambda) * 0.5_dp * (u_old + u_new)

        a0 = a
        a  = a0 + dt*j_avg
        v  = v  + dt*a0 + 0.5_dp*dt*dt*j_avg
        x  = x  + dt*v  + 0.5_dp*dt*dt*a0 + (1.0_dp/6.0_dp)*dt**3*j_avg

        u_old = u_new
        msd(it) = msd(it) + x*x
      end do
    end do
  end subroutine run_persistent_ou

  !---------------- PASSIVE WHITE: D_pass fixed by FDT-like mapping -------------------
  subroutine run_passive_white(msd, Ntraj, Nsteps, dt, lambda, gamma_act, D)
    real(dp), intent(inout) :: msd(0:Nsteps)
    integer,  intent(in)    :: Ntraj, Nsteps
    real(dp), intent(in)    :: dt, lambda, gamma_act, D

    real(dp) :: x, v, a, a0, u_mid, j, g, sig
    integer  :: n, it

    sig = sqrt(2.0_dp*D/dt)

    do n=1,Ntraj
      x=0.0_dp; v=0.0_dp; a=0.0_dp
      msd(0)=msd(0)+x*x
      do it=1,Nsteps
        call gaussian(g)
        u_mid = sig * g
        j = (gamma_act/lambda) * u_mid

        a0 = a
        a  = a0 + dt*j
        v  = v  + dt*a0 + 0.5_dp*dt*dt*j
        x  = x  + dt*v  + 0.5_dp*dt*dt*a0 + (1.0_dp/6.0_dp)*dt**3*j

        msd(it) = msd(it) + x*x
      end do
    end do
  end subroutine run_passive_white

  !---------------- Theory (Eq. 23) for persistent OU --------------------------------
  subroutine write_theory_eq23(fname, Nsteps, dt, taup, v0, gamma_act, lambda, normY)
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: Nsteps
    real(dp),         intent(in) :: dt, taup, v0, gamma_act, lambda, normY
    integer  :: it, uout
    real(dp) :: t, s, val
    open(newunit=uout, file=fname, status="replace", action="write")
    write(uout,*) "# s=t/taup   MSD_theory_norm (Eq.23)"
    do it=0,Nsteps
      t = it*dt;  s = t/taup
      val = (gamma_act*gamma_act)*(v0*v0)*(taup**6)/(lambda*lambda) * &
            ( s**5/10.0_dp - s**4/4.0_dp + s**3/3.0_dp - 2.0_dp + exp(-s)*(s*s + 2.0_dp*s + 2.0_dp) )
      write(uout,'(2(1X,ES20.12))') s, val/normY
    end do
    close(uout)
  end subroutine write_theory_eq23

  !---------------- Gaussian via Box–Muller (uses grnd() from mt.f90) ----------------
  subroutine gaussian(z)
    real(dp), intent(out) :: z
    real(dp) :: x1, x2, w
    real(dp), external :: grnd
    w = 2.0_dp
    do while (w > 1.0_dp .or. w == 0.0_dp)
      x1 = 1.0_dp - 2.0_dp*grnd()
      x2 = 1.0_dp - 2.0_dp*grnd()
      w  = x1*x1 + x2*x2
    end do
    z = sqrt(-2.0_dp*log(w)/w) * x1
  end subroutine gaussian

end program pure_jerky_all_three

! RNG implementation (must provide sgrnd, grnd):
include 'mt.f90'

