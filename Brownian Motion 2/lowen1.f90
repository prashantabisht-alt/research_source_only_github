!===============================================================

! Goal: Generate u(t) from AOUP SDE and verify Eqs. (6)–(8)
!variance is a measure of the strength of propulsion noise btw
!===============================================================
program aoup_eq678
  implicit none
  
  integer,  parameter :: Nsteps = 2000000      ! total steps (in time seried of u(t) it has 2 million u's in   timesteps 0,t,2dt..kmax*dt))
  real*8,   parameter :: dt     = 0.001d0      ! time step
  real*8,   parameter :: taup   = 0.50d0       ! persistence time τ_p
  real*8,   parameter :: v0     = 1.00d0       ! stationary std of u (variance df stationary u is vo^2)
  real*8,   parameter :: gamma  = 0.8d0        ! friction (for M(t)=γ^2 C(t))
  integer,  parameter :: Kfrac  = 20           ! max lag = N_eff/Kfrac for 20 tmax is 200 times taup
  integer,  parameter :: seed   = 1002         

  
  logical,  parameter :: stationary_init = .true.  ! .true.: u0~N(0,v0^2); .false.: u0=0 + burn-in
  real*8,   parameter :: burnin_tau = 5.0d0        ! burn-in length in units of τ_p (if stationary_init=.false.)

  ! ---------- Derived / bookkeeping ----------
  integer :: ibegin, N_eff, Kmax
  real*8  :: a, b, D_theory

  ! ---------- Arrays ----------
  real*8, allocatable :: u(:)
  real*8, allocatable :: C_est(:), C_theory(:)
  real*8, allocatable :: M_est(:), M_theory(:)

  ! ---------- Work vars ----------
  integer :: n, k
  real*8  :: s, um, var, tmp, tlag, sumC, D_est
  character(len=32) :: mode

  ! ---------- RNG ----------
  call sgrnd(seed)

  ! ---------- Exact OU coefficients ----------
  a = dexp(-dt/taup)
  b = v0 * dsqrt(1.0d0 - a*a)
  D_theory = v0*v0*taup

  ! ---------- Allocate ----------
  allocate(u(Nsteps))

  ! ---------- Initialize u(1) ----------
  if (stationary_init) then
     call gaussian(s)         ! N(0,1)
     u(1) = v0 * s            ! u0 ~ N(0, v0^2) -> instant stationarity
     ibegin = 1               ! include all data in stats
     mode = 'stationary_init'
  else
     u(1) = 0.0d0             ! non-stationary start; discard early segment in stats
     ibegin = max(1, int(burnin_tau * taup / dt))
     mode = 'zero_init+burnin'
  end if

  ! ---------- Generate u(t) trajectory ----------
  do n = 2, Nsteps
     call gaussian(s)
     u(n) = a*u(n-1) + b*s
  end do

  ! ---------- Effective sample size for stats ----------
  N_eff = Nsteps - ibegin + 1
  if (N_eff < 1000) then
     print *, 'Warning: N_eff is small. Increase Nsteps or reduce burn-in.'
  end if

  ! Use Kmax based on stationary slice length (not full record)
  Kmax =  N_eff / Kfrac

  allocate(C_est(0:Kmax), C_theory(0:Kmax))
  allocate(M_est(0:Kmax), M_theory(0:Kmax))

  ! ---------- Mean & variance over stationary slice ----------
  um = 0.0d0
  do n = ibegin, Nsteps
     um = um + u(n)
  end do
  um = um / dble(N_eff)

  var = 0.0d0
  do n = ibegin, Nsteps
     tmp = u(n) - um
     var = var + tmp*tmp
  end do
  var = var / dble(N_eff)

  ! ---------- Autocorrelation C[k] over stationary slice ----------
  do k = 0, Kmax
     C_est(k) = 0.0d0
     ! n+k must remain inside [ibegin, Nsteps]
     do n = ibegin, Nsteps - k
        C_est(k) = C_est(k) + (u(n)-um)*(u(n+k)-um)
     end do
     C_est(k) = C_est(k) / dble(Nsteps - k - ibegin + 1)
     C_theory(k) = v0*v0 * (a**k)      ! Eq. (6): v0^2 * exp(-k*dt/taup)
  end do

  ! ---------- Memory kernel M[k] and integral check ----------
  do k = 0, Kmax
     M_est(k)    = gamma*gamma * C_est(k)        ! Eq. (7)
     M_theory(k) = gamma*gamma * C_theory(k)
  end do

  ! Integral of C(τ): trapezoid over the lags we computed
  sumC = 0.0d0
  do k = 0, Kmax
     if (k == 0 .or. k == Kmax) then
        sumC = sumC + 0.5d0 * C_est(k)
     else
        sumC = sumC + C_est(k)
     end if
  end do
  D_est = dt * sumC                          ! ≈ v0^2 * taup (Eq. area under C)

  ! ---------- Write outputs ----------
  open(unit=11, file='u_stats.txt', status='replace')
  write(11,'(A)')        '# AOUP u(t) statistics'
  write(11,'(A,1X,A)')   '# init_mode',          trim(mode)
  write(11,'(A,1X,I0)')  '# Nsteps',             Nsteps
  write(11,'(A,1X,ES20.12)') '# dt',             dt
  write(11,'(A,1X,ES20.12)') '# taup',           taup
  write(11,'(A,1X,ES20.12)') '# v0',             v0
  write(11,'(A,1X,ES20.12)') '# gamma',          gamma
  write(11,'(A,1X,I0)')  '# ibegin',             ibegin
  write(11,'(A,1X,I0)')  '# N_eff',              N_eff
  write(11,'(A,1X,ES20.12)') '# mean_est',       um
  write(11,'(A,1X,ES20.12)') '# var_est',        var
  write(11,'(A,1X,ES20.12)') '# v0^2_theory',    v0*v0
  write(11,'(A,1X,ES20.12)') '# D_est(∫C dt)',   D_est
  write(11,'(A,1X,ES20.12)') '# D_theory',       D_theory
  close(11)

  open(unit=12, file='u_autocorr.txt', status='replace')
  write(12,*) '# tau   C_est   C_theory'
  do k = 0, Kmax
     tlag = dble(k) * dt
     write(12,'(3(1X,ES20.12))') tlag, C_est(k), C_theory(k)
  end do
  close(12)

  open(unit=13, file='M_kernel.txt', status='replace')
  write(13,*) '# tau   M_est   M_theory'
  do k = 0, Kmax
     tlag = dble(k) * dt
     write(13,'(3(1X,ES20.12))') tlag, M_est(k), M_theory(k)
  end do
  close(13)

  print *, 'Done. Files: u_stats.txt, u_autocorr.txt, M_kernel.txt'
  print *, 'Check: var≈v0^2,  C(τ)≈v0^2*exp(-τ/τp),  ∫C dτ≈v0^2 τp.'

  deallocate(u, C_est, C_theory, M_est, M_theory)
end program aoup_eq678

! ------------------------------------------------------------
! RNG & Gaussian (your style)
include 'mt.f90'

! Box–Muller transform: returns s ~ N(0,1) using uniform grnd()
subroutine gaussian(s)
  implicit none
  real*8 :: x1, x2, w, s, grnd
  w = 2.0d0
  do while (w > 1.0d0 .or. w == 0.0d0)
     x1 = 1.0d0 - 2.0d0*grnd()
     x2 = 1.0d0 - 2.0d0*grnd()
     w  = x1*x1 + x2*x2
  end do
  s = dsqrt(-2.0d0 * dlog(w) / w) * x1
end subroutine gaussian
