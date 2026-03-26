!===============================================================
! Goal: Generate u(t) from AOUP SDE and verify Eqs. (6)–(8)
! - variance is a measure of the strength of propulsion noise btw
!===============================================================
program aoup_eq678
  implicit none
  
  ! ---------- Simulation parameters ----------
  integer,  parameter :: Nsteps = 2000000      ! total steps (like 0,dt,2dt..kmax*dt))
  real*8,   parameter :: dt     = 0.001d0      ! time step
  real*8,   parameter :: taup   = 1.0d0       ! persistence time τ_p(fot persistent limit increase taup but it will include kmax hnece more running time )
  real*8,   parameter :: v0     = 1.00d0       ! stationary std of u
  !real*8,   parameter :: taup   = 0.010d0    ! for white noise limit
  real*8,   parameter :: gamma  = 0.8d0        ! friction (for M(t)=γ^2 C(t))
  integer,  parameter :: seed   = 100          ! RNG seed
  !real*8, parameter :: Dfix = 0.5d0
  !real*8, parameter :: v0   = dsqrt(Dfix/taup)
  integer, parameter :: Lcorr = 50             ! Kmax ≈ 12 τp
  
  logical,  parameter :: stationary_init = .true.   ! .true.: u0~N(0,v0^2)
  real*8,   parameter :: burnin_tau = 5.0d0         ! burn-in length in τ_p units
  !integer,  parameter :: Lcorr = 12                 ! correlation window ~12 τp
  integer :: uac, umk, ust
  ! ---------- Derived / bookkeeping ----------
  integer :: ibegin, N_eff, Kmax
  real*8  :: a, b, D_theory, tail

  ! ---------- Arrays ----------
  real*8, allocatable :: u(:)
  real*8, allocatable :: C_est(:), C_theory(:)
  real*8, allocatable :: M_est(:), M_theory(:)

  ! ---------- Work vars ----------
  integer :: n, k, k_tau
  real*8  :: s, um, var, tmp, tlag, sumC, D_est
  real*8  :: tau_test, C_est_test, C_theory_test
  character(len=32) :: mode

  ! ---------- RNG ----------
  call sgrnd(seed)

  ! ---------- Exact OU coefficients ----------
  a = dexp(-dt/taup)
  b = v0 * dsqrt(1.0d0 - a*a)
  D_theory = v0*v0*taup

  ! ---------- Allocate ----------
  allocate(u(Nsteps))

  ! ---------- Initialize u(1) ----
  if (stationary_init) then
     call gaussian(s)
     u(1) = v0 * s            ! u0 ~ N(0,v0^2)  stationary start
     ibegin = 1
     mode = 'stationary_init'
  else
     u(1) = 0.0d0
     ibegin = max(1, int(burnin_tau * taup / dt))
     mode = 'zero_init+burnin'
  end if

  ! ---------- Generate u(t) trajectory ----------
  do n = 2, Nsteps
     call gaussian(s)
     u(n) = a*u(n-1) + b*s
  end do

  ! ---------- Effective sample size ----------
  N_eff = Nsteps - ibegin + 1
  if (N_eff < 1000) print *, 'Warning: N_eff is small.'

  ! ---------- Choose Kmax based on physics (~12 τp) ----------
  Kmax = ceiling( Lcorr * taup / dt ) ! max lag ~12 τp 
  Kmax = min(Kmax, N_eff-1)           ! (minimal safety: keep within usable data)

  allocate(C_est(0:Kmax), C_theory(0:Kmax))
  allocate(M_est(0:Kmax), M_theory(0:Kmax))

  ! ---------- Mean & variance ----------
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

  ! ---------- Autocorrelation C[k] ----------
  do k = 0, Kmax
     C_est(k) = 0.0d0
     do n = ibegin, Nsteps - k
        C_est(k) = C_est(k) + (u(n)-um)*(u(n+k)-um)
     end do
     C_est(k) = C_est(k) / dble(Nsteps - k - ibegin + 1)
     C_theory(k) = v0*v0 * (a**k)
  end do

  ! ---------- Write autocorrelation to file ----------
  
  open(newunit=uac, file="u_autocorr.txt", status="replace", action="write")
  write(uac,*) "# tau   C_est   C_theory"
  do k = 0, Kmax
     tlag = dble(k) * dt
     write(uac,'(3(1X,ES20.12))') tlag, C_est(k), C_theory(k)
  end do
  close(uac)

  ! ---------- Memory kernel M[k] ----
  do k = 0, Kmax
     M_est(k)    = gamma*gamma * C_est(k)
     M_theory(k) = gamma*gamma * C_theory(k)
  end do

  ! ---------- Write memory kernel to file ----------
  
  open(newunit=umk, file="M_kernel.txt", status="replace", action="write")
  write(umk,*) "# tau   M_est   M_theory"
  do k = 0, Kmax
     tlag = dble(k) * dt
     write(umk,'(3(1X,ES20.12))') tlag, M_est(k), M_theory(k)
  end do
  close(umk)

  ! ---------- Integral of C(τ) ----------
  sumC = 0.0d0
  do k = 0, Kmax
     if (k == 0 .or. k == Kmax) then
        sumC = sumC + 0.5d0 * C_est(k)
     else
        sumC = sumC + C_est(k)
     end if
  end do
  D_est = dt * sumC

  ! ---------- Analytic tail correction ----------
  tail = v0*v0 * taup * dexp( - dble(Kmax)*dt / taup )
  D_est = D_est + tail

  ! ---------- Quick checks in terminal ----------
  print *, '-------------------------------------------'
  print *, 'Quick checks (theory vs estimates):'
  print *, 'Variance:  var_est=', var, '   v0^2=', v0*v0

  k_tau = nint(taup/dt)        ! nearest lag to τp
  if (k_tau > Kmax) k_tau = Kmax    ! safety if Kmax is small
  tau_test = k_tau*dt
  C_est_test = C_est(k_tau)
  C_theory_test = v0*v0*exp(-tau_test/taup)
  print *, 'Correlation at τ≈τp:'
  print *, '   C_est=', C_est_test, '   C_theory=', C_theory_test

  print *, 'Integral ∫C dτ:  D_est=', D_est, '   D_theory=', D_theory
  print *, 'Relative error in ∫C dτ =', abs(D_est-D_theory)/D_theory
  print *, '-------------------------------------------'

  ! ---------- Write run stats to file ----------
  
  open(newunit=ust, file="u_stats.txt", status="replace", action="write")
  write(ust,'(A)') '# AOUP run stats'
  write(ust,'(A,1X,A)')   '# init_mode',          trim(mode)
  write(ust,'(A,1X,I0)')  '# Nsteps',             Nsteps
  write(ust,'(A,1X,I0)')  '# Kmax',               Kmax
  write(ust,'(A,1X,ES20.12)') '# dt',             dt
  write(ust,'(A,1X,ES20.12)') '# taup',           taup
  write(ust,'(A,1X,ES20.12)') '# v0',             v0
  write(ust,'(A,1X,ES20.12)') '# gamma',          gamma
  write(ust,'(A,1X,ES20.12)') '# var_est',        var
  write(ust,'(A,1X,ES20.12)') '# D_est(∫C dt)',   D_est
  write(ust,'(A,1X,ES20.12)') '# D_theory',       D_theory
  close(ust)

  ! ---------- Optional transient variance ----------
  if (.not. stationary_init) then
     open(unit=14, file='var_transient.txt', status='replace')
     write(14,*) '# time   var_est   theory_var'
     var = 0.0d0
     do n = 1, Nsteps
        var = var + u(n)**2
        if (mod(n,1000) == 0) then
           write(14,'(3(1X,ES20.12))') n*dt, var/n, v0*v0*(1.0d0 - exp(-2.0d0*n*dt/taup))
        end if
     end do
     close(14)
     print *, 'Transient variance written to var_transient.txt'
  end if

  deallocate(u, C_est, C_theory, M_est, M_theory)
end program aoup_eq678

! ------------------------------------------------------------
! RNG & Gaussian
include 'mt.f90'

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


!btw check aoup_eq678_6 and aoup_eq678_7 for gnuplot scripts!!!!!!!
