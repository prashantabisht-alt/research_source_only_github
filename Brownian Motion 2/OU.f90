!=====================================================================
!  1) Legacy RNG included as standalone program units (no scoping mix)
!=====================================================================
include 'mt.f90'     ! <-- your Mersenne Twister file here (no PROGRAM inside)

!=====================================================================
!  2) Ornstein–Uhlenbeck (OU) — Single long trajectory verification
!     - Exact discrete update (no dt bias)
!     - Stationary start
!     - Outputs:
!         ou_series.txt
!         ou_stats.txt
!         ou_hist_vs_gauss.txt
!         ou_cov_eq6.txt          (unnormalized Eq. 6)
!         ou_acf_eq6_norm.txt     (normalized ACF vs exp(-tau/tau_p))
!=====================================================================
program SimulationOU_WithIncludedMT
  implicit none

  ! --- declare legacy RNG routines provided by mt.f90 ---
  double precision, external :: grnd   ! uniform(0,1)
  external :: sgrnd                    ! seed

  ! -----------------------------
  ! Physical Parameters (OU)
  ! -----------------------------
  real*8, parameter :: tau_p = 1.0d0   ! persistence time  τp
  real*8, parameter :: v0    = 0.5d0   ! stationary std of u,  Var[u]=v0^2

  ! -----------------------------
  ! Simulation Control
  ! -----------------------------
  integer, parameter :: nt       = 300000   ! number of time steps
  real*8,  parameter :: dt       = 1.0d-3   ! time step
  integer, parameter :: nlag     = 400      ! max lag for cov/ACF
  integer, parameter :: num_bins = 200      ! histogram bins (PDF)

  ! -----------------------------
  ! Data Arrays
  ! -----------------------------
  real*8, allocatable :: u(:)         ! OU time series
  real*8, allocatable :: Cemp(:)      ! empirical covariance (0..nlag)
  real*8, allocatable :: ACFn(:)      ! normalized ACF (0..nlag)
  integer,  dimension(num_bins) :: histogram

  ! -----------------------------
  ! Working Variables
  ! -----------------------------
  real*8 :: a_fac, b_fac, g             ! exact OU coeffs & gaussian sample
  real*8 :: PI
  real*8 :: mean_u, var_u
  real*8 :: tval, Cthe
  real*8 :: u_min, u_max, bin_width, u_center
  real*8 :: P_sim, P_theory, norm_factor
  integer :: i, k, seed, bin_index, kept, N_eff

  ! -----------------------------
  ! Initialization
  ! -----------------------------
  seed = 20250812
  call sgrnd(seed)                    ! seed legacy RNG
  PI = dacos(-1.0d0)

  a_fac = dexp(-dt/tau_p)             ! a = exp(-dt/τp)
  b_fac = v0 * dsqrt(max(0.0d0, 1.0d0 - a_fac*a_fac))  ! b = v0*sqrt(1-a^2)

  allocate(u(nt), Cemp(0:nlag), ACFn(0:nlag))
  histogram = 0

  print*, "=== OU Simulation (Exact AR(1) update) with included mt.f90 ==="
  print*, "tau_p =", tau_p, "  v0 =", v0, "  dt =", dt
  print*, "a =", a_fac, "  b =", b_fac

  ! -----------------------------------------
  ! Main Simulation Loop (Single Trajectory)
  ! -----------------------------------------
  call gaussian01(g)                   ! g ~ N(0,1)
  u(1) = v0 * g                        ! stationary start: N(0, v0^2)

  do i = 2, nt
    call gaussian01(g)                 ! η ~ N(0,1)
    u(i) = a_fac * u(i-1) + b_fac * g  ! exact OU update
  end do

  ! -----------------------------------------
  ! Basic Stats: mean and variance
  ! -----------------------------------------
  mean_u = sum(u)/dble(nt)
  var_u  = sum( (u - mean_u)**2 )/dble(nt)

  open(unit=20, file="ou_stats.txt", status="replace")
  write(20,'(A,F10.6)') 'dt     = ', dt
  write(20,'(A,F10.6)') 'tau_p  = ', tau_p
  write(20,'(A,F10.6)') 'v0     = ', v0
  write(20,'(A,ES14.6)') 'mean(u)= ', mean_u
  write(20,'(A,ES14.6)') 'var(u) = ', var_u
  write(20,'(A,ES14.6)') 'theory v0^2 = ', v0*v0
  close(20)

  ! -----------------------------------------
  ! Output 1: Time Series (n, t, u)
  ! -----------------------------------------
  open(unit=10, file="ou_series.txt", status="replace")
  write(10,*) "# n     t(0.00)       u(t)"
  do i = 1, nt
    tval = dble(i-1)*dt
    write(10,'(I8,1X,F10.2,1X,ES20.10)') i, tval, u(i)
  end do
  close(10)

  ! -----------------------------------------
  ! Output 2: Stationary PDF of u  (Histogram vs Gaussian)
  ! -----------------------------------------
  u_min = -3.0d0*v0
  u_max =  3.0d0*v0
  bin_width = (u_max - u_min) / dble(num_bins)
  histogram = 0
  kept = 0

  do i = 1, nt
    if (u(i) > u_min .and. u(i) < u_max) then
      bin_index = floor( (u(i) - u_min)/bin_width ) + 1
      bin_index = max(1, min(num_bins, bin_index))
      histogram(bin_index) = histogram(bin_index) + 1
    end if
    kept = kept + 1
  end do
  N_eff = kept

  open(unit=40, file="ou_hist_vs_gauss.txt", status="replace")
  write(40,*) "# u_center(0.00)      P_sim(u)                 P_theory(u)"
  norm_factor = 1.0d0 / ( v0*dsqrt(2.0d0*PI) )   ! N(0, v0^2)
  do bin_index = 1, num_bins
    u_center = u_min + (dble(bin_index) - 0.5d0)*bin_width
    P_sim    = dble(histogram(bin_index)) / ( dble(N_eff)*bin_width )
    P_theory = norm_factor * dexp( -0.5d0*(u_center/v0)**2 )
    write(40,'(F10.2,1X,ES20.10,1X,ES20.10)') u_center, P_sim, P_theory
  end do
  close(40)

  ! -----------------------------------------
  ! Output 3: Eq. (6) — UNNORMALIZED Covariance
  !           C_emp(k) = < (u_t - mean)*(u_{t+k} - mean) >
  !           Theory   = v0^2 * exp(-tau/tau_p)
  ! -----------------------------------------
  do k = 0, nlag
    Cemp(k) = sum( (u(1:nt-k) - mean_u) * (u(1+k:nt) - mean_u) ) / dble(nt - k)
  end do

  open(unit=30, file="ou_cov_eq6.txt", status="replace")
  write(30,*) "# k     tau(0.00)      C_emp(tau)              v0^2*exp(-tau/tau_p)"
  do k = 0, nlag
    tval = dble(k)*dt
    Cthe = (v0*v0) * dexp( - tval / tau_p )
    write(30,'(I8,1X,F10.2,1X,ES20.10,1X,ES20.10)') k, tval, Cemp(k), Cthe
  end do
  close(30)

  ! -----------------------------------------
  ! Output 4: Normalized ACF  ACFn(k) = Cemp(k) / Cemp(0)
  !           Theory          = exp(-tau/tau_p)
  ! -----------------------------------------
  if (Cemp(0) > 0.0d0) then
    ACFn = Cemp / Cemp(0)
  else
    ACFn = 0.0d0
    print*, "WARNING: zero variance? ACFn set to 0."
  end if

  open(unit=31, file="ou_acf_eq6_norm.txt", status="replace")
  write(31,*) "# k     tau(0.00)      ACF_emp_norm            exp(-tau/tau_p)"
  do k = 0, nlag
    tval = dble(k)*dt
    write(31,'(I8,1X,F10.2,1X,ES20.10,1X,ES20.10)') k, tval, ACFn(k), dexp( - tval / tau_p )
  end do
  close(31)

  print*, "Done. Files:"
  print*, "  ou_series.txt"
  print*, "  ou_stats.txt"
  print*, "  ou_hist_vs_gauss.txt"
  print*, "  ou_cov_eq6.txt"
  print*, "  ou_acf_eq6_norm.txt"

contains
  !===============================================================
  ! gaussian01 : N(0,1) via Polar Box–Muller using legacy grnd()
  !   - Uses grnd() :: double precision uniform(0,1)
  !   - Fully contained: no global state here
  !===============================================================
  subroutine gaussian01(z)
    real*8, intent(out) :: z
    real*8 :: u1, u2, s, r

    s = 2.0d0
    do while (s >= 1.0d0 .or. s == 0.0d0)
      u1 = 2.0d0 * dble(grnd()) - 1.0d0   ! ensure REAL*8
      u2 = 2.0d0 * dble(grnd()) - 1.0d0
      s  = u1*u1 + u2*u2
    end do

    r = dsqrt(-2.0d0 * dlog(s) / s)
    z = u1 * r
  end subroutine gaussian01

end program SimulationOU_WithIncludedMT
