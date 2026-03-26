program AOUP_EQ6_EQ7
  implicit none

  ! ===== user parameters =====
  integer,          parameter :: Nsamp = 1000000
  integer,          parameter :: Kfrac = 10            ! safety cap for max lag
  real*8,           parameter :: taup  = 1.0d0         ! persistence time τ_p
  real*8,           parameter :: v0    = 1.0d0         ! stationary std of u(t)
  real*8,           parameter :: dt    = 0.001d0       ! timestep
  real*8,           parameter :: Tburn = 5.0d0         ! burn-in time (~5 τ_p)
  real*8,           parameter :: gamma = 0.8d0         ! friction (scales M)
  integer,          parameter :: seed  = 1002          ! RNG seed

  character(*),     parameter :: f_u  = 'u_samples.txt'
  character(*),     parameter :: f_ru = 'eq6_ru.txt'
  character(*),     parameter :: f_M  = 'eq7_M.txt'

  ! ===== externals from mt.f90 and our gaussian =====
  external :: sgrnd, gaussian
  real*8   :: grnd
  external :: grnd

  ! ===== work vars =====
  integer :: n, k, Nburn, Kmax
  real*8  :: phi, sigma, u, t, eta
  real*8  :: mean_u, var_u, lag_t
  real*8  :: area_est, area_th_pos
  real*8, allocatable :: us(:), u0(:), Ru(:), Ru_th(:), M(:), M_th(:)

  ! --- seed and exact OU coefficients ---
  call sgrnd(seed)
  Nburn = int(Tburn/dt); if (Nburn < 0) Nburn = 0
  phi   = dexp(-dt/taup)                             ! memory per step
  sigma = v0 * dsqrt(dmax1(0.0d0, 1.0d0 - phi*phi))  ! innovation scale

  ! choose max lag up to ~20 τ_p, but not exceeding Nsamp/Kfrac
  Kmax = int(20.0d0 * taup / dt)
  if (Kmax <= 0) Kmax = 1
  Kmax = min(Kmax, Nsamp / Kfrac)

  ! --- sample u(t) with exact OU AR(1) ---
  allocate(us(Nsamp))
  call gaussian(eta)             ! init from stationary
  u = v0 * eta

  do n = 1, Nburn
    call gaussian(eta)
    u = phi*u + sigma*eta
  end do

  open(11, file=f_u, status='replace')
  write(11,'(a)') '# t    u'
  t = 0.0d0
  do n = 1, Nsamp
    call gaussian(eta)
    u = phi*u + sigma*eta
    us(n) = u
    write(11,'(f16.8,1x,es20.10)') t, u
    t = t + dt
  end do
  close(11)

  ! --- stats and de-mean for autocorrelation ---
  mean_u = sum(us) / dble(Nsamp)
  var_u  = sum( (us - mean_u)**2 ) / dble(Nsamp)

  allocate(u0(Nsamp))
  u0 = us - mean_u    ! remove DC bias before Ru

  write(*,*) '--- AOUP sample stats (demeaned for R_u) ---'
  write(*,'(a,1x,es12.5)') 'phi   =', phi
  write(*,'(a,1x,es12.5)') 'sigma =', sigma
  write(*,'(a,1x,es12.5)') 'mean(u)=', mean_u
  write(*,'(a,1x,es12.5)') 'var(u) =', var_u, ' (target ~ ', v0*v0, ')'
  write(*,'(a,i0)')        'Kmax   =', Kmax

  ! --- empirical Ru(t) and M(t) vs theory (Eqs. 6 & 7) ---
  allocate(Ru(0:Kmax), Ru_th(0:Kmax), M(0:Kmax), M_th(0:Kmax))
  do k = 0, Kmax
    Ru(k) = 0.0d0
    do n = 1, Nsamp - k
      Ru(k) = Ru(k) + u0(n) * u0(n+k)
    end do
    Ru(k)    = Ru(k) / dble(Nsamp - k)
    Ru_th(k) = v0*v0 * dexp( - (dble(k)*dt)/taup )     ! Eq. (6)
    M(k)     = gamma*gamma * Ru(k)
    M_th(k)  = gamma*gamma * Ru_th(k)                  ! Eq. (7)
  end do

  ! area over t≥0 (theory: ∫_0^∞ M dt = γ^2 v0^2 τ_p)
  area_est    = dt * sum(M)
  area_th_pos = gamma*gamma * v0*v0 * taup
  write(*,'(a,1x,es12.5)') '∫_0^∞ M dt (est) ~', area_est
  write(*,'(a,1x,es12.5)') '∫_0^∞ M dt (th ) =', area_th_pos

  ! --- write txt files (space-separated) ---
  open(12, file=f_ru, status='replace')
  write(12,'(a)') '# t    Ru_empirical    Ru_theory'
  do k = 0, Kmax
    lag_t = dble(k) * dt
    write(12,'(f12.6,1x,es20.10,1x,es20.10)') lag_t, Ru(k), Ru_th(k)
  end do
  close(12)

  open(13, file=f_M, status='replace')
  write(13,'(a)') '# t    M_empirical    M_theory'
  do k = 0, Kmax
    lag_t = dble(k) * dt
    write(13,'(f12.6,1x,es20.10,1x,es20.10)') lag_t, M(k), M_th(k)
  end do
  close(13)

  deallocate(us, u0, Ru, Ru_th, M, M_th)

  write(*,*) 'Wrote:'
  write(*,*) '  ', trim(f_u)
  write(*,*) '  ', trim(f_ru)
  write(*,*) '  ', trim(f_M)
  write(*,*) 'Plot lin-lin first; then log-y. Mask negatives on log-y tail.'

end program AOUP_EQ6_EQ7

! keep the MT include OUTSIDE the program (avoids clashes with implicit)
include 'mt.f90'

! === your polar Box–Muller (unchanged) ===
subroutine gaussian(s)
  implicit none
  real*8, intent(out) :: s
  real*8 :: x1, x2, w
  real*8 :: grnd
  external :: grnd

  w = 2.0d0
  do while (w .gt. 1.0d0 .or. w .eq. 0.0d0)
    x1 = 1.0d0 - 2.0d0*grnd()
    x2 = 1.0d0 - 2.0d0*grnd()
    w  = x1*x1 + x2*x2
  end do

  s = dsqrt(-2.0d0 * dlog(w) / w) * x1
end subroutine gaussian
