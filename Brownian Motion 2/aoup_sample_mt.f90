program AOUP_sample
  implicit none

  ! -------- user params --------
  integer,          parameter :: Nsamp = 1000000   ! number of recorded samples
  double precision, parameter :: taup  = 1.0d0     ! persistence time τ_p
  double precision, parameter :: v0    = 1.0d0     ! stationary std of u(t)
  double precision, parameter :: dt    = 0.001d0   ! timestep  (≲ τ_p/50)
  double precision, parameter :: Tburn = 5.0d0     ! burn-in (~ 5 τ_p)
  integer,          parameter :: seed  = 1002      ! RNG seed
  character(*),     parameter :: outfile = 'u_samples.txt'
  ! --------------------------------

  ! declare external RNG entry points provided by mt.f90
  external :: sgrnd
  double precision :: grnd
  external :: grnd

  ! work vars
  integer :: n, Nburn, unit
  double precision :: phi, sigma, u, t, eta, mean_u, var_u, sum_u, sum_u2

  ! seed MT
  call sgrnd(seed)

  ! exact OU coefficients
  Nburn = int(Tburn/dt); if (Nburn < 0) Nburn = 0
  phi   = dexp(-dt/taup)
  sigma = v0 * dsqrt(dmax1(0.0d0, 1.0d0 - phi*phi))

  ! init from stationary: u0 ~ N(0, v0^2) to shorten burn-in
  call gaussian(eta)     ! eta ~ N(0,1) using grnd()
  u = v0 * eta

  ! burn-in
  do n = 1, Nburn
    call gaussian(eta)
    u = phi*u + sigma*eta
  end do

  ! record samples to file
  open(newunit=unit, file=outfile, status='replace', action='write')
  write(unit,'(a)') '# time    u'
  t = 0.0d0
  sum_u  = 0.0d0
  sum_u2 = 0.0d0

  do n = 1, Nsamp
    call gaussian(eta)
    u = phi*u + sigma*eta
    write(unit,'(f16.8,1x,es20.10)') t, u
    t = t + dt
    sum_u  = sum_u  + u
    sum_u2 = sum_u2 + u*u
  end do
  close(unit)

  mean_u = sum_u / dble(Nsamp)
  var_u  = (sum_u2 / dble(Nsamp)) - mean_u*mean_u

  write(*,*) '--- AOUP sampling (exact OU AR(1)) ---'
  write(*,'(a,1x,es12.5)') 'phi   =', phi
  write(*,'(a,1x,es12.5)') 'sigma =', sigma
  write(*,'(a,1x,es12.5)') 'mean(u)=', mean_u
  write(*,'(a,1x,es12.5)') 'var(u) =', var_u, ' (target ~ ', v0*v0, ')'
  write(*,*) 'Wrote ', Nsamp, ' samples to ', trim(outfile)

contains
  ! Box–Muller using your mt.f90 uniform grnd()
  subroutine gaussian(z)
    implicit none
    double precision, intent(out) :: z
    double precision :: u1, u2, r
    ! grnd() is visible from the host (declared above)

    u1 = grnd(); if (u1 <= 1.0d-300) u1 = 1.0d-300
    u2 = grnd()
    r  = dsqrt(-2.0d0*dlog(u1))
    z  = r * dcos(6.2831853071795864769d0 * u2)   ! 2π = 6.28...
  end subroutine gaussian
end program AOUP_sample

! <<< IMPORTANT: include mt.f90 OUTSIDE the program >>>
include 'mt.f90'
! btw check aoup_eq678_6 and aoup_eq678_7 for gnuplot scripts!!!!!!!
