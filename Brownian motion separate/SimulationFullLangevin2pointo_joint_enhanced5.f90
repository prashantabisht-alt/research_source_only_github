!==========================================================
! GOAL: Under-damped free Langevin (V=0)
!  - P(v): Maxwell–Boltzmann
!  - P(x,t): Gaussian with exact finite-time variance (cold start)
!  - P(x,v,t): bivariate Gaussian; write theory grid and perpendicular axes.
!==========================================================
program SimulationfullLangevintwopointo
  implicit none

  !-----Physical Parameters-----
  real*8, parameter :: m=1.0d0, gamma=0.5d0, kB_T=1.0d0

  !-----Simulation Control Parameters----
  integer, parameter :: N_particles = 1000000 !its 10^6 (million)
  integer, parameter :: num_bins    = 100  ! did shorten it to 100 from 200 to make more smooth so more particles in one bin more smooth
  real*8, parameter  :: T_max       = 5.0d0
  real*8, parameter  :: dt          = 0.01d0

  !-----Data Storage arrays-----
  real*8, dimension(N_particles) :: velocities
  real*8, dimension(N_particles) :: positions
  integer,  dimension(num_bins)  :: histogram
  integer,  dimension(num_bins, num_bins) :: hist2D

  ! --- Working Variables ---
  real*8 :: D, v_old, friction_term, noise_term, random_kick
  real*8, external :: grnd
  real*8 :: PI
  real*8 :: v_min, v_max, x_min, x_max
  real*8 :: bin_width
  real*8 :: v_center, x_center
  real*8 :: P_sim, P_theory, norm_factor
  real*8 :: variance_pos_theory
  integer :: i, j, seed, bin_index
  integer :: num_steps, i_step
  real*8 :: dx, dv, P_joint
  integer :: ix, iv

  !------initialise------
  seed=6666
  call sgrnd(seed)

  D = kB_T/gamma
  velocities = 0.0d0
  positions  = 0.0d0
  PI = dacos(-1.0d0)

  !--------------------------------------
  ! Main simulation loop (Euler–Maruyama)
  !--------------------------------------
  print*,"Simulating particle velocities and positions..."
  num_steps = nint(T_max/dt)

  do i_step = 1, num_steps
     do i = 1, N_particles
        call gaussian(random_kick)
        v_old = velocities(i)
        friction_term = (gamma / m) * v_old * dt
        noise_term    = gamma * dsqrt(2.0d0 * D * dt) / m * random_kick
        velocities(i) = v_old - friction_term + noise_term
        positions(i)  = positions(i) + velocities(i) * dt
     end do
  end do
  print*, "Simulation finished."

  !--------------------------------------------------------------------------
  ! Velocity histogram + theory (Maxwell–Boltzmann)
  !--------------------------------------------------------------------------
  print*, "Building velocity histogram..."
  v_min = -5.0d0
  v_max =  5.0d0
  bin_width = (v_max - v_min) / dble(num_bins)
  histogram = 0

  do i = 1, N_particles
     if (velocities(i) > v_min .and. velocities(i) < v_max) then
        bin_index = int( (velocities(i) - v_min) / bin_width ) + 1
        if (bin_index < 1) bin_index = 1
        if (bin_index > num_bins) bin_index = num_bins
        histogram(bin_index) = histogram(bin_index) + 1
     end if
  end do

  print*, "Writing velocity output file..."
  open(unit=12, file="velocity_distributions5.txt", status='replace', action='write')
  write(12,*) "# v_center, P_simulated, P_theoretical"
  do j = 1, num_bins
     v_center   = v_min + (dble(j) - 0.5d0) * bin_width
     P_sim      = dble(histogram(j)) / (dble(N_particles) * bin_width)
     norm_factor= dsqrt(m / (2.0d0 * PI * kB_T))
     P_theory   = norm_factor * dexp(-m * v_center**2 / (2.0d0 * kB_T))
     write(12, '(3E20.8)') v_center, P_sim, P_theory
  end do
  close(12)

  !==========================================================================
  ! Position histogram + theory (exact finite-time variance, cold start)
  !==========================================================================
  print*, "Building position histogram..."
  x_min = -50.0d0
  x_max =  50.0d0
  bin_width = (x_max - x_min) / dble(num_bins)
  histogram = 0

  do i = 1, N_particles
     if (positions(i) > x_min .and. positions(i) < x_max) then
        bin_index = int( (positions(i) - x_min) / bin_width ) + 1
        if (bin_index < 1) bin_index = 1
        if (bin_index > num_bins) bin_index = num_bins
        histogram(bin_index) = histogram(bin_index) + 1
     end if
  end do

  print*, "Writing position output file..."
  open(unit=13, file="position_distributions5.txt", status='replace', action='write')
  write(13,*) "# x_center, P_simulated, P_theoretical (exact cold-start variance)"
  variance_pos_theory = 2.0d0*(kB_T/gamma) * ( T_max                               &
                        - 2.0d0*(m/gamma)*(1.0d0 - dexp(-T_max*gamma/m))            &
                        + 0.5d0*(m/gamma)*(1.0d0 - dexp(-2.0d0*T_max*gamma/m)) )

  do j = 1, num_bins
     x_center   = x_min + (dble(j) - 0.5d0) * bin_width
     P_sim      = dble(histogram(j)) / (dble(N_particles) * bin_width)
     norm_factor= 1.0d0 / dsqrt(2.0d0 * PI * variance_pos_theory)
     P_theory   = norm_factor * dexp(-x_center**2 / (2.0d0 * variance_pos_theory))
     write(13, '(3E20.8)') x_center, P_sim, P_theory
  end do
  close(13)

  !==========================================================================
  ! Joint distribution P(x,v) @ T_max: histogram grid
  !==========================================================================
  print*, "Building joint histogram of positions and velocities..."
  dv = (v_max - v_min) / dble(num_bins)
  dx = (x_max - x_min) / dble(num_bins)
  hist2D = 0

  do i = 1, N_particles
     if (positions(i) > x_min .and. positions(i) < x_max .and. &
         velocities(i) > v_min .and. velocities(i) < v_max) then
        ix = int( (positions(i) - x_min) / dx ) + 1
        iv = int( (velocities(i) - v_min) / dv ) + 1
        if (ix < 1) ix = 1
        if (ix > num_bins) ix = num_bins
        if (iv < 1) iv = 1
        if (iv > num_bins) iv = num_bins
        hist2D(ix, iv) = hist2D(ix, iv) + 1
     end if
  end do

  print*, "Writing joint distribution output file..."
  open(unit=14, file="joint_distribution5.txt", status='replace', action='write')
  write(14,*) "# x_center, v_center, P_joint(x,v) @ T_max"
  do ix = 1, num_bins
     x_center = x_min + (dble(ix) - 0.5d0) * dx
     do iv = 1, num_bins
        v_center = v_min + (dble(iv) - 0.5d0) * dv
        P_joint  = dble(hist2D(ix, iv)) / (dble(N_particles) * dx * dv)
        write(14, '(3E20.8)') x_center, v_center, P_joint
     end do
     write(14,*)
  end do
  close(14)

  ! Raw phase-space points (optional)
  open(unit=15, file="phase_space_points5.txt", status='replace', action='write')
  write(15,*) "# x, v (final time)"
  do i = 1, N_particles
     write(15,'(2E20.8)') positions(i), velocities(i)
  end do
  close(15)

  !==== THEORY GRID + DIAGNOSTICS + AXIS LINES (theta) ====
  call write_theory_and_diagnostics( T_max, m, gamma, kB_T, &
       positions, velocities, N_particles, &
       x_min, x_max, v_min, v_max, num_bins, &
       dx, dv )

  print*, "Done. Files:"
  print*, "  velocity_distributions5.txt"
  print*, "  position_distributions5.txt"
  print*, "  joint_distribution5.txt, joint_theory5.txt"
  print*, "  axes_theta_theory5.txt, axes_theta_sim5.txt"
  print*, "  diagnostics_Tmax5.txt"
end program SimulationfullLangevintwopointo

!========================================================
! Mersenne Twister RNG and Box–Muller Gaussian
!========================================================
include 'mt.f90'

! Unit Gaussian via polar Box–Muller: N(0,1)
subroutine gaussian(s)
  implicit none
  real*8 :: x1, x2, w, s
  real*8, external :: grnd
  do
     x1 = 1.0d0 - 2.0d0*grnd()
     x2 = 1.0d0 - 2.0d0*grnd()
     w  = x1*x1 + x2*x2
     if (w <= 1.0d0 .and. w > 1.0d-300) exit
  end do
  s = dsqrt(-2.0d0 * dlog(w)/w) * x1
end subroutine gaussian

!========================================================
! Theory overlay + diagnostics + perpendicular axis lines (using theta)
!  - joint_theory5.txt: theory P(x,v) on same grid as hist
!  - diagnostics_Tmax5.txt: theory vs sample stats
!  - axes_theta_theory5.txt, axes_theta_sim5.txt: two perpendicular lines each
!========================================================
subroutine write_theory_and_diagnostics(tnow, m, gamma, kBT, &
                                        x, v, Np, &
                                        x_min, x_max, v_min, v_max, num_bins, &
                                        dx, dv)
  implicit none
  ! inputs
  integer, intent(in) :: Np, num_bins
  real*8,  intent(in) :: tnow, m, gamma, kBT
  real*8,  intent(in) :: x(Np), v(Np)
  real*8,  intent(in) :: x_min, x_max, v_min, v_max, dx, dv

  ! locals
  integer :: i, ix, iv, u_the, u_inf, u_ax_th, u_ax_si
  real*8 :: PI, tau, a, D
  real*8 :: sigv2_th, sigx2_th, cov_th, detS, inv11, inv12, inv22
  real*8 :: x_c, v_c, Pth, rho_th, theta_th
  real*8 :: xmean, vmean, xvar, vvar, cov
  real*8 :: rho_sim, theta_sim, tiny
  real*8 :: norm_check_the
  real*8 :: L, cth, sth, csi, ssi
  character(len=64) :: fthy, finfo

  PI   = dacos(-1.0d0)
  tau  = m / gamma
  a    = dexp( - tnow / tau )
  D    = kBT / gamma
  tiny = 1.0d-300

  ! --- exact (cold-start) theory moments at t = tnow ---
  sigv2_th = (kBT/m) * (1.0d0 - a*a)
  cov_th   = D * (1.0d0 - a)**2
  sigx2_th = 2.0d0*D * ( tnow - 2.0d0*tau*(1.0d0 - a) + 0.5d0*tau*(1.0d0 - a*a) )

  detS  = max( sigx2_th*sigv2_th - cov_th*cov_th, tiny )
  inv11 =  sigv2_th / detS
  inv12 = -cov_th   / detS
  inv22 =  sigx2_th / detS
  rho_th   = cov_th / dsqrt(sigx2_th*sigv2_th)
  theta_th = 0.5d0 * datan2( 2.0d0*cov_th, sigx2_th - sigv2_th )  ! radians

  ! --- write theoretical joint PDF on the same (x,v) grid ---
  fthy = "joint_theory5.txt"
  open(newunit=u_the, file=fthy, status='replace', action='write')
  write(u_the,*) "# x_center, v_center, P_theory(x,v) @ t=", tnow
  norm_check_the = 0.0d0

  do ix = 1, num_bins
     x_c = x_min + (dble(ix) - 0.5d0)*dx
     do iv = 1, num_bins
        v_c = v_min + (dble(iv) - 0.5d0)*dv
        Pth = dexp( -0.5d0*( inv11*x_c*x_c + 2.0d0*inv12*x_c*v_c + inv22*v_c*v_c ) ) &
              / ( 2.0d0*PI*dsqrt(detS) )
        write(u_the,'(3E20.8)') x_c, v_c, Pth
        norm_check_the = norm_check_the + Pth*dx*dv
     end do
     write(u_the,*)
  end do
  close(u_the)

  ! --- sample stats from particles (means, vars, cov) ---
  xmean = 0.0d0; vmean = 0.0d0
  do i = 1, Np
     xmean = xmean + x(i)
     vmean = vmean + v(i)
  end do
  xmean = xmean / dble(Np)
  vmean = vmean / dble(Np)

  xvar = 0.0d0; vvar = 0.0d0; cov = 0.0d0
  do i = 1, Np
     xvar = xvar + (x(i) - xmean)*(x(i) - xmean)
     vvar = vvar + (v(i) - vmean)*(v(i) - vmean)
     cov  = cov  + (x(i) - xmean)*(v(i) - vmean)
  end do
  xvar = xvar / dble(Np)
  vvar = vvar / dble(Np)
  cov  = cov  / dble(Np)

  rho_sim   = cov / max( dsqrt(xvar*vvar), tiny )
  theta_sim = 0.5d0 * datan2( 2.0d0*cov, xvar - vvar )

  ! --- print + save diagnostics ---
  finfo = "diagnostics_Tmax5.txt"
  open(newunit=u_inf, file=finfo, status='replace', action='write')
  write(u_inf,'(A,F10.4)') "t = ", tnow
  write(u_inf,*) "---- THEORY (cold start) ----"
  write(u_inf,'(A,1PE16.8)') "sigma_v^2(th) = ", sigv2_th
  write(u_inf,'(A,1PE16.8)') "sigma_x^2(th) = ", sigx2_th
  write(u_inf,'(A,1PE16.8)') "Cov[x,v](th)  = ", cov_th
  write(u_inf,'(A,1PE16.8)') "rho(th)       = ", rho_th
  write(u_inf,'(A,1PE16.8)') "theta(th)[deg]= ", theta_th*180.0d0/PI

  write(u_inf,*) "---- SIMULATION (sample) ----"
  write(u_inf,'(A,1PE16.8)') "x_mean(sim)   = ", xmean
  write(u_inf,'(A,1PE16.8)') "v_mean(sim)   = ", vmean
  write(u_inf,'(A,1PE16.8)') "sigma_x^2(sim)= ", xvar
  write(u_inf,'(A,1PE16.8)') "sigma_v^2(sim)= ", vvar
  write(u_inf,'(A,1PE16.8)') "Cov[x,v](sim) = ", cov
  write(u_inf,'(A,1PE16.8)') "rho(sim)      = ", rho_sim
  write(u_inf,'(A,1PE16.8)') "theta(sim)[deg]= ", theta_sim*180.0d0/PI

  write(u_inf,*) "---- NORMALIZATION CHECK ----"
  write(u_inf,'(A,1PE16.8)') "∫ P_theory dx dv  ≈ ", norm_check_the
  close(u_inf)

  print *, "Diagnostics @ t=", tnow
  print *, "  theory:  sigv2=", sigv2_th, "  sigx2=", sigx2_th, "  Cov=", cov_th
  print *, "           rho=", rho_th, "  theta[deg]=", theta_th*180.0d0/PI
  print *, "  sample:  sigv2=", vvar,     "  sigx2=", xvar,     "  Cov=", cov
  print *, "           rho=", rho_sim, "  theta[deg]=", theta_sim*180.0d0/PI
  print *, "  ∫P_theory ≈ ", norm_check_the

  ! --- two perpendicular lines (through origin) using just theta ---
  L   = 0.5d0 * dsqrt( (x_max - x_min)**2 + (v_max - v_min)**2 )

  ! THEORY axes
  cth = dcos(theta_th);  sth = dsin(theta_th)
  open(newunit=u_ax_th, file="axes_theta_theory5.txt", status='replace', action='write')
  write(u_ax_th,*) "# two perpendicular lines at theta_th; z=0 for splot"
  ! along theta
  write(u_ax_th,'(3E20.8)') -L*cth, -L*sth, 0.0d0
  write(u_ax_th,'(3E20.8)')  L*cth,  L*sth, 0.0d0
  write(u_ax_th,*)
  ! perpendicular
  write(u_ax_th,'(3E20.8)') -L*sth,  L*cth, 0.0d0
  write(u_ax_th,'(3E20.8)')  L*sth, -L*cth, 0.0d0
  close(u_ax_th)

  ! SIM axes
  csi = dcos(theta_sim); ssi = dsin(theta_sim)
  open(newunit=u_ax_si, file="axes_theta_sim5.txt", status='replace', action='write')
  write(u_ax_si,*) "# two perpendicular lines at theta_sim; z=0 for splot"
  write(u_ax_si,'(3E20.8)') -L*csi, -L*ssi, 0.0d0
  write(u_ax_si,'(3E20.8)')  L*csi,  L*ssi, 0.0d0
  write(u_ax_si,*)
  write(u_ax_si,'(3E20.8)') -L*ssi,  L*csi, 0.0d0
  write(u_ax_si,'(3E20.8)')  L*ssi, -L*csi, 0.0d0
  close(u_ax_si)

  print *, "Wrote joint_theory5.txt, diagnostics_Tmax5.txt,"
  print *, "      axes_theta_theory5.txt, axes_theta_sim5.txt"
end subroutine write_theory_and_diagnostics