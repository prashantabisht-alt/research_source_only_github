program ThirdOrderLambdaLangevin
    implicit none
    ! -------------------------------
    ! Physical parameters
    real*8, parameter :: lambda = 1.0d0, kB_T = 1.0d0  ! Damping coefficient, thermal energy
    ! Simulation control
    integer, parameter :: N_particles = 1000000  ! Number of trajectories
    integer, parameter :: Tmax = 10000           ! Increased time steps (t=10)
    integer, parameter :: nbins = 500            ! Increased bins
    real*8, parameter :: dt = 0.001d0            ! Reduced time step
    ! Arrays
    real*8, dimension(N_particles) :: positions, velocities, accelerations, jerks
    integer, dimension(nbins) :: histogram
    real*8, allocatable :: msd(:)
    ! Variables
    real*8 :: sigma, PI, variance_pos_theory, random_kick, a_old, j_old, noise_term
    real*8 :: x_min, x_max, bin_width, x_center, P_sim, P_theory, norm_factor
    integer :: i, j, bin_index, timestep, seed
    real*8 :: grnd

    ! Initialize
    sigma = sqrt(2.0d0 * kB_T / lambda)  ! Noise strength from FDT
    PI = acos(-1.0d0)
    seed = 1002
    call sgrnd(seed)
    positions = 0.0d0
    velocities = 0.0d0
    accelerations = 0.0d0
    call gaussian(random_kick)
    jerks = random_kick * sigma  ! Random initial jerk
    allocate(msd(Tmax))
    msd = 0.0d0

    ! Main simulation loop
    do j = 1, N_particles
      if (mod(j, 100000) == 0) print *, "Trajectory:", j
      do timestep = 1, Tmax
        call gaussian(random_kick)
        j_old = (sigma * random_kick) / lambda  ! λ j = ξ(t)
        a_old = accelerations(j)
        accelerations(j) = a_old + j_old * dt
        velocities(j) = velocities(j) + accelerations(j) * dt
        positions(j) = positions(j) + velocities(j) * dt
        msd(timestep) = msd(timestep) + positions(j)**2
      end do
    end do

    ! Normalize MSD
    msd = msd / dble(N_particles)

    ! Histogram of position
    x_min = -50.0d0
    x_max = 50.0d0
    bin_width = (x_max - x_min) / dble(nbins)
    histogram = 0
    do i = 1, N_particles
      if (positions(i) > x_min .and. positions(i) < x_max) then
        bin_index = int((positions(i) - x_min) / bin_width) + 1
        histogram(bin_index) = histogram(bin_index) + 1
      end if
    end do

    ! Output position histogram
    open(unit=43, file="pdf_third_order_lambda.txt")
    variance_pos_theory = 2.0d0 * (kB_T / lambda) * (Tmax * dt)  ! Approximate variance
    norm_factor = 1.0d0 / sqrt(2.0d0 * PI * variance_pos_theory)
    do bin_index = 1, nbins
      x_center = x_min + (dble(bin_index) - 0.5d0) * bin_width
      P_sim = dble(histogram(bin_index)) / (dble(N_particles) * bin_width)
      P_theory = norm_factor * exp(-x_center**2 / (2.0d0 * variance_pos_theory))
      write(43, *) x_center, P_sim, P_theory
    end do
    close(43)

    ! Output MSD
    open(unit=44, file="msd_third_order_lambda.txt")
    do i = 1, Tmax
      write(44, *) i * dt, msd(i)
    end do
    close(44)

    deallocate(msd)
    print *, 'Simulation complete. Data written to pdf_third_order_lambda.txt and msd_third_order_lambda.txt'

  contains
    subroutine gaussian(s)
      implicit none
      real*8, intent(out) :: s
      real*8 :: x1, x2, w, grnd
      w = 2.0d0
      do while (w >= 1.0d0 .or. w == 0.0d0)
        x1 = 2.0d0 * grnd() - 1.0d0
        x2 = 2.0d0 * grnd() - 1.0d0
        w = x1 * x1 + x2 * x2
      end do
      s = dsqrt(-2.0d0 * dlog(w) / w) * x1
    end subroutine gaussian
  end program ThirdOrderLambdaLangevin

  include 'mt.f90'