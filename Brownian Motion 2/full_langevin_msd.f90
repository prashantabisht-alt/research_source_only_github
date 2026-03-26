program full_langevin_msd
    implicit none
    integer, parameter :: N_particles = 100000  ! Particles
    integer, parameter :: Tmax = 5000           ! Time steps (t=50)
    integer :: i, j, timestep, seed, pot_type, noise_type
    real*8 :: dt, m, gamma, kBT, k, a, D
    real*8 :: v, x, v_old, x_old, friction_term, potential_force_term, noise_term
    real*8 :: gaussianvariable, eta, lambda
    real*8, allocatable :: msd_free_gauss(:), msd_free_exp(:)
    real*8, allocatable :: msd_harm_gauss(:), msd_harm_exp(:)
    real*8, allocatable :: msd_quart_gauss(:), msd_quart_exp(:)
    integer, parameter :: nbins = 400          ! Histogram bins (optional)
    real*8 :: xmin, xmax, bin_width, vmin, vmax
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0        ! Time step
    m = 1.0d0          ! Mass
    gamma = 0.5d0      ! Friction coefficient
    kBT = 1.0d0        ! Thermal energy
    D = kBT / gamma    ! Diffusion coefficient
    k = 0.2d0          ! Harmonic spring constant
    a = 0.1d0          ! Quartic coefficient
    lambda = 1.0d0     ! Exponential noise parameter
    seed = 6666        ! Random seed
    xmin = -15.0d0     ! Position range (optional)
    xmax = 15.0d0
    vmin = -5.0d0      ! Velocity range (optional)
    vmax = 5.0d0
    bin_width = (xmax - xmin) / dble(nbins)

    ! Allocate arrays
    allocate(msd_free_gauss(Tmax), msd_free_exp(Tmax))
    allocate(msd_harm_gauss(Tmax), msd_harm_exp(Tmax))
    allocate(msd_quart_gauss(Tmax), msd_quart_exp(Tmax))
    msd_free_gauss = 0.0d0
    msd_free_exp = 0.0d0
    msd_harm_gauss = 0.0d0
    msd_harm_exp = 0.0d0
    msd_quart_gauss = 0.0d0
    msd_quart_exp = 0.0d0

    call sgrnd(seed)

    ! Simulation loop for different potentials and noise types
    do pot_type = 0, 2  ! 0: Free, 1: Harmonic, 2: Quartic
      do noise_type = 1, 2  ! 1: Gaussian, 2: Exponential
        do j = 1, N_particles
          if (mod(j, 10000) == 0) print *, "Pot Type:", pot_type, "Noise Type:", noise_type, "Particle:", j
          x = 0.0d0
          v = 0.0d0
          do timestep = 1, Tmax
            v_old = v
            x_old = x
            if (noise_type == 1) then
              call gaussian(gaussianvariable)
              noise_term = gamma * sqrt(2.0d0 * D * dt) / m * gaussianvariable
            else
              call exponential(eta, lambda)
              noise_term = gamma * sqrt(2.0d0 * D * dt) / m * (eta / sqrt(2.0d0 / lambda**2))
            end if
            friction_term = (gamma / m) * v_old * dt
            if (pot_type == 1) potential_force_term = (k * x_old / m) * dt  ! Harmonic
            if (pot_type == 2) potential_force_term = (4.0d0 * a * x_old**3 / m) * dt  ! Quartic
            if (pot_type == 0) potential_force_term = 0.0d0  ! Free case
            v = v_old - friction_term - potential_force_term + noise_term
            x = x_old + v * dt
            if (pot_type == 0 .and. noise_type == 1) msd_free_gauss(timestep) = msd_free_gauss(timestep) + x**2
            if (pot_type == 0 .and. noise_type == 2) msd_free_exp(timestep) = msd_free_exp(timestep) + x**2
            if (pot_type == 1 .and. noise_type == 1) msd_harm_gauss(timestep) = msd_harm_gauss(timestep) + x**2
            if (pot_type == 1 .and. noise_type == 2) msd_harm_exp(timestep) = msd_harm_exp(timestep) + x**2
            if (pot_type == 2 .and. noise_type == 1) msd_quart_gauss(timestep) = msd_quart_gauss(timestep) + x**2
            if (pot_type == 2 .and. noise_type == 2) msd_quart_exp(timestep) = msd_quart_exp(timestep) + x**2
          end do
        end do
      end do
    end do

    ! Normalize MSD
    msd_free_gauss = msd_free_gauss / dble(N_particles)
    msd_free_exp = msd_free_exp / dble(N_particles)
    msd_harm_gauss = msd_harm_gauss / dble(N_particles)
    msd_harm_exp = msd_harm_exp / dble(N_particles)
    msd_quart_gauss = msd_quart_gauss / dble(N_particles)
    msd_quart_exp = msd_quart_exp / dble(N_particles)

    ! Write MSD to files
    open(unit=44, file="msd_free_gauss.txt")
    do i = 1, Tmax
      write(44, *) i * dt, msd_free_gauss(i)
    end do
    close(44)

    open(unit=45, file="msd_free_exp.txt")
    do i = 1, Tmax
      write(45, *) i * dt, msd_free_exp(i)
    end do
    close(45)

    open(unit=46, file="msd_harm_gauss.txt")
    do i = 1, Tmax
      write(46, *) i * dt, msd_harm_gauss(i)
    end do
    close(46)

    open(unit=47, file="msd_harm_exp.txt")
    do i = 1, Tmax
      write(47, *) i * dt, msd_harm_exp(i)
    end do
    close(47)

    open(unit=48, file="msd_quart_gauss.txt")
    do i = 1, Tmax
      write(48, *) i * dt, msd_quart_gauss(i)
    end do
    close(48)

    open(unit=49, file="msd_quart_exp.txt")
    do i = 1, Tmax
      write(49, *) i * dt, msd_quart_exp(i)
    end do
    close(49)

    print *, "Simulation complete. Outputs: msd_*.txt"

  contains

    subroutine gaussian(s)
      implicit none
      real*8 :: x1, x2, w, s, grnd
      w = 2.0d0
      do while (w > 1.0d0 .or. w == 0.0d0)
        x1 = 1.0d0 - 2.0d0 * grnd()
        x2 = 1.0d0 - 2.0d0 * grnd()
        w = x1 * x1 + x2 * x2
      end do
      s = dsqrt(-2.0d0 * dlog(w) / w) * x1
    end subroutine gaussian

    subroutine exponential(s, lambda)
      implicit none
      real*8, intent(out) :: s
      real*8, intent(in) :: lambda
      real*8 :: u, grnd
      u = grnd()
      if (u < 0.5d0) then
        s = -dlog(grnd()) / lambda
      else
        s = dlog(grnd()) / lambda
      endif
    end subroutine exponential

  end program full_langevin_msd

  include 'mt.f90'