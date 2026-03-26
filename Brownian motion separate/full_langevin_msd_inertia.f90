program full_langevin_msd_inertia_only
    implicit none
    integer, parameter :: N_particles = 100000  ! Particles
    integer, parameter :: Tmax = 5000           ! Time steps (t=50)
    integer :: i, j, timestep, seed, noise_type
    real*8 :: dt, m, gamma, kBT, D
    real*8 :: v, x, v_old, x_old, friction_term, noise_term
    real*8 :: gaussianvariable, eta, lambda
    real*8, allocatable :: msd_gauss(:), msd_exp(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0        ! Time step
    m = 1.0d0          ! Mass
    gamma = 0.5d0      ! Friction coefficient
    kBT = 1.0d0        ! Thermal energy
    D = kBT / gamma    ! Diffusion coefficient (theoretical)
    lambda = 1.0d0     ! Exponential noise parameter
    seed = 6666        ! Random seed

    ! Allocate arrays
    allocate(msd_gauss(Tmax), msd_exp(Tmax))
    msd_gauss = 0.0d0
    msd_exp = 0.0d0

    call sgrnd(seed)

    ! Simulation loop for free case with both noise types
    do noise_type = 1, 2  ! 1: Gaussian, 2: Exponential
      do j = 1, N_particles
        if (mod(j, 10000) == 0) print *, "Noise Type:", noise_type, "Particle:", j
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
          v = v_old - friction_term + noise_term  ! No potential force
          x = x_old + v * dt
          if (noise_type == 1) msd_gauss(timestep) = msd_gauss(timestep) + x**2
          if (noise_type == 2) msd_exp(timestep) = msd_exp(timestep) + x**2
        end do
      end do
    end do

    ! Normalize MSD
    msd_gauss = msd_gauss / dble(N_particles)
    msd_exp = msd_exp / dble(N_particles)

    ! Write MSD to files
    open(unit=44, file="msd_gauss_inertia_only.txt")
    do i = 1, Tmax
      write(44, *) i * dt, msd_gauss(i)
    end do
    close(44)

    open(unit=45, file="msd_exp_inertia_only.txt")
    do i = 1, Tmax
      write(45, *) i * dt, msd_exp(i)
    end do
    close(45)

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

  end program full_langevin_msd_inertia_only

  include 'mt.f90'