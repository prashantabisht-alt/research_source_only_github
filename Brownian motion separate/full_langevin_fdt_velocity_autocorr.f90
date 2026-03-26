program full_langevin_fdt_velocity_autocorr
    implicit none
    integer, parameter :: N_particles = 100000  ! Particles
    integer, parameter :: Tmax = 5000           ! Time steps (t=50)
    integer :: i, j, timestep, seed, noise_type
    real*8 :: dt, m, gamma, kBT, D, tau
    real*8 :: v, x, v_old, x_old, friction_term, noise_term, v0
    real*8 :: gaussianvariable, eta, lambda
    real*8, allocatable :: msd_gauss(:), msd_exp(:)
    real*8, allocatable :: autocorr_gauss(:), autocorr_exp(:)
    real*8 :: grnd

    ! ---- Parameters ----
    dt = 0.01d0        ! Time step
    m = 1.0d0          ! Mass
    gamma = 0.5d0      ! Friction coefficient
    kBT = 1.0d0        ! Thermal energy
    D = kBT / gamma    ! Diffusion coefficient
    lambda = 1.0d0     ! Exponential noise parameter
    tau = m / gamma    ! Relaxation time
    seed = 6666        ! Random seed

    ! ---- Allocate arrays ----
    allocate(msd_gauss(Tmax), msd_exp(Tmax))
    allocate(autocorr_gauss(Tmax), autocorr_exp(Tmax))
    msd_gauss = 0.0d0
    msd_exp = 0.0d0
    autocorr_gauss = 0.0d0
    autocorr_exp = 0.0d0

    call sgrnd(seed)

    ! ---- Simulation loop ----
    do noise_type = 1, 2  ! 1: Gaussian, 2: Exponential
        do j = 1, N_particles
            if (mod(j, 10000) == 0) print *, "Noise Type:", noise_type, "Particle:", j

            ! Initial velocity from Maxwell-Boltzmann
            call gaussian(gaussianvariable)
            v = gaussianvariable * sqrt(kBT / m)
            x = 0.0d0
            v0 = v  ! Store initial velocity for VACF

            do timestep = 1, Tmax
                v_old = v
                x_old = x

                ! ---- Noise term ----
                if (noise_type == 1) then
                    call gaussian(gaussianvariable)
                    noise_term = gamma * sqrt(2.0d0 * D * dt) / m * gaussianvariable
                else
                    call exponential(eta, lambda)
                    noise_term = gamma * sqrt(2.0d0 * D * dt) / m * (eta / sqrt(2.0d0 / lambda**2))
                end if

                ! ---- Friction ----
                friction_term = (gamma / m) * v_old * dt

                ! ---- Velocity & position update ----
                v = v_old - friction_term + noise_term
                x = x_old + v * dt

                ! ---- Accumulate MSD & VACF ----
                if (noise_type == 1) then
                    msd_gauss(timestep) = msd_gauss(timestep) + x**2
                    autocorr_gauss(timestep) = autocorr_gauss(timestep) + v * v0
                else
                    msd_exp(timestep) = msd_exp(timestep) + x**2
                    autocorr_exp(timestep) = autocorr_exp(timestep) + v * v0
                end if
            end do
        end do
    end do

    ! ---- Normalize ----
    msd_gauss = msd_gauss / dble(N_particles)
    msd_exp = msd_exp / dble(N_particles)
    autocorr_gauss = autocorr_gauss / dble(N_particles)
    autocorr_exp = autocorr_exp / dble(N_particles)

    ! ---- Write outputs ----
    open(unit=44, file="msd_gauss_fdt_full.txt")
    do i = 1, Tmax
        write(44, *) i * dt, msd_gauss(i)
    end do
    close(44)

    open(unit=45, file="msd_exp_fdt_full.txt")
    do i = 1, Tmax
        write(45, *) i * dt, msd_exp(i)
    end do
    close(45)

    open(unit=46, file="autocorr_gauss_fdt_full.txt")
    do i = 1, Tmax
        write(46, *) i * dt, autocorr_gauss(i)
    end do
    close(46)

    open(unit=47, file="autocorr_exp_fdt_full.txt")
    do i = 1, Tmax
        write(47, *) i * dt, autocorr_exp(i)
    end do
    close(47)

    print *, "Simulation complete. Outputs: msd_*.txt, autocorr_*.txt"

  contains

    ! ---- Gaussian random generator (Box-Muller) ----
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

    ! ---- Symmetric exponential (Laplace) random generator ----
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

end program full_langevin_fdt_velocity_autocorr

include 'mt.f90'