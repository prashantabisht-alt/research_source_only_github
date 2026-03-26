program overdamped_langevin_drive_compare
    implicit none
    integer, parameter :: N_particles = 100000, Tmax = 1000
    integer :: i, j, timestep, seed, run_case
    real*8 :: dt, gamma, kBT, D
    real*8 :: x, noise_term, gaussianvariable
    real*8 :: A, omega, t
    real*8, allocatable :: msd_undriven(:), meanpos_undriven(:)
    real*8, allocatable :: msd_driven(:), meanpos_driven(:)
    real*8 :: grnd

    ! ---- Parameters ----
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    A = 1.0d0       ! driving amplitude for driven case
    omega = 1.0d0   ! driving frequency
    seed = 1234

    allocate(msd_undriven(Tmax), meanpos_undriven(Tmax))
    allocate(msd_driven(Tmax), meanpos_driven(Tmax))
    msd_undriven = 0.0d0
    meanpos_undriven = 0.0d0
    msd_driven = 0.0d0
    meanpos_driven = 0.0d0

    call sgrnd(seed)

    ! =================================================
    ! First run: Undriven case (A = 0)
    ! =================================================
    do j = 1, N_particles
        x = 0.0d0
        do timestep = 1, Tmax
            t = timestep * dt
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            x = x + noise_term
            msd_undriven(timestep) = msd_undriven(timestep) + x**2
            meanpos_undriven(timestep) = meanpos_undriven(timestep) + x
        end do
    end do

    ! =================================================
    ! Second run: Driven case (A > 0)
    ! =================================================
    do j = 1, N_particles
        x = 0.0d0
        do timestep = 1, Tmax
            t = timestep * dt
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            x = x + (A * sin(omega * t) / gamma) * dt + noise_term
            msd_driven(timestep) = msd_driven(timestep) + x**2
            meanpos_driven(timestep) = meanpos_driven(timestep) + x
        end do
    end do

    ! Normalize
    msd_undriven = msd_undriven / dble(N_particles)
    meanpos_undriven = meanpos_undriven / dble(N_particles)
    msd_driven = msd_driven / dble(N_particles)
    meanpos_driven = meanpos_driven / dble(N_particles)

    ! Write outputs
    open(unit=10, file="mean_pos_compare.txt")
    do i = 1, Tmax
        write(10, *) i * dt, meanpos_undriven(i), meanpos_driven(i)
    end do
    close(10)

    open(unit=11, file="msd_compare.txt")
    do i = 1, Tmax
        write(11, *) i * dt, msd_undriven(i), msd_driven(i)
    end do
    close(11)

    print *, "Simulation complete. Outputs: mean_pos_compare.txt, msd_compare.txt"

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

end program overdamped_langevin_drive_compare

include 'mt.f90'