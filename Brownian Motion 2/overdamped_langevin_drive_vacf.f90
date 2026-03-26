program overdamped_langevin_drive_vacf
    implicit none
    integer, parameter :: N_particles = 100000, Tmax = 1000
    integer :: i, j, timestep, seed
    real*8 :: dt, gamma, kBT, D
    real*8 :: x, v, v_old, noise_term, gaussianvariable
    real*8 :: A, omega, t, v0
    real*8, allocatable :: vacf_undriven(:), vacf_driven(:)
    real*8 :: grnd

    ! ---- Parameters ----
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    A = 1.0d0
    omega = 1.0d0
    seed = 5678

    allocate(vacf_undriven(Tmax), vacf_driven(Tmax))
    vacf_undriven = 0.0d0
    vacf_driven = 0.0d0

    call sgrnd(seed)

    ! =================================================
    ! Case 1: Undriven (A = 0)
    ! =================================================
    do j = 1, N_particles
        v = sqrt(kBT) * 0.0d0   ! start with v(0) = 0 for overdamped case
        x = 0.0d0
        v0 = v
        do timestep = 1, Tmax
            t = timestep * dt
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            v_old = v
            v = ( -gamma * v_old * dt + noise_term ) / dt   ! Overdamped v ~ Δx/dt
            x = x + v * dt
            vacf_undriven(timestep) = vacf_undriven(timestep) + v * v0
        end do
    end do

    ! =================================================
    ! Case 2: Driven (A > 0)
    ! =================================================
    do j = 1, N_particles
        v = 0.0d0
        x = 0.0d0
        v0 = v
        do timestep = 1, Tmax
            t = timestep * dt
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            v_old = v
            v = ( A * sin(omega * t) / gamma * dt + noise_term - gamma * v_old * dt ) / dt
            x = x + v * dt
            vacf_driven(timestep) = vacf_driven(timestep) + v * v0
        end do
    end do

    ! Normalize
    vacf_undriven = vacf_undriven / dble(N_particles)
    vacf_driven = vacf_driven / dble(N_particles)

    ! Output
    open(unit=10, file="vacf_compare.txt")
    do i = 1, Tmax
        write(10, *) i * dt, vacf_undriven(i), vacf_driven(i)
    end do
    close(10)

    print *, "Simulation complete. Output: vacf_compare.txt"

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

end program overdamped_langevin_drive_vacf

include 'mt.f90'
