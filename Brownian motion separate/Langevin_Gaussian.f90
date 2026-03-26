! ===================================================================
!         Full Langevin Simulation with GAUSSIAN Noise
!
! GOAL: Show that particle velocities thermalize to the
!       Maxwell-Boltzmann distribution.
! ===================================================================
program Langevin_Gaussian
    implicit none

    !-----Physical Parameters-----
    real(8), parameter :: m = 1.0d0, gamma = 0.5d0, kB_T = 1.0d0

    !-----Simulation Control Paramters----
    integer, parameter :: N_particles = 100000
    integer, parameter :: num_bins = 200
    real(8), parameter :: T_max = 50.0d0
    real(8), parameter :: dt = 0.01d0

    !-----Data Storage and Variables-----
    real(8), dimension(N_particles) :: velocities
    integer, dimension(num_bins) :: histogram
    real(8) :: D, v_old, friction_term, noise_term, random_kick
    real(8) :: PI, v_min, v_max, bin_width, v_center, P_sim, P_theory, norm_factor
    integer :: i, j, seed, bin_index, num_steps, i_step

    !------Initialization------
    seed = 6666
    call sgrnd(seed)
    D = kB_T / gamma ! Einstein relation
    velocities = 0.0d0

    !--------------------------------------
    ! Main simulation loop (Euler-Maruyama)
    !--------------------------------------
    print *, "Simulating with GAUSSIAN noise..."
    num_steps = nint(T_max / dt)

    do i_step = 1, num_steps
        do i = 1, N_particles
            call gaussian(random_kick)
            v_old = velocities(i)
            friction_term = (gamma / m) * v_old * dt
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * random_kick
            velocities(i) = v_old - friction_term + noise_term
        end do
    end do
    print *, "Simulation finished."

    !--------------------------------------
    ! Data Analysis: Build Histogram
    !--------------------------------------
    v_min = -5.0d0
    v_max = 5.0d0
    bin_width = (v_max - v_min) / dble(num_bins)
    histogram = 0

    do i = 1, N_particles
        if (velocities(i) > v_min .and. velocities(i) < v_max) then
            bin_index = floor((velocities(i) - v_min) / bin_width) + 1
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    !--------------------------------------
    ! Output Results to File
    !--------------------------------------
    open(unit=29, file="dist_gaussian.txt")
    write(29, *) "# v_center, P_sim_gaussian, P_theoretical"
    PI = acos(-1.0d0)

    do j = 1, num_bins
        v_center = v_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)
        norm_factor = sqrt(m / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-m * v_center**2 / (2.0d0 * kB_T))
        write(29, '(3E20.8)') v_center, P_sim, P_theory
    end do

    close(29)
    print *, "Done. Check dist_gaussian.txt"

contains
    

    subroutine gaussian(s)
        real(8), intent(out) :: s
        real(8) :: x1, x2, w, grnd
        w = 2.0d0
        do while (w >= 1.0d0 .or. w == 0.0d0)
            x1 = 2.0d0 * grnd() - 1.0d0
            x2 = 2.0d0 * grnd() - 1.0d0
            w = x1 * x1 + x2 * x2
        end do
        s = dsqrt(-2.0d0 * dlog(w) / w) * x1
    end subroutine gaussian

end program Langevin_Gaussian
