! ===================================================================
!         Full Langevin Simulation with EXPONENTIAL Noise
!
! GOAL: Test if the stationary velocity distribution is universal,
!       even with non-Gaussian thermal noise.
! ===================================================================
program Langevin_Exponential
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
    real(8) :: D, v_old, friction_term, noise_term, random_kick, lambda, noise_scaling
    real(8) :: v_min, v_max, bin_width, v_center, P_sim
    integer :: i, j, seed, bin_index, num_steps, i_step

    !------Initialization------
    seed = 6666
    lambda = 1.0d0 ! Parameter for exponential distribution
    call sgrnd(seed)
    D = kB_T / gamma ! Einstein relation
    velocities = 0.0d0

    !--------------------------------------
    ! Main simulation loop (Euler-Maruyama)
    !--------------------------------------
    print *, "Simulating with EXPONENTIAL noise..."
    num_steps = nint(T_max / dt)

    do i_step = 1, num_steps
        do i = 1, N_particles
            call exponential(random_kick, lambda)
            v_old = velocities(i)
            friction_term = (gamma / m) * v_old * dt

            ! The noise term is scaled to have the same variance as the Gaussian case.
            ! Var(Laplace) = 2/lambda^2. Var(Gaussian) = 1.
            ! We must divide by sqrt(Var(Laplace)) to normalize it.
            noise_scaling = sqrt(2.0d0) / lambda
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * (random_kick / noise_scaling)

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
    open(unit=30, file="dist_exponential.txt")
    write(30, *) "# v_center, P_sim_exponential"

    do j = 1, num_bins
        v_center = v_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)
        write(30, '(2E20.8)') v_center, P_sim
    end do

    close(30)
    print *, "Done. Check dist_exponential.txt"

contains
    

    ! Symmetric exponential (Laplace) generator
    subroutine exponential(s, lambda)
        implicit none
        real(8), intent(out) :: s
        real(8), intent(in) :: lambda
        real(8) :: u, grnd
        u = grnd()
        if (u < 0.5d0) then
            s = log(2.0d0 * u) / lambda
        else
            s = -log(2.0d0 * (1.0d0 - u)) / lambda
        endif
    end subroutine exponential

end program Langevin_Exponential
