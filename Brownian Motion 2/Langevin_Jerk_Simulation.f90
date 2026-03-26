! ===================================================================
!         Simulation of a "Jerky" Langevin Particle
!         (m*a_dot + gamma*a = xi)
! ===================================================================
program SimulationJerkyLangevin
    implicit none
    !-----Physical Parameters-----
    real*8, parameter :: m = 1.0d0, gamma = 0.5d0, kBT = 1.0d0
    !-----Simulation Control Paramters----
    integer, parameter :: N_particles = 100000
    integer, parameter :: num_bins = 200
    real*8, parameter :: T_max = 50.0d0, dt = 0.01d0
    !-----Data Storage Arrays-----
    real*8, dimension(N_particles) :: positions, velocities, accelerations
    integer, dimension(num_bins)   :: histogram
    !--- Working Variables ---
    real*8 :: noise_term, random_kick, a_old, bin_width
    real*8 :: a_center, v_center, x_center, P_sim
    integer :: i, j, seed, bin_index, num_steps, i_step

    !------Initialization------
    seed = 2025
    call sgrnd(seed)
    positions = 0.0d0; velocities = 0.0d0; accelerations = 0.0d0

    !-----------------------------------------
    ! Main Simulation Loop (Jerky Langevin Dynamics)
    !-----------------------------------------
    print *, "Simulating 'jerky' particle dynamics..."
    num_steps = nint(T_max / dt)

    do i_step = 1, num_steps
        do i = 1, N_particles
            call gaussian(random_kick)

            a_old      = accelerations(i)
            noise_term = sqrt(2.0d0*gamma*kBT*dt)/m * random_kick
            accelerations(i) = a_old - (gamma/m)*a_old*dt + noise_term

            velocities(i) = velocities(i) + accelerations(i)*dt
            positions(i)  = positions(i)  + velocities(i)*dt
        end do
    end do
    print *, "Simulation finished."

    ! =================================================================
    ! ACCELERATION histogram
    ! =================================================================
    print*, "Building histogram of final accelerations..."
    histogram = 0
    bin_width = 10.0d0 / dble(num_bins)           ! << moved outside loop
    do i = 1, N_particles
        if (accelerations(i) > -5.d0 .and. accelerations(i) < 5.d0) then
            bin_index = int(floor((accelerations(i)+5.d0)/bin_width)) + 1
            bin_index = min(max(1,bin_index), num_bins)               ! clamp
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    open(unit=91, file="acceleration_dist_jerky.txt")
    write(91,*) "# a_center  P_sim"
    do j = 1, num_bins
        a_center = -5.d0 + (dble(j)-0.5d0)*bin_width
        P_sim    = dble(histogram(j))/(dble(N_particles)*bin_width)
        write(91,'(2E20.8)') a_center, P_sim
    end do
    close(91)

    ! =================================================================
    ! VELOCITY histogram
    ! =================================================================
    print*, "Building histogram of final velocities..."
    histogram = 0
    bin_width = 40.0d0 / dble(num_bins)           ! << moved outside loop
    do i = 1, N_particles
        if (velocities(i) > -20.d0 .and. velocities(i) < 20.d0) then
            bin_index = int(floor((velocities(i)+20.d0)/bin_width)) + 1
            bin_index = min(max(1,bin_index), num_bins)
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    open(unit=92, file="velocity_dist_jerky.txt")
    write(92,*) "# v_center  P_sim"
    do j = 1, num_bins
        v_center = -20.d0 + (dble(j)-0.5d0)*bin_width
        P_sim    = dble(histogram(j))/(dble(N_particles)*bin_width)
        write(92,'(2E20.8)') v_center, P_sim
    end do
    close(92)

    ! =================================================================
    ! POSITION histogram
    ! =================================================================
    print*, "Building histogram of final positions..."
    histogram = 0
    bin_width = 400.0d0 / dble(num_bins)           ! << moved outside loop
    do i = 1, N_particles
        if (positions(i) > -200.d0 .and. positions(i) < 200.d0) then
            bin_index = int(floor((positions(i)+200.d0)/bin_width)) + 1
            bin_index = min(max(1,bin_index), num_bins)
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    open(unit=93, file="position_dist_jerky.txt")
    write(93,*) "# x_center  P_sim"
    do j = 1, num_bins
        x_center = -200.d0 + (dble(j)-0.5d0)*bin_width
        P_sim    = dble(histogram(j))/(dble(N_particles)*bin_width)
        write(93,'(2E20.8)') x_center, P_sim
    end do
    close(93)

    print *, "Done. Check output files."

contains
    ! Gaussian N(0,1) via Box–Muller
    subroutine gaussian(s)
        real*8, intent(out) :: s
        real*8 :: x1, x2, w, grnd
        do
            x1 = 2.d0*grnd() - 1.d0
            x2 = 2.d0*grnd() - 1.d0
            w  = x1*x1 + x2*x2
            if (w > 1.d-18 .and. w < 1.d0) exit
        end do
        s = sqrt(-2.d0*log(w)/w) * x1
    end subroutine gaussian
end program SimulationJerkyLangevin

include 'mt.f90'   ! RNG source
