program VerifyFokkerPlanck
    implicit none

    !======================================================================
    ! STEP 2: DEFINE PARAMETERS AND VARIABLES
    !======================================================================

    ! Simulation Controls
    integer, parameter :: N_particles = 100000
    integer, parameter :: num_bins = 200
    real*8, parameter :: T_max = 50.0d0
    real*8, parameter :: dt = 0.01d0

    ! Physical Parameters
    real*8, parameter :: m = 1.0d0
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: kB_T = 1.0d0

    ! Data Arrays
    real*8, dimension(N_particles) :: velocities
    integer, dimension(num_bins) :: histogram

    ! Other Variables
    real*8 :: D, time, v_old, friction_term, noise_term
    real*8 :: random_kick, grnd, PI
    real*8 :: v_min, v_max, bin_width, v_center
    real*8 :: P_sim, P_theory, norm_factor
    integer :: i, j, seed, bin_index

    !======================================================================
    ! STEP 3: INITIALIZATION
    !======================================================================

    seed = 1234
    call sgrnd(seed)

    ! Calculate Diffusion Constant from Einstein Relation
    D = kB_T / gamma

    ! Initialize all particles at rest
    do i = 1, N_particles
       velocities(i) = 0.0d0
    end do

    !======================================================================
    ! STEP 4: MAIN SIMULATION LOOP (LANGEVIN DYNAMICS)
    !======================================================================
    print*, "Simulating..."
    do time = dt, T_max, dt
        do i = 1, N_particles
            ! Get a normally distributed random number
            call gaussian(random_kick)

            ! Discretized Langevin Equation for velocity v
            v_old = velocities(i)
            friction_term = (gamma / m) * v_old * dt
            noise_term = sqrt(2.0d0 * D * dt) / m * random_kick

            velocities(i) = v_old - friction_term + noise_term
        end do
    end do
    print*, "Simulation finished."

    !======================================================================
    ! STEP 5: DATA ANALYSIS (BUILD HISTOGRAM)
    !======================================================================
    print*, "Building histogram..."
    ! Define Histogram Bounds
    v_min = -5.0d0
    v_max = 5.0d0
    bin_width = (v_max - v_min) / dble(num_bins)

    ! Initialize histogram array
    histogram = 0

    ! Bin the data
    do i = 1, N_particles
        if (velocities(i) > v_min .and. velocities(i) < v_max) then
            bin_index = floor((velocities(i) - v_min) / bin_width) + 1
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    !======================================================================
    ! STEP 6: OUTPUT FOR VISUALIZATION
    !======================================================================
    print*, "Writing output file..."
    open(unit=42, file="velocity_distribution.txt")
    write(42,*) "# v_center, P_simulated, P_theoretical" ! Header

    PI = acos(-1.0d0)

    do j = 1, num_bins
        ! Velocity at the center of the bin
        v_center = v_min + (dble(j) - 0.5d0) * bin_width

        ! Simulated probability density (normalized by area)
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)

        ! Theoretical Maxwell-Boltzmann probability density
        norm_factor = sqrt(m / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-m * v_center**2 / (2.0d0 * kB_T))

        ! Write the three columns to the file
        write(42, '(3E20.8)') v_center, P_sim, P_theory
    end do

    close(42)
    print*, "Done. Check velocity_distribution.txt"

end program VerifyFokkerPlanck

!=============================================================================================
! This subroutine generates a unit gaussian random variable using the polar Box-Muller transform
[cite_start]! with distribution 1/sqrt(2 pi) exp (-x^2/2) [cite: 770]
!=============================================================================================
subroutine gaussian(s)
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while(w > 1.0d0)
        x1 = 1.0d0 - 2.0d0*grnd()
        x2 = 1.0d0 - 2.0d0*grnd()
        w = x1*x1 + x2*x2
    enddo
    s = dsqrt(-2.0d0 * dlog(w)/w) * x1
end subroutine gaussian