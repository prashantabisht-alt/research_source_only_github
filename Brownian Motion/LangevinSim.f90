! =============================================================================
! PROGRAM: VerifyFokkerPlanck
! PURPOSE:
! This program numerically demonstrates the connection between the microscopic
! Langevin equation and the macroscopic Fokker-Planck equation. It simulates
! a large number of particles whose velocities evolve according to friction
! and random thermal kicks. It then shows that the final equilibrium
! probability distribution of these velocities matches the theoretical
! Maxwell-Boltzmann distribution, which is the solution to the
! Fokker-Planck equation at equilibrium.
! =============================================================================
program VerifyFokkerPlanck
    implicit none

    !--------------------------------------------------------------------------
    ! Variable Declarations
    !--------------------------------------------------------------------------

    ! --- Simulation Control Parameters ---
    integer, parameter :: N_particles = 100000 ! Number of particles to simulate.
    integer, parameter :: num_bins = 200       ! Number of bins for the final histogram.
    real*8, parameter :: T_max = 50.0d0       ! Total simulation time.
    real*8, parameter :: dt = 0.01d0          ! The time step for numerical integration.

    ! --- Physical Parameters ---
    real*8, parameter :: m = 1.0d0            ! Mass of each particle.
    real*8, parameter :: gamma = 0.5d0        ! The friction coefficient (γ).
    real*8, parameter :: kB_T = 1.0d0         ! Thermal energy (kB * T).

    ! --- Data Storage Arrays ---
    real*8, dimension(N_particles) :: velocities ! Stores the velocity of each particle.
    integer, dimension(num_bins) :: histogram    ! Stores the binned counts of final velocities.

    ! --- Working Variables ---
    real*8 :: D                     ! The diffusion constant, derived from the Einstein relation.
    real*8 :: v_old                 ! Temporary storage for a particle's velocity.
    real*8 :: friction_term         ! The velocity change due to friction in one time step.
    real*8 :: noise_term            ! The velocity change due to random noise in one time step.
    real*8 :: random_kick           ! A random number from a standard Gaussian distribution.
    real*8 :: grnd                  ! Function for the raw random number generator.
    real*8 :: PI                    ! The constant Pi.
    real*8 :: v_min, v_max          ! The velocity range for the histogram.
    real*8 :: bin_width             ! The width of a single histogram bin.
    real*8 :: v_center              ! The velocity at the center of a histogram bin.
    real*8 :: P_sim, P_theory       ! Simulated and theoretical probability densities.
    real*8 :: norm_factor           ! Normalization constant for the Maxwell-Boltzmann distribution.
    integer :: i, j, seed, bin_index ! Loop counters and other integers.
    integer :: num_steps, i_step    ! Integer loop variables for time evolution.

    !--------------------------------------------------------------------------
    ! Initialization
    !--------------------------------------------------------------------------

    ! Seed the random number generator for reproducible results.
    seed = 1234
    call sgrnd(seed)

    ! Calculate the diffusion constant D using the Einstein Relation.
    D = kB_T / gamma

    ! Initialize all particles at rest (velocity = 0). This allows us to
    ! clearly see the system evolve towards thermal equilibrium.
    velocities = 0.0d0

    !--------------------------------------------------------------------------
    ! Main Simulation Loop (Langevin Dynamics)
    !--------------------------------------------------------------------------
    print*, "Simulating particle velocities..."
    ! Calculate the total number of integer steps required.
    num_steps = nint(T_max / dt)

    ! Loop forward in time for a fixed number of steps.
    do i_step = 1, num_steps
        ! In each time step, update the velocity of every particle.
        do i = 1, N_particles
            ! Get a normally distributed random number
            call gaussian(random_kick)

            ! This block implements the discretized Langevin equation:
            ! v_new = v_old - (gamma/m)*v_old*dt + noise
            v_old = velocities(i)

            ! 1. Calculate the deterministic change due to friction.
            friction_term = (gamma / m) * v_old * dt

            ! 2. Calculate the stochastic change due to random thermal kicks.
            !    The scaling factors ensure the noise has the correct variance
            !    consistent with the fluctuation-dissipation theorem.
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * random_kick

            ! 3. Update the particle's velocity for the next time step.
            velocities(i) = v_old - friction_term + noise_term
        end do
    end do
    print*, "Simulation finished."

    !--------------------------------------------------------------------------
    ! Data Analysis: Build a Histogram of Final Velocities
    !--------------------------------------------------------------------------
    print*, "Building histogram of final velocities..."
    ! Define the velocity range and bin width for our histogram.
    v_min = -5.0d0
    v_max = 5.0d0
    bin_width = (v_max - v_min) / dble(num_bins)

    ! Initialize all bin counts to zero.
    histogram = 0

    ! Loop through all particles and place their final velocities into bins.
    do i = 1, N_particles
        if (velocities(i) > v_min .and. velocities(i) < v_max) then
            ! Calculate the correct bin index for the particle's velocity.
            bin_index = floor((velocities(i) - v_min) / bin_width) + 1
            ! Increment the counter for that bin.
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    !--------------------------------------------------------------------------
    ! Output: Write Results to a File for Plotting
    !--------------------------------------------------------------------------
    print*, "Writing output file..."
    open(unit=42, file="velocity_distribution.txt")
    ! Write a header for clarity in the output file.
    write(42,*) "# v_center, P_simulated, P_theoretical"

    ! Define the constant PI for use in the theoretical formula.
    PI = acos(-1.0d0)

    ! Loop through each bin of the histogram to calculate and write the results.
    do j = 1, num_bins
        ! Find the velocity at the center of the current bin.
        v_center = v_min + (dble(j) - 0.5d0) * bin_width

        ! Normalize the histogram count to get a probability density.
        ! This is done by dividing by the total number of particles and
        ! the width of the bin, so the total area under the curve is 1.
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)

        ! Calculate the theoretical Maxwell-Boltzmann probability density
        ! for the velocity at the center of the bin.
        norm_factor = sqrt(m / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-m * v_center**2 / (2.0d0 * kB_T))

        ! Write the data for this bin as one line in the output file.
        write(42, '(3E20.8)') v_center, P_sim, P_theory
    end do

    close(42)
    print*, "Done. Check velocity_distribution.txt"

end program VerifyFokkerPlanck

! This line pastes the source code for the Mersenne Twister random
! number generator here, so the linker can find it.
include 'mt.f90'

!=============================================================================================
! SUBROUTINE: gaussian
! PURPOSE:
! Generates a random variable 's' from a standard Gaussian distribution.
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
