!========================================================================
!GOAL:Our main imp gaal is to show collection of particles evolving by full langevin equation (not neglecting
!inertia term but WITH A POTENTIAL V(x))for velocity will end up with velocities that follow the Maxwell Boltzam Distribution
!and positions that follow the Boltzmann Distribution P(x) ~ exp(-V(x)/k_B T).
!===========================================================================
program SimulationWithPotential
    implicit none 
    !-------------------
    !lets declare all variables acc to new fortran manuals 
    !-------------------
    !-----Physical Parameters-----
    real*8, parameter::m=1.0d0,gamma=0.5d0,kB_T=1.0d0
    real*8, parameter :: k = 0.2d0  ! *** ADDED: Spring constant for the harmonic potential V(x)=0.5*k*x^2
    !d0 is bcz real*8 using 8 bytes its high precison (to match data type )
    
    !-----Simulation Control Paramters----
    integer, parameter :: N_particles = 100000 ! Total Number of particles to simulate.
    integer, parameter :: num_bins = 200       ! Number of bins for the final histogram.
    real*8, parameter :: T_max = 50.0d0       ! Total simulation time.
    real*8, parameter :: dt = 0.01d0          ! The time step for numerical integration.

    !-----Data Storage arrays-----
    real*8, dimension(N_particles)::velocities !this stores velocity for each particle
    real*8, dimension(N_particles)::positions  ! *** newbie: Stores position for each particle
    integer,dimension(num_bins)::histogram !this stores the binned counts for final velocties
    
    ! --- Working Variables ---
    real*8 :: D                     ! the diffusion constant, derived from the Einstein relation.
    real*8 :: v_old                 ! Temporary storage for a particle's velocity.
    real*8 :: friction_term         ! the velocity change due to friction in one time step.
    real*8 :: noise_term            ! The velocity change due to random noise in one time step.
    real*8 :: potential_force_term  ! *** ADDED: The velocity change due to the potential.
    real*8 :: random_kick           ! A random number from a standard Gaussian distribution.
    real*8 :: grnd                  ! Function for the raw random number generator.
    real*8 :: PI                    ! The constant Pi.
    real*8 :: v_min, v_max          ! The velocity range for the histogram.
    real*8 :: x_min, x_max          ! *** newbie: The position range for the histogram.
    real*8 :: bin_width             ! The width of a single histogram bin.
    real*8 :: v_center, x_center    ! *** little modified: Center of velocity and position bins.
    real*8 :: P_sim, P_theory       ! Simulated and theoretical probability densities.
    real*8 :: norm_factor           ! Normalization constant for the probability distribution.
    integer :: i, j, seed, bin_index ! Loop counters and other integers.
    integer :: num_steps, i_step    ! Integer loop variables for time evolution.

!------lets initialise------
    seed=6666 !if i include saqme seed it will give same pre determined seq of random numbers (helps in reproduce result)
    call sgrnd(seed) !this subrutine inclded in mt.f90 file (mersenne twister algorithm)#seed generate random

    D=kB_T/gamma !(thats einstein relation)D stands for diffusion constant 
    velocities=0.0d0
    positions=0.0d0  ! *** ADDED: Initialize all positions to zero.

!-----------------------------------------
! Main simulation loop (Langevin Dynamics)#Algorithm-Euler-maruyama Method
!--------------------------------------
    print*,"Simulating particle velocities is happening so chill"
    num_steps = nint(T_max/dt)!nint is inbuild standard fortran function Nearest INTeger bcz resut of division is real 5000.0 but mr fortran need an integer to tell how many times to run

! do outer loop - advancing time & inner loop nested in time loop for each single step it updates every particle one by one 
    do i_step = 1, num_steps
        ! In each time step, update the velocity of every particle.
        do i = 1, N_particles
            ! Get a normally distributed random number
            call gaussian(random_kick)!unique no. from gaussian distribution bcz white nause is gaussian (gaussian random variable)

            v_old = velocities(i)!temporary copies of particle velocities before change

            ! 1. Calculate the deterministic change due to friction.(thats euler part)
            friction_term = (gamma / m) * v_old * dt

            ! 2. Calculate the stochastic change due to random thermal kicks.
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * random_kick

            ! 3. *** ADDED: Calculate the deterministic change due to the potential's force. ***
            !    The force is F = -k*x, so the change in velocity is (F/m)*dt.
            potential_force_term = (k * positions(i) / m) * dt

            ! 4. Update the particle's velocity for the next time step with all three terms.
            velocities(i) = v_old - friction_term - potential_force_term + noise_term

            ! 5. Update the particle's position using the new velocity.
            positions(i) = positions(i) + velocities(i) * dt
        end do
    end do
    print*, "Simulation finished."
    
    !--------------------------------------------------------------------------
    ! Data Analysis: Build a Histogram of Final Velocities
    !--------------------------------------------------------------------------
    print*, "Building histogram of final velocities..."
    v_min = -5.0d0
    v_max = 5.0d0
    bin_width = (v_max - v_min) / dble(num_bins)
    histogram = 0
    do i = 1, N_particles
        if (velocities(i) > v_min .and. velocities(i) < v_max) then !(make sure it falls under -5 to 5)
            bin_index = floor((velocities(i) - v_min) / bin_width) + 1
            histogram(bin_index) = histogram(bin_index) + 1! adding 1 bcz fortan are numbered satrting from 1.
        endif
    end do
    
    !--------------------------------------------------------------------------
    ! Output: Writing Results to a File for Plotting (VELOCITY)
    !--------------------------------------------------------------------------
    print*, "Writing velocity output file..."
    open(unit=29, file="velocity_distribution_potential.txt")
    write(29,*) "# v_center, P_simulated, P_theoretical"
    PI = acos(-1.0d0)
    do j = 1, num_bins
        v_center = v_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)
        norm_factor = sqrt(m / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-m * v_center**2 / (2.0d0 * kB_T))
        write(29, '(3E20.8)') v_center, P_sim, P_theory
    end do
    close(29)

    !==========================================================================
    ! *** MODIFIED SECTION: Data Analysis and Output for Final Positions ***
    !==========================================================================

    print*, "Building histogram of final positions..."
    ! The position range will be smaller now because the potential traps the particles.
    x_min = -10.0d0
    x_max = 10.0d0
    bin_width = (x_max - x_min) / dble(num_bins)
    histogram = 0
    do i = 1, N_particles
        if (positions(i) > x_min .and. positions(i) < x_max) then
            bin_index = floor((positions(i) - x_min) / bin_width) + 1
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    print*, "Writing position output file..."
    open(unit=30, file="position_distribution_potential1.txt")
    write(30,*) "# x_center, P_simulated, P_theoretical"

    ! *** MODIFIED: The theoretical distribution is now the Boltzmann distribution. ***
    ! For V(x) = 0.5*k*x^2, this is a Gaussian with variance = kB_T / k.
    do j = 1, num_bins
        x_center = x_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)

        ! Calculate the theoretical Boltzmann probability density for position.
        norm_factor = sqrt(k / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-0.5d0 * k * x_center**2 / kB_T)

        write(30, '(3E20.8)') x_center, P_sim, P_theory
    end do
    close(30)

    print*, "Done. Check velocity_distribution_potential.txt and position_distribution_potential1.txt"

end program SimulationWithPotential

! This  pastes the source code for the Mersenne Twister random
! number generator.
include 'mt.f90'

!=============================================================================================
!This subroutine generates a unit gaussian random variable using the polar Box-Muller transform
!with distribution 1/sqrt(2 pi) exp (-x^2/2)
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
