!==========================================================
!GOAL:Our main imp gaal is to show collection of particles evolving by full langevin equation (not neglecting
!inertia term but V=0)for velocity will end up with velocities that follow the Maxwell Boltzam Distribution (I'm 
!using David Tong Lectures on Stochastic Process for formulaes)which is essentially gaussian which 
!solution of Fokker Planck Equation when you take derivative  of velocity Probability Distribution
! zero wrt time.This program will show like connection b/w microscopic Langevin Eq and macroscopic 
!Fokker Planck Eq.
!
!MODIFIED GOAL: This version now ALSO calculates the final position distribution P(x,t)
!and compares it to the theoretical diffusive Gaussian.
!===========================================================================
program SimulationfullLangevintwopointo
    implicit none 
    !-------------------
    !lets declare all variables acc to new fortran manuals 
    !-------------------
    !-----Physical Parameters-----
    real*8, parameter :: m=1.0d0, gamma=0.5d0, kB_T=1.0d0
    !d0 is bcz real*8 using 8 bytes its high precison (to match data type )
    
    !-----Simulation Control Paramters----
    integer, parameter :: N_particles = 1000000 ! Total Number of particles to simulate.
    integer, parameter :: num_bins = 100       ! Number of bins for the final histogram.
    real*8, parameter :: T_max = 20.0d0       ! Total simulation time.
    real*8, parameter :: dt = 0.01d0          ! The time step for numerical integration.

    !-----Data Storage arrays-----
    real*8, dimension(N_particles) :: velocities !this stores velocity for each particle
    real*8, dimension(N_particles) :: positions  ! *** newbie: Stores position for each particle
    integer, dimension(num_bins) :: histogram !this stores the binned counts for final velocties
    ! *** ADDED (minimal): 2D histogram for joint P(x,v)
    integer, dimension(num_bins, num_bins) :: hist2D

    ! --- Working Variables ---
    real*8 :: D                     ! the diffusion constant, derived from the Einstein relation.
    real*8 :: v_old                 ! Temporary storage for a particle's velocity.
    real*8 :: friction_term         ! the velocity change due to friction in one time step.
    real*8 :: noise_term            ! The velocity change due to random noise in one time step.
    real*8 :: random_kick           ! A random number from a standard Gaussian distribution.
    real*8 :: grnd                  ! Function for the raw random number generator.
    real*8 :: PI                    ! The constant Pi.
    real*8 :: v_min, v_max          ! The velocity range for the histogram.
    real*8 :: x_min, x_max          ! *** newbie: The position range for the histogram.
    real*8 :: bin_width             ! The width of a single histogram bin.
    real*8 :: v_center, x_center    ! *** little modified: Center of velocity and position bins.
    real*8 :: P_sim, P_theory       ! Simulated and theoretical probability densities.
    real*8 :: norm_factor           ! Normalization constant for the probability distribution.
    real*8 :: variance_pos_theory   ! *** ADDED: Theoretical variance for position.
    integer :: i, j, seed, bin_index ! Loop counters and other integers.
    integer :: num_steps, i_step 
    real*8 :: tau
  ! Integer loop variables for time evolution.

    ! *** ADDED  for joint P(x,v):
    real*8 :: dx, dv, P_joint
    integer :: ix, iv

!------lets initialise------
    seed=6666
    call sgrnd(seed)

    D=kB_T/gamma
    velocities=0.0d0
    positions=0.0d0

!-----------------------------------------
! Main simulation loop (Langevin Dynamics)#Algorithm-Euler-maruyama Method
!--------------------------------------
    print*,"Simulating particle velocities is happening so chill"
    num_steps = nint(T_max/dt)

    do i_step = 1, num_steps
        do i = 1, N_particles
            call gaussian(random_kick)

            v_old = velocities(i)
            friction_term = (gamma / m) * v_old * dt
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * random_kick
            velocities(i) = v_old - friction_term + noise_term
            positions(i)  = positions(i) + velocities(i) * dt
        end do
    end do
    print*, "Simulation finished."
    
    !--------------------------------------------------------------------------
    ! Velocity histogram
    !--------------------------------------------------------------------------
    print*, "Building histogram of final velocities..."
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
    
    print*, "Writing velocity output file..."
    open(unit=72, file="velocity_distributions20.txt")
    write(72,*) "# v_center, P_simulated, P_theoretical"
    PI = acos(-1.0d0)

    do j = 1, num_bins
        v_center = v_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)
        norm_factor = sqrt(m / (2.0d0 * PI * kB_T))
        P_theory = norm_factor * exp(-m * v_center**2 / (2.0d0 * kB_T))
        write(72, '(3E20.8)') v_center, P_sim, P_theory
    end do
    close(72)

    !==========================================================================
    ! Position histogram
    !==========================================================================
    print*, "Building histogram of final positions..."
    x_min = -50.0d0
    x_max = 50.0d0
    bin_width = (x_max - x_min) / dble(num_bins)
    histogram = 0

    do i = 1, N_particles
        if (positions(i) > x_min .and. positions(i) < x_max) then
            bin_index = floor((positions(i) - x_min) / bin_width) + 1
            histogram(bin_index) = histogram(bin_index) + 1
        endif
    end do

    print*, "Writing position output file..."
    open(unit=90, file="position_distributions20.txt")
    write(90,*) "# x_center, P_simulated, P_theoretical"
     tau = m/gamma
     variance_pos_theory = 2.0d0*D * ( T_max- 2.0d0*tau*(1.0d0 - exp(-T_max/tau))+ 0.5d0*tau*(1.0d0 - exp(-2.0d0*T_max/tau)) )


     do j = 1, num_bins
        x_center = x_min + (dble(j) - 0.5d0) * bin_width
        P_sim = dble(histogram(j)) / (dble(N_particles) * bin_width)
        norm_factor = 1.0d0 / sqrt(2.0d0 * PI * variance_pos_theory)
        P_theory = norm_factor * exp(-x_center**2 / (2.0d0 * variance_pos_theory))
        write(90, '(3E20.8)') x_center, P_sim, P_theory
     end do
    close(90)

    !==========================================================================
    ! Joint distribution P(x,v)
    !==========================================================================
    print*, "Building joint histogram of positions and velocities..."
    dv = (v_max - v_min) / dble(num_bins)
    dx = (x_max - x_min) / dble(num_bins)
    hist2D = 0

    do i = 1, N_particles
        if (positions(i) > x_min .and. positions(i) < x_max .and. &
            velocities(i) > v_min .and. velocities(i) < v_max) then
            ix = floor((positions(i) - x_min) / dx) + 1
            iv = floor((velocities(i) - v_min) / dv) + 1
            hist2D(ix, iv) = hist2D(ix, iv) + 1
        end if
    end do

    print*, "Writing joint distribution output file..."
    open(unit=110, file="joint_distribution20.txt")
    write(110,*) "# x_center, v_center, P_joint(x,v)"
    do ix = 1, num_bins
        x_center = x_min + (dble(ix) - 0.5d0) * dx
        do iv = 1, num_bins
            v_center = v_min + (dble(iv) - 0.5d0) * dv
            P_joint = dble(hist2D(ix, iv)) / (dble(N_particles) * dx * dv)
            write(110, '(3E20.8)') x_center, v_center, P_joint
        end do
        ! >>> ADDED: blank line to separate scans for gnuplot pm3d <<<
        write(110,*)
    end do
    close(110)     
    open(unit=150, file="phase_space_points.txt")
    write(150,*) "# x, v (final time)"
    do i = 1, N_particles
      write(150,'(2E20.8)') positions(i), velocities(i)
    end do
    close(150)


    print*, "Done. Check velocity_distributions20.txt, position_distributions20.txt, joint_distribution20.txt"

end program SimulationfullLangevintwopointo

include 'mt.f90'

!========================================================
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

