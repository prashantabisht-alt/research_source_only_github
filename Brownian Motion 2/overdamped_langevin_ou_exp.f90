program overdamped_langevin_ou_exp
    implicit none
    integer, parameter :: N_particles = 100000, Tmax = 1000
    integer :: j, timestep, seed, bin_index
    real*8 :: dt, gamma, kBT, D, tau_c, sigma, lambda
    real*8 :: x, eta, noise_term, expvariable
    real*8 :: xmin, xmax, bin_width
    integer, parameter :: nbins = 200
    real*8, allocatable :: hist(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    tau_c = 0.1d0   ! correlation time
    lambda = 1.0d0  ! exponential rate parameter
    sigma = sqrt(2.0d0 * D)  ! variance scaling
    xmin = -50.0d0
    xmax = 50.0d0
    bin_width = (xmax - xmin) / dble(nbins)
    seed = 777

    allocate(hist(nbins))
    hist = 0.0d0

    call sgrnd(seed)

    ! Loop over particles
    do j = 1, N_particles
        x = 0.0d0
        eta = 0.0d0
        do timestep = 1, Tmax
            ! Generate symmetric exponential (Laplace) white noise
            call exponential_laplace(expvariable, lambda)

            ! OU update with Laplace noise
            eta = eta * exp(-dt/tau_c) + sigma * sqrt(1.0d0 - exp(-2.0d0*dt/tau_c)) * expvariable

            ! Overdamped position update
            noise_term = eta * sqrt(dt)
            x = x + noise_term
        end do

        if (x >= xmin .and. x < xmax) then
            bin_index = int((x - xmin) / bin_width) + 1
            hist(bin_index) = hist(bin_index) + 1.0d0
        end if
    end do

    ! Normalize
    hist = hist / (sum(hist) * bin_width)

    ! Output
    open(unit=10, file="pdf_ou_exp.txt")
    do bin_index = 1, nbins
        write(10,*) xmin + (bin_index - 0.5d0) * bin_width, hist(bin_index)
    end do
    close(10)

    print *, "Done: pdf_ou_exp.txt"

contains

    subroutine exponential_laplace(s, lambda)
        implicit none
        real*8, intent(out) :: s
        real*8, intent(in) :: lambda
        real*8 :: u, grnd
        u = grnd()
        if (u < 0.5d0) then
            s = -dlog(2.0d0*u) / lambda
        else
            s =  dlog(2.0d0*(1.0d0-u)) / lambda
        end if
    end subroutine exponential_laplace

end program overdamped_langevin_ou_exp

include 'mt.f90'
