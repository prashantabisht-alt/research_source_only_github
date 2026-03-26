program overdamped_langevin_ou
    implicit none
    integer, parameter :: N_particles = 100000, Tmax = 1000
    integer :: i, j, timestep, seed
    real*8 :: dt, gamma, kBT, D, tau_c, sigma
    real*8 :: x, noise_term, eta
    real*8 :: gaussianvariable
    real*8 :: xmin, xmax, bin_width
    integer, parameter :: nbins = 200
    integer :: bin_index
    real*8, allocatable :: hist(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    tau_c = 0.1d0  ! correlation time
    sigma = sqrt(2.0d0 * D)  ! match white noise variance
    xmin = -50.0d0
    xmax = 50.0d0
    bin_width = (xmax - xmin) / dble(nbins)
    seed = 42

    allocate(hist(nbins))
    hist = 0.0d0

    call sgrnd(seed)

    ! Loop over particles
    do j = 1, N_particles
        x = 0.0d0
        eta = 0.0d0  ! start with zero OU noise
        do timestep = 1, Tmax
            ! Update OU noise
            call gaussian(gaussianvariable)
            eta = eta * exp(-dt/tau_c) + sigma * sqrt(1.0d0 - exp(-2.0d0*dt/tau_c)) * gaussianvariable

            ! Overdamped position update
            noise_term = eta * sqrt(dt)  ! scale appropriately
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
    open(unit=10, file="pdf_ou.txt")
    do bin_index = 1, nbins
        write(10,*) xmin + (bin_index - 0.5d0) * bin_width, hist(bin_index)
    end do
    close(10)

    print *, "Done: pdf_ou.txt"

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

end program overdamped_langevin_ou

include 'mt.f90'
