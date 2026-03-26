program ou_noise_all_cases
    implicit none
    integer, parameter :: N_particles = 100000, Tmax = 1000, nbins = 200
    integer :: j, timestep, seed, bin_index
    integer :: noise_type_idx, tau_idx
    real*8, parameter :: tau_list(3) = (/ 0.0d0, 0.1d0, 1.0d0 /)
    character(len=10) :: noise_type(2)
    character(len=50) :: filename
    real*8 :: dt, gamma, kBT, D, tau_c, sigma, lambda
    real*8 :: x, eta, noise_term, rndvar
    real*8 :: xmin, xmax, bin_width
    real*8, allocatable :: hist(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    lambda = 1.0d0      ! Laplace rate
    sigma = sqrt(2.0d0 * D)
    xmin = -50.0d0
    xmax = 50.0d0
    bin_width = (xmax - xmin) / dble(nbins)
    seed = 12345

    noise_type(1) = 'gauss'
    noise_type(2) = 'exp'

    allocate(hist(nbins))

    call sgrnd(seed)

    do noise_type_idx = 1, 2
        do tau_idx = 1, 3
            tau_c = tau_list(tau_idx)
            hist = 0.0d0

            ! --- Simulate ---
            do j = 1, N_particles
                x = 0.0d0
                eta = 0.0d0
                do timestep = 1, Tmax
                    ! Generate noise
                    if (noise_type_idx == 1) then
                        call gaussian(rndvar)
                    else
                        call laplace(rndvar, lambda)
                    end if

                    ! OU update
                    if (tau_c == 0.0d0) then
                        ! White noise case
                        eta = sigma * rndvar
                    else
                        eta = eta * exp(-dt/tau_c) + sigma * sqrt(1.0d0 - exp(-2.0d0*dt/tau_c)) * rndvar
                    end if

                    ! Position update (overdamped)
                    noise_term = eta * sqrt(dt)
                    x = x + noise_term
                end do

                ! Histogram
                if (x >= xmin .and. x < xmax) then
                    bin_index = int((x - xmin) / bin_width) + 1
                    hist(bin_index) = hist(bin_index) + 1.0d0
                end if
            end do

            ! Normalise to PDF
            hist = hist / (sum(hist) * bin_width)

            ! Output file
            if (tau_c == 0.0d0) then
                write(filename,'("pdf_",a,"_tau0.txt")') trim(noise_type(noise_type_idx))
            else if (tau_c == 0.1d0) then
                write(filename,'("pdf_",a,"_tau0.1.txt")') trim(noise_type(noise_type_idx))
            else
                write(filename,'("pdf_",a,"_tau1.0.txt")') trim(noise_type(noise_type_idx))
            end if

            open(unit=10, file=trim(filename))
            do bin_index = 1, nbins
                write(10,*) xmin + (bin_index - 0.5d0) * bin_width, hist(bin_index)
            end do
            close(10)

            print *, "Wrote:", trim(filename)

        end do
    end do

    print *, "All cases complete."

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

    subroutine laplace(s, lambda)
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
    end subroutine laplace

end program ou_noise_all_cases

include 'mt.f90'
