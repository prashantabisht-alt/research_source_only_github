program overdamped_langevin_drive_pdf
    implicit none
    integer, parameter :: N_particles = 1000000, Tmax = 1000, nbins = 200
    integer :: j, timestep, seed, bin_index
    real*8 :: dt, gamma, kBT, D
    real*8 :: x, noise_term, gaussianvariable
    real*8 :: A, omega, t
    real*8, allocatable :: hist_undriven(:), hist_driven(:)
    real*8 :: xmin, xmax, bin_width
    real*8 :: grnd

    ! Parameters
    dt = 0.01d0
    gamma = 1.0d0
    kBT = 1.0d0
    D = kBT / gamma
    A = 1.0d0       ! driving amplitude for driven case
    omega = 1.0d0   ! driving frequency
    seed = 2025

    xmin = -50.0d0
    xmax = 50.0d0
    bin_width = (xmax - xmin) / dble(nbins)

    allocate(hist_undriven(nbins), hist_driven(nbins))
    hist_undriven = 0.0d0
    hist_driven = 0.0d0

    call sgrnd(seed)

    ! ---- Undriven case ----
    do j = 1, N_particles
        x = 0.0d0
        do timestep = 1, Tmax
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            x = x + noise_term
        end do
        if (x >= xmin .and. x < xmax) then
            bin_index = int((x - xmin) / bin_width) + 1
            hist_undriven(bin_index) = hist_undriven(bin_index) + 1.0d0
        end if
    end do

    ! ---- Driven case ----
    do j = 1, N_particles
        x = 0.0d0
        do timestep = 1, Tmax
            t = timestep * dt
            call gaussian(gaussianvariable)
            noise_term = sqrt(2.0d0 * D * dt) * gaussianvariable
            x = x + (A * sin(omega * t) / gamma) * dt + noise_term
        end do
        if (x >= xmin .and. x < xmax) then
            bin_index = int((x - xmin) / bin_width) + 1
            hist_driven(bin_index) = hist_driven(bin_index) + 1.0d0
        end if
    end do

    ! ---- Normalise to PDF ----
    hist_undriven = hist_undriven / (sum(hist_undriven) * bin_width)
    hist_driven = hist_driven / (sum(hist_driven) * bin_width)

    ! ---- Write outputs ----
    open(unit=10, file="pdf_undriven.txt")
    do bin_index = 1, nbins
        t = xmin + (bin_index - 0.5d0) * bin_width
        write(10,*) t, hist_undriven(bin_index)
    end do
    close(10)

    open(unit=11, file="pdf_driven.txt")
    do bin_index = 1, nbins
        t = xmin + (bin_index - 0.5d0) * bin_width
        write(11,*) t, hist_driven(bin_index)
    end do
    close(11)

    print *, "PDFs written: pdf_undriven.txt, pdf_driven.txt"

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

end program overdamped_langevin_drive_pdf

include 'mt.f90'
