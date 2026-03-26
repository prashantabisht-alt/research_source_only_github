program overdamped_ldt
  implicit none
  integer, parameter :: ntraj = 1000000, Tmax = 1000, nbins = 200
  real*8 :: dt, diff, lambda, ttotal, bin_width, xmin, xmax
  real*8, allocatable :: hist_gauss(:), hist_laplace(:)
  integer :: seed

  ! Parameters
  dt   = 0.01d0
  diff = 1.0d0
  xmin = -50.0d0
  xmax =  50.0d0
  lambda = 1.0d0 / dsqrt(2.0d0 * diff * dt)  ! for Laplace scaling
  seed = 1234
  ttotal = Tmax * dt

  allocate(hist_gauss(nbins), hist_laplace(nbins))
  hist_gauss = 0.0d0
  hist_laplace = 0.0d0
  bin_width = (xmax - xmin) / nbins

  call sgrnd(seed)

  print*, "Running Gaussian noise simulation..."
  call run_simulation('gaussian', hist_gauss, lambda, dt, diff, nbins, xmin, xmax, bin_width)

  print*, "Running Laplace noise simulation..."
  call run_simulation('laplace', hist_laplace, lambda, dt, diff, nbins, xmin, xmax, bin_width)

  call write_pdf(hist_gauss, nbins, xmin, xmax, "pdf_gaussian.txt", bin_width)
  call write_pdf(hist_laplace, nbins, xmin, xmax, "pdf_laplace.txt", bin_width)

  call write_ldf(hist_gauss, nbins, xmin, xmax, "ldf_gaussian.txt", bin_width, ttotal)
  call write_ldf(hist_laplace, nbins, xmin, xmax, "ldf_laplace.txt", bin_width, ttotal)

  print*, "Done. Output files written."

contains

  subroutine run_simulation(noise_type, hist, lambda, dt, diff, nbins, xmin, xmax, bin_width)
    character(len=*), intent(in) :: noise_type
    real*8, intent(out) :: hist(:)
    real*8, intent(in) :: lambda, dt, diff, xmin, xmax, bin_width
    integer, intent(in) :: nbins
    real*8 :: position, eta, gv
    integer :: i, j, bin_index

    do j = 1, ntraj
      position = 0.0d0
      do i = 1, Tmax
        if (trim(noise_type) == 'gaussian') then
          call gaussian(gv)
          eta = gv * dsqrt(2.0d0 * diff * dt)
        else
          call laplace(eta, lambda)
          eta = eta * dsqrt(2.0d0 * diff * dt) / dsqrt(2.0d0 / lambda**2)
        end if
        position = position + eta
      end do

      if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin) / bin_width) + 1
        hist(bin_index) = hist(bin_index) + 1.0d0
      end if
    end do
  end subroutine run_simulation

  subroutine write_pdf(hist, nbins, xmin, xmax, filename, bin_width)
    integer, intent(in) :: nbins
    real*8, intent(in) :: hist(nbins), xmin, xmax, bin_width
    character(len=*), intent(in) :: filename
    real*8 :: norm, x
    integer :: i
    norm = sum(hist) * bin_width
    open(unit=10, file=filename)
    do i = 1, nbins
      if (hist(i) > 0.0d0) then
        x = xmin + (i - 0.5d0) * bin_width
        write(10,*) x, hist(i) / norm
      end if
    end do
    close(10)
  end subroutine write_pdf

  subroutine write_ldf(hist, nbins, xmin, xmax, filename, bin_width, ttotal)
    integer, intent(in) :: nbins
    real*8, intent(in) :: hist(nbins), xmin, xmax, bin_width, ttotal
    character(len=*), intent(in) :: filename
    real*8 :: norm, x, p
    integer :: i
    norm = sum(hist) * bin_width
    open(unit=11, file=filename)
    do i = 1, nbins
      p = hist(i) / norm
      if (p > 0.0d0) then
        x = xmin + (i - 0.5d0) * bin_width
        write(11,*) x, -dlog(p) / ttotal
      end if
    end do
    close(11)
  end subroutine write_ldf

  subroutine gaussian(s)
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
      x1 = 2.0d0 * grnd() - 1.0d0
      x2 = 2.0d0 * grnd() - 1.0d0
      w = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0 * dlog(w) / w) * x1
  end subroutine gaussian

  subroutine laplace(s, lambda)
    real*8, intent(out) :: s
    real*8, intent(in) :: lambda
    real*8 :: u, grnd
    u = grnd() - 0.5d0
    s = -sign(1.0d0, u) * dlog(1.0d0 - 2.0d0 * abs(u)) / lambda
  end subroutine laplace

end program overdamped_ldt

include 'mt.f90'
