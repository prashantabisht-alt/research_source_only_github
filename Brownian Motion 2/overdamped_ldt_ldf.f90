program overdamped_ldt_ldf
  implicit none
  integer, parameter :: ntraj = 1000000
  integer, parameter :: Tmax = 1000
  integer, parameter :: nbins = 200
  real*8 :: dt, diff, lambda, total_time
  real*8, allocatable :: hist_gauss(:), hist_exp(:)
  real*8 :: xmin, xmax, bin_width
  integer :: seed

  ! Parameters
  dt = 0.01d0
  diff = 1.0d0
  lambda = 1.0d0 / dsqrt(2.0d0 * diff * dt)
  xmin = -50.0d0
  xmax =  50.0d0
  seed = 12345
  total_time = Tmax * dt

  allocate(hist_gauss(nbins), hist_exp(nbins))
  hist_gauss = 0.0d0
  hist_exp = 0.0d0
  bin_width = (xmax - xmin)/nbins

  call sgrnd(seed)

  print *, "Running Gaussian noise simulation..."
  call run_simulation('gaussian', hist_gauss, lambda, dt, diff, nbins, xmin, xmax, bin_width)

  print *, "Running Laplace noise simulation..."
  call run_simulation('laplace', hist_exp, lambda, dt, diff, nbins, xmin, xmax, bin_width)

  call write_pdf(hist_gauss, nbins, xmin, xmax, "pdf_gaussian.txt")
  call write_pdf(hist_exp, nbins, xmin, xmax, "pdf_exponential.txt")

  call write_ldf(hist_gauss, nbins, xmin, xmax, total_time, "ldf_gaussian.txt")
  call write_ldf(hist_exp, nbins, xmin, xmax, total_time, "ldf_exponential.txt")

  print *, "Simulation complete."

contains

  subroutine run_simulation(noise_type, hist, lambda, dt, diff, nbins, xmin, xmax, bin_width)
    character(len=*), intent(in) :: noise_type
    integer, intent(in) :: nbins
    real*8, intent(out) :: hist(nbins)
    real*8, intent(in) :: lambda, dt, diff, xmin, xmax, bin_width
    real*8 :: position, eta, gaussianvariable
    integer :: j, timestep, bin_index

    hist = 0.0d0
    do j = 1, ntraj
      position = 0.0d0
      do timestep = 1, Tmax
        if (trim(noise_type) == 'gaussian') then
          call gaussian(gaussianvariable)
          eta = gaussianvariable * dsqrt(2.0d0 * diff * dt)
        else if (trim(noise_type) == 'laplace') then
          call laplace(eta, lambda)
          eta = eta * dsqrt(2.0d0 * diff * dt) / dsqrt(2.0d0 / lambda**2)
        end if
        position = position + eta
      end do

      if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin)/bin_width) + 1
        if (bin_index >= 1 .and. bin_index <= nbins) then
          hist(bin_index) = hist(bin_index) + 1.0d0
        end if
      end if
    end do
  end subroutine run_simulation

  subroutine write_pdf(hist, nbins, xmin, xmax, filename)
    integer, intent(in) :: nbins
    real*8, intent(in) :: hist(nbins), xmin, xmax
    character(len=*), intent(in) :: filename
    real*8 :: bin_width, norm, x
    integer :: i
    open(unit=10, file=trim(filename))
    bin_width = (xmax - xmin)/nbins
    norm = sum(hist) * bin_width
    if (norm == 0.0d0) norm = 1.0d0
    do i = 1, nbins
      x = xmin + (i - 0.5d0) * bin_width
      write(10,'(2F15.8)') x, hist(i)/norm
    end do
    close(10)
  end subroutine write_pdf

  subroutine write_ldf(hist, nbins, xmin, xmax, total_time, filename)
    integer, intent(in) :: nbins
    real*8, intent(in) :: hist(nbins), xmin, xmax, total_time
    character(len=*), intent(in) :: filename
    real*8 :: bin_width, norm, x, p
    integer :: i
    open(unit=20, file=trim(filename))
    bin_width = (xmax - xmin)/nbins
    norm = sum(hist) * bin_width
    if (norm == 0.0d0) norm = 1.0d0
    do i = 1, nbins
      x = xmin + (i - 0.5d0) * bin_width
      p = hist(i)/norm
      if (p > 0.0d0) then
        write(20,'(2F15.8)') x, -log(p)/total_time
      else
        write(20,'(2F15.8)') x, 0.0d0
      end if
    end do
    close(20)
  end subroutine write_ldf

  subroutine gaussian(s)
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
      x1 = 2.0d0*grnd() - 1.0d0
      x2 = 2.0d0*grnd() - 1.0d0
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

end program overdamped_ldt_ldf

include 'mt.f90'
