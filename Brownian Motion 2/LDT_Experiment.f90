program LDT_Experiment
  implicit none

  ! Parameters
  integer, parameter :: ntraj = 1000000, Tmax = 50, nbins = 1000
  real*8 :: dt, diff, lambda

  ! Variables
  integer :: seed
  real*8, allocatable :: position_hist(:)
  real*8 :: xmin, xmax

  ! --- Initialization ---
  dt = 0.001d0
  diff = 1.0d0
  seed = 1001
  xmin = -5.0d0
  xmax = 5.0d0
  lambda = 1.0d0 / dsqrt(2.0d0 * diff * dt) ! Correctly matched lambda

  call sgrnd(seed)
  allocate(position_hist(nbins))

  ! --- Run Simulations ---
  print*, "Running GAUSSIAN simulation..."
  call run_simulation('gaussian', position_hist, lambda)
  call write_pdf(position_hist, nbins, xmin, xmax, "pdf_gaussian.txt")
  print*, "Gaussian data written."

  print*, "Running EXPONENTIAL simulation..."
  call run_simulation('exponential', position_hist, lambda)
  call write_pdf(position_hist, nbins, xmin, xmax, "pdf_exponential.txt")
  print*, "Exponential data written."

  deallocate(position_hist)

contains

  subroutine run_simulation(noise_type, hist, lambda_in)
    character(len=*), intent(in) :: noise_type
    real*8, intent(out) :: hist(:)
    real*8, intent(in) :: lambda_in
    real*8 :: position, eta, gaussianvariable, bin_width
    integer :: j, timestep, bin_index

    hist = 0.0d0 ! Reset histogram for each run
    bin_width = (xmax - xmin)/nbins

    do j = 1, ntraj
      position = 0.0d0
      do timestep = 1, Tmax
        if (trim(noise_type) == 'gaussian') then
          call gaussian(gaussianvariable)
          eta = gaussianvariable * dsqrt(2.0d0 * diff * dt)
        else if (trim(noise_type) == 'exponential') then
          call exponential(eta, lambda_in)
        end if
        position = position + eta
      end do

      if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin)/bin_width) + 1
        hist(bin_index) = hist(bin_index) + 1.0d0
      end if
    end do
  end subroutine run_simulation

  subroutine write_pdf(hist, nbins, xmin, xmax, filename)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nbins
    real*8, intent(in) :: hist(nbins), xmin, xmax
    real*8 :: bin_width, norm, x
    integer :: i
    open(unit=43, file=trim(filename))
    bin_width = (xmax - xmin)/nbins
    norm = sum(hist) * bin_width
    if (norm == 0.0d0) norm = 1.0d0
    do i = 1, nbins
      x = xmin + (i - 0.5d0) * bin_width
      write(43,*) x, hist(i) / norm
    end do
    close(43)
  end subroutine write_pdf

  subroutine gaussian(s)
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while(w > 1.0d0 .or. w == 0.0d0)
      x1 = 2.0d0*grnd() - 1.0d0; x2 = 2.0d0*grnd() - 1.0d0
      w = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0 * dlog(w) / w) * x1
  end subroutine gaussian

  subroutine exponential(s, lambda)
    real*8, intent(out) :: s; real*8, intent(in) :: lambda
    real*8 :: u, grnd
    do
      u = grnd()
      if (u > 0.d0) exit
    end do
    s = -dlog(u) / lambda
  end subroutine exponential

  

end program LDT_Experiment
