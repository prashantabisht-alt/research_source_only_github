! ===================================================================
!       LDT Task 1: Simulation with EXPONENTIAL Noise
! ===================================================================
program exponential_simulation
  implicit none

  ! --- Parameters ---
  integer, parameter :: ntraj = 1000000 ! Number of trajectories
  integer, parameter :: Tmax = 100      ! Number of steps in each trajectory
  integer, parameter :: nbins = 500     ! Number of bins for the histogram

  ! --- Variables ---
  integer :: i, j, timestep, seed, bin_index
  real*8 :: position, eta
  real*8, allocatable :: position_hist(:)
  real*8 :: xmin, xmax, bin_width
  real*8 :: lambda
  character(len=30), parameter :: output_filename = 'pdf_exponentials.txt'

  ! --- Initial values ---
  lambda = 1.0d0
  seed = 1002
  call sgrnd(seed)

  ! Histogram range and binning
  ! For exponential noise, the distribution is asymmetric. We shift the range.
  xmin = -10.0d0
  xmax = 150.0d0
  bin_width = (xmax - xmin) / nbins
  allocate(position_hist(nbins))
  position_hist = 0.0d0

  ! --- Main simulation loop ---
  do j = 1, ntraj
    position = 0.0d0
    do timestep = 1, Tmax
      ! Get a random kick from an Exponential distribution
      call exponential(eta, lambda)
      position = position + eta
    end do

    ! Bin the final position
    if (position >= xmin .and. position < xmax) then
      bin_index = int((position - xmin) / bin_width) + 1
      position_hist(bin_index) = position_hist(bin_index) + 1
    end if
  end do

  ! Normalize and write the probability distribution function (PDF)
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj, output_filename)
  deallocate(position_hist)
  print *, 'Exponential simulation complete. Data written to ', trim(output_filename)

contains

  

  ! Exponential generator using inversion method
  subroutine exponential(s, lambda)
    implicit none
    real*8, intent(out) :: s
    real*8, intent(in)  :: lambda
    real*8 :: u, grnd
    do
      u = grnd()
      if (u > 0.d0) exit
    end do
    s = -dlog(u) / lambda
  end subroutine exponential

  ! Write normalized PDF to file
  subroutine write_pdf(hist, nbins, xmin, xmax, ntraj, fname)
    implicit none
    integer, intent(in) :: nbins, ntraj
    real*8, intent(in) :: hist(nbins), xmin, xmax
    character(len=*), intent(in) :: fname
    real*8 :: bin_width, norm_factor, x
    integer :: i
    open(unit=43, file=fname)
    bin_width = (xmax - xmin) / nbins
    norm_factor = sum(hist) * bin_width
    if (norm_factor > 0.0d0) then
        do i = 1, nbins
          x = xmin + (i - 0.5d0) * bin_width
          if (hist(i) > 0) then
             write(43, *) x, hist(i) / norm_factor
          endif
        end do
    endif
    close(43)
  end subroutine write_pdf

end program exponential_simulation
 