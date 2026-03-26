program Brownian_motion
  implicit none

  integer, parameter :: ntraj = 1000000
  integer, parameter :: Tmax = 50
  integer :: i, j, timestep, seed, bin_index
  real*8 :: dt, diff, eta, position, gaussianvariable
  real*8, allocatable :: position_hist(:)
  integer, parameter :: nbins = 1000
  real*8 :: xmin, xmax, bin_width, x
  real*8 :: lambda

  lambda = 1.0d0  ! or whatever value you want

  ! Parameters
  dt = 0.001d0
  diff = 1.0d0
  seed = 1002

  ! Histogram range and binning
  xmin = -200.0d0
  xmax = 200.0d0
  bin_width = (xmax - xmin)/nbins
  allocate(position_hist(nbins))
  position_hist = 0.0d0

  call sgrnd(seed)

  ! Main simulation loop
  do j = 1, ntraj
    print*, "Trajectory:", j
    position = 0.0d0

    do timestep = 1, Tmax
      call gaussian(gaussianvariable)
      eta = gaussianvariable * dsqrt(2.0d0 * diff * dt)
      !call exponential(eta, lambda)
      position = position + eta
    end do

    ! Bin the final position
    if (position .ge. xmin .and. position .lt. xmax) then
      bin_index = int((position - xmin)/bin_width) + 1
      position_hist(bin_index) = position_hist(bin_index) + 1.0d0
    end if
  end do

  ! Normalize and write the PDF
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj)

  deallocate(position_hist)

end program Brownian_motion

include 'mt.f90'

!=============================================================================================
! Gaussian generator using Box-Muller transform
subroutine gaussian(s)
  implicit none
  real*8 :: x1, x2, w, s, grnd

  w = 2.0d0
  do while(w .gt. 1.0d0 .or. w == 0.0d0)
    x1 = 1.0d0 - 2.0d0 * grnd()
    x2 = 1.0d0 - 2.0d0 * grnd()
    w = x1 * x1 + x2 * x2
  end do

  s = dsqrt(-2.0d0 * dlog(w) / w) * x1
end subroutine gaussian
!=============================================================================================
!=============================================================================================
! Exponential generator 
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
!=============================================================================================

!=============================================================================================
! Write normalized PDF to file
subroutine write_pdf(hist, nbins, xmin, xmax, ntraj)
  implicit none
  integer, intent(in) :: nbins, ntraj
  real*8, intent(in) :: hist(nbins), xmin, xmax
  real*8 :: bin_width, norm, x
  integer :: i
  open(unit=43, file="pdf_new.txt")

  bin_width = (xmax - xmin)/nbins
  norm = sum(hist) * bin_width

  do i = 1, nbins
    x = xmin + (i - 0.5d0) * bin_width
    write(43,*) x, hist(i) / norm
  end do

  close(43)
end subroutine write_pdf
!=============================================================================================