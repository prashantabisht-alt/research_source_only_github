program Brownian_motion_gaussian
  implicit none
  integer, parameter :: ntraj = 1000000  ! Trajectories for tail statistics
  integer, parameter :: Tmax = 1000      ! Time steps
  integer :: i, j, timestep, seed, bin_index
  real*8 :: dt, diff, eta, position, gaussianvariable, kBT, gamma
  real*8, allocatable :: position_hist(:), msd(:)
  integer, parameter :: nbins = 1000     ! Histogram bins
  real*8 :: xmin, xmax, bin_width
  real*8 :: grnd

  ! Parameters
  dt = 0.001d0       ! Time step
  diff = 1.0d0       ! Diffusion coefficient
  kBT = 1.0d0        ! Thermal energy
  gamma = kBT / diff ! Friction coefficient (Einstein relation)
  seed = 1002        ! Random seed
  xmin = -50.0d0     ! Histogram range
  xmax = 50.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! Allocate arrays
  allocate(position_hist(nbins), msd(Tmax))
  position_hist = 0.0d0
  msd = 0.0d0

  ! Initialize random number generator
  call sgrnd(seed)

  ! Simulation loop: Gaussian noise
  do j = 1, ntraj
    if (mod(j, 100000) == 0) print *, "Gaussian Trajectory:", j
    position = 0.0d0
    do timestep = 1, Tmax
      call gaussian(gaussianvariable)
      eta = gaussianvariable * dsqrt(2.0d0 * diff / dt)  ! Gaussian noise scaled for FDT
      position = position + eta * dt
      msd(timestep) = msd(timestep) + position**2
    end do
    ! Bin position
    if (position >= xmin .and. position < xmax) then
      bin_index = int((position - xmin) / bin_width) + 1
      position_hist(bin_index) = position_hist(bin_index) + 1.0d0
    end if
  end do

  ! Normalize MSD
  msd = msd / dble(ntraj)

  ! Write MSD to file
  open(unit=44, file="msd_gaussian.txt")
  do i = 1, Tmax
    write(44, *) i * dt, msd(i)
  end do
  close(44)

  ! Write normalized PDF
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj, "pdf_gaussian.txt")

  ! Deallocate arrays.
  deallocate(position_hist, msd)

  print *, "Gaussian simulation complete. Outputs: pdf_gaussian.txt, msd_gaussian.txt"

end program Brownian_motion_gaussian

! Gaussian generator
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

! Write normalized PDF to file
subroutine write_pdf(hist, nbins, xmin, xmax, ntraj, filename)
  implicit none
  integer, intent(in) :: nbins, ntraj
  real*8, intent(in) :: hist(nbins), xmin, xmax
  character(len=*), intent(in) :: filename
  real*8 :: bin_width, norm, x
  integer :: i
  open(unit=43, file=filename)
  bin_width = (xmax - xmin) / dble(nbins)
  norm = sum(hist) * bin_width
  do i = 1, nbins
    x = xmin + (i - 0.5d0) * bin_width
    write(43, *) x, hist(i) / norm
  end do
  close(43)
end subroutine write_pdf

include 'mt.f90'