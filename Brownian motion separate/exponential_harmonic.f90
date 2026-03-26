program Brownian_motion_exponential_harmonic
  implicit none
  integer, parameter :: ntraj = 1000000  ! Trajectories
  integer, parameter :: Tmax = 10000     ! Time steps (t=10)
  integer :: i, j, timestep, seed, bin_index
  real*8 :: dt, diff, eta, position, kBT, gamma, lambda, k
  real*8, allocatable :: position_hist(:), msd(:)
  integer, parameter :: nbins = 2000     ! Histogram bins
  real*8 :: xmin, xmax, bin_width
  real*8 :: grnd

  ! Parameters
  dt = 0.001d0       ! Time step
  diff = 1.0d0       ! Diffusion coefficient
  kBT = 1.0d0        ! Thermal energy
  gamma = kBT / diff ! Friction coefficient
  k = 1.0d0          ! Spring constant
  lambda = 1.0d0     ! Exponential noise parameter
  seed = 1002        ! Random seed
  xmin = -10.0d0     ! Histogram range
  xmax = 10.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! Allocate arrays
  allocate(position_hist(nbins), msd(Tmax))
  position_hist = 0.0d0
  msd = 0.0d0

  ! Initialize random number generator
  call sgrnd(seed)

  ! Simulation loop: Exponential noise with harmonic potential
  do j = 1, ntraj
    if (mod(j, 100000) == 0) print *, "Exponential Trajectory:", j
    position = 0.0d0
    do timestep = 1, Tmax
      call exponential(eta, lambda)
      eta = eta * dsqrt(2.0d0 * diff / dt) / dsqrt(2.0d0 / lambda**2)  ! Scaled
      position = position + (-k / gamma * position + eta / gamma) * dt
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
  open(unit=44, file="msd_exp_harmonic.txt")
  do i = 1, Tmax
    write(44, *) i * dt, msd(i)
  end do
  close(44)

  ! Write normalized PDF
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj, "pdf_exp_harmonic.txt")

  ! Deallocate arrays
  deallocate(position_hist, msd)

  print *, "Exponential simulation complete. Outputs: pdf_exp_harmonic.txt, msd_exp_harmonic.txt"

end program Brownian_motion_exponential_harmonic

! Symmetric exponential (Laplace) generator
subroutine exponential(s, lambda)
  implicit none
  real*8, intent(out) :: s
  real*8, intent(in) :: lambda
  real*8 :: u, grnd
  u = grnd()
  if (u < 0.5d0) then
    s = -dlog(grnd()) / lambda  ! Negative half
  else
    s = dlog(grnd()) / lambda   ! Positive half
  end if
end subroutine exponential

! Write normalized PDF to file
subroutine write_pdf(hist, nbins, xmin, xmax, ntraj, filename)
  implicit none
  integer, intent(in) :: nbins, ntraj
  real*8, intent(in) :: hist(nbins), xmin, xmax
  character(len=*) , intent(in) :: filename
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