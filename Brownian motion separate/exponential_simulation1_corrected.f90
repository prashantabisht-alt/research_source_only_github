program Brownian_motion_exponential
  implicit none
  integer, parameter :: ntraj = 1000000
  integer, parameter :: Tmax = 1000


  integer, parameter :: nbins = 1000
  integer :: j, timestep, seed, bin_index
  real*8 :: dt, diff, eta, position, kBT, gamma, lambda
  real*8, allocatable :: position_hist(:)
  real*8 :: xmin, xmax, bin_width
  real*8 :: grnd

  ! --- Parameters ---
  dt    = 0.001d0       ! Time step
  diff  = 1.0d0         ! Diffusion coefficient
  kBT   = 1.0d0         ! Thermal energy
  gamma = kBT / diff    ! Einstein relation
  lambda = 1.0d0        ! Laplace distribution parameter
  seed  = 1002
  xmin  = -200.0d0
  xmax  = 200.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! --- Allocate histogram ---
  allocate(position_hist(nbins))
  position_hist = 0.0d0

  ! --- Initialize RNG ---
  call sgrnd(seed)

  ! --- Simulation loop ---
  do j = 1, ntraj
    if (mod(j, 100000) == 0) print *, "Trajectory:", j
    position = 0.0d0
    do timestep = 1, Tmax
      call exponential(eta, lambda)
      ! Scale Laplace noise to have variance = 2*diff/dt
      eta = eta * dsqrt(2.0d0 * diff / dt) / dsqrt(2.0d0 / lambda**2)
      position = position + eta * dt
    end do
    ! Bin final position
    if (position >= xmin .and. position < xmax) then
      bin_index = int((position - xmin) / bin_width) + 1
      if (bin_index > 0 .and. bin_index <= nbins) then
        position_hist(bin_index) = position_hist(bin_index) + 1.0d0
      end if
    end if
  end do

  ! --- Write normalized PDF (2 columns) ---
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj, "exponential_corrected1.txt")

  deallocate(position_hist)

  print *, "Simulation complete: exponential_corrected1.txt"
end program Brownian_motion_exponential

! --- Symmetric exponential (Laplace) generator ---
subroutine exponential(s, lambda)
  implicit none
  real*8, intent(out) :: s
  real*8, intent(in)  :: lambda
  real*8 :: u, grnd
  u = grnd()
  if (u < 0.5d0) then
    s = -dlog(grnd()) / lambda
  else
    s =  dlog(grnd()) / lambda
  end if
end subroutine exponential

! --- Write normalized PDF (2 columns) ---
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
