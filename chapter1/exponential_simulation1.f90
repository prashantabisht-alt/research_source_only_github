program Brownian_motion_exponential
  implicit none
  integer, parameter :: ntraj = 10000000  ! Number of trajectories
  integer, parameter :: Tmax = 200   ! Time steps
  integer, parameter :: nbins = 1000     ! Histogram bins
  integer :: j, timestep, seed, bin_index
  real*8 :: dt, diff, eta, position, kBT, gamma, lambda, mean_shift
  real*8, allocatable :: position_hist(:)
  real*8 :: xmin, xmax, bin_width
  real*8 :: grnd

  ! --- Parameters ---
  dt     = 0.001d0
  diff   = 1.0d0
  kBT    = 1.0d0
  gamma  = kBT / diff
  lambda = 1.0d0
  seed   = 1002
  xmin   = -10.0d0
  xmax   = 10.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  mean_shift = 1.0d0 / lambda   

  
  allocate(position_hist(nbins))
  position_hist = 0.0d0

  
  call sgrnd(seed)

  
  do j = 1, ntraj
    if (mod(j, 100000) == 0) print *, "Trajectory:", j
    position = 0.0d0
    do timestep = 1, Tmax
      call exponential(eta, lambda)

      
      eta = eta - mean_shift

      ! Scale to match Gaussian variance (2D/dt)
      eta = eta * dsqrt(2.0d0 * diff / dt) / dsqrt(1.0d0 / lambda**2)

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

  ! --- Write normalized PDF ---
  call write_pdf(position_hist, nbins, xmin, xmax, ntraj, "exponential_shifted.txt")

  deallocate(position_hist)

  print *, "Exponential simulation complete. Output: exponential_shifted.txt"

end program Brownian_motion_exponential



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


! --- Write normalized PDF ---
subroutine write_pdf(hist, nbins, xmin, xmax, ntraj, filename)
  implicit none
  integer, intent(in) :: nbins, ntraj
  real*8, intent(in)  :: hist(nbins), xmin, xmax
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



