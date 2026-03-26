program Brownian_motion_laplace
  implicit none
  integer, parameter :: ntraj = 10000000
  integer, parameter :: Tmax  = 200
  integer, parameter :: nbins = 1000
  integer :: j, timestep, seed, bin_index
  real*8 :: dt, diff, position, eta, b
  real*8 :: xmin, xmax, bin_width
  real*8, allocatable :: position_hist(:)
  real*8 :: grnd

  ! --- Parameters ---
  dt    = 0.001d0
  diff  = 1.0d0            ! D
  seed  = 1002
  xmin  = -10.0d0
  xmax  =  10.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! --- Laplace scale b chosen so Var(eta) = 2D/dt  (since Var(Laplace(0,b)) = 2 b^2)
  !     => 2 b^2 = 2*diff/dt  ->  b = sqrt(diff/dt)
  b = dsqrt(diff/dt)

  ! --- Allocate histogram ---
  allocate(position_hist(nbins)); position_hist = 0.0d0

  ! --- Initialize RNG ---
  call sgrnd(seed)

  ! --- Simulation loop ---
  do j = 1, ntraj
     if (mod(j,100000)==0) print*, "Trajectory:", j
     position = 0.0d0
     do timestep = 1, Tmax
        call laplace(eta, b)      ! eta ~ Laplace(0, b), zero mean, symmetric
        position = position + eta * dt
     end do

     ! Bin final position x(T)
     if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin) / bin_width) + 1
        if (bin_index > 0 .and. bin_index <= nbins) position_hist(bin_index) = position_hist(bin_index) + 1.0d0
     end if
  end do

  call write_pdf(position_hist, nbins, xmin, xmax, "laplace.txt")

  deallocate(position_hist)
  print*, "Laplace simulation complete. Output: laplace.txt"
end program Brownian_motion_laplace


! ---- Laplace(0,b) sampler via inverse-CDF (U ~ Uniform(0,1)) ----
subroutine laplace(s, b)
  implicit none
  real*8, intent(out) :: s
  real*8, intent(in)  :: b
  real*8 :: u, grnd
  do
     u = grnd()
     if (u>0.d0 .and. u<1.d0) exit
  end do
  if (u <= 0.5d0) then
     s =  b * dlog(2.d0*u)          ! negative side
  else
     s = -b * dlog(2.d0*(1.d0-u))   ! positive side
  end if
end subroutine laplace


! ---- Write normalized PDF from histogram ----
subroutine write_pdf(hist, nbins, xmin, xmax, filename)
  implicit none
  integer, intent(in) :: nbins
  real*8,  intent(in) :: hist(nbins), xmin, xmax
  character(len=*), intent(in) :: filename
  real*8 :: bin_width, norm, x
  integer :: i
  open(unit=43, file=filename)
  bin_width = (xmax - xmin) / dble(nbins)
  norm = sum(hist) * bin_width
  do i = 1, nbins
     x = xmin + (i - 0.5d0) * bin_width
     write(43,*) x, hist(i) / norm
  end do
  close(43)
end subroutine write_pdf

include 'mt.f90'
