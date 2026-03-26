program Langevin_In_Periodic_Improved
  implicit none

  !-----Physical Parameters-----
  real(8), parameter :: gamma = 1.0d0, kBT = 0.5d0
  real(8), parameter :: V0 = 1.0d0, k = 2.0d0

  !-----Simulation Control Parameters----
  integer, parameter :: N_particles = 100000
  integer, parameter :: num_bins = 500
  real(8), parameter :: T_max = 500.0d0    ! Extended time
  real(8), parameter :: dt = 0.01d0

  !-----Data Storage and Variables-----
  real(8), dimension(N_particles) :: positions
  real(8), allocatable :: position_hist_gauss(:), position_hist_exp(:), msd(:)
  real(8) :: diff, random_kick, lambda, force
  real(8) :: x_min, x_max, bin_width
  integer :: i, seed, bin_index, num_steps, i_step

  !------Initialization------
  seed = 2026
  lambda = 1.0d0
  diff = kBT / gamma
  x_min = -15.0d0    ! Wider range
  x_max = 15.0d0
  bin_width = (x_max - x_min) / dble(num_bins)
  num_steps = nint(T_max / dt)

  allocate(position_hist_gauss(num_bins), position_hist_exp(num_bins), msd(num_steps))
  position_hist_gauss = 0.0d0
  position_hist_exp = 0.0d0
  msd = 0.0d0

  call sgrnd(seed)

  !======================================
  ! 1. SIMULATION WITH GAUSSIAN NOISE
  !======================================
  print *, "Simulating with GAUSSIAN noise in a periodic potential..."
  positions = 0.0d0

  do i_step = 1, num_steps
    do i = 1, N_particles
      force = -V0 * k * sin(k * positions(i))
      call gaussian(random_kick)
      positions(i) = positions(i) + (force / gamma) * dt + (sqrt(2.0d0 * diff * dt) / gamma) * random_kick
      msd(i_step) = msd(i_step) + positions(i)**2  ! Track MSD
    end do
  end do

  position_hist_gauss = 0
  do i = 1, N_particles
    if (positions(i) > x_min .and. positions(i) < x_max) then
      bin_index = floor((positions(i) - x_min) / bin_width) + 1
      position_hist_gauss(bin_index) = position_hist_gauss(bin_index) + 1
    endif
  end do
  call write_pdf(position_hist_gauss, num_bins, x_min, x_max, N_particles, "pos_dist_periodic_gaussian.txt")
  msd = msd / dble(N_particles)  ! Normalize MSD
  call write_msd(msd, num_steps, dt, "msd_periodic_gaussian.txt")
  print *, "Gaussian simulation finished."

  !======================================
  ! 2. SIMULATION WITH EXPONENTIAL NOISE
  !======================================
  print *, "Simulating with EXPONENTIAL noise in a periodic potential..."
  call sgrnd(seed)
  positions = 0.0d0
  msd = 0.0d0

  do i_step = 1, num_steps
    do i = 1, N_particles
      force = -V0 * k * sin(k * positions(i))
      call exponential(random_kick, lambda)
      positions(i) = positions(i) + (force / gamma) * dt + (sqrt(2.0d0 * diff * dt) / gamma) * (random_kick / sqrt(2.0d0))
      msd(i_step) = msd(i_step) + positions(i)**2
    end do
  end do

  position_hist_exp = 0
  do i = 1, N_particles
    if (positions(i) > x_min .and. positions(i) < x_max) then
      bin_index = floor((positions(i) - x_min) / bin_width) + 1
      position_hist_exp(bin_index) = position_hist_exp(bin_index) + 1
    endif
  end do
  call write_pdf(position_hist_exp, num_bins, x_min, x_max, N_particles, "pos_dist_periodic_exponential.txt")
  msd = msd / dble(N_particles)
  call write_msd(msd, num_steps, dt, "msd_periodic_exponential.txt")
  print *, "Exponential simulation finished."

  deallocate(position_hist_gauss, position_hist_exp, msd)

contains

  subroutine gaussian(s)
    real(8), intent(out) :: s
    real(8) :: x1, x2, w, grnd
    w = 2.0d0
    do while (w >= 1.0d0 .or. w == 0.0d0)
      x1 = 2.0d0 * grnd() - 1.0d0
      x2 = 2.0d0 * grnd() - 1.0d0
      w = x1 * x1 + x2 * x2
    end do
    s = dsqrt(-2.0d0 * dlog(w) / w) * x1
  end subroutine gaussian

  subroutine exponential(s, lambda)
    implicit none
    real(8), intent(out) :: s
    real(8), intent(in) :: lambda
    real(8) :: u, grnd
    u = grnd()
    if (u < 0.5d0) then
      s = log(2.0d0 * u) / lambda
    else
      s = -log(2.0d0 * (1.0d0 - u)) / lambda
    endif
  end subroutine exponential

  subroutine write_pdf(hist, nbins, xmin, xmax, ntraj, filename)
    implicit none
    integer, intent(in) :: nbins, ntraj
    real*8, intent(in) :: hist(nbins), xmin, xmax
    character(len=*), intent(in) :: filename
    real*8 :: bin_width, norm, x
    integer :: i
    open(unit=34, file=filename, status='replace')
    write(34, *) "# x_center, P_sim"
    bin_width = (xmax - xmin) / dble(nbins)
    norm = sum(hist) * bin_width
    if (norm > 0.0d0) then
      do i = 1, nbins
        x = xmin + (dble(i) - 0.5d0) * bin_width
        write(34, '(2E20.8)') x, hist(i) / norm
      end do
    endif
    close(34)
  end subroutine write_pdf

  subroutine write_msd(msd, nsteps, dt, filename)
    implicit none
    integer, intent(in) :: nsteps
    real*8, intent(in) :: msd(nsteps), dt
    character(len=*), intent(in) :: filename
    integer :: i
    open(unit=35, file=filename, status='replace')
    write(35, *) "# t, MSD"
    do i = 1, nsteps
      write(35, '(2E20.8)') i * dt, msd(i)
    end do
    close(35)
  end subroutine write_msd

end program Langevin_In_Periodic_Improved

include 'mt.f90'