program overdamped_multi_potential
  implicit none
  integer, parameter :: ntraj = 1000000  ! Trajectories
  integer, parameter :: Tmax = 1000    ! Time steps (t=10)
  integer :: i, j, timestep, seed, bin_index, pot_type
  real*8 :: dt, diff, eta, position, gaussianvariable, kBT, gamma, k, a
  real*8, allocatable :: position_hist_gauss_harm(:), position_hist_exp_harm(:)
  real*8, allocatable :: position_hist_gauss_quart(:), position_hist_exp_quart(:)
  real*8, allocatable :: msd_harm(:), msd_quart(:)
  integer, parameter :: nbins = 2000     ! Histogram bins
  real*8 :: xmin, xmax, bin_width, lambda
  real*8 :: grnd

  ! Parameters
  dt = 0.001d0       ! Time step
  diff = 1.0d0       ! Diffusion coefficient
  kBT = 1.0d0        ! Thermal energy
  gamma = kBT / diff ! Friction coefficient
  k = 1.0d0          ! Spring constant for harmonic
  a = 0.1d0          ! Coefficient for quartic
  lambda = 1.0d0     ! Exponential noise parameter
  seed = 1002        ! Random seed
  xmin = -10.0d0     ! Histogram range (adjusted for confinement)
  xmax = 10.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! Allocate arrays
  allocate(position_hist_gauss_harm(nbins), position_hist_exp_harm(nbins))
  allocate(position_hist_gauss_quart(nbins), position_hist_exp_quart(nbins))
  allocate(msd_harm(Tmax), msd_quart(Tmax))
  position_hist_gauss_harm = 0.0d0
  position_hist_exp_harm = 0.0d0
  position_hist_gauss_quart = 0.0d0
  position_hist_exp_quart = 0.0d0
  msd_harm = 0.0d0
  msd_quart = 0.0d0

  call sgrnd(seed)

  ! Simulation loop for harmonic potential
  do pot_type = 1, 2  ! 1 for harmonic, 2 for quartic
    do j = 1, ntraj
      if (mod(j, 100000) == 0) then
        if (pot_type == 1) print *, "Harmonic Gaussian Trajectory:", j
        if (pot_type == 2) print *, "Quartic Gaussian Trajectory:", j
      endif
      position = 0.0d0
      do timestep = 1, Tmax
        call gaussian(gaussianvariable)
        eta = gaussianvariable * dsqrt(2.0d0 * diff / dt)
        if (pot_type == 1) position = position + (-k / gamma * position + eta / gamma) * dt
        if (pot_type == 2) position = position + (-4.0d0 * a * position**3 / gamma + eta / gamma) * dt
        if (pot_type == 1) msd_harm(timestep) = msd_harm(timestep) + position**2
        if (pot_type == 2) msd_quart(timestep) = msd_quart(timestep) + position**2
      end do
      if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin) / bin_width) + 1
        if (pot_type == 1) position_hist_gauss_harm(bin_index) = position_hist_gauss_harm(bin_index) + 1.0d0
        if (pot_type == 2) position_hist_gauss_quart(bin_index) = position_hist_gauss_quart(bin_index) + 1.0d0
      end if
    end do

    call sgrnd(seed)  ! Reset seed for consistency

    do j = 1, ntraj
      if (mod(j, 100000) == 0) then
        if (pot_type == 1) print *, "Harmonic Exponential Trajectory:", j
        if (pot_type == 2) print *, "Quartic Exponential Trajectory:", j
      endif
      position = 0.0d0
      do timestep = 1, Tmax
        call exponential(eta, lambda)
        eta = eta * dsqrt(2.0d0 * diff / dt) / dsqrt(2.0d0 / lambda**2)
        if (pot_type == 1) position = position + (-k / gamma * position + eta / gamma) * dt
        if (pot_type == 2) position = position + (-4.0d0 * a * position**3 / gamma + eta / gamma) * dt
        if (pot_type == 1) msd_harm(timestep) = msd_harm(timestep) + position**2
        if (pot_type == 2) msd_quart(timestep) = msd_quart(timestep) + position**2
      end do
      if (position >= xmin .and. position < xmax) then
        bin_index = int((position - xmin) / bin_width) + 1
        if (pot_type == 1) position_hist_exp_harm(bin_index) = position_hist_exp_harm(bin_index) + 1.0d0
        if (pot_type == 2) position_hist_exp_quart(bin_index) = position_hist_exp_quart(bin_index) + 1.0d0
      end if
    end do
  end do

  ! Normalize MSD
  msd_harm = msd_harm / dble(ntraj)
  msd_quart = msd_quart / dble(ntraj)

  ! Write MSD to files
  open(unit=44, file="msd_harmonic.txt")
  do i = 1, Tmax
    write(44, *) i * dt, msd_harm(i)
  end do
  close(44)

  open(unit=45, file="msd_quartic.txt")
  do i = 1, Tmax
    write(45, *) i * dt, msd_quart(i)
  end do
  close(45)

  ! Write normalized PDFs
  call write_pdf(position_hist_gauss_harm, nbins, xmin, xmax, ntraj, "pdf_gauss_harmonic.txt")
  call write_pdf(position_hist_exp_harm, nbins, xmin, xmax, ntraj, "pdf_exp_harmonic.txt")
  call write_pdf(position_hist_gauss_quart, nbins, xmin, xmax, ntraj, "pdf_gauss_quartic.txt")
  call write_pdf(position_hist_exp_quart, nbins, xmin, xmax, ntraj, "pdf_exp_quartic.txt")

  ! Deallocate arrays
  deallocate(position_hist_gauss_harm, position_hist_exp_harm)
  deallocate(position_hist_gauss_quart, position_hist_exp_quart)
  deallocate(msd_harm, msd_quart)

  print *, "Simulation complete. Outputs: pdf_*.txt, msd_*.txt"

contains

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

  subroutine exponential(s, lambda)
    implicit none
    real*8, intent(out) :: s
    real*8, intent(in) :: lambda
    real*8 :: u, grnd
    u = grnd()
    if (u < 0.5d0) then
      s = -dlog(grnd()) / lambda
    else
      s = dlog(grnd()) / lambda
    endif
  end subroutine exponential

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

end program overdamped_multi_potential

include 'mt.f90'