program full_langevin_multi_potential
  implicit none
  integer, parameter :: N_particles = 100000  ! Particles
  integer, parameter :: Tmax = 5000           ! Time steps (t=50)
  integer :: i, j, timestep, seed, bin_index, pot_type, noise_type
  real*8 :: dt, m, gamma, kBT, k, a, D
  real*8 :: v, x, v_old, x_old, friction_term, potential_force_term, noise_term
  real*8 :: gaussianvariable, eta, lambda
  real*8, allocatable :: position_hist_harm_gauss(:), position_hist_harm_exp(:)
  real*8, allocatable :: position_hist_quart_gauss(:), position_hist_quart_exp(:)
  real*8, allocatable :: msd_harm(:), msd_quart(:)
  integer, parameter :: nbins = 400          ! Histogram bins
  real*8 :: xmin, xmax, bin_width
  real*8 :: grnd

  ! Parameters
  dt = 0.01d0        ! Time step
  m = 1.0d0          ! Mass
  gamma = 0.5d0      ! Friction coefficient
  kBT = 1.0d0        ! Thermal energy
  D = kBT / gamma    ! Diffusion coefficient
  k = 0.2d0          ! Harmonic spring constant
  a = 0.1d0          ! Quartic coefficient
  lambda = 1.0d0     ! Exponential noise parameter
  seed = 6666        ! Random seed
  xmin = -15.0d0     ! Histogram range
  xmax = 15.0d0
  bin_width = (xmax - xmin) / dble(nbins)

  ! Allocate arrays
  allocate(position_hist_harm_gauss(nbins), position_hist_harm_exp(nbins))
  allocate(position_hist_quart_gauss(nbins), position_hist_quart_exp(nbins))
  allocate(msd_harm(Tmax), msd_quart(Tmax))
  position_hist_harm_gauss = 0.0d0
  position_hist_harm_exp = 0.0d0
  position_hist_quart_gauss = 0.0d0
  position_hist_quart_exp = 0.0d0
  msd_harm = 0.0d0
  msd_quart = 0.0d0

  call sgrnd(seed)

  ! Simulation loop for both potentials and noise types
  do pot_type = 1, 2  ! 1: Harmonic, 2: Quartic
    do noise_type = 1, 2  ! 1: Gaussian, 2: Exponential
      do j = 1, N_particles
        if (mod(j, 10000) == 0) print *, "Pot Type:", pot_type, "Noise Type:", noise_type, "Particle:", j
        x = 0.0d0  ! Initial position
        v = 0.0d0  ! Initial velocity
        do timestep = 1, Tmax
          v_old = v
          x_old = x
          if (noise_type == 1) then
            call gaussian(gaussianvariable)
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * gaussianvariable
          else
            call exponential(eta, lambda)
            noise_term = gamma * sqrt(2.0d0 * D * dt) / m * (eta / sqrt(2.0d0 / lambda**2))
          end if
          friction_term = (gamma / m) * v_old * dt
          if (pot_type == 1) potential_force_term = (k * x_old / m) * dt  ! Harmonic
          if (pot_type == 2) potential_force_term = (4.0d0 * a * x_old**3 / m) * dt  ! Quartic
          v = v_old - friction_term - potential_force_term + noise_term
          x = x_old + v * dt
          if (pot_type == 1) msd_harm(timestep) = msd_harm(timestep) + x**2
          if (pot_type == 2) msd_quart(timestep) = msd_quart(timestep) + x**2
        end do
        if (x >= xmin .and. x < xmax) then
          bin_index = int((x - xmin) / bin_width) + 1
          if (pot_type == 1 .and. noise_type == 1) position_hist_harm_gauss(bin_index) = position_hist_harm_gauss(bin_index) + 1.0d0
          if (pot_type == 1 .and. noise_type == 2) position_hist_harm_exp(bin_index) = position_hist_harm_exp(bin_index) + 1.0d0
          if (pot_type == 2 .and. noise_type == 1) position_hist_quart_gauss(bin_index) = position_hist_quart_gauss(bin_index) + 1.0d0
          if (pot_type == 2 .and. noise_type == 2) position_hist_quart_exp(bin_index) = position_hist_quart_exp(bin_index) + 1.0d0
        end if
      end do
    end do
  end do

  ! Normalize MSD
  msd_harm = msd_harm / dble(N_particles * 2)  ! Average over particles and noise types
  msd_quart = msd_quart / dble(N_particles * 2)

  ! Write MSD to files
  open(unit=44, file="msd_harm_full.txt")
  do i = 1, Tmax
    write(44, *) i * dt, msd_harm(i)
  end do
  close(44)

  open(unit=45, file="msd_quart_full.txt")
  do i = 1, Tmax
    write(45, *) i * dt, msd_quart(i)
  end do
  close(45)

  ! Write normalized PDFs
  call write_pdf(position_hist_harm_gauss, nbins, xmin, xmax, N_particles, "pdf_harm_gauss_full.txt")
  call write_pdf(position_hist_harm_exp, nbins, xmin, xmax, N_particles, "pdf_harm_exp_full.txt")
  call write_pdf(position_hist_quart_gauss, nbins, xmin, xmax, N_particles, "pdf_quart_gauss_full.txt")
  call write_pdf(position_hist_quart_exp, nbins, xmin, xmax, N_particles, "pdf_quart_exp_full.txt")

  ! Deallocate arrays
  deallocate(position_hist_harm_gauss, position_hist_harm_exp)
  deallocate(position_hist_quart_gauss, position_hist_quart_exp)
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

end program full_langevin_multi_potential

include 'mt.f90'