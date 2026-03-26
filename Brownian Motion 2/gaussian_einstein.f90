program Brownian_motion_gaussian_einstein
    implicit none
    integer, parameter :: ntraj = 1000000  ! Trajectories
    integer, parameter :: Tmax = 100
    integer :: i, j, timestep, seed
    real*8 :: dt, diff, eta, position, gaussianvariable, kBT, gamma
    real*8, allocatable :: msd(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.001d0       ! Time step
    diff = 1.0d0       ! Diffusion coefficient (to be verified)
    kBT = 1.0d0        ! Thermal energy
    gamma = kBT / diff ! Friction coefficient (initial guess)
    seed = 1002        ! Random seed

    ! Allocate arrays
    allocate(msd(Tmax))
    msd = 0.0d0

    ! Initialize random number generator
    call sgrnd(seed)

    ! Simulation loop: Gaussian noise without potential
    do j = 1, ntraj
      if (mod(j, 100000) == 0) print *, "Gaussian Trajectory:", j
      position = 0.0d0
      do timestep = 1, Tmax
        call gaussian(gaussianvariable)
        eta = gaussianvariable * dsqrt(2.0d0 * diff / dt)  ! Noise scaled
        position = position + (eta / gamma) * dt
        msd(timestep) = msd(timestep) + position**2
      end do
    end do

    ! Normalize MSD
    msd = msd / dble(ntraj)

    ! Write MSD to file
    open(unit=44, file="msd_gauss_einstein.txt")
    do i = 1, Tmax
      write(44, *) i * dt, msd(i)
    end do
    close(44)

    ! Deallocate arrays
    deallocate(msd)

    print *, "Gaussian simulation complete. Output: msd_gauss_einstein.txt"

  end program Brownian_motion_gaussian_einstein

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

  include 'mt.f90'