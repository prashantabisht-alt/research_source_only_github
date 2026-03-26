program Brownian_motion_exponential_einstein
    implicit none
    integer, parameter :: ntraj = 1000000  ! Trajectories
    integer, parameter :: Tmax = 100
    integer :: i, j, timestep, seed
    real*8 :: dt, diff, eta, position, kBT, gamma, lambda
    real*8, allocatable :: msd(:)
    real*8 :: grnd

    ! Parameters
    dt = 0.001d0       ! Time step
    diff = 1.0d0       ! Diffusion coefficient (to be verified)
    kBT = 1.0d0        ! Thermal energy
    gamma = kBT / diff ! Friction coefficient (initial guess)
    lambda = 1.0d0     ! Exponential noise parameter
    seed = 1002        ! Random seed

    ! Allocate arrays
    allocate(msd(Tmax))
    msd = 0.0d0

    ! Initialize random number generator
    call sgrnd(seed)

    ! Simulation loop: Exponential noise without potential
    do j = 1, ntraj
      if (mod(j, 100000) == 0) print *, "Exponential Trajectory:", j
      position = 0.0d0
      do timestep = 1, Tmax
        call exponential(eta, lambda)
        eta = eta * dsqrt(2.0d0 * diff / dt) / dsqrt(2.0d0 / lambda**2)
        position = position + (eta / gamma) * dt
        msd(timestep) = msd(timestep) + position**2
      end do
    end do

    ! Normalize MSD
    msd = msd / dble(ntraj)

    ! Write MSD to file
    open(unit=44, file="msd_exp_einstein.txt")
    do i = 1, Tmax
      write(44, *) i * dt, msd(i)
    end do
    close(44)

    ! Deallocate arrays
    deallocate(msd)

    print *, "Exponential simulation complete. Output: msd_exp_einstein.txt"

  end program Brownian_motion_exponential_einstein

  ! Symmetric exponential (Laplace) generator
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
    end if
  end subroutine exponential

  include 'mt.f90'