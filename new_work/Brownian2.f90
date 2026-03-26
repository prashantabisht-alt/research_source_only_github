program Brownian_motion_ensemble_CSRW
  implicit none

  !------------ User params -------------
  integer, parameter :: NumTraj = 10000
  integer, parameter :: Tmax    = 100
  real*8,  parameter :: dt      = 1.0d-3
  real*8,  parameter :: D       = 1.0d0
  integer,  parameter :: Seed    = 1002
  !--------------------------------------

  integer :: j, timestep
  real*8  :: g, eta, t
  real*8, allocatable :: x(:)
  real*8  :: msd
  character(len=*), parameter :: traj_file = "Brownian_trajectories.txt"
  character(len=*), parameter :: msd_file  = "Brownian_msd.txt"

  call sgrnd(Seed)

  allocate(x(NumTraj))
  x(:) = 0.0d0

  open(unit=10, file=traj_file, status="replace")  ! columns: t, traj_id, x
  open(unit=11, file=msd_file,  status="replace")  ! columns: t, <x^2>

  do timestep = 1, Tmax
    t = timestep * dt

    ! Update all trajectories (Euler–Maruyama)
    do j = 1, NumTraj
      call gaussian(g)                             ! N(0,1)
      eta = g * dsqrt(2.0d0 * D / dt)              ! noise acceleration term
      x(j) = x(j) + eta*dt                         ! = x + g*sqrt(2 D dt)
      write(10,*) t, j, x(j)
    end do

    ! Ensemble MSD at this time
    msd = sum(x*x) / real(NumTraj, kind=8)
    write(11,*) t, msd

    ! light progress
    if (mod(timestep, 10000) == 0) print '(A, F7.3)', "t = ", t
  end do

  close(10); close(11)
  deallocate(x)

  print *, "Wrote:", trim(traj_file), "and", trim(msd_file)
end program Brownian_motion_ensemble_CSRW


!=============================================================================================
! Mersenne Twister RNG (provide mt.f90 in the same folder)
include 'mt.f90'
!=============================================================================================
! Unit Gaussian via polar Box–Muller: N(0,1)
subroutine gaussian(s)
  implicit none
  real*8, intent(out) :: s
  real*8 :: x1, x2, w, grnd

  w = 2.0d0
  do while (w > 1.0d0 .or. w == 0.0d0)
    x1 = 1.0d0 - 2.0d0*grnd()
    x2 = 1.0d0 - 2.0d0*grnd()
    w  = x1*x1 + x2*x2
  end do
  s = dsqrt(-2.0d0 * dlog(w) / w) * x1
end subroutine gaussian
