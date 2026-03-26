program ctrw_single_particle
  implicit none

  integer, parameter :: L = 1000
  integer, parameter :: nsteps = 1000000
  integer :: step, pos, dir, seed
  integer :: x_unwrap

  real*8  :: t, dt, lambda, r1, r2
  real*8  :: grnd
  external :: grnd

  ! Output
  open(unit=42, file="ctrw_single_particle.txt", status="replace")
  open(unit=43, file="waiting_times.txt",       status="replace")

  ! Initialization
  pos    = L/2
  x_unwrap = 0 
  t      = 0.0d0
  lambda = 1.0d0
  seed   = 2025
  call sgrnd(seed)

  ! Simulation loop
  do step = 1, nsteps
     ! Draw waiting time from exponential distribution
     r1 = grnd()
     if (r1 <= 0.0d0) r1 = 1.0d-308   ! safety against log(0)
     dt = -log(r1)/lambda
     t  = t + dt
     write(43,*) dt

     ! Choose jump direction: +1 or -1 with equal probability
     r2 = grnd()
     if (r2 < 0.5d0) then
        dir = 1
     else
        dir = -1
     end if
     x_unwrap = x_unwrap + dir

     ! Update position with periodic boundary conditions
     pos = mod(pos + dir + L, L)

     ! Write time and position to file
     write(42,*) t, pos, x_unwrap
  end do

  close(42)
  close(43)
end program ctrw_single_particle

! RNG (Mersenne Twister) and seeding
include 'mt.f90'

