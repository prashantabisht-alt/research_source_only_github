program Brownian_motion_single_CSRW
  implicit none

  integer :: timestep, Tmax, seed
  real*8  :: g, eta, D, dt, position

  ! RNG seed
  seed = 1002
  call sgrnd(seed)

  ! Parameters
  dt   = 1.0d-3
  D    = 1.0d0
  Tmax = 100000

  position = 0.0d0

  open(unit=42, file="Brownian_single.txt", status="replace")

  do timestep = 1, Tmax
     call gaussian(g)                          ! Gaussian random N(0,1)
     eta = g * dsqrt(2.0d0 * D / dt)           ! noise term
     position = position + eta*dt              ! update position
     write(42,*) timestep*dt, position
  end do

  close(42)

  print *, "Trajectory saved to Brownian_single.txt"
end program Brownian_motion_single_CSRW

include 'mt.f90'

!=============================================================================================
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
