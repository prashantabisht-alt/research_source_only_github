!===========================================================
! Dump two single-trajectory paths for ABP in 2D
!   Short-time  (Tshort << 1/DR)
!   Long-time   (Tlong  >> 1/DR)
! Files: traj_short.txt  (n, x, y)
!        traj_long.txt   (n, x, y)
! Compile: gfortran -O3 abp_traj_dump.f90 -o abp_traj
!          (or: gfortran -O3 abp_traj_dump.f90 mt.f90 -o abp_traj)
!===========================================================
program abp_traj_dump
  implicit none
  real*8, parameter :: v0   = 1.0d0
  real*8, parameter :: DR   = 1.0d-2       ! tau_R = 100
  real*8, parameter :: dt   = 2.0d-2

  real*8, parameter :: Tshort = 5.0d0      ! 0.05 * tau_R  (t << 1/DR)
  real*8, parameter :: Tlong  = 1000.0d0   ! 10 * tau_R    (t >> 1/DR)

  integer :: nS, nL, n, stride_long, uS, uL
  real*8  :: x, y, phi, g, sq2DRdt, pi

  external grnd
  integer :: seed
  seed = 20250903
  call sgrnd(seed)

  pi = dacos(-1.0d0)
  sq2DRdt = dsqrt(2.0d0*DR*dt)

  ! ---- SHORT ----
  nS = nint(Tshort/dt)
  x=0d0; y=0d0; phi=0d0
  open(newunit=uS, file="traj_short.txt", status="replace")
  write(uS,'(I10,1X,2(ES20.12,1X))') 0, x, y
  do n=1,nS
     call gaussian(g)
     phi = phi + sq2DRdt*g
     if (phi >  pi) phi = phi - 2d0*pi
     if (phi <= -pi) phi = phi + 2d0*pi
     x = x + v0*dcos(phi)*dt
     y = y + v0*dsin(phi)*dt
     write(uS,'(I10,1X,2(ES20.12,1X))') n, x, y
  end do
  close(uS)

  ! ---- LONG (downsampled writes) ----
  nL = nint(Tlong/dt)
  stride_long = 10                       ! write every 10th step
  x=0d0; y=0d0; phi=0d0
  open(newunit=uL, file="traj_long.txt", status="replace")
  write(uL,'(I10,1X,2(ES20.12,1X))') 0, x, y
  do n=1,nL
     call gaussian(g)
     phi = phi + sq2DRdt*g
     if (phi >  pi) phi = phi - 2d0*pi
     if (phi <= -pi) phi = phi + 2d0*pi
     x = x + v0*dcos(phi)*dt
     y = y + v0*dsin(phi)*dt
     if (mod(n,stride_long)==0) write(uL,'(I10,1X,2(ES20.12,1X))') n, x, y
  end do
  close(uL)

  print *,"Wrote: traj_short.txt , traj_long.txt"

contains
  subroutine gaussian(s)
    real*8 :: s, x1, x2, w, grnd
    w = 2d0
    do while (w > 1d0 .or. w == 0d0)
      x1 = 1d0 - 2d0*grnd()
      x2 = 1d0 - 2d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    s = dsqrt(-2d0*dlog(w)/w) * x1
  end subroutine gaussian
end program abp_traj_dump

! If you keep your RNG in a separate file, include it:
include 'mt.f90'
