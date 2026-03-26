program ABP_Fig1b_PDF
  implicit none

  !-------------------- Physical parameters --------------------
  real*8, parameter :: v0   = 1.0d0
  real*8, parameter :: DR   = 1.0d-1      ! tau_R = 10; choose tsnap << 10

  !-------------------- Snapshot & numerics --------------------
  real*8, parameter :: tsnap = 1.0d0      ! snapshot time (early: tsnap/tau_R = 0.1)
  real*8, parameter :: dt    = 5.0d-3
  integer,  parameter :: ntraj = 200000   ! raise for smoother heatmap

  !-------------------- 2D histogram window & resolution -------
  integer,  parameter :: Nx = 320, Ny = 280
  real*8,  parameter :: cy = 5.0d0        ! y-window half-width in σy units
  real*8,  parameter :: a1_max = 3.0d0    ! x-window: a1=(v0 t - x)/(v0 DR t^2) in [0, a1_max]

  !-------------------- Derived quantities ---------------------
  integer :: nsteps, i, j, ix, iy, fid_xy, fid_mat, seed
  real*8 :: t, phi, x, y, z, sq2DRdt
  real*8 :: sigy, xmin, xmax, ymin, ymax, dx, dy

  !-------------------- 2D histogram (counts) ------------------
  integer, allocatable :: H(:,:)

  !-------------------- RNG init --------------------------------
  seed = 1234567
  call sgrnd(seed)

  !-------------------- Time grid ------------------------------
  nsteps   = nint(tsnap/dt)
  sq2DRdt  = dsqrt(2.0d0*DR*dt)

  !-------------------- Theory-based window --------------------
  ! σy^2 = (2/3) v0^2 DR t^3  (short-time result used for y-window)
  sigy = v0 * dsqrt( (2.0d0/3.0d0) * DR * tsnap**3 )
  xmax = v0 * tsnap
  xmin = xmax - a1_max * (v0*DR*tsnap**2)
  ymin = -cy * sigy
  ymax =  cy * sigy

  dx = (xmax - xmin)/dble(Nx)
  dy = (ymax - ymin)/dble(Ny)

  allocate(H(Ny, Nx)); H = 0

  !-------------------- Open outputs ---------------------------
  open(newunit=fid_xy,  file="xy_t1.00.txt",     status="replace")
  open(newunit=fid_mat, file="hist2d_t1.00.dat", status="replace")

  !-------------------- Monte Carlo for ABP --------------------
  do i = 1, ntraj
    t   = 0.0d0
    phi = 0.0d0          ! φ(0)=0 picks +x as initial direction (breaks x/y symmetry)
    x   = 0.0d0
    y   = 0.0d0

    do j = 1, nsteps
      call gaussian(z)                 ! Z ~ N(0,1)
      phi = phi + sq2DRdt * z          ! dφ = sqrt(2DR dt) Z

      ! wrap φ into (-π, π] for hygiene (cos/sin use φ mod 2π)
      if (phi >  dacos(-1.0d0)) phi = phi - 2.0d0*dacos(-1.0d0)
      if (phi <=-dacos(-1.0d0)) phi = phi + 2.0d0*dacos(-1.0d0)

      x = x + v0*dcos(phi)*dt          ! dx = v0 cosφ dt
      y = y + v0*dsin(phi)*dt          ! dy = v0 sinφ dt
    end do

    ! raw point cloud (for scatter if needed)
    write(fid_xy,'(2(1X,ES20.12))') x, y

    ! 2D bin (count)
    ix = int( (x - xmin)/dx ) + 1
    iy = int( (y - ymin)/dy ) + 1
    if (ix>=1 .and. ix<=Nx .and. iy>=1 .and. iy<=Ny) H(iy,ix) = H(iy,ix) + 1
  end do

  !-------------------- Write normalized 2D PDF ----------------
  call write_pdf2d_matrix(fid_mat, H, Nx, Ny, dx, dy, ntraj)

  close(fid_xy); close(fid_mat)
  deallocate(H)

  print *, "Done. Files:"
  print *, "  xy_t1.00.txt       (raw cloud)"
  print *, "  hist2d_t1.00.dat   (Ny x Nx matrix of P(x,y,t*))"

contains

  !-------------------- N(0,1) via polar Box–Muller -----------
  subroutine gaussian(s)
    implicit none
    real*8 :: s, x1, x2, w, grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
      x1 = 1.0d0 - 2.0d0*grnd()
      x2 = 1.0d0 - 2.0d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0*dlog(w)/w) * x1
  end subroutine gaussian

  !-------------------- Write matrix as probability density ----
  subroutine write_pdf2d_matrix(fid, H, Nx, Ny, dx, dy, ntraj)
    implicit none
    integer, intent(in) :: fid, Nx, Ny, ntraj
    integer, intent(in) :: H(Ny, Nx)
    real*8, intent(in)  :: dx, dy
    integer :: j, k
    real*8 :: density
    do j = 1, Ny
      do k = 1, Nx
        density = dble(H(j,k)) / ( dble(ntraj) * dx * dy )
        write(fid,'(1X,ES16.8)', advance='no') density
      end do
      write(fid,*)
    end do
  end subroutine write_pdf2d_matrix

end program ABP_Fig1b_PDF

!---------------- RNG backend (e.g., Mersenne Twister) --------
include 'mt.f90'
