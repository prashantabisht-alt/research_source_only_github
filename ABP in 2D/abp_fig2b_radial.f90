!============================================================================
! ABP Fig. 2(b): Short-time radial scaling + THEORY (Brownian-functional MC)
!
! Scaling (paper, Eqs. (12)-(13)):
!   z = (v0^2 t^2 - r^2) / (2 DR v0^2 t^3)  >= 0,  r^2 = x^2 + y^2
!   P(r^2,t) = [1/(2 DR v0^2 t^3)] * f_r(z),   with
!   Laplace{ f_r }(λ) = [ sqrt(2λ) / sinh(sqrt(2λ)) ]^{1/2}
!   and z = a1 - a2^2 with a1=∫_0^1 B(s)^2 ds, a2=∫_0^1 B(s) ds.
!
! This code:
!   1) Simulates ABP at t ∈ {0.25, 0.50, 1.00} (v0=1, DR=0.1) and writes
!         fr_tXX.txt : columns (z_center, f_r(z))   [normalized density].
!   2) Computes the universal THEORY curve f_r(z) by Monte Carlo sampling
!      of Brownian paths on [0,1], outputs:
!         fr_theory.txt : columns (z, f_r(z))       [normalized density].
!
! Notes:
!   - For short-time validity, we DO NOT wrap φ; φ follows unbounded BM.
!   - Uses your RNG backend mt.f90 providing  sgrnd(seed), grnd().
!============================================================================
program abp_fig2b_with_theory
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)

  ! -------- Physics (ABP) --------
  real(dp), parameter :: v0 = 1.0_dp
  real(dp), parameter :: DR = 1.0e-1_dp     ! τR = 10

  ! -------- Simulation controls --------
  integer,  parameter :: Ntraj = 300000     ! per snapshot
  real(dp), parameter :: dt    = 5.0e-3_dp
  integer                       :: nsteps

  ! Snapshots (short-time)
  integer,  parameter :: Ns = 3
  real(dp), parameter :: ts(Ns) = (/ 0.25_dp, 0.50_dp, 1.00_dp /)

  ! Histogram for z (sim and theory)
  integer,  parameter :: Nz = 360
  real(dp), parameter :: z_max = 4.0_dp
  real(dp)            :: dz

  ! -------- THEORY (Brownian functional MC) --------
  integer,  parameter :: N_FR = 600000      ! # of Brownian paths
  integer,  parameter :: M_FR = 600         ! time steps on [0,1] for B(s)
  real(dp), parameter :: z_max_th = 4.0_dp
  integer,  parameter :: Nz_th = 360

  ! Work vars
  integer :: s, i, j, bin, uid, seed
  character(len=32) :: fout
  real(dp) :: tstar, phi, x, y, zeta, sq2DRdt
  real(dp) :: r2, v0t, vr2, z
  integer, allocatable :: Hz(:), Hth(:)

  ! RNG
  seed = 24681357
  call sgrnd(seed)

  !------------------- simulate & bin f_r(z) -------------------
  allocate(Hz(Nz))
  dz = z_max / real(Nz,dp)

  do s = 1, Ns
    tstar   = ts(s)
    nsteps  = nint(tstar/dt)
    sq2DRdt = sqrt(2.0_dp*DR*dt)
    Hz      = 0

    do i = 1, Ntraj
      ! ICs
      phi = 0.0_dp
      x   = 0.0_dp
      y   = 0.0_dp

      ! Evolve: unwrapped φ (Gaussian increments), midpoint for pos
      do j = 1, nsteps
        call gaussian(zeta)                     ! zeta ~ N(0,1)
        phi = phi + sq2DRdt * zeta
        x = x + v0*cos(phi - 0.5_dp*sq2DRdt*zeta)*dt
        y = y + v0*sin(phi - 0.5_dp*sq2DRdt*zeta)*dt
      end do

      ! Radial scaling variable z >= 0
      r2  = x*x + y*y
      v0t = v0*tstar
      vr2 = v0t*v0t - r2
      z   = vr2 / ( 2.0_dp*DR*v0*v0*tstar*tstar*tstar )

      if (z >= 0.0_dp) then
        bin = int( z / dz ) + 1
        if (bin>=1 .and. bin<=Nz) Hz(bin) = Hz(bin) + 1
      end if
    end do

    ! Dump normalized density f_r(z) from sim
    write(fout,'("fr_t",F4.2,".txt")') tstar
    call dump_histogram(trim(adjustl(fout)), Hz, Nz, 0.0_dp, z_max, Ntraj)
  end do
  deallocate(Hz)

  !------------------- THEORY: MC of Brownian functional -------------------
  allocate(Hth(Nz_th))
  Hth = 0
  call theory_fr_mc(Hth, Nz_th, z_max_th, N_FR, M_FR)

  call dump_histogram("fr_theory.txt", Hth, Nz_th, 0.0_dp, z_max_th, N_FR)
  deallocate(Hth)

  print *, "Done:  fr_t*.txt (sim)   and   fr_theory.txt (short-time theory)."

contains
  !---------------------- Gaussian N(0,1) via polar Box–Muller ----------------------
  subroutine gaussian(s)
    implicit none
    real(dp) :: s, x1, x2, w, grnd
    w = 2.0_dp
    do while (w > 1.0_dp .or. w == 0.0_dp)
      x1 = 1.0_dp - 2.0_dp*grnd()
      x2 = 1.0_dp - 2.0_dp*grnd()
      w  = x1*x1 + x2*x2
    end do
    s = sqrt(-2.0_dp*log(w)/w) * x1
  end subroutine gaussian

  !---------------------- Dump histogram as normalized density ----------------------
  subroutine dump_histogram(fname, H, nbin, xmin, xmax, Nsamples)
    implicit none
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: H(:), nbin, Nsamples
    real(dp),         intent(in) :: xmin, xmax
    integer :: j, uid
    real(dp) :: binw, xcenter
    binw = (xmax - xmin) / real(nbin,dp)
    open(newunit=uid, file=fname, status="replace")
    do j = 1, nbin
      xcenter = xmin + (real(j,dp)-0.5_dp)*binw
      write(uid,'(2(1X,ES20.12))') xcenter, real(H(j),dp)/( real(Nsamples,dp)*binw )
    end do
    close(uid)
    print *, "Wrote ", trim(fname)
  end subroutine dump_histogram

  !---------------------- THEORY: MC for f_r(z) via Brownian functional -------------
  ! Draw B(s) on [0,1] with M steps; compute z = ∫B^2 ds  -  (∫B ds)^2
  ! Repeat N_FR times; bin into [0, z_max_th].
  subroutine theory_fr_mc(Hth, Nz, zmax, N_FR, M)
    implicit none
    integer,  intent(in)    :: Nz, N_FR, M
    real(dp), intent(in)    :: zmax
    integer,  intent(inout) :: Hth(Nz)

    real(dp) :: du, B, dB, a1, a2, z, binw
    integer  :: n, k, bin
    binw = zmax / real(Nz,dp)
    du   = 1.0_dp / real(M,dp)

    do n = 1, N_FR
      B  = 0.0_dp
      a1 = 0.0_dp
      a2 = 0.0_dp
      do k = 1, M
        dB = sqrt(du) * randn_std()
        B  = B + dB
        a1 = a1 + B*B*du
        a2 = a2 + B*du
      end do
      z = a1 - a2*a2
      if (z >= 0.0_dp .and. z < zmax) then
        bin = int( z / binw ) + 1
        if (bin < 1) bin = 1
        if (bin > Nz) bin = Nz
        Hth(bin) = Hth(bin) + 1
      end if
    end do
  end subroutine theory_fr_mc

  !---------------------- N(0,1) using same RNG core ----------------------
  real(dp) function randn_std()
    implicit none
    real(dp) :: u1, u2
    call random_number(u1)  ! uses the same seed state as mt.f90 through random_seed?
    call random_number(u2)
    if (u1 <= 1.0e-16_dp) u1 = 1.0e-16_dp
    randn_std = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn_std

end program abp_fig2b_with_theory

! RNG backend
include 'mt.f90'

