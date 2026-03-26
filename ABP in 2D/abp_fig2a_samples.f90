!============================================================================
! ABP Fig. 2(a) marginals + analytic theory curve f_x(a1) (Eq. 11)
!  - Simulates ABP at t ∈ {0.25, 0.50, 1.00} with v0=1, DR=0.1 (τR=10)
!  - Writes PDFs:
!      * y_scaled_tXX.txt : PDF of y/σy  (σy^2 = (2/3) v0^2 DR t^3)
!      * a1_edge_tXX.txt  : PDF of a1 = (v0 t - x)/(v0 DR t^2)  (≥ 0)
!  - Writes analytic theory: fx_theory.txt for overlay of f_x(a1) (Eq. 11)
!
! Files are normalized densities (integrate ≈ 1 over their domain).
! Requires: your RNG backend (mt.f90) providing sgrnd(seed), grnd().
!================================
============================================
program abp_fig2a_with_theory
  implicit none
  ! ---------------- Physics ----------------
  integer, parameter :: dp = selected_real_kind(15,300)
  real(dp), parameter :: v0 = 1.0_dp
  real(dp), parameter :: DR = 1.0e-1_dp      ! tau_R = 10

  ! ---------------- Simulation controls ----------------
  integer,  parameter :: Ntraj = 300000      ! per snapshot (raise for smoother)
  real(dp), parameter :: dt    = 1.0e-3_dp
  integer                      :: nsteps

  ! Snapshots
  integer,  parameter :: Ns = 3
  real(dp), parameter :: ts(Ns) = (/ 0.25_dp, 0.50_dp, 1.00_dp /)

  ! 1D hist settings
  integer,  parameter :: Nb_a1 = 300, Nb_y = 300
  real(dp), parameter :: a1_max = 4.0_dp       ! x-edge left tail coverage
  real(dp), parameter :: yscaled_span = 4.0_dp ! y/σ_y in [-4,4]

  ! Theory curve settings for f_x(a1) (Eq. 11)
  integer,  parameter :: KMAX = 300        ! terms in the series (converges fast)
  integer,  parameter :: Na   = 1000      ! number of a1 points
  real(dp), parameter :: a1_min_th = 1.0e-3_dp
  real(dp), parameter :: a1_max_th = 4.0_dp

  ! ---------------- Work arrays ----------------
  integer :: s, i, j, bina1, biny, uid
  character(len=32) :: fa1, fy
  real(dp) :: tstar, sigy, sq2DRdt, phi, x, y, z, a1, yscaled
  real(dp) :: binw_a1, binw_y, center, integ
  integer, allocatable :: Ha1(:), Hy(:)

  ! RNG
  integer :: seed
  seed = 9876543
  call sgrnd(seed)

  allocate(Ha1(Nb_a1), Hy(Nb_y))

  !==================== main snapshots loop ====================
  do s = 1, Ns
    tstar   = ts(s)
    nsteps  = nint(tstar/dt)
    sq2DRdt = sqrt(2.0_dp*DR*dt)
    sigy    = v0*sqrt((2.0_dp/3.0_dp)*DR*tstar**3)

    Ha1 = 0; Hy = 0
    binw_a1 = a1_max / real(Nb_a1,dp)
    binw_y  = (2.0_dp*yscaled_span) / real(Nb_y,dp)

    do i = 1, Ntraj
      ! ICs
      phi = 0.0_dp; x = 0.0_dp; y = 0.0_dp

      ! evolve to tstar
      do j = 1, nsteps
        call gaussian(z)
        phi = phi + sq2DRdt * z
        
        x = x + v0*cos(phi)*dt
        y = y + v0*sin(phi)*dt
      end do

      ! ---- y Gaussian scaling ----
      yscaled = y / sigy
      biny = int( (yscaled + yscaled_span) / binw_y ) + 1
      if (biny>=1 .and. biny<=Nb_y) Hy(biny) = Hy(biny) + 1

      ! ---- x edge scaling ----
      a1 = (v0*tstar - x) / (v0*DR*tstar**2)
      if (a1 >= 0.0_dp) then
        bina1 = int( a1 / binw_a1 ) + 1
        if (bina1>=1 .and. bina1<=Nb_a1) Ha1(bina1) = Ha1(bina1) + 1
      end if
    end do

    ! ---- write y_scaled density ----
    write(fy,'("y_scaled_t",F4.2,".txt")') tstar
    open(newunit=uid, file=trim(adjustl(fy)), status="replace")
    integ = 0.0_dp
    do j = 1, Nb_y
      center = -yscaled_span + ( (real(j,dp)-0.5_dp)*binw_y )
      call write_pair(uid, center, real(Hy(j),dp)/( real(Ntraj,dp)*binw_y ))
      integ = integ + ( real(Hy(j),dp)/( real(Ntraj,dp) ) )
    end do
    close(uid)
    print *, "Wrote ", trim(adjustl(fy)), " (sum≈", integ, ")"

    ! ---- write a1 density ----
    write(fa1,'("a1_edge_t",F4.2,".txt")') tstar
    open(newunit=uid, file=trim(adjustl(fa1)), status="replace")
    integ = 0.0_dp
    do j = 1, Nb_a1
      center = (real(j,dp)-0.5_dp) * binw_a1
      call write_pair(uid, center, real(Ha1(j),dp)/( real(Ntraj,dp)*binw_a1 ))
      integ = integ + ( real(Ha1(j),dp)/( real(Ntraj,dp) ) )
    end do
    close(uid)
    print *, "Wrote ", trim(adjustl(fa1)), " (sum≈", integ, ")"
  end do
  deallocate(Ha1, Hy)

  !==================== analytic theory curve f_x(a1) ====================
  call write_fx_theory("fx_theory.txt", a1_min_th, a1_max_th, Na, KMAX)

  print *, "All done: y_scaled_*, a1_edge_*, fx_theory.txt"

contains

  !---------------------- write helper ----------------------
  subroutine write_pair(uid, x, p)
    implicit none
    integer, intent(in) :: uid
    real(dp), intent(in) :: x, p
    write(uid,'(2(1X,ES20.12))') x, p
  end subroutine write_pair

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

  !---------------------- Theory f_x(a1) writer (Eq. 11) ----------------------
  subroutine write_fx_theory(fname, a1min, a1max, Na, KMAX)
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: Na, KMAX
    real(dp), intent(in) :: a1min, a1max
    integer :: i, uid
    real(dp) :: a1, da, fx
    da = (a1max - a1min)/real(Na,dp)
    open(newunit=uid, file=fname, status="replace")
    do i = 0, Na
      a1 = a1min + real(i,dp)*da
      if (a1 <= 0.0_dp) then
        fx = 0.0_dp
      else
        fx = fx_series_Eq11(a1, KMAX)   ! series sum (without prefactor)
        fx = fx / ( 2.0_dp*sqrt(acos(-1.0_dp)*a1**3) )   ! 1/(2√(π a1^3))
      end if
      write(uid,'(2(1X,ES20.12))') a1, fx
    end do
    close(uid)
    print *, "Wrote fx_theory.txt (Eq. 11)"
  end subroutine write_fx_theory

  !---------------------- Eq. (11) series core (without prefactor) ----------------------
  ! fx(a1) = 1/(2√(π a1^3)) * Σ_{k=0..∞} [ (-1)^k/(4k+1) * ( (2k choose k)/2^{2k} ) * exp(-(4k+1)^2/(8 a1)) ]
  function fx_series_Eq11(a1, KMAX) result(sumfx)
    implicit none
    real(dp), intent(in) :: a1
    integer,  intent(in) :: KMAX
    integer :: k
    real(dp) :: sumfx, term, cbin
    sumfx = 0.0_dp
    do k = 0, KMAX
      cbin = central_binom_over_2pow2k(k)     ! (2k choose k) / 2^{2k}
      term = ((-1.0_dp)**k) / real(4*k+1,dp) * cbin * exp( - real(4*k+1,dp)**2 / (8.0_dp*a1) )
      sumfx = sumfx + term
    end do
  end function fx_series_Eq11

  !---------------------- (2k choose k) / 2^{2k}, CORRECT product ----------------------
  function central_binom_over_2pow2k(k) result(val)
    implicit none
    integer, intent(in) :: k
    integer :: m
    real(dp) :: val
    ! (2k choose k) / 2^{2k} = Π_{m=1..k} (2m-1)/(2m)  with value 1 at k=0
    val = 1.0_dp
    do m = 1, k
      val = val * real(2*m-1,dp) / real(2*m,dp)
    end do
  end function central_binom_over_2pow2k

end program abp_fig2a_with_theory

!---------------- RNG backend (Mersenne Twister etc.) ----------------
include 'mt.f90'

