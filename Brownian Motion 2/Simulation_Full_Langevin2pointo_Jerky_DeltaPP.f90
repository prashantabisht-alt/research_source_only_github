!======================================================================
!  Simulation_Full_Langevin2pointo_Jerky_DeltaPP (Euler–Maruyama style)
!  Force-level FDT:  <xi_x(t) xi_x(t')> = 2*gamma*kBT * delta''(t-t')
!  Model: m a_dot + gamma a = xi_x(t),  v_dot = a,  x_dot = v
!  Discretization:
!    dW_n ~ N(0, dt)
!    eta_n = dW_n / dt  (white)
!    eta_dot_n ≈ (eta_n - eta_{n-1})/dt = (dW_n - dW_{n-1})/dt^2
!    xi_n = 2*gamma*kBT * eta_dot_n
!
!  Update per particle:
!    a <- a + dt*( - (gamma/m) a + xi_n / m )
!    v <- v + dt * a
!    x <- x + dt * v
!
!  Outputs:
!    jerkpp_stats.txt          : t, means/vars (no theory lines)
!    acceleration_dist_pp.txt  : a_center, P_sim
!    velocity_dist_pp.txt      : v_center, P_sim
!    position_dist_pp.txt      : x_center, P_sim
!======================================================================
program Simulation_Full_Langevin2pointo_Jerky_DeltaPP
  implicit none
  ! ---------------- user params ----------------
  integer,  parameter :: Np = 200000
  real*8,   parameter :: T_max = 50.0d0
  real*8,   parameter :: dt    = 0.001d0   ! smaller dt recommended for delta'' forcing
  integer,  parameter :: save_every = 500

  ! Physics
  real*8, parameter :: m=1.0d0, gamma=0.5d0, kBT=1.0d0

  ! Hist settings
  integer, parameter :: nbins=200
  real*8, parameter :: a_min=-6d0,   a_max=6d0
  real*8, parameter :: v_min=-80d0,  v_max=80d0
  real*8, parameter :: x_min=-2500d0, x_max=2500d0

  ! State
  real*8, dimension(Np) :: x, v, a
  ! Store previous Wiener increment per particle for delta'' construction
  real*8, dimension(Np) :: dW_prev

  ! Histograms
  integer, dimension(nbins) :: hist_a, hist_v, hist_x

  ! Work vars
  integer :: i, n, nsteps, bin, seed
  real*8  :: t
  real*8  :: mean_a, var_a, mean_v, var_v, mean_x, var_x

  ! noise construction
  real*8  :: z, dW, xi
  real*8  :: two_gkBT_over_dt2

  ! histogram bin widths + work
  real*8  :: a_bw, v_bw, x_bw
  real*8  :: ac, vc, xc, Psim

  ! IO
  integer :: fa, fv, fx, fs

  ! RNG
  real*8 :: grnd

  ! ---------------- setup ----------------
  two_gkBT_over_dt2 = 2d0*gamma*kBT/(dt*dt)

  seed=7771
  call sgrnd(seed)

  x=0d0; v=0d0; a=0d0
  dW_prev = 0d0
  nsteps = nint(T_max/dt)

  open(newunit=fs, file="jerkpp_stats.txt")
  write(fs,*) "# t  mean_a var_a  mean_v var_v  mean_x var_x"

  ! ---------------- time loop ----------------
  do n=1, nsteps
    t = dble(n)*dt

    do i=1, Np
      call gaussian(z)          ! z ~ N(0,1)
      dW = sqrt(dt)*z           ! Wiener increment dW_n ~ N(0, dt)

      ! xi_n = 2*gamma*kBT * (dW_n - dW_{n-1}) / dt^2
      xi = two_gkBT_over_dt2 * (dW - dW_prev(i))

      ! a update: m a_dot + gamma a = xi => a_{n+1} = a_n + dt*( -gamma/m a_n + xi/m )
      a(i) = a(i) + dt*( - (gamma/m)*a(i) + xi/m )

      ! chained integrals (use "new" a)
      v(i) = v(i) + dt*a(i)
      x(i) = x(i) + dt*v(i)

      ! store for next step
      dW_prev(i) = dW
    end do

    if (mod(n, save_every) == 0) then
      call moments(a, mean_a, var_a)
      call moments(v, mean_v, var_v)
      call moments(x, mean_x, var_x)
      write(fs,'(7E20.10)') t, mean_a, var_a, mean_v, var_v, mean_x, var_x
    end if
  end do
  close(fs)

  ! ----------- PDFs at final time -----------
  hist_a=0; hist_v=0; hist_x=0
  a_bw=(a_max-a_min)/dble(nbins)
  v_bw=(v_max-v_min)/dble(nbins)
  x_bw=(x_max-x_min)/dble(nbins)

  do i=1, Np
    if (a(i)>=a_min .and. a(i)<a_max) then
      bin = int((a(i)-a_min)/a_bw)+1
      bin = max(1, min(nbins, bin))
      hist_a(bin)=hist_a(bin)+1
    end if
    if (v(i)>=v_min .and. v(i)<v_max) then
      bin = int((v(i)-v_min)/v_bw)+1
      bin = max(1, min(nbins, bin))
      hist_v(bin)=hist_v(bin)+1
    end if
    if (x(i)>=x_min .and. x(i)<x_max) then
      bin = int((x(i)-x_min)/x_bw)+1
      bin = max(1, min(nbins, bin))
      hist_x(bin)=hist_x(bin)+1
    end if
  end do

  open(newunit=fa, file="acceleration_dist_pp.txt")
  write(fa,*) "# a_center   P_sim"
  do bin=1, nbins
    ac  = a_min + (dble(bin)-0.5d0)*a_bw
    Psim = dble(hist_a(bin))/(dble(Np)*a_bw)
    write(fa,'(2E20.10)') ac, Psim
  end do
  close(fa)

  open(newunit=fv, file="velocity_dist_pp.txt")
  write(fv,*) "# v_center   P_sim"
  do bin=1, nbins
    vc  = v_min + (dble(bin)-0.5d0)*v_bw
    Psim = dble(hist_v(bin))/(dble(Np)*v_bw)
    write(fv,'(2E20.10)') vc, Psim
  end do
  close(fv)

  open(newunit=fx, file="position_dist_pp.txt")
  write(fx,*) "# x_center   P_sim"
  do bin=1, nbins
    xc  = x_min + (dble(bin)-0.5d0)*x_bw
    Psim = dble(hist_x(bin))/(dble(Np)*x_bw)
    write(fx,'(2E20.10)') xc, Psim
  end do
  close(fx)

  print *, "Done (delta'' forcing). Files:"
  print *, "  jerkpp_stats.txt"
  print *, "  acceleration_dist_pp.txt"
  print *, "  velocity_dist_pp.txt"
  print *, "  position_dist_pp.txt"

contains

  subroutine moments(arr, mean, var)
    real*8, intent(in)  :: arr(:)
    real*8, intent(out) :: mean, var
    real*8 :: s1, s2
    integer :: L, k
    L = size(arr); s1=0d0; s2=0d0
    do k=1, L
      s1 = s1 + arr(k)
      s2 = s2 + arr(k)*arr(k)
    end do
    mean = s1/dble(L)
    var  = max(0d0, s2/dble(L) - mean*mean)
  end subroutine moments

  ! Gaussian N(0,1) via polar Box–Muller (same style as before)
  subroutine gaussian(s)
    real*8, intent(out) :: s
    real*8 :: u1, u2, r2, fac, grnd
    do
      u1 = 2d0*grnd() - 1d0
      u2 = 2d0*grnd() - 1d0
      r2 = u1*u1 + u2*u2
      if (r2 > 1d-18 .and. r2 < 1d0) exit
    end do
    fac = sqrt(-2d0*log(r2)/r2)
    s = u1*fac
  end subroutine gaussian

end program Simulation_Full_Langevin2pointo_Jerky_DeltaPP

include 'mt.f90'
