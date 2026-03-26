!======================================================================
!  Simulation_Full_Langevin2pointo_Jerky (Euler–Maruyama)
!  m a_dot + gamma a = xi(t),   <xi xi> = 2 gamma kBT delta
!  a = x_ddot, v = x_dot
!
!  Outputs:
!   jerk_stats.txt          : t, means/vars + theory
!   acceleration_dist.txt   : a_center, P_sim, P_theory
!   velocity_dist.txt       : v_center, P_sim, P_theory
!   position_dist.txt       : x_center, P_sim, P_theory
!======================================================================
program Simulation_Full_Langevin2pointo_Jerky
  implicit none
  integer,  parameter :: ntraj = 200000
  real*8,   parameter :: T_max = 50.0d0, dt = 0.005d0
  integer,  parameter :: save_every = 100
  real*8,   parameter :: m = 1.0d0, gamma = 0.5d0, kBT = 1.0d0

  ! Histogram settings
  integer, parameter :: nbins = 200
  real*8, parameter :: a_min = -5d0,   a_max = 5d0
  real*8, parameter :: v_min = -60d0,  v_max = 60d0
  real*8, parameter :: x_min = -1600d0, x_max = 1600d0

  ! State arrays
  real*8, dimension(ntraj) :: x, v, a
  integer, dimension(nbins) :: hist_a, hist_v, hist_x

  ! Time + physics vars
  integer :: i, n, nsteps, bin, seed
  real*8 :: tau_a, noise
  real*8 :: mean_a, var_a, mean_v, var_v, mean_x, var_x
  real*8 :: th_var_a, th_var_v, th_var_x
  real*8 :: t, PI

  ! Histogram bin widths + work vars
  real*8 :: a_bw, v_bw, x_bw
  real*8 :: ac, vc, xc, Psim, Pth

  ! RNG vars
  real*8 :: z1, grnd
  integer :: fa, fv, fx, fs

  ! Setup constants
  PI     = acos(-1.0d0)
  tau_a  = m / gamma
  noise  = sqrt(2d0 * gamma * kBT * dt) / m
  seed   = 901
  call sgrnd(seed)

  ! Init
  x = 0d0; v = 0d0; a = 0d0
  nsteps = nint(T_max / dt)

  open(newunit=fs, file="jerk_stats.txt")
  write(fs,*) "# t  mean_a var_a mean_v var_v mean_x var_x th_var_a th_var_v th_var_x"

  ! ---------------- Time loop ----------------
  do n = 1, nsteps
    t = dble(n) * dt

    ! Euler–Maruyama update
    do i = 1, ntraj
      call gaussian(z1)
      a(i) = a(i) - (gamma/m) * a(i) * dt + noise * z1
      v(i) = v(i) + a(i) * dt
      x(i) = x(i) + v(i) * dt
    end do

    ! Save stats
    if (mod(n, save_every) == 0) then
      call moments(a, mean_a, var_a)
      call moments(v, mean_v, var_v)
      call moments(x, mean_x, var_x)

      ! Theory (OU acceleration → integrated for v, x)
      th_var_a = (kBT/m) * (1d0 - exp(-2d0 * t / tau_a))
      th_var_v = (2d0 * kBT / m) * (tau_a * t - tau_a**2 * (1d0 - exp(-t/tau_a)))
      th_var_x = (tau_a * kBT / (3d0*m)) * (2d0*t**3 - 3d0*tau_a*t**2 + 6d0*tau_a**3 &
                   - 6d0*tau_a**2 * (t+tau_a) * exp(-t/tau_a))

      write(fs,'(10E20.10)') t, mean_a, var_a, mean_v, var_v, mean_x, var_x, &
                             th_var_a, th_var_v, th_var_x
    end if
  end do
  close(fs)

  ! ---------------- Histograms ----------------
  hist_a = 0; hist_v = 0; hist_x = 0
  a_bw = (a_max - a_min) / dble(nbins)
  v_bw = (v_max - v_min) / dble(nbins)
  x_bw = (x_max - x_min) / dble(nbins)

  do i = 1, ntraj
    if (a(i) >= a_min .and. a(i) < a_max) then
      bin = int((a(i) - a_min) / a_bw) + 1
      hist_a(bin) = hist_a(bin) + 1
    end if
    if (v(i) >= v_min .and. v(i) < v_max) then
      bin = int((v(i) - v_min) / v_bw) + 1
      hist_v(bin) = hist_v(bin) + 1
    end if
    if (x(i) >= x_min .and. x(i) < x_max) then
      bin = int((x(i) - x_min) / x_bw) + 1
      hist_x(bin) = hist_x(bin) + 1
    end if
  end do

  ! --- Acceleration PDF ---
  open(newunit=fa, file="acceleration_dist.txt")
  write(fa,*) "# a_center  P_sim  P_theory"
  do bin = 1, nbins
    ac   = a_min + (dble(bin) - 0.5d0) * a_bw
    Psim = dble(hist_a(bin)) / (dble(ntraj) * a_bw)
    Pth  = sqrt(m/(2d0*PI*kBT)) * exp(-m*ac*ac/(2d0*kBT))
    write(fa,'(3E20.10)') ac, Psim, Pth
  end do
  close(fa)

  ! --- Velocity PDF ---
  open(newunit=fv, file="velocity_dist.txt")
  write(fv,*) "# v_center  P_sim  P_theory"
  do bin = 1, nbins
    vc   = v_min + (dble(bin) - 0.5d0) * v_bw
    Psim = dble(hist_v(bin)) / (dble(ntraj) * v_bw)
    ! Gaussian with variance = th_var_v at final time
    Pth  = sqrt(1d0/(2d0*PI*th_var_v)) * exp(-vc*vc/(2d0*th_var_v))
    write(fv,'(3E20.10)') vc, Psim, Pth
  end do
  close(fv)

  ! --- Position PDF ---
  open(newunit=fx, file="position_dist.txt")
  write(fx,*) "# x_center  P_sim  P_theory"
  do bin = 1, nbins
    xc   = x_min + (dble(bin) - 0.5d0) * x_bw
    Psim = dble(hist_x(bin)) / (dble(ntraj) * x_bw)
    Pth  = sqrt(1d0/(2d0*PI*th_var_x)) * exp(-xc*xc/(2d0*th_var_x))
    write(fx,'(3E20.10)') xc, Psim, Pth
  end do
  close(fx)

  print *, "Done: jerk_stats.txt, acceleration_dist.txt, velocity_dist.txt, position_dist.txt"

contains

  subroutine moments(arr, mean, var)
    real*8, intent(in)  :: arr(:)
    real*8, intent(out) :: mean, var
    real*8 :: s1, s2
    integer :: L, k
    L = size(arr)
    s1 = 0d0; s2 = 0d0
    do k = 1, L
      s1 = s1 + arr(k)
      s2 = s2 + arr(k)*arr(k)
    end do
    mean = s1/dble(L)
    var  = max(0d0, s2/dble(L) - mean*mean)
  end subroutine moments

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
    s   = u1*fac
  end subroutine gaussian

end program Simulation_Full_Langevin2pointo_Jerky

include 'mt.f90'
