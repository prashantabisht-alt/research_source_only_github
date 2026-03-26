!======================================================================
!  Simulation_Full_Langevin2pointo_Jerky_MixedFric_Gaussian
!  m x‴ + gamma_a x″ + gamma_v x′ = ξ(t)
!  <ξ(t) ξ(t')> = 2 gamma_a kBT δ(t-t')
!
!  State: a = x″, v = x′
!  da = [-(gamma_a/m) a  - (gamma_v/m) v] dt + (sqrt(2 gamma_a kBT)/m) dW
!  dv = a dt
!  dx = v dt
!
!  Outputs:
!   mixed_stats_gaussian.txt
!   acceleration_dist_mixed_gaussian.txt
!   velocity_dist_mixed_gaussian.txt
!   position_dist_mixed_gaussian.txt
!======================================================================
program Simulation_Full_Langevin2pointo_Jerky_MixedFric_Gaussian
  implicit none
  integer,  parameter :: Np = 200000
  real*8,   parameter :: T_max = 50.0d0, dt = 0.005d0
  integer,  parameter :: save_every = 100
  real*8,   parameter :: m=1.0d0, gamma_a=0.5d0, gamma_v=0.5d0, kBT=1.0d0

  integer, parameter :: nbins=200
  real*8, parameter :: a_min=-5d0,   a_max=5d0
  real*8, parameter :: v_min=-20d0,  v_max=20d0
  real*8, parameter :: x_min=-200d0, x_max=200d0

  real*8, dimension(Np) :: x, v, a
  integer, dimension(nbins) :: hist_a, hist_v, hist_x

  integer :: i, n, nsteps, bin, seed
  real*8  :: t, PI, a_bw, v_bw, x_bw
  real*8  :: mean_a, var_a, mean_v, var_v, mean_x, var_x
  real*8  :: z, noise
  real*8  :: ac, vc, xc, Psim

  integer :: fa,fv,fx,fs
  real*8 :: grnd

  PI = acos(-1.0d0)
  noise = sqrt(2d0*gamma_a*kBT*dt)/m

  seed=2026
  call sgrnd(seed)

  x=0d0; v=0d0; a=0d0
  nsteps = nint(T_max/dt)

  open(newunit=fs, file="mixed_stats_gaussian.txt")
  write(fs,*) "# t  mean_a var_a  mean_v var_v  mean_x var_x"

  do n=1, nsteps
    t = dble(n)*dt
    do i=1, Np
      call gaussian(z)
      a(i) = a(i) - (gamma_a/m)*a(i)*dt - (gamma_v/m)*v(i)*dt + noise*z
      v(i) = v(i) + a(i)*dt
      x(i) = x(i) + v(i)*dt
    end do

    if (mod(n, save_every) == 0) then
      call moments(a, mean_a, var_a)
      call moments(v, mean_v, var_v)
      call moments(x, mean_x, var_x)
      write(fs,'(7E20.10)') t, mean_a, var_a, mean_v, var_v, mean_x, var_x
    end if
  end do
  close(fs)

  hist_a=0; hist_v=0; hist_x=0
  a_bw=(a_max-a_min)/dble(nbins)
  v_bw=(v_max-v_min)/dble(nbins)
  x_bw=(x_max-x_min)/dble(nbins)

  do i=1, Np
    if (a(i)>=a_min .and. a(i)<a_max) then
      bin = int((a(i)-a_min)/a_bw)+1
      hist_a(bin)=hist_a(bin)+1
    end if
    if (v(i)>=v_min .and. v(i)<v_max) then
      bin = int((v(i)-v_min)/v_bw)+1
      hist_v(bin)=hist_v(bin)+1
    end if
    if (x(i)>=x_min .and. x(i)<x_max) then
      bin = int((x(i)-x_min)/x_bw)+1
      hist_x(bin)=hist_x(bin)+1
    end if
  end do

  open(newunit=fa, file="acceleration_dist_mixed_gaussian.txt")
  write(fa,*) "# a_center   P_sim"
  do bin=1, nbins
    ac  = a_min + (dble(bin)-0.5d0)*a_bw
    Psim = dble(hist_a(bin))/(dble(Np)*a_bw)
    write(fa,'(2E20.10)') ac, Psim
  end do
  close(fa)

  open(newunit=fv, file="velocity_dist_mixed_gaussian.txt")
  write(fv,*) "# v_center   P_sim"
  do bin=1, nbins
    vc  = v_min + (dble(bin)-0.5d0)*v_bw
    Psim = dble(hist_v(bin))/(dble(Np)*v_bw)
    write(fv,'(2E20.10)') vc, Psim
  end do
  close(fv)

  open(newunit=fx, file="position_dist_mixed_gaussian.txt")
  write(fx,*) "# x_center   P_sim"
  do bin=1, nbins
    xc  = x_min + (dble(bin)-0.5d0)*x_bw
    Psim = dble(hist_x(bin))/(dble(Np)*x_bw)
    write(fx,'(2E20.10)') xc, Psim
  end do
  close(fx)

  print *, "Done (Gaussian noise)."

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

end program Simulation_Full_Langevin2pointo_Jerky_MixedFric_Gaussian

include 'mt.f90'
