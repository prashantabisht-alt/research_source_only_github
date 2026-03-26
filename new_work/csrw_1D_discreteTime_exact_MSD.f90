program csrw_1D_discreteTime_exact_MSD
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)

  ! --------- knobs (match your CTRW) ---------
  integer,  parameter :: ntraj      = 50000
  integer,  parameter :: NT         = 120
  logical,  parameter :: use_log_grid = .true.
  real(dp), parameter :: D      = 1.0_dp
  real(dp), parameter :: tmax   = 1.0e4_dp
  integer,  parameter :: seed0  = 2025

  ! --------- state ---------
  integer :: traj, k, seed
  real(dp) :: tprev, tcur, dt, x, z
  real(dp), allocatable :: tg(:), sumx(:), sumx2(:)
  real(dp) :: grnd
  external  :: grnd

  ! --------- init ---------
  seed = seed0
  call sgrnd(seed)

  allocate(tg(NT), sumx(NT), sumx2(NT))
  sumx  = 0.0_dp
  sumx2 = 0.0_dp

  call build_time_grid(tg, tmax, use_log_grid)

  ! --------- ensemble loop ---------
  do traj = 1, ntraj
     x = 0.0_dp
     tprev = 0.0_dp

     do k = 1, NT
        tcur = tg(k)
        dt   = max(0.0_dp, tcur - tprev)
        if (dt > 0.0_dp) then
           call gaussian(z)                           ! N(0,1)
           x = x + sqrt(2.0_dp*D*dt) * z             ! discrete-time Gaussian step
        end if
        sumx(k)  = sumx(k)  + x
        sumx2(k) = sumx2(k) + x*x
        tprev = tcur
     end do
  end do

  call write_msd("msd_ensemble_csrw.txt", tg, sumx, sumx2, ntraj)

  deallocate(tg, sumx, sumx2)

contains
  subroutine build_time_grid(tg, tmax, loggrid)
    real(dp), intent(out) :: tg(:)
    real(dp), intent(in)  :: tmax
    logical, intent(in)   :: loggrid
    integer :: i, NTloc
    real(dp) :: t0, t1, r
    NTloc = size(tg)
    if (loggrid) then
       t0 = max(1.0e-3_dp, tmax/1.0e6_dp)  ! avoid zero on log scale
       t1 = tmax
       r  = (log(t1) - log(t0)) / real(NTloc-1,dp)
       do i=1,NTloc
          tg(i) = exp( log(t0) + real(i-1,dp)*r )
       end do
    else
       do i=1,NTloc
          tg(i) = (real(i-1,dp)/real(NTloc-1,dp)) * tmax
       end do
    end if
  end subroutine build_time_grid

  subroutine gaussian(s)
    real(dp), intent(out) :: s
    real(dp) :: x1,x2,w
    w = 2.0_dp
    do while (w >= 1.0_dp .or. w == 0.0_dp)
       x1 = 1.0_dp - 2.0_dp*grnd()
       x2 = 1.0_dp - 2.0_dp*grnd()
       w  = x1*x1 + x2*x2
    end do
    s = sqrt(-2.0_dp*log(w)/w) * x1
  end subroutine gaussian

  subroutine write_msd(fname, tg, sumx, sumx2, ntraj)
    character(*), intent(in) :: fname
    real(dp),     intent(in) :: tg(:), sumx(:), sumx2(:)
    integer,      intent(in) :: ntraj
    integer :: u,i
    real(dp) :: meanx, msd
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# t    MSD(t)"
    do i=1,size(tg)
       meanx = sumx(i)/real(ntraj,dp)
       msd   = sumx2(i)/real(ntraj,dp) - meanx*meanx
       write(u,'(2es20.10)') tg(i), msd
    end do
    close(u)
  end subroutine write_msd
end program csrw_1D_discreteTime_exact_MSD

! RNG (Mersenne Twister)
include 'mt.f90'
