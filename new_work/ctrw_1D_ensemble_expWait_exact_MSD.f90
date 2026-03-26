program ctrw_1D_ensemble_expWait_exact_MSD
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)

  ! --------- knobs ---------
  integer,  parameter :: ntraj      = 50000      ! number of trajectories
  integer,  parameter :: nsteps_max = 1000000    ! safety guard
  integer,  parameter :: NT         = 120        ! # time-grid points
  logical,  parameter :: use_log_grid = .true.   ! log or linear grid
  real(dp), parameter :: lambda = 1.0_dp         ! Exp rate
  real(dp), parameter :: D      = 1.0_dp         ! diffusion coefficient
  real(dp), parameter :: tmax   = 1.0e4_dp       ! max time for MSD
  integer,  parameter :: seed0  = 2025           ! RNG seed

  ! --------- state ---------
  integer :: traj, step, k, seed
  real(dp) :: t, dt_cand, rem, tk, piece, x, z, r1
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
     t = 0.0_dp; x = 0.0_dp; k = 1
     if (.not.use_log_grid) then
        if (tg(1) == 0.0_dp) then
           sumx(1)  = sumx(1)  + x
           sumx2(1) = sumx2(1) + x*x
           k = 2
        end if
     end if

     do step = 1, nsteps_max
        if (k > NT) exit

        ! propose next waiting time (inverse CDF for exponential)
        r1 = grnd(); if (r1 <= 0.0_dp) r1 = 1.0e-308_dp
        dt_cand = -log(r1)/lambda
        rem = dt_cand

        ! consume this candidate, splitting at grid times as needed
        do
           if (k > NT) exit
           tk = tg(k)

           if (t + rem < tk) then
              call gaussian(z)
              x = x + sqrt(2.0_dp*D*rem)*z
              t = t + rem
              exit
           else
              piece = tk - t
              if (piece > 0.0_dp) then
                 call gaussian(z)
                 x = x + sqrt(2.0_dp*D*piece)*z
                 t = tk
              else
                 t = tk
              end if

              ! record snapshot at tg(k)
              sumx(k)  = sumx(k)  + x
              sumx2(k) = sumx2(k) + x*x
              k = k + 1

              rem = rem - piece
              if (rem <= 0.0_dp) exit
           end if
        end do
     end do
  end do

  ! --------- write MSD(t) ---------
  call write_msd("msd_ensemble.txt", tg, sumx, sumx2, ntraj)

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
end program ctrw_1D_ensemble_expWait_exact_MSD

! RNG (Mersenne Twister)
include 'mt.f90'
