!===============================================================================
! 1D Discrete-Space, Discrete-Time Random Walk (simplified ensemble)
! - Steps: +1 with prob p, -1 with prob (1-p)
! - Outputs (consistent with earlier DSDT code):
!   * traj1.txt                        : n, X_n (trajectory #1)
!   * msd_dsdt_pXX_nstepsNNNN.txt      : MSD vs n (MC + exact)
!   * xCounts_dsdt_pXX_nstepsNNNN.txt  : raw counts at final time
!   * xPDF_dsdt_pXX_nstepsNNNN.txt     : normalized pmf at final time (MC)
!   * xPDF_exact_pXX_nstepsNNNN.txt    : exact pmf at final time
!===============================================================================
program rw_1D_DSDT_simple_named
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)
  ! ---- knobs ----
  integer,  parameter :: ntraj  = 100000
  integer,  parameter :: nsteps = 10000
  real(dp), parameter :: p      = 0.50_dp     ! prob(step=+1)
  integer,  parameter :: seed0  = 2025
  real(dp), parameter :: ksigma = 6.0_dp      ! histogram half-width in sigma

  ! ---- state ----
  integer :: i, n, x, dir, seed
  real(dp) :: r
  real(dp), allocatable :: sumx(:), sumx2(:)
  integer :: xmin, xmax, nbins, idx
  integer, allocatable :: hcount(:)

  real(dp) :: mu, sigma, meanx, msd, meanx_th, var_th
  integer  :: u_traj, u_msd, u_counts, u_pmf, u_pmfex

  real(dp) :: grnd
  external  :: grnd

  character(len=64) :: tag

  ! ---- init ----
  seed = seed0
  call sgrnd(seed)

  allocate(sumx(0:nsteps), sumx2(0:nsteps)); sumx = 0.0_dp; sumx2 = 0.0_dp

  ! final-time mean, std for range (theory)
  mu    = (2.0_dp*p - 1.0_dp) * real(nsteps,dp)
  sigma = sqrt( max(0.0_dp, 4.0_dp*p*(1.0_dp-p) * real(nsteps,dp)) )
  xmin  = int( floor(mu - ksigma*sigma) ) - 2
  xmax  = int( ceiling(mu + ksigma*sigma) ) + 2
  if (xmin > xmax) then; xmin = -2; xmax = 2; end if
  nbins = xmax - xmin + 1
  allocate(hcount(nbins)); hcount = 0

  ! filename tag: pXX_nstepsNNNN
  write(tag,'("p",I0.2,"_nsteps",I0)') int(nint(p*100.0_dp)), nsteps

  ! open outputs
  open(newunit=u_traj, file="traj1.txt", status="replace", action="write")
  write(u_traj,*) "# n   X_n"

  open(newunit=u_msd, file="msd_dsdt_"//trim(tag)//".txt", status="replace", action="write")
  write(u_msd,*) "# n    <X>_MC     MSD_MC      <X>_exact     MSD_exact"

  ! ---- ensemble ----
  do i = 1, ntraj
     x = 0
     ! n=0 contributions
     sumx(0)  = sumx(0)  + 0.0_dp
     sumx2(0) = sumx2(0) + 0.0_dp
     if (i == 1) write(u_traj,'(i10,1x,i10)') 0, x

     do n = 1, nsteps
        r = grnd()
        if (r < p) then
           dir = +1
        else
           dir = -1
        end if
        x = x + dir

        sumx(n)  = sumx(n)  + real(x,dp)
        sumx2(n) = sumx2(n) + real(x*x,dp)

        if (i == 1) write(u_traj,'(i10,1x,i10)') n, x
     end do

     ! final-time histogram
     idx = x - xmin + 1
     if (idx >= 1 .and. idx <= nbins) hcount(idx) = hcount(idx) + 1
  end do
  close(u_traj)

  ! ---- write MSD time series (with exact overlays) ----
  do n = 0, nsteps
     meanx    = sumx(n)/real(ntraj,dp)
     msd      = sumx2(n)/real(ntraj,dp) - meanx*meanx
     meanx_th = (2.0_dp*p - 1.0_dp)*real(n,dp)
     var_th   = 4.0_dp*p*(1.0_dp-p)*real(n,dp)
     write(u_msd,'(i8,1x,4es20.10)') n, meanx, msd, meanx_th, var_th
  end do
  close(u_msd)

  ! ---- write final-time counts ----
  open(newunit=u_counts, file="xCounts_dsdt_"//trim(tag)//".txt", status="replace", action="write")
  write(u_counts,*) "# x   counts"
  do idx = 1, nbins
     write(u_counts,'(i12,1x,i12)') (xmin + (idx-1)), hcount(idx)
  end do
  close(u_counts)

  ! ---- write final-time PMF (MC) ----
  open(newunit=u_pmf, file="xPDF_dsdt_"//trim(tag)//".txt", status="replace", action="write")
  write(u_pmf,*) "# x   P_MC(x, nsteps)"
  if (sum(hcount) > 0) then
     do idx = 1, nbins
        write(u_pmf,'(i12,1x,es20.10)') (xmin + (idx-1)), &
             real(hcount(idx),dp) / real(sum(hcount),dp)
     end do
  end if
  close(u_pmf)

  ! ---- write final-time PMF (exact, on lattice) ----
  open(newunit=u_pmfex, file="xPDF_exact_"//trim(tag)//".txt", status="replace", action="write")
  write(u_pmfex,*) "# x   P_exact(x, nsteps)"
  do idx = 1, nbins
     call write_exact_line(u_pmfex, xmin+(idx-1), nsteps, p)
  end do
  close(u_pmfex)

  deallocate(sumx, sumx2, hcount)

contains
  ! write exact P(X_n = x) with parity/range check
  subroutine write_exact_line(u, x, n, p)
    integer, intent(in) :: u, x, n
    real(dp), intent(in):: p
    real(dp) :: px
    integer  :: k
    if (abs(x) > n .or. mod(n+x,2) /= 0) then
       px = 0.0_dp
    else
       k  = (n + x)/2
       ! compute binomial pmf stably via logs
       px = exp( log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - &
                 log_gamma(real(n-k+1,dp)) + real(k,dp)*log(p) + &
                 real(n-k,dp)*log(1.0_dp-p) )
    end if
    write(u,'(i12,1x,es20.10)') x, px
  end subroutine write_exact_line
end program rw_1D_DSDT_simple_named

! --- RNG (Mersenne Twister) ---
include 'mt.f90'

