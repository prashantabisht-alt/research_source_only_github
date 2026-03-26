!===============================================================
! alpha_from_msd_binned.f90
! Smooth α(t) = d ln(MSD)/d ln t from a 2-col file (t/τI, MSD_norm)
! Steps:
!   1) ignore comments/#, require t>0, MSD>0
!   2) log-bin MSD (geometric mean) with BINS_PER_DECADE
!   3) sliding-window least-squares slope on log-log (window grows with t)
! Outputs: alpha_<infile>.txt (t/τI, alpha)
!===============================================================
program alpha_from_msd_binned
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  character(len=*), parameter :: infile  = "msd_active_taupratio_1e3two.txt"
  character(len=*), parameter :: outfile = "alpha_active_taupratio_1e3two.txt"

  integer,  parameter :: BINS_PER_DECADE = 40   ! 30–50 works well
  integer,  parameter :: MIN_PER_BIN     = 10   ! discard sparse bins
  integer,  parameter :: WIN_MIN         = 5    ! min regression half-width (bins)
  real(dp), parameter :: WIN_GROW        = 0.03_dp ! how fast window grows with index

  ! raw arrays
  integer :: i, nmax, n, ios, uin, uout
  real(dp), allocatable :: t(:), y(:)
  character(len=512) :: line
  real(dp) :: tt, yy

  ! binning arrays (accumulate in log space for geometric means)
  integer :: nbins, k, cnt
  real(dp) :: log10t_min, log10t_max, dlog10
  real(dp), allocatable :: sum_logt(:), sum_logy(:); integer, allocatable :: count(:)
  real(dp), allocatable :: tbin(:), ybin(:)  ! binned series (linear domain)
  real(dp), allocatable :: lx(:), ly(:)      ! logs of binned series
  integer :: m, j, j1, j2, win

  ! alpha
  real(dp), allocatable :: alpha(:)

  ! ---------- pass 1: count data ----------
  nmax = 0
  open(newunit=uin, file=infile, status="old", action="read", iostat=ios)
  if (ios /= 0) stop "Cannot open input file"
  do
    read(uin,'(A)', iostat=ios) line; if (ios /= 0) exit
    if (len_trim(line)==0) cycle
    if (line(1:1)=="#") cycle
    nmax = nmax + 1
  end do
  close(uin)
  if (nmax < 10) stop "Not enough data"

  ! ---------- read data ----------
  allocate(t(nmax), y(nmax)); n=0
  open(newunit=uin, file=infile, status="old", action="read")
  do
    read(uin,'(A)', iostat=ios) line; if (ios /= 0) exit
    if (len_trim(line)==0) cycle
    if (line(1:1)=="#") cycle
    read(line,*,iostat=ios) tt, yy
    if (ios==0 .and. tt>0.0_dp .and. yy>0.0_dp) then
      n = n + 1; t(n)=tt; y(n)=yy
    end if
  end do
  close(uin)
  if (n<10) stop "Not enough positive points"

  ! ---------- define log-time bins ----------
  log10t_min = log10(t(1)); log10t_max = log10(t(n))
  dlog10     = 1.0_dp / real(BINS_PER_DECADE,dp)
  nbins      = int((log10t_max - log10t_min)/dlog10) + 1

  allocate(sum_logt(nbins), sum_logy(nbins), count(nbins))
  sum_logt=0.0_dp; sum_logy=0.0_dp; count=0

  do i=1,n
    k = int( (log10(t(i)) - log10t_min) / dlog10 ) + 1
    if (k<1 .or. k>nbins) cycle
    sum_logt(k) = sum_logt(k) + log(t(i))
    sum_logy(k) = sum_logy(k) + log(y(i))
    count(k)    = count(k) + 1
  end do

  ! ---------- compact binned series (geometric means) ----------
  allocate(tbin(nbins), ybin(nbins), lx(nbins), ly(nbins))
  m = 0
  do k=1,nbins
    if (count(k) >= MIN_PER_BIN) then
      m = m + 1
      lx(m)  = sum_logt(k) / real(count(k),dp)   ! log t̄
      ly(m)  = sum_logy(k) / real(count(k),dp)   ! log ȳ
      tbin(m)= exp(lx(m)); ybin(m)=exp(ly(m))
    end if
  end do
  if (m < 10) stop "Too few bins after filtering"

  ! ---------- sliding-window least-squares slope ----------
  allocate(alpha(m))
  do j=1,m
    ! window half-width grows with j, but at least WIN_MIN
    win = max( WIN_MIN, int(WIN_GROW*real(j,dp)) )
    j1 = max(1, j - win); j2 = min(m, j + win)
    call local_slope(lx(j1:j2), ly(j1:j2), alpha(j))
  end do

  ! ---------- write ----------
  open(newunit=uout, file=outfile, status="replace", action="write")
  write(uout,*) "# t/τI    alpha(t)  (log-binned + local regression)"
  do j=1,m
    write(uout,'(2(1X,ES20.12))') tbin(j), alpha(j)
  end do
  close(uout)
  print *, "Wrote:", trim(outfile)

contains
  subroutine local_slope(x, y, a)
    real(dp), intent(in) :: x(:), y(:)   ! logs
    real(dp), intent(out):: a
    integer :: nloc, k
    real(dp) :: sx, sy, sxx, sxy, w
    nloc = size(x); if (nloc<2) then; a=0.0_dp; return; end if
    sx=0.0_dp; sy=0.0_dp; sxx=0.0_dp; sxy=0.0_dp
    do k=1,nloc
      w = 1.0_dp                          ! (could add weights if desired)
      sx  = sx  + w*x(k)
      sy  = sy  + w*y(k)
      sxx = sxx + w*x(k)*x(k)
      sxy = sxy + w*x(k)*y(k)
    end do
    a = (real(nloc,dp)*sxy - sx*sy) / max(1.0e-30_dp, (real(nloc,dp)*sxx - sx*sx))
  end subroutine local_slope
end program alpha_from_msd_binned
