program alpha_from_msd
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  character(len=*), parameter :: infile  = "msd_active_taupratio_1e-3two.txt"
  character(len=*), parameter :: outfile = "alpha_active_taupratio_1e-3two.txt"


  integer,  parameter :: SMOOTH_WIN = 2     ! 0 = no smoothing, 2 ≈ mild

  ! --- storage ---
  integer :: i, n, nmax, ios, uin, uout, win, j, j1, j2, cnt
  real(dp), allocatable :: s(:), y(:), lxs(:), lys(:), a(:), as(:)
  logical,  allocatable :: good(:)
  character(len=512) :: line
  real(dp) :: xs, ys

  ! -------- first pass: count data lines --------
  nmax = 0
  open(newunit=uin, file=infile, status="old", action="read", iostat=ios)
  if (ios /= 0) then
    print *, "ERROR: cannot open ", trim(infile)
    stop 1
  end if
  do
    read(uin,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == "#") cycle
    nmax = nmax + 1
  end do
  close(uin)
  if (nmax < 3) then
    print *, "ERROR: not enough data points in ", trim(infile)
    stop 1
  end if

  ! -------- read data --------
  allocate(s(nmax), y(nmax))
  n = 0
  open(newunit=uin, file=infile, status="old", action="read")
  do
    read(uin,'(A)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == "#") cycle
    read(line,*, iostat=ios) xs, ys
    if (ios == 0) then
      n = n + 1
      s(n) = xs
      y(n) = ys
    end if
  end do
  close(uin)

  ! -------- prepare logs & validity --------
  allocate(lxs(n), lys(n), a(n), as(n), good(n))
  do i = 1, n
    good(i) = (s(i) > 0.0_dp) .and. (y(i) > 0.0_dp)
    if (good(i)) then
      lxs(i) = log(s(i))
      lys(i) = log(y(i))
    else
      lxs(i) = 0.0_dp
      lys(i) = 0.0_dp
    end if
  end do

  ! -------- finite-difference derivative on log-log --------
  do i = 1, n
    if (.not. good(i)) then
      a(i) = 0.0_dp
    else if (i == 1) then
      if (good(2)) then
        a(i) = (lys(2) - lys(1)) / (lxs(2) - lxs(1))
      else
        a(i) = 0.0_dp
      end if
    else if (i == n) then
      if (good(n-1)) then
        a(i) = (lys(n) - lys(n-1)) / (lxs(n) - lxs(n-1))
      else
        a(i) = 0.0_dp
      end if
    else
      if (good(i-1) .and. good(i+1)) then
        a(i) = (lys(i+1) - lys(i-1)) / (lxs(i+1) - lxs(i-1))
      else if (good(i+1)) then
        a(i) = (lys(i+1) - lys(i)) / (lxs(i+1) - lxs(i))
      else if (good(i-1)) then
        a(i) = (lys(i) - lys(i-1)) / (lxs(i) - lxs(i-1))
      else
        a(i) = 0.0_dp
      end if
    end if
  end do

  ! -------- optional smoothing --------
  win = SMOOTH_WIN
  if (win > 0) then
    do i = 1, n
      j1  = max(1, i - win)
      j2  = min(n, i + win)
      cnt = 0
      as(i) = 0.0_dp
      do j = j1, j2
        if (good(j)) then
          as(i) = as(i) + a(j)
          cnt   = cnt + 1
        end if
      end do
      if (cnt > 0) then
        as(i) = as(i) / real(cnt, dp)
      else
        as(i) = a(i)
      end if
    end do
  else
    as = a
  end if

  ! -------- write output --------
  open(newunit=uout, file=outfile, status="replace", action="write")
  write(uout,*) "# s=t/tau_p    alpha(s) = d ln(MSD)/d ln s"
  do i = 1, n
    if (good(i)) write(uout,'(2(1X,ES20.12))') s(i), as(i)
  end do
  close(uout)

  print *, "Wrote: ", trim(outfile)
end program alpha_from_msd

