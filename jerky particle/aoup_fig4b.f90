!======================================================================
! Fig. 4(b): alpha(s) = d ln(MSD)/d ln s for inertial-jerky AOUP
! Processes two MSD files from Fig. 4(a):
!   - msd_active_taupratio_1e-3two.txt  (τp = τI/1000 → 6→5→3)
!   - msd_active_taupratio_1e3two.txt   (τp = 1000 τI → 6→4→3)
!
! Writes:
!   - alpha_taupratio_1e-3two.txt   (s, α(s))
!   - alpha_taupratio_1e3two.txt    (s, α(s))
!
! Notes:
!  - s here is t/τI (x-axis in Fig. 4).
!  - Input MSD must be normalized already (as you wrote in your sim).
!======================================================================
program alpha_fig4b
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: SMOOTH_WIN = 2     ! 0 = off, 2 ~ mild smoothing
  character(len=*), parameter :: in_fast  = "msd_active_taupratio_1e-3two.txt"
  character(len=*), parameter :: in_slow  = "msd_active_taupratio_1e3two.txt"
  character(len=*), parameter :: out_fast = "alpha_taupratio_1e-3two.txt"
  character(len=*), parameter :: out_slow = "alpha_taupratio_1e3two.txt"

  call process_alpha(in_fast, out_fast, SMOOTH_WIN)
  call process_alpha(in_slow, out_slow, SMOOTH_WIN)

  print *, "Wrote:"
  print *, "  ", trim(out_fast)
  print *, "  ", trim(out_slow)

contains

  subroutine process_alpha(infile, outfile, smooth_win)
    implicit none
    character(len=*), intent(in) :: infile, outfile
    integer,          intent(in) :: smooth_win

    integer :: i, n, nmax, ios, uin, uout, win, j, j1, j2, cnt
    real(dp), allocatable :: s(:), y(:), lxs(:), lys(:), a(:), as(:)
    logical,  allocatable :: good(:)
    character(len=512)    :: line
    real(dp), parameter   :: eps = 1.0e-14_dp
    real(dp) :: xs, ys

    ! ---------- pass 1: count numeric rows (ignore # / blanks) ----------
    nmax = 0
    open(newunit=uin, file=infile, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      print *, "ERROR: cannot open ", trim(infile); stop 1
    end if
    do
      read(uin,'(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == "#")    cycle
      nmax = nmax + 1
    end do
    close(uin)
    if (nmax < 5) then
      print *, "ERROR: need >= 5 data rows in ", trim(infile); stop 1
    end if

    allocate(s(nmax), y(nmax))
    n = 0

    ! ---------------- pass 2: read numeric data ----------------
    open(newunit=uin, file=infile, status="old", action="read")
    do
      read(uin,'(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == "#")    cycle
      read(line,*, iostat=ios) xs, ys
      if (ios == 0) then
        n = n + 1
        s(n) = xs
        y(n) = ys
      end if
    end do
    close(uin)

    allocate(lxs(n), lys(n), a(n), as(n), good(n))

    ! -------- logs & valid flags --------
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

    ! -------- alpha = d ln(y) / d ln(s) --------
    do i = 1, n
      if (.not. good(i)) then
        a(i) = 0.0_dp

      else if (i == 1) then
        if (good(2) .and. abs(lxs(2)-lxs(1)) > eps) then
          a(i) = (lys(2) - lys(1)) / (lxs(2) - lxs(1))
        else
          a(i) = 0.0_dp
        end if

      else if (i == n) then
        if (good(n-1) .and. abs(lxs(n)-lxs(n-1)) > eps) then
          a(i) = (lys(n) - lys(n-1)) / (lxs(n) - lxs(n-1))
        else
          a(i) = 0.0_dp
        end if

      else
        if (good(i-1) .and. good(i+1) .and. abs(lxs(i+1)-lxs(i-1)) > eps) then
          a(i) = (lys(i+1) - lys(i-1)) / (lxs(i+1) - lxs(i-1))
        else if (good(i+1) .and. abs(lxs(i+1)-lxs(i)) > eps) then
          a(i) = (lys(i+1) - lys(i)) / (lxs(i+1) - lxs(i))
        else if (good(i-1) .and. abs(lxs(i)-lxs(i-1)) > eps) then
          a(i) = (lys(i) - lys(i-1)) / (lxs(i) - lxs(i-1))
        else
          a(i) = 0.0_dp
        end if
      end if
    end do

    ! -------- optional smoothing --------
    win = smooth_win
    if (win > 0) then
      do i = 1, n
        j1 = max(1, i - win); j2 = min(n, i + win)
        cnt = 0; as(i) = 0.0_dp
        do j = j1, j2
          if (good(j)) then
            as(i) = as(i) + a(j); cnt = cnt + 1
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

    ! -------- write trimmed output (skip 2-edge points) --------
    open(newunit=uout, file=outfile, status="replace", action="write")
    write(uout,*) "# s=t/tau_I    alpha(s) = d ln(MSD)/d ln s"
    do i = 3, n-2
      if (good(i)) write(uout,'(2(1X,ES20.12))') s(i), as(i)
    end do
    close(uout)

    deallocate(s, y, lxs, lys, a, as, good)
  end subroutine process_alpha

end program alpha_fig4b
