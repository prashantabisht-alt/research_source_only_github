!===========================================================================
! Fig. 1 (only): Monte Carlo P(x,t) for RTP on infinite line
! OUTPUT per D: one file "fig1_D<val>.txt" with columns: x  t0.2  t0.4 ... t2.0
!===========================================================================

program rtp_fig1_only_mc
  implicit none
  ! ----------------- knobs -----------------
  integer,  parameter :: Ntraj    = 1000000
  integer,  parameter :: nbins    = 300
  integer,  parameter :: ntimes   = 10           ! 0.2,0.4,...,2.0
  real*8,   parameter :: dt_step  = 0.2d0
  real*8,   parameter :: x_min    = -2.2d0
  real*8,   parameter :: x_max    =  2.2d0
  integer,  parameter :: nD       = 3
  ! -----------------------------------------

  integer :: k, Didx, seed
  real*8  :: Dlist(nD), tgrid(ntimes)
  real*8  :: binw
  integer :: hist(nbins, ntimes)

  ! trajectory state
  real*8  :: x, tcur, sigma, tobs, dt, flip_wait, next_flip, g
  integer :: idx

  ! RNG from mt.f90
  double precision grnd
  external grnd

  ! ---- build time grid ----
  do k = 1, ntimes
     tgrid(k) = dt_step * dble(k)   ! 0.2, 0.4, ..., 2.0
  end do

  ! ---- D values ----
  Dlist = (/ 0.03d0, 0.10d0, 0.20d0 /)

  ! ---- histogram geometry ----
  binw = (x_max - x_min) / dble(nbins)

  ! ---- seed RNG ----
  seed = 31415926
  call sgrnd(seed)

  ! ---- loop over D ----
  do Didx = 1, nD
     call simulate_and_dump(Dlist(Didx), tgrid, ntimes, binw, hist)
  end do

contains

  subroutine simulate_and_dump(D, tgrid, ntimes, binw, hist)
    implicit none
    real*8, intent(in) :: D, tgrid(ntimes), binw
    integer, intent(in) :: ntimes
    integer, intent(inout) :: hist(nbins, ntimes)
    integer :: n, k, b, unitout
    real*8  :: x, tcur, sigma, tobs, dt, flip_wait, next_flip, g
    integer :: idx
    real*8  :: x_center, r, normchk
    character(len=64) :: fname

    hist = 0

    print *,"[D=", D, "] Simulating", Ntraj, "trajectories ..."

    do n = 1, Ntraj
       ! --- init trajectory ---
       x    = 0.0d0
       tcur = 0.0d0
       if (grnd() < 0.5d0) then
         sigma = -1.0d0
       else
         sigma =  1.0d0
       end if
       call exponential(flip_wait, 1.0d0)
       next_flip = tcur + flip_wait

       ! --- loop over observation times ---
       do k = 1, ntimes
          tobs = tgrid(k)

          ! march to tobs across any flips
          do
            if (next_flip < tobs) then
              dt = next_flip - tcur
              if (dt > 0.0d0) then
                call gaussian(g)
                x = x + sigma*dt + dsqrt(2.0d0*D*dt)*g
              end if
              tcur = next_flip
              sigma = -sigma
              call exponential(flip_wait, 1.0d0)
              next_flip = tcur + flip_wait
            else
              dt = tobs - tcur
              if (dt > 0.0d0) then
                call gaussian(g)
                x = x + sigma*dt + dsqrt(2.0d0*D*dt)*g
              end if
              tcur = tobs
              exit
            end if
          end do

          ! bin at t = tgrid(k)
          if (x > x_min .and. x < x_max) then
            idx = floor((x - x_min)/binw) + 1
            if (idx >= 1 .and. idx <= nbins) hist(idx,k) = hist(idx,k) + 1
          end if
       end do

       if (mod(n, 50000) == 0) print *,"  ...", n, "trajectories"
    end do

    ! ---- dump ONE file per D with all times as columns ----
    call make_D_filename(fname, D)
    print *, "Writing:", trim(fname)
    open(newunit=unitout, file=trim(fname), status="replace", action="write")

    ! header
    write(unitout,'(A)', advance="no") "# x"
    do k = 1, ntimes
      write(unitout,'(A)', advance="no") "   t"//trim(adjustl(fmt_tlabel(tgrid(k))))
    end do
    write(unitout,*)

    ! rows: x_center, then P(x,t) for each time
    do b = 1, nbins
      x_center = x_min + (dble(b) - 0.5d0)*binw
      write(unitout,'(1X,ES20.12)', advance="no") x_center
      do k = 1, ntimes
        r = dble(hist(b,k)) / (dble(Ntraj)*binw)
        write(unitout,'(1X,ES20.12)', advance="no") r
      end do
      write(unitout,*)
    end do
    close(unitout)

    ! quick normalization check at last time column
    normchk = 0.0d0
    do b = 1, nbins
      r = dble(hist(b,ntimes)) / (dble(Ntraj)*binw)
      normchk = normchk + r*binw
    end do
    print *,"[D=",D,"] integral (t=",tgrid(ntimes),") ~", normchk

  end subroutine simulate_and_dump

  ! --- predictable filename: "fig1_D0.03.txt", "fig1_D0.10.txt", "fig1_D0.20.txt"
  subroutine make_D_filename(fname, D)
    implicit none
    character(len=*), intent(out) :: fname
    real*8, intent(in) :: D
    integer :: frac
    character(len=6) :: tag
    frac = nint(100.0d0*D)        ! 3, 10, 20
    write(tag,'(I2.2)') frac      ! "03","10","20"
    fname = "fig1_D0."//trim(tag)//".txt"
  end subroutine make_D_filename

  ! --- make time labels "0.2","0.4",...,"2.0" reliably
  function fmt_tlabel(t) result(lbl)
    implicit none
    real*8, intent(in) :: t
    character(len=16) :: lbl
    integer :: tenths
    tenths = nint(10.0d0 * t)      ! 2,4,...,20
    write(lbl,'(F4.1)') dble(tenths)/10.0d0
  end function fmt_tlabel

  ! --- exponential waiting time (rate lambda) ---
  subroutine exponential(s, lambda)
    implicit none
    real*8, intent(out) :: s
    real*8, intent(in)  :: lambda
    real*8 :: u
    do
      u = grnd()
      if (u > 0.d0) exit
    end do
    s = -dlog(u) / lambda
  end subroutine exponential

  ! --- standard normal via Box–Muller (polar form) ---
  subroutine gaussian(z)
    implicit none
    real*8, intent(out) :: z
    real*8 :: x1, x2, w
    double precision grnd
    external grnd
    w = 2.0d0
    do while(w > 1.0d0 .or. w == 0.0d0)
      x1 = 1.0d0 - 2.0d0*grnd()
      x2 = 1.0d0 - 2.0d0*grnd()
      w = x1*x1 + x2*x2
    end do
    z = dsqrt(-2.0d0 * dlog(w)/w) * x1
  end subroutine gaussian

end program rtp_fig1_only_mc

! Mersenne Twister RNG (provides sgrnd, grnd)
include 'mt.f90'
