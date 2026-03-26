!===========================================================================
! GOAL: Reproduce Fig. 1 for a 1D Run-and-Tumble Particle (infinite line)
! MODEL (dimensionless): xdot = sigma(t) + sqrt(2D) * eta(t),
!   where sigma(t) = +/- 1 flips at rate 1 (Poisson), eta is white Gaussian.
! INITIAL CONDITION: x(0)=0, sigma(0)=± with prob 1/2
!
! WHAT WE OUTPUT (per D):
!   P_D<val>_t<val>.txt         : x, P(x,t) for t=0.2,0.4,...,2.0
!   P0_vs_t_D<val>.txt          : t, P(0,t) estimate (center bin)
!   curvature0_vs_t_D<val>.txt  : t, discrete curvature at x=0
!   msd_vs_t_D<val>.txt         : t, <x^2>(t) (MC)
!

!===========================================================================

program rtp_infinite_mc_fig1
  implicit none
  ! ----------------- knobs you can tweak -----------------
  integer,  parameter :: Ntraj    = 10000000      ! trajectories per D
  integer,  parameter :: nbins    = 300         ! histogram bins over x-range
  integer,  parameter :: ntimes   = 10          ! 0.2, 0.4, ..., 2.0
  real*8,   parameter :: dx_times = 0.2d0
  real*8,   parameter :: x_min    = -2.8d0
  real*8,   parameter :: x_max    =  2.8d0
  ! -------------------------------------------------------

  integer :: i, k, b, seed, Didx
  real*8  :: Dlist(3)
  real*8  :: tgrid(ntimes), binw, x_center
  integer :: hist(nbins, ntimes)
  real*8  :: sumx2(ntimes)
  real*8  :: P0(ntimes), curv0(ntimes)
  real*8  :: normchk, t, x, sigma
  real*8  :: next_flip, tobs, dt
  real*8  :: g, r
  real*8 :: grnd
  integer :: idx, idx0, unitout
  character(len=64) :: fname

  ! ---- init time grid ----
  do k = 1, ntimes
     tgrid(k) = dx_times * dble(k)
  end do

  ! ---- D values (dimensionless) ----
  Dlist = (/ 0.03d0, 0.1d0, 0.2d0 /)

  ! ---- bins ----
  binw = (x_max - x_min) / dble(nbins)
  idx0 = floor((0.0d0 - x_min) / binw) + 1   ! center bin index (for P(0,t))

  ! ---- RNG ----
  seed = 31415926
  call sgrnd(seed)

  do Didx = 1, 3
     call simulate_one_D(Dlist(Didx))
  end do

contains

  subroutine simulate_one_D(D)
    implicit none
    real*8, intent(in) :: D
    integer :: n, k, b
    real*8  :: tcur, x, g, flip_wait, sigma
    real*8  :: next_flip, tobs, dt, Pval_left, Pval_mid, Pval_right
    integer :: idx, unitout

    hist   = 0
    sumx2  = 0.0d0
    P0     = 0.0d0
    curv0  = 0.0d0

    print *,"[D=", D, "] Simulating", Ntraj, "trajectories ..."

    do n = 1, Ntraj
       ! --- initial condition ---
       x    = 0.0d0
       tcur = 0.0d0
       if (grnd() < 0.5d0) then
          sigma = -1.0d0
       else
          sigma =  1.0d0
       end if
       call exponential(flip_wait, 1.0d0)  ! rate = 1
       next_flip = tcur + flip_wait

       do k = 1, ntimes
          tobs = tgrid(k)

          ! march from tcur to tobs, possibly through several flips
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

          ! bin x at observation time tgrid(k)
          if (x > x_min .and. x < x_max) then
             idx = floor((x - x_min)/binw) + 1
             if (idx >= 1 .and. idx <= nbins) hist(idx, k) = hist(idx, k) + 1
          end if

          ! accumulate MSD at this t
          sumx2(k) = sumx2(k) + x*x
       end do

       if (mod(n, 50000) == 0) print *,"  ... finished", n, "trajectories"
    end do

    ! ----- write PDFs per time, plus P(0,t) and curvature at center -----
    do k = 1, ntimes
       write(fname,'("P_D",f0.2,"_t",f0.1,".txt")') D, tgrid(k)
       call sanitize_filename(fname)
       open(newunit=unitout, file=trim(fname), status="replace", action="write")
       write(unitout,*) "# x   P(x,t)   (D=",D,", t=",tgrid(k),")"
       normchk = 0.0d0
       do b = 1, nbins
          x_center = x_min + (dble(b) - 0.5d0)*binw
          ! normalize histogram to PDF
          r = dble(hist(b,k)) / (dble(Ntraj)*binw)
          write(unitout,'(2(1X,ES20.12))') x_center, r
          normchk = normchk + r*binw
       end do
       close(unitout)

       ! center diagnostics using neighboring bins (second difference)
       if (idx0>=2 .and. idx0<=nbins-1) then
          Pval_left  = dble(hist(idx0-1,k)) / (dble(Ntraj)*binw)
          Pval_mid   = dble(hist(idx0  ,k)) / (dble(Ntraj)*binw)
          Pval_right = dble(hist(idx0+1,k)) / (dble(Ntraj)*binw)
          P0(k)      = Pval_mid
          curv0(k)   = (Pval_right - 2.0d0*Pval_mid + Pval_left) / (binw*binw)
       else
          P0(k)      = dble(hist(idx0,k)) / (dble(Ntraj)*binw)
          curv0(k)   = 0.0d0
       end if

       ! quick normalization check
       if (k==ntimes) print *,"[D=",D,"] norm check (t=",tgrid(k),") ~", normchk
    end do

    ! write P(0,t)
    write(fname,'("P0_vs_t_D",f0.2,".txt")') D
    call sanitize_filename(fname)
    open(newunit=unitout, file=trim(fname), status="replace", action="write")
    write(unitout,*) "# t   P(0,t)   (D=",D,")"
    do k = 1, ntimes
       write(unitout,'(2(1X,ES20.12))') tgrid(k), P0(k)
    end do
    close(unitout)

    ! write curvature at center
    write(fname,'("curvature0_vs_t_D",f0.2,".txt")') D
    call sanitize_filename(fname)
    open(newunit=unitout, file=trim(fname), status="replace", action="write")
    write(unitout,*) "# t   curvature@0   (D=",D,")"
    do k = 1, ntimes
       write(unitout,'(2(1X,ES20.12))') tgrid(k), curv0(k)
    end do
    close(unitout)

    ! write MSD
    write(fname,'("msd_vs_t_D",f0.2,".txt")') D
    !call sanitize_filename(fname)
    open(newunit=unitout, file=trim(fname), status="replace", action="write")
    write(unitout,*) "# t   <x^2>(t)   (D=",D,")"
    do k = 1, ntimes
       write(unitout,'(2(1X,ES20.12))') tgrid(k), sumx2(k)/dble(Ntraj)
    end do
    close(unitout)

  end subroutine simulate_one_D

  ! --- helper: exponential waiting time with rate lambda ---
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

  ! --- helper: replace '.' with 'p' in filenames to keep things tidy ---
  subroutine sanitize_filename(str)
    character(len=*), intent(inout) :: str
    integer :: i
    do i = 1, len_trim(str)
      if (str(i:i)=='.') str(i:i) = 'p'
      if (str(i:i)==' ') str(i:i) = '_'
    end do
  end subroutine sanitize_filename

  ! ---------- your Gaussian (Box–Muller), using mt.f90::grnd() ----------
  subroutine gaussian(s)
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while(w > 1.0d0 .or. w == 0.0d0)
        x1 = 1.0d0 - 2.0d0*grnd()
        x2 = 1.0d0 - 2.0d0*grnd()
        w = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0 * dlog(w)/w) * x1
  end subroutine gaussian

end program rtp_infinite_mc_fig1

! your RNG
include 'mt.f90'
