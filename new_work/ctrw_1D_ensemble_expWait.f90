program ctrw_1D_ensemble_expWait_exactT
  implicit none
  ! ============================================================
  ! Continuous-Time Random Walk (CTRW), 1D, ENSEMBLE version
  ! - Waiting times Δt ~ Exp(λ)
  ! - Jumps Δx ~ Normal(0, 2 D Δt)
  ! - For each particle, evolve until EXACTLY t = tmax by trimming
  !   the final interval to leftover time dt_left = tmax - t.
  !
  ! Output: normalized PDF of positions P(x, tmax) and raw counts.
  ! Gaussian generator uses the same interface as your Brownian code:
  !   subroutine gaussian(s)  -> returns one N(0,1) in 's'
  ! ============================================================

  ! ---------- Parameters (tune for your laptop) ----------
  integer,  parameter :: ntraj      = 100000      ! number of particles
  integer,  parameter :: nbins      = 200         ! histogram bins in x
  integer,  parameter :: nsteps_max = 1000000     ! safety cap per trajectory
  real*8,   parameter :: lambda     = 1.0d0       ! waiting-time rate
  real*8,   parameter :: D          = 1.0d0       ! diffusion coefficient
  real*8,   parameter :: tmax       = 10000.0d0      ! observation time (snapshot)
  real*8,   parameter :: xmin       = -800.0d0    ! histogram range
  real*8,   parameter :: xmax       =  800.0d0
  integer,  parameter :: seed0      = 2025

  ! ---------- State ----------
  integer :: traj, step, bin_index, seed
  real*8  :: t, dt_cand, dt_left, binw
  real*8  :: x, jump, r1, z          ! z = standard normal from gaussian()
  real*8, allocatable :: x_hist(:)
  real*8  :: grnd                    ! from mt.f90
  external :: grnd

  ! ---------- File tag ----------
  character(len=*), parameter :: base_tag = "ctrw_1D_ensemble_expWait_exactT_lambda1p0_D1p0_tmax10"

  ! ---------- Init ----------
  seed = seed0
  call sgrnd(seed)                   ! seed the Mersenne Twister RNG
  binw = (xmax - xmin)/nbins
  allocate(x_hist(nbins)); x_hist = 0.0d0

  ! ============================================================
  ! ENSEMBLE LOOP: simulate ntraj independent particles
  ! ============================================================
  do traj = 1, ntraj
     t = 0.0d0
     x = 0.0d0

     ! Keep proposing exponential waits until we hit tmax.
     do step = 1, nsteps_max
        ! Draw candidate waiting time Δt_cand ~ Exp(λ) via inverse CDF:
        r1 = grnd(); if (r1 <= 0.0d0) r1 = 1.0d-308
        dt_cand = -log(r1) / lambda

        if (t + dt_cand >= tmax) then
           ! Would overshoot: trim to leftover so we land EXACTLY at tmax
           dt_left = tmax - t
           if (dt_left > 0.0d0) then
              call gaussian(z)                               ! z ~ N(0,1)
              jump = sqrt(2.0d0 * D * dt_left) * z           ! Δx ~ N(0, 2D dt_left)
              x    = x + jump
           end if
           t = tmax
           exit                                              ! finished this particle
        else
           ! Safe to take the full candidate interval
           t = t + dt_cand
           call gaussian(z)                                  ! z ~ N(0,1)
           jump = sqrt(2.0d0 * D * dt_cand) * z              ! Δx ~ N(0, 2D dt_cand)
           x    = x + jump
        end if
     end do
     ! (nsteps_max is just a guard; with sane params we exit earlier.)

     ! Bin the final position x at EXACT time t = tmax
     if (x >= xmin .and. x < xmax) then
        bin_index = int((x - xmin)/binw) + 1
        if (bin_index < 1)    bin_index = 1
        if (bin_index > nbins) bin_index = nbins
        x_hist(bin_index) = x_hist(bin_index) + 1.0d0
     end if
  end do

  ! ---------- Outputs ----------
  call write_pdf   (x_hist, nbins, xmin, xmax, "xPDF_"   //base_tag//".txt")
  call write_counts(x_hist, nbins, xmin, xmax, "xCounts_"//base_tag//".txt")

  deallocate(x_hist)

contains
  ! ============================================================
  ! Gaussian RNG (Marsaglia polar / Box–Muller style as in your code)
  ! Returns ONE N(0,1) sample in 's'
  ! ============================================================
  subroutine gaussian(s)
    real*8, intent(out) :: s
    real*8 :: x1, x2, w
    w = 2.0d0
    do while (w >= 1.0d0 .or. w == 0.0d0)
       x1 = 1.0d0 - 2.0d0 * grnd()
       x2 = 1.0d0 - 2.0d0 * grnd()
       w  = x1*x1 + x2*x2
    end do
    s = sqrt(-2.0d0 * log(w) / w) * x1
  end subroutine gaussian

  ! ============================================================
  ! Write normalized PDF: columns = (bin_center, p(x))
  ! Area ≈ 1 over [xmin, xmax)
  ! ============================================================
  subroutine write_pdf(hist, nb, xmin, xmax, fname)
    real*8, intent(in) :: hist(:), xmin, xmax
    integer, intent(in):: nb
    character(*), intent(in) :: fname
    real*8 :: bw, norm, xc
    integer :: i, u
    bw   = (xmax - xmin)/nb
    norm = sum(hist)*bw
    if (norm <= 0.0d0) return
    open(newunit=u, file=fname, status="replace")
    do i = 1, nb
      xc = xmin + (i - 0.5d0)*bw
      write(u,*) xc, hist(i)/norm
    end do
    close(u)
  end subroutine write_pdf

  ! ============================================================
  ! Write raw counts: columns = (bin_center, counts)
  ! Useful for error bars (Poisson) or re-normalization
  ! ============================================================
  subroutine write_counts(hist, nb, xmin, xmax, fname)
    real*8, intent(in) :: hist(:), xmin, xmax
    integer, intent(in):: nb
    character(*), intent(in) :: fname
    real*8 :: bw, xc
    integer :: i, u
    bw = (xmax - xmin)/nb
    open(newunit=u, file=fname, status="replace")
    do i = 1, nb
      xc = xmin + (i - 0.5d0)*bw
      write(u,*) xc, hist(i)
    end do
    close(u)
  end subroutine write_counts

end program ctrw_1D_ensemble_expWait_exactT

! RNG (Mersenne Twister): provides grnd() and sgrnd()
include 'mt.f90'
