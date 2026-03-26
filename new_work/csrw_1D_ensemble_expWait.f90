program csrw_1D_ensemble_exactT
  implicit none
  ! ============================================================
  ! CSRW (discrete-time Gaussian steps), 1D, ENSEMBLE version
  ! Exact snapshot at t = tmax. Since x(tmax) ~ N(0, 2 D tmax),
  ! we can either simulate increments or draw once per particle.
  ! Here we draw ONCE per trajectory for exactness & speed.
  ! Output: normalized PDF P(x, tmax) and raw counts (for error bars).
  ! ============================================================

  ! ---------- Parameters (match CTRW) ----------
  integer,  parameter :: ntraj  = 100000      ! number of particles
  integer,  parameter :: nbins  = 200         ! histogram bins in x
  real*8,   parameter :: D      = 1.0d0       ! diffusion coefficient
  real*8,   parameter :: tmax   = 10000.0d0   ! observation time
  real*8,   parameter :: xmin   = -800.0d0    ! histogram range
  real*8,   parameter :: xmax   =  800.0d0
  integer,  parameter :: seed0  = 2025

  ! ---------- State ----------
  integer :: traj, bin_index, seed
  real*8  :: binw, x, z
  real*8, allocatable :: x_hist(:)
  real*8  :: grnd
  external :: grnd

  ! ---------- File tag (distinct from CTRW) ----------
  character(len=*), parameter :: base_tag = "csrw_1D_ensemble_exactT_lambdaNA_D1p0_tmax10"

  ! ---------- Init ----------
  seed = seed0
  call sgrnd(seed)
  binw = (xmax - xmin)/nbins
  allocate(x_hist(nbins)); x_hist = 0.0d0

  ! ============================================================
  ! ENSEMBLE LOOP: x(tmax) ~ N(0, 2 D tmax)
  ! ============================================================
  do traj = 1, ntraj
     call gaussian(z)                           ! z ~ N(0,1)
     x = sqrt(2.0d0*D*tmax) * z                 ! exact snapshot
     if (x >= xmin .and. x < xmax) then
        bin_index = int((x - xmin)/binw) + 1
        if (bin_index < 1)     bin_index = 1
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
  ! Gaussian RNG (Marsaglia polar / Box–Muller)
  ! ============================================================
  subroutine gaussian(s)
    real*8, intent(out) :: s
    real*8 :: x1, x2, w
    w = 2.0d0
    do while (w >= 1.0d0 .or. w == 0.0d0)
       x1 = 1.0d0 - 2.0d0*grnd()
       x2 = 1.0d0 - 2.0d0*grnd()
       w  = x1*x1 + x2*x2
    end do
    s = sqrt(-2.0d0 * log(w) / w) * x1
  end subroutine gaussian

  ! ============================================================
  ! Normalized PDF: (bin_center, p(x))
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
  ! Raw counts: (bin_center, counts) — for error bars
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

end program csrw_1D_ensemble_exactT

! RNG (Mersenne Twister): provides grnd() and sgrnd()
include 'mt.f90'
