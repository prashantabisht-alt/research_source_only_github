program ctrw_2D_contSpace_contTime_gaussJump_expWait
  implicit none

  ! --------- parameters ----------
  integer,  parameter :: nsteps  = 10000
  integer,  parameter :: nbins   = 200
  real*8,   parameter :: lambda  = 1.0d0
  real*8,   parameter :: D       = 1.0d0
  real*8,   parameter :: dtmin   = 0.0d0
  real*8,   parameter :: dtmax   = 10.0d0
  integer,  parameter :: seed0   = 2025

  ! --------- state ----------
  integer :: step, seed, bin_index
  real*8  :: t, dt, bin_width
  real*8  :: x, y, jx, jy
  real*8, allocatable :: dt_hist(:)
  real*8  :: r1
  real*8  :: grnd
  external :: grnd

  ! --------- files ----------
  character(len=*), parameter :: base_tag = "ctrw_2D_contSpace_contTime_gaussJump_expWait_lambda1p0_D1p0"
  integer, parameter :: U_TRAJ=42, U_DTS=43

  ! --------- open outputs ----------
  open(unit=U_TRAJ, file="trajectory_"//base_tag//".txt", status="replace")
  open(unit=U_DTS,  file="waitingTimes_"//base_tag//".txt", status="replace")

  ! --------- init ----------
  t    = 0.0d0
  x    = 0.0d0
  y    = 0.0d0
  seed = seed0
  call sgrnd(seed)

  bin_width = (dtmax - dtmin)/nbins
  allocate(dt_hist(nbins)); dt_hist = 0.0d0

  ! --------- simulate ----------
  do step = 1, nsteps
     ! waiting time from exponential distribution
     r1 = grnd()
     if (r1 <= 0.0d0) r1 = 1.0d-308
     dt = -log(r1) / lambda
     t  = t + dt
     write(U_DTS,*) dt

     ! add Δt to histogram
     if (dt >= dtmin .and. dt < dtmax) then
        bin_index = int((dt - dtmin)/bin_width) + 1
        if (bin_index < 1) bin_index = 1
        if (bin_index > nbins) bin_index = nbins
        dt_hist(bin_index) = dt_hist(bin_index) + 1.0d0
     end if

     ! Gaussian jumps in x and y
     jx = sqrt(2.0d0 * D * dt) * gaussian()
     jy = sqrt(2.0d0 * D * dt) * gaussian()
     x  = x + jx
     y  = y + jy

     ! log trajectory
     write(U_TRAJ,*) t, x, y
  end do

  close(U_TRAJ); close(U_DTS)

  ! --------- write histogram outputs ----------
  call write_pdf(dt_hist, nbins, dtmin, dtmax, "waitPDF_"//base_tag//".txt")
  call write_counts(dt_hist, nbins, dtmin, dtmax, "waitCounts_"//base_tag//".txt")
  call print_summary(dt_hist, nbins, dtmin, dtmax, bin_width)

  deallocate(dt_hist)

contains
  ! Gaussian RNG using Box–Muller
  real*8 function gaussian()
    real*8 :: x1, x2, w
    w = 2.0d0
    do while(w >= 1.0d0 .or. w == 0.0d0)
       x1 = 1.0d0 - 2.0d0 * grnd()
       x2 = 1.0d0 - 2.0d0 * grnd()
       w = x1*x1 + x2*x2
    end do
    gaussian = sqrt(-2.0d0 * log(w) / w) * x1
  end function gaussian

  subroutine write_pdf(hist, nbins, xmin, xmax, fname)
    real*8, intent(in) :: hist(nbins), xmin, xmax
    integer, intent(in) :: nbins
    character(*), intent(in) :: fname
    real*8 :: bw, norm, xc
    integer :: i
    bw   = (xmax - xmin)/nbins
    norm = sum(hist)*bw
    if (norm <= 0.0d0) return
    open(unit=69, file=fname, status="replace")
    do i = 1, nbins
       xc = xmin + (i - 0.5d0)*bw
       write(69,*) xc, hist(i)/norm
    end do
    close(69)
  end subroutine write_pdf

  subroutine write_counts(hist, nbins, xmin, xmax, fname)
    real*8, intent(in) :: hist(nbins), xmin, xmax
    integer, intent(in) :: nbins
    character(*), intent(in) :: fname
    real*8 :: bw, xc
    integer :: i
    bw = (xmax - xmin)/nbins
    open(unit=68, file=fname, status="replace")
    do i = 1, nbins
       xc = xmin + (i - 0.5d0)*bw
       write(68,*) xc, hist(i)
    end do
    close(68)
  end subroutine write_counts

  subroutine print_summary(hist, nbins, xmin, xmax, bw)
    real*8, intent(in) :: hist(nbins), xmin, xmax, bw
    integer, intent(in) :: nbins
    real*8 :: xc, sum_dt, total, mean_dt, lambda_hat
    integer :: i
    total = sum(hist)
    if (total > 0.0d0) then
       sum_dt = 0.0d0
       do i = 1, nbins
          xc = xmin + (i - 0.5d0)*bw
          sum_dt = sum_dt + xc * hist(i)
       end do
       mean_dt = (sum_dt * bw) / (total * bw)
       lambda_hat = 1.0d0 / mean_dt
       write(*,'(A,F10.6)') "Estimated lambda_hat = ", lambda_hat
       write(*,'(A,F10.6)') "Mean dt = ", mean_dt
       write(*,'(A,F10.6)') "Fraction outside histogram = ", 1.0d0 - (total / nsteps)
    end if
  end subroutine print_summary

end program ctrw_2D_contSpace_contTime_gaussJump_expWait

include 'mt.f90'
