program ctrw_contSpace_contTime_gaussJump_expWait
  implicit none

  ! --------- parameters ----------
  integer,  parameter :: nsteps  = 1000000  ! number of steps to simulate (1 million)nsteps waiting times in total.
  integer,  parameter :: nbins   = 200
  real*8,   parameter :: lambda  = 1.0d0      ! waiting time rate / whatever you want
  real*8,   parameter :: D       = 1.0d0      ! diffusion coefficient
  real*8,   parameter :: dtmin   = 0.0d0
  real*8,   parameter :: dtmax   = 10.0d0
  integer,  parameter :: seed   = 2025     ! SAME seed for reproducibility

  ! --------- state ----------
  integer :: step, bin_index !bin_index to track which bin the waiting time falls into
  real*8  :: t, dt, bin_width
  real*8  :: position, jump
  real*8, allocatable :: dt_hist(:) !array to hold histogram of waiting times. we'll allocate later and size would be nbins
  real*8  :: r1
  real*8  :: grnd
  external :: grnd !just to declare its external doesnt belong here. it belongs in mt.f90

  ! --------- files ----------
  character(len=*), parameter :: base_tag = "ctrw_contSpace_contTime_gaussJump_expWait_lambda1p0_D1p0"!base tag is to idntify files like ehich program gave this lol.
  integer, parameter :: U_TRAJ=42, U_DTS=43 !just to identify the files

  ! --------- open outputs ----------
  open(unit=U_TRAJ, file="trajectory_"//base_tag//".txt", status="replace")
  open(unit=U_DTS,  file="waitingTimes_"//base_tag//".txt", status="replace")

  ! --------- initialisation ----------
  t        = 0.0d0
  position = 0.0d0
  
  call sgrnd(seed)

  bin_width = (dtmax - dtmin)/nbins
  allocate(dt_hist(nbins)); dt_hist = 0.0d0 !allocate allocates array memory dynamially and histogram array and set all values to 0

  ! --------- simulate ----------
  do step = 1, nsteps
     ! waiting time from exponential distribution
     r1 = grnd() 
     if (r1 <= 0.0d0) r1 = 1.0d-308 ! safety to avoid log(0)
     ! inverse transform sampling to get exponential random variable with rate lambda
     dt = -log(r1) / lambda
     t  = t + dt
     write(U_DTS,*) dt ! log waiting time to file

     ! add Δt to histogram if in range
     if (dt >= dtmin .and. dt < dtmax) then !not less tahn equal to dtmax to avoid double counting on the edge
        ! find which bin this dt falls into
        bin_index = int((dt - dtmin)/bin_width) + 1
        if (bin_index < 1) bin_index = 1 !safety clamp thats it .
        if (bin_index > nbins) bin_index = nbins
        dt_hist(bin_index) = dt_hist(bin_index) + 1.0d0
     end if

     ! Gaussian jump length ~ N(0, 2D*dt)
     jump = sqrt(2.0d0 * D * dt) * gaussian() !variance depend on dt which is waiting time from exponential distribution
     position = position + jump

     ! log trajectory
     write(U_TRAJ,*) t, position
  end do

  close(U_TRAJ); close(U_DTS)

  ! --------- write histogram outputs ----------
  call write_pdf(dt_hist, nbins, dtmin, dtmax, "waitPDF_"//base_tag//".txt")
  call write_counts(dt_hist, nbins, dtmin, dtmax, "waitCounts_"//base_tag//".txt")

  ! --------- summary ----------
  call print_summary(dt_hist, nbins, dtmin, dtmax, bin_width)

  deallocate(dt_hist)

contains
  ! Gaussian RNG using Box–Muller (from Brownian_motion)
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

  ! write normalized PDF to file
  subroutine write_pdf(hist, nbins, xmin, xmax, fname)
    real*8, intent(in) :: hist(nbins), xmin, xmax
    integer, intent(in) :: nbins
    character(*), intent(in) :: fname
    real*8 :: bw, norm, xc
    integer :: i, u
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

  ! write raw counts to file
  subroutine write_counts(hist, nbins, xmin, xmax, fname)
    real*8, intent(in) :: hist(nbins), xmin, xmax
    integer, intent(in) :: nbins
    character(*), intent(in) :: fname
    real*8 :: bw, xc
    integer :: i, u
    bw = (xmax - xmin)/nbins
    open(unit=68, file=fname, status="replace")
    do i = 1, nbins
       xc = xmin + (i - 0.5d0)*bw
       write(68,*) xc, hist(i)
    end do
    close(68)
  end subroutine write_counts

  ! print quick summary
  subroutine print_summary(hist, nbins, xmin, xmax, bw)
    real*8, intent(in) :: hist(nbins), xmin, xmax, bw
    integer, intent(in) :: nbins
    real*8 :: xc, sum_dt, total, mean_dt, lambda_hat !sum_dt is sum of xc*count over all bins, total is total counts in histogram, mean_dt is mean waiting time from histogram, lambda_hat is estimated lambda from mean_dt
    integer :: i
    total = sum(hist)
    if (total > 0.0d0) then !if no in-range samples, skip everything to avoid division by zero.
       sum_dt = 0.0d0
       do i = 1, nbins
          xc = xmin + (i - 0.5d0)*bw
          sum_dt = sum_dt + xc * hist(i)
       end do
       mean_dt = (sum_dt * bw) / (total * bw)  ! bw cancels actually, but kept for clarity numerator is xc* probabilty
       lambda_hat = 1.0d0 / mean_dt
       write(*,'(A,F10.6)') "Estimated lambda_hat = ", lambda_hat
       write(*,'(A,F10.6)') "Mean dt = ", mean_dt
       write(*,'(A,F10.6)') "Fraction outside histogram = ", 1.0d0 - (total / nsteps)
    end if
  end subroutine print_summary

end program ctrw_contSpace_contTime_gaussJump_expWait

! Mersenne Twister RNG
include 'mt.f90'
