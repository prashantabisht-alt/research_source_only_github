!===============================================================================
! 1D Discrete-Space, Discrete-Time Random Walk (Ensemble, with exact overlays)
!===============================================================================
program rw_1D_DSDT_ensemble
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)

  ! ---------------- knobs ----------------
  integer,  parameter :: ntraj   = 100000
  integer,  parameter :: nsteps  = 10000
  real(dp), parameter :: p_step_right = 0.5_dp
  integer,  parameter :: seed0   = 2025

  real(dp), parameter :: ksigma = 6.0_dp

  logical,  parameter :: OUTPUT_CHECKPOINTS = .true.
  integer,  parameter :: Kout = 5
  integer,  parameter :: n_out(Kout) = (/ 100, 500, 1000, 2000, 5000 /)

  ! ---------------- state ----------------
  integer :: i, n, seed, dir
  integer :: x
  real(dp) :: r

  real(dp), allocatable :: sumx(:), sumx2(:)

  real(dp) :: mu, sigma
  integer  :: xmin, xmax, nbins, idx
  integer, allocatable :: hcount(:)

  integer, allocatable :: hcount_out(:,:)
  integer, allocatable :: step_index(:)

  real(dp) :: grnd
  external  :: grnd

  character(len=64) :: tag

  ! ---------------- init ----------------
  seed = seed0
  call sgrnd(seed)

  allocate(sumx(0:nsteps), sumx2(0:nsteps))
  sumx  = 0.0_dp
  sumx2 = 0.0_dp

  mu    = (2.0_dp*p_step_right - 1.0_dp) * real(nsteps,dp)
  sigma = sqrt( max(0.0_dp, 4.0_dp*p_step_right*(1.0_dp - p_step_right) * real(nsteps,dp)) )

  xmin  = int( floor(mu - ksigma*sigma) ) - 2
  xmax  = int( ceiling(mu + ksigma*sigma) ) + 2
  if (xmin > xmax) then
     xmin = -2; xmax = 2
  end if

  nbins = xmax - xmin + 1
  allocate(hcount(nbins)); hcount = 0

  if (OUTPUT_CHECKPOINTS) then
     allocate(hcount_out(Kout, nbins)); hcount_out = 0
     allocate(step_index(0:nsteps)); step_index = -1
  end if

  if (OUTPUT_CHECKPOINTS) then
     do i = 1, Kout
        if (n_out(i) >= 0 .and. n_out(i) <= nsteps) step_index(n_out(i)) = i
     end do
  end if

  write(tag,'("p",I0.2,"_nsteps",I0)') int(nint(p_step_right*100.0_dp)), nsteps

  ! ---------------- ensemble loop ----------------
  do i = 1, ntraj
     x = 0

     sumx(0)  = sumx(0)  + 0.0_dp
     sumx2(0) = sumx2(0) + 0.0_dp
     if (OUTPUT_CHECKPOINTS) then
        if (step_index(0) /= -1) then
           idx = 0 - xmin + 1
           if (idx < 1) idx = 1
           if (idx > nbins) idx = nbins
           hcount_out(step_index(0), idx) = hcount_out(step_index(0), idx) + 1
        end if
     end if

     do n = 1, nsteps
        r = grnd()
        if (r < p_step_right) then
           dir = +1
        else
           dir = -1
        end if
        x = x + dir

        sumx(n)  = sumx(n)  + real(x, dp)
        sumx2(n) = sumx2(n) + real(x*x, dp)

        if (OUTPUT_CHECKPOINTS) then
           if (step_index(n) /= -1) then
              idx = x - xmin + 1
              if (idx < 1) idx = 1
              if (idx > nbins) idx = nbins
              hcount_out(step_index(n), idx) = hcount_out(step_index(n), idx) + 1
           end if
        end if
     end do

     idx = x - xmin + 1
     if (idx < 1)   idx = 1
     if (idx > nbins) idx = nbins
     hcount(idx) = hcount(idx) + 1
  end do

  call write_msd_centered("msd_dsdt_"//trim(tag)//".txt", nsteps, ntraj, p_step_right, sumx, sumx2)

  call write_counts_int("xCounts_dsdt_"//trim(tag)//".txt", xmin, nbins, hcount)
  call write_pdf_int    ("xPDF_dsdt_"   //trim(tag)//".txt", xmin, nbins, hcount, 1.0_dp)
  call write_pdf_exact  ("xPDF_exact_"  //trim(tag)//".txt", nsteps, p_step_right, xmin, xmax)

  if (OUTPUT_CHECKPOINTS) then
     do i = 1, Kout
        call write_counts_int("xCounts_dsdt_"//trim(tag)//"_n"//i4(n_out(i))//".txt", xmin, nbins, hcount_out(i,:))
        call write_pdf_int   ("xPDF_dsdt_"   //trim(tag)//"_n"//i4(n_out(i))//".txt", xmin, nbins, hcount_out(i,:), 1.0_dp)
        call write_pdf_exact ("xPDF_exact_"  //trim(tag)//"_n"//i4(n_out(i))//".txt", n_out(i), p_step_right, xmin, xmax)
     end do
  end if

  deallocate(sumx, sumx2, hcount)
  if (OUTPUT_CHECKPOINTS) then
     deallocate(hcount_out, step_index)
  end if

contains
  pure function i4(n) result(s)
    integer, intent(in) :: n
    character(len=8) :: s
    write(s,'(I0.4)') n
  end function i4

  subroutine write_msd_centered(fname, nsteps, ntraj, p, sumx, sumx2)
    character(*), intent(in) :: fname
    integer,      intent(in) :: nsteps, ntraj
    real(dp),     intent(in) :: p
    real(dp),     intent(in) :: sumx(0:nsteps), sumx2(0:nsteps)
    integer :: u, n
    real(dp) :: meanx, msd, meanx_exact, var_exact
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# n   <x>(n)   MSD_MC(n)   <x>_exact(n)   MSD_exact(n)"
    do n = 0, nsteps
       meanx = sumx(n) / real(ntraj,dp)
       msd   = sumx2(n)/ real(ntraj,dp) - meanx*meanx
       meanx_exact = (2.0_dp*p - 1.0_dp) * real(n,dp)
       var_exact   = 4.0_dp*p*(1.0_dp - p) * real(n,dp)
       write(u,'(i10,1x,es20.10,1x,es20.10,1x,es20.10,1x,es20.10)') n, meanx, msd, meanx_exact, var_exact
    end do
    close(u)
  end subroutine write_msd_centered

  subroutine write_counts_int(fname, xmin, nbins, h)
    character(*), intent(in) :: fname
    integer,      intent(in) :: xmin, nbins
    integer,      intent(in) :: h(nbins)
    integer :: u, i, x
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# x   counts"
    do i = 1, nbins
       x = xmin + (i-1)
       write(u,'(i12,1x,i12)') x, h(i)
    end do
    close(u)
  end subroutine write_counts_int

  subroutine write_pdf_int(fname, xmin, nbins, h, bw)
    character(*), intent(in) :: fname
    integer,      intent(in) :: xmin, nbins
    integer,      intent(in) :: h(nbins)
    real(dp),     intent(in) :: bw
    integer :: u, i, x, total
    real(dp) :: norm
    total = 0
    do i=1,nbins
       total = total + h(i)
    end do
    if (total <= 0) return
    norm = real(total,dp) * bw
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# x   P_MC(x,n)"
    do i = 1, nbins
      x = xmin + (i-1)
      write(u,'(i12,1x,es20.10)') x, real(h(i),dp)/norm
    end do
    close(u)
  end subroutine write_pdf_int

  subroutine write_pdf_exact(fname, n, p, xmin, xmax)
    character(*), intent(in) :: fname
    integer,      intent(in) :: n, xmin, xmax
    real(dp),     intent(in) :: p
    integer :: u, x
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# x   P_exact(x,n)"
    do x = xmin, xmax
       write(u,'(i12,1x,es20.10)') x, p_exact_x_given_n_p(x, n, p)
    end do
    close(u)
  end subroutine write_pdf_exact

  pure function p_exact_x_given_n_p(x, n, p) result(px)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: p
    real(dp) :: px
    integer :: k
    if (abs(x) > n .or. mod(n+x,2) /= 0) then
       px = 0.0_dp
    else
       k  = (n + x)/2
       ! LOG_GAMMA is intrinsic; just call it (no declarations).
       px = exp( log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - &
                 log_gamma(real(n-k+1,dp)) + real(k,dp)*log(p) + &
                 real(n-k,dp)*log(1.0_dp-p) )
    end if
  end function p_exact_x_given_n_p

end program rw_1D_DSDT_ensemble

! ======= RNG (Mersenne Twister): provides grnd() and sgrnd() =======
include 'mt.f90'

