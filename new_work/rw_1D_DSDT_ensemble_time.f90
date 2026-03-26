!===============================================================================
! 1D Discrete-Space, Discrete-Time Random Walk (Ensemble) — with physical time t
! - Each tick takes dt_tick units of physical time → t = n * dt_tick
! - Step length = 1, step probabilities: P(+1)=p_step_right, P(-1)=1-p
!
! Outputs (tagged with p and nsteps; files also carry t in names/headers):
!   1) msd_dsdt_<tag>.txt
!      # t   <x>(t)   MSD_MC(t)   <x>_exact(t)   MSD_exact(t)
!   2) xCounts_dsdt_<tag>_t<ttag>.txt   (x, counts)   at final t=N*dt_tick
!   3) xPDF_dsdt_<tag>_t<ttag>.txt      (x, P_MC(x,t))
!   4) xPDF_exact_<tag>_t<ttag>.txt     (x, P_exact(x,t))
!   5) (optional) same for checkpoints t_out(:)
!
! Compile:  gfortran -O3 rw_1D_DSDT_ensemble_time.f90 -o rw_walk_t
!===============================================================================
program rw_1D_DSDT_ensemble_time
  implicit none
  integer,  parameter :: dp = selected_real_kind(15,300)

  ! ---------------- knobs ----------------
  integer,  parameter :: ntraj   = 100000          ! ensemble size
  integer,  parameter :: nsteps  = 10000           ! number of ticks
  real(dp), parameter :: dt_tick = 1.0_dp          ! physical time per tick
  real(dp), parameter :: p_step_right = 0.5_dp     ! bias: P(+1)=p
  integer,  parameter :: seed0   = 2025            ! RNG seed
  real(dp), parameter :: ksigma  = 6.0_dp          ! histogram half-span in sigmas

  ! Optional: output distributions at selected physical times
  logical,  parameter :: OUTPUT_CHECKPOINTS = .true.
  integer,  parameter :: Kout = 5
  real(dp), parameter :: t_out(Kout) = (/ 100.0_dp, 500.0_dp, 1000.0_dp, 2000.0_dp, 5000.0_dp /)

  ! ---------------- state ----------------
  integer :: i, n, seed, dir
  integer :: x
  real(dp) :: r

  real(dp), allocatable :: sumx(:), sumx2(:)       ! accumulate over trajectories for each n (0..nsteps)
  real(dp), allocatable :: tgrid(:)                ! t(n) = n * dt_tick

  ! histogram setup using final n = nsteps
  real(dp) :: mu_steps, sigma_steps
  integer  :: xmin, xmax, nbins, idx
  integer,  allocatable :: hcount(:)               ! final-step counts

  ! checkpoint handling
  integer,  allocatable :: n_out(:)                ! nearest steps for each t_out
  integer,  allocatable :: step_index(:)           ! map 0..nsteps -> {1..Kout} or -1
  integer,  allocatable :: hcount_out(:,:)         ! (Kout, nbins)

  ! RNG
  real(dp) :: grnd
  external  :: grnd

  ! tags
  character(len=64) :: tag
  character(len=16) :: tfin_tag

  ! ---------------- init ----------------
  seed = seed0
  call sgrnd(seed)

  allocate(sumx(0:nsteps), sumx2(0:nsteps), tgrid(0:nsteps))
  sumx  = 0.0_dp
  sumx2 = 0.0_dp
  do n = 0, nsteps
     tgrid(n) = real(n,dp) * dt_tick
  end do

  ! mean & sigma in "step units" for final n; used only for histogram span
  mu_steps    = (2.0_dp*p_step_right - 1.0_dp) * real(nsteps,dp)
  sigma_steps = sqrt( max(0.0_dp, 4.0_dp*p_step_right*(1.0_dp - p_step_right) * real(nsteps,dp)) )

  xmin  = int( floor(mu_steps - ksigma*sigma_steps) ) - 2
  xmax  = int( ceiling(mu_steps + ksigma*sigma_steps) ) + 2
  if (xmin > xmax) then
     xmin = -2; xmax = 2
  end if
  nbins = xmax - xmin + 1
  allocate(hcount(nbins)); hcount = 0

  ! checkpoint mapping t_out -> nearest step indices (bounded)
  if (OUTPUT_CHECKPOINTS) then
     allocate(n_out(Kout))
     allocate(step_index(0:nsteps)); step_index = -1
     allocate(hcount_out(Kout, nbins)); hcount_out = 0
     do i = 1, Kout
        n_out(i) = max(0, min(nsteps, int( nint( t_out(i) / dt_tick ) )))
        step_index(n_out(i)) = i
     end do
  end if

  ! filename tag
  write(tag,'("p",I0.2,"_nsteps",I0)') int(nint(p_step_right*100.0_dp)), nsteps
  tfin_tag = ttag( tgrid(nsteps) )

  ! ---------------- ensemble loop ----------------
  do i = 1, ntraj
     x = 0

     ! n=0
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

     ! final histogram at t = nsteps*dt_tick
     idx = x - xmin + 1
     if (idx < 1)   idx = 1
     if (idx > nbins) idx = nbins
     hcount(idx) = hcount(idx) + 1
  end do

  ! ---------------- write MSD (with exact vs t) ----------------
  call write_msd_time("msd_dsdt_"//trim(tag)//".txt", nsteps, ntraj, dt_tick, p_step_right, tgrid, sumx, sumx2)

  ! ---------------- write final PDFs/Counts (labeled by t) -----
  call write_counts_int("xCounts_dsdt_"//trim(tag)//"_t"//trim(tfin_tag)//".txt", xmin, nbins, hcount)
  call write_pdf_int   ("xPDF_dsdt_"   //trim(tag)//"_t"//trim(tfin_tag)//".txt", xmin, nbins, hcount, 1.0_dp)
  call write_pdf_exact ("xPDF_exact_"  //trim(tag)//"_t"//trim(tfin_tag)//".txt", nsteps, p_step_right, xmin, xmax)

  ! ---------------- checkpoints (optional) ---------------------
  if (OUTPUT_CHECKPOINTS) then
     do i = 1, Kout
        call write_counts_int("xCounts_dsdt_"//trim(tag)//"_t"//trim(ttag(t_out(i)))//".txt", xmin, nbins, hcount_out(i,:))
        call write_pdf_int   ("xPDF_dsdt_"   //trim(tag)//"_t"//trim(ttag(t_out(i)))//".txt", xmin, nbins, hcount_out(i,:), 1.0_dp)
        call write_pdf_exact ("xPDF_exact_"  //trim(tag)//"_t"//trim(ttag(t_out(i)))//".txt", n_out(i), p_step_right, xmin, xmax)
     end do
  end if

  deallocate(sumx, sumx2, tgrid, hcount)
  if (OUTPUT_CHECKPOINTS) then
     deallocate(n_out, step_index, hcount_out)
  end if

contains
  !---------------- format helper: "0500p000" for t=500.000 ----------
  pure function ttag(t) result(s)
    real(dp), intent(in) :: t
    character(len=16) :: s, tmp
    integer :: j
    write(tmp,'(F0.3)') t           ! e.g., "500.000"
    s = tmp
    do j=1,len_trim(s)
       if (s(j:j) == ' ') s(j:j) = '0'
    end do
    ! replace '.' with 'p' to be filename-safe and consistent with your style
    do j=1,len_trim(s)
       if (s(j:j) == '.') s(j:j) = 'p'
    end do
    ! zero-pad to width 7-ish if you like; current string is fine
  end function ttag

  !---------------- MSD writer versus t ------------------------------
  subroutine write_msd_time(fname, nsteps, ntraj, dt, p, tgrid, sumx, sumx2)
    character(*), intent(in) :: fname
    integer,      intent(in) :: nsteps, ntraj
    real(dp),     intent(in) :: dt, p
    real(dp),     intent(in) :: tgrid(0:nsteps), sumx(0:nsteps), sumx2(0:nsteps)
    integer :: u, n
    real(dp) :: t, meanx, msd, meanx_exact, var_exact
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# t   <x>(t)   MSD_MC(t)   <x>_exact(t)   MSD_exact(t)"
    do n = 0, nsteps
       t = tgrid(n)
       meanx = sumx(n) / real(ntraj,dp)
       msd   = sumx2(n)/ real(ntraj,dp) - meanx*meanx
       ! exact formulas: replace n = t/dt
       meanx_exact = (2.0_dp*p - 1.0_dp) * (t/dt)
       var_exact   = 4.0_dp*p*(1.0_dp - p) * (t/dt)
       write(u,'(es20.10,1x,es20.10,1x,es20.10,1x,es20.10,1x,es20.10)') t, meanx, msd, meanx_exact, var_exact
    end do
    close(u)
  end subroutine write_msd_time

  !---------------- integer-lattice counts: (x, count) --------------
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

  !---------------- normalized PDF (bin width = 1) -------------------
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
    write(u,'(a)') "# x   P_MC(x,t)"
    do i = 1, nbins
       x = xmin + (i-1)
       write(u,'(i12,1x,es20.10)') x, real(h(i),dp)/norm
    end do
    close(u)
  end subroutine write_pdf_int

  !---------------- exact PDF for given step n (used at time t=n*dt) -
  subroutine write_pdf_exact(fname, n, p, xmin, xmax)
    character(*), intent(in) :: fname
    integer,      intent(in) :: n, xmin, xmax
    real(dp),     intent(in) :: p
    integer :: u, x
    open(newunit=u, file=fname, status="replace")
    write(u,'(a)') "# x   P_exact(x,t)   # t = n*dt, n="//trim(adjustl(i2s(n)))
    do x = xmin, xmax
       write(u,'(i12,1x,es20.10)') x, p_exact_x_given_n_p(x, n, p)
    end do
    close(u)
  end subroutine write_pdf_exact

  !---------------- exact mass function P(x,n) (biased) -------------
  pure function p_exact_x_given_n_p(x, n, p) result(px)
    integer, intent(in) :: x, n
    real(dp), intent(in) :: p
    real(dp) :: px
    integer :: k
    if (abs(x) > n .or. mod(n+x,2) /= 0) then
       px = 0.0_dp
    else
       k  = (n + x)/2
       ! LOG_GAMMA is intrinsic; call directly
       px = exp( log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - &
                 log_gamma(real(n-k+1,dp)) + real(k,dp)*log(p) + &
                 real(n-k,dp)*log(1.0_dp-p) )
    end if
  end function p_exact_x_given_n_p

  !---------------- tiny helpers ------------------------------------
  pure function i2s(n) result(s)
    integer, intent(in) :: n
    character(len=16) :: s
    write(s,'(I0)') n
  end function i2s

end program rw_1D_DSDT_ensemble_time

! ======= RNG (Mersenne Twister) ====================================
include 'mt.f90'
