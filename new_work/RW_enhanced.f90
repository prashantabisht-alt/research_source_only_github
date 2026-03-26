program ctrw_single_particle_enhanced
  implicit none

  ! --------- params ----------
  integer,  parameter :: L = 1000
  integer,  parameter :: nsteps = 1000000
  integer,  parameter :: nbins  = 200
  real*8,   parameter :: lambda = 1.0d0
  real*8,   parameter :: dtmin  = 0.0d0
  real*8,   parameter :: dtmax  = 10.0d0

  ! --------- state ----------
  integer :: step, pos, dir, seed, bin_index
  real*8  :: t, dt, r1, r2, bin_width
  real*8, allocatable :: dt_hist(:)
  real*8  :: grnd
  external :: grnd

  ! --------- files ----------
  integer, parameter :: U_TRAJ=42, U_DTS=43

  ! --------- open outputs ----------
  open(unit=U_TRAJ, file="ctrw_single_particle.txt", status="replace")
  open(unit=U_DTS,  file="waiting_times.txt",       status="replace")

  ! --------- init ----------
  pos  = L/2
  t    = 0.0d0
  seed = 2025
  call sgrnd(seed)

  bin_width = (dtmax - dtmin)/nbins
  allocate(dt_hist(nbins)); dt_hist = 0.0d0

  ! --------- simulate ----------
  do step = 1, nsteps
     ! waiting time ~ Exp(lambda)
     r1 = grnd()
     if (r1 <= 0.0d0) r1 = 1.0d-308          ! safety
     dt = -log(r1)/lambda
     t  = t + dt
     write(U_DTS,*) dt

     ! bin Δt
     if (dt >= dtmin .and. dt < dtmax) then
        bin_index = int((dt - dtmin)/bin_width) + 1
        if (bin_index < 1) bin_index = 1
        if (bin_index > nbins) bin_index = nbins
        dt_hist(bin_index) = dt_hist(bin_index) + 1.0d0
     end if

     ! ±1 jump
     r2 = grnd()
     if (r2 < 0.5d0) then
        dir = 1
     else
        dir = -1
     end if

     ! periodic boundary
     pos = mod(pos + dir + L, L)

     ! log trajectory (time, position)
     write(U_TRAJ,*) t, pos
  end do

  close(U_TRAJ); close(U_DTS)

  ! --------- write histogram outputs ----------
  call write_pdf(dt_hist, nbins, dtmin, dtmax, "waiting_time_hist_pdf.txt")
  call write_counts(dt_hist, nbins, dtmin, dtmax, "waiting_time_hist_counts.txt")

  deallocate(dt_hist)

end program ctrw_single_particle_enhanced


! ========= RNG (Mersenne Twister) =========
include 'mt.f90'


! ========= write normalized PDF (area=1 over [dtmin,dtmax)) =========
subroutine write_pdf(hist, nbins, dtmin, dtmax, fname)
  implicit none
  integer, intent(in) :: nbins
  real*8, intent(in)  :: hist(nbins), dtmin, dtmax
  character(*), intent(in) :: fname
  real*8 :: bw, norm, x
  integer :: i, u

  bw   = (dtmax - dtmin)/nbins
  norm = sum(hist)*bw
  if (norm <= 0.0d0) return

  u = 69
  open(unit=u, file=fname, status="replace")
  do i = 1, nbins
     x = dtmin + (i - 0.5d0)*bw
     write(u,*) x, hist(i)/norm
  end do
  close(u)
end subroutine write_pdf


! ========= write raw counts (bin_center, count) =========
subroutine write_counts(hist, nbins, dtmin, dtmax, fname)
  implicit none
  integer, intent(in) :: nbins
  real*8, intent(in)  :: hist(nbins), dtmin, dtmax
  character(*), intent(in) :: fname
  real*8 :: bw, x
  integer :: i, u

  bw = (dtmax - dtmin)/nbins
  u = 68
  open(unit=u, file=fname, status="replace")
  do i = 1, nbins
     x = dtmin + (i - 0.5d0)*bw
     write(u,*) x, hist(i)
  end do
  close(u)
end subroutine write_counts

