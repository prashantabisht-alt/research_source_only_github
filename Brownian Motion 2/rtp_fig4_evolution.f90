!============================================================================
! RTP Fig. 4: Time evolution P(x,t) in a finite interval (ell=1) with
! stick-until-flip walls, D=0.1, IC: x(0)=0, sigma(0)=±1 equiprobable.
! Writes one file per snapshot time: fig4_D0p10_t<tp>ptxt
! Columns: x   P_MC(x,t)   P_ss(x)   (so you can overlay directly)
! Also writes: fig4_curvature_vs_t_D0p10ptxt  (t, d2P(0,t))
!              fig4_L1_to_ss_D0p10ptxt       (t, L1 distance to P_ss)
!============================================================================
program rtp_fig4_evolution
  implicit none
  ! ---------------- user knobs ----------------
  integer,  parameter :: Ntraj  = 400000      ! raise to 1e6 for super-smooth lines
  real*8,   parameter :: D      = 0.10d0
  real*8,   parameter :: ell    = 1.0d0
  real*8,   parameter :: Tfin   = 2.00d0      ! last snapshot (matches paper)
  real*8,   parameter :: dt     = 1.0d-3      ! reduce (e.g. 5e-4) if needed at walls
  integer,  parameter :: nbins  = 401         ! odd → center bin at x=0
  ! snapshots (t=0.1:0.1:2.0)
  integer,  parameter :: nsnaps = 20
  ! -------------------------------------------

  real*8  :: tgrid(nsnaps)
  integer :: k, steps_total, steps_to_k, n, s, idx, idx0
  real*8  :: binw, xmin, xmax, x, sigma, pflip, w, next_t
  integer, allocatable :: hist(:,:)      ! (nbins, nsnaps)
  real*8,  allocatable :: Pss(:), xgrid(:)
  character(len=64) :: fname
  integer :: uout, ucurv, uerr

  ! RNG
  integer :: seed
  real*8, external :: grnd
  seed = 24681357
  call sgrnd(seed)

  ! time grid
  do k=1, nsnaps
     tgrid(k) = 0.1d0 * dble(k)
  end do
  steps_total = nint(Tfin/dt)

  ! bins / grid
  xmin = -ell; xmax = ell
  binw = (xmax - xmin) / dble(nbins)
  allocate(xgrid(nbins), Pss(nbins))
  do k=1, nbins
    xgrid(k) = xmin + (dble(k)-0.5d0)*binw
    Pss(k)   = Pss_exact(D, ell, xgrid(k))   ! analytic steady state (Eq. 16)
  end do
  idx0 = int((0.0d0 - xmin)/binw) + 1
  if (mod(nbins,2)==0 .or. idx0<2 .or. idx0>nbins-1) then
     print *,"Choose odd nbins; center bin must be interior."; stop
  end if

  allocate(hist(nbins, nsnaps)); hist = 0

  ! ==== Monte Carlo ensemble ====
  pflip = 1.0d0 - dexp(-dt)    ! unit flip rate

  do n = 1, Ntraj
     ! IC
     x = 0.0d0
     if (grnd() < 0.5d0) then; sigma = -1.0d0; else; sigma = 1.0d0; end if

     s = 0
     do k = 1, nsnaps
        steps_to_k = nint(tgrid(k)/dt) - s
        do while (steps_to_k > 0)
          ! tumble?
          if (grnd() < pflip) sigma = -sigma

          ! ballistic piece with "stick" outward
          if (x >= ell-1d-12 .and. sigma>0.0d0) then
             x = ell
          else if (x <= -ell+1d-12 .and. sigma<0.0d0) then
             x = -ell
          else
             x = x + sigma*dt
             if (x > ell)  x = ell
             if (x < -ell) x = -ell
          end if

          ! Brownian piece: inward-only when sitting at wall
          call gaussian(w)
          w = dsqrt(2.0d0*D*dt) * w
          if (x >= ell-1d-12 .and. w > 0.0d0)  w = 0.0d0
          if (x <= -ell+1d-12 .and. w < 0.0d0) w = 0.0d0
          x = x + w
          if (x > ell)  x = ell
          if (x < -ell) x = -ell

          s = s + 1
          steps_to_k = steps_to_k - 1
        end do

        ! bin at snapshot k
        if (x > xmin .and. x < xmax) then
           idx = int((x - xmin)/binw) + 1
           if (idx >= 1 .and. idx <= nbins) hist(idx,k) = hist(idx,k) + 1
        end if
     end do

     if (mod(n,100000)==0) print *,"  ...", n, "trajectories"
  end do

  ! ==== write one file per time: x, P_MC(x,t), P_ss(x) ====
  do k = 1, nsnaps
     write(fname,'("fig4_D0p10_t",f0.1,".txt")') tgrid(k)
     call sanitize_filename(fname)
     open(newunit=uout, file=trim(fname), status="replace", action="write")
     write(uout,*) "# x   P_MC(x,t)   P_ss(x)   (D=0.1, ell=1, t=",tgrid(k),")"
     call dump_snapshot(uout, nbins, xgrid, hist(:,k), binw, Pss, Ntraj)
     close(uout)
  end do
  print *,"Wrote fig4_D0p10_t*.ptxt snapshots."

  ! ==== extra diagnostics: curvature at 0 vs t, and L1 distance to P_ss ====
  open(newunit=ucurv, file="fig4_curvature_vs_t_D0p10ptxt", status="replace")
  open(newunit=uerr,  file="fig4_L1_to_ss_D0p10ptxt",     status="replace")
  write(ucurv,*) "# t   d2P(0,t)   (5-pt stencil)"
  write(uerr, *) "# t   L1 = ∑|P_MC-P_ss| dx"

  do k = 1, nsnaps
     call write_curvature_line(ucurv, tgrid(k), hist(:,k), idx0, binw, Ntraj)
     call write_L1_line(uerr,    tgrid(k), hist(:,k), Pss, binw, Ntraj)
  end do
  close(ucurv); close(uerr)

contains

  subroutine dump_snapshot(uout, nbins, xgrid, histk, binw, Pss, Ntraj)
    implicit none
    integer, intent(in) :: uout, nbins, histk(nbins), Ntraj
    real*8, intent(in)  :: xgrid(nbins), binw, Pss(nbins)
    integer :: b
    real*8  :: Pmc
    do b = 1, nbins
       Pmc = dble(histk(b)) / (dble(Ntraj)*binw)
       write(uout,'(3(1X,ES20.12))') xgrid(b), Pmc, Pss(b)
    end do
  end subroutine dump_snapshot

  subroutine write_curvature_line(u, t, histk, idx0, h, Ntraj)
    implicit none
    integer, intent(in) :: u, histk(:), idx0, Ntraj
    real*8, intent(in)  :: t, h
    real*8 :: Pm2,Pm1,P0,Pp1,Pp2, d2P
    Pm2 = dble(histk(idx0-2)) / (dble(Ntraj)*h)
    Pm1 = dble(histk(idx0-1)) / (dble(Ntraj)*h)
    P0  = dble(histk(idx0  )) / (dble(Ntraj)*h)
    Pp1 = dble(histk(idx0+1)) / (dble(Ntraj)*h)
    Pp2 = dble(histk(idx0+2)) / (dble(Ntraj)*h)
    d2P = (-Pm2 + 16d0*Pm1 - 30d0*P0 + 16d0*Pp1 - Pp2)/(12d0*h*h)
    write(u,'(2(1X,ES20.12))') t, d2P
  end subroutine write_curvature_line

  subroutine write_L1_line(u, t, histk, Pss, h, Ntraj)
    implicit none
    integer, intent(in) :: u, histk(:), Ntraj
    real*8, intent(in)  :: t, Pss(:), h
    integer :: b, nb
    real*8  :: Pmc, L1
    nb = size(histk)
    L1 = 0.0d0
    do b=1, nb
      Pmc = dble(histk(b)) / (dble(Ntraj)*h)
      L1  = L1 + dabs(Pmc - Pss(b)) * h
    end do
    write(u,'(2(1X,ES20.12))') t, L1
  end subroutine write_L1_line

  ! ---- analytic steady state P_ss(x) (paper Eq. 16) ----
  function Pss_exact(D, ell, x) result(P)
    implicit none
    real*8, intent(in) :: D, ell, x
    real*8 :: P, mu, normA, denom
    mu    = dsqrt(2.0d0*D + 1.0d0) / D
    normA = 1.0d0 / ( dtanh(mu*ell)/dsqrt(2.0d0*D+1.0d0) + 2.0d0*ell )
    denom = 2.0d0*D*dcosh(mu*ell) + 1.0d0
    P = normA * dcosh(mu*x) / denom
  end function Pss_exact

  ! ---- RNG helpers ----
  subroutine gaussian(z)
    real*8, intent(out) :: z
    real*8 :: x1,x2,w
    real*8, external :: grnd
    w = 2.0d0
    do while (w>1.0d0 .or. w==0.0d0)
      x1 = 1.0d0 - 2.0d0*grnd()
      x2 = 1.0d0 - 2.0d0*grnd()
      w  = x1*x1 + x2*x2
    end do
    z = dsqrt(-2.0d0*dlog(w)/w) * x1
  end subroutine gaussian

  subroutine sanitize_filename(str)
    character(len=*), intent(inout) :: str
    integer :: i, L
    L = len_trim(str)
    do i=1,L
      if (str(i:i)=='.') str(i:i)='p'
      if (str(i:i)==' ') str(i:i)='_'
    end do
    if (L>=4) then
      if (str(L-3:L)==".txt") str(L-3:L)="ptxt"
    end if
  end subroutine sanitize_filename

end program rtp_fig4_evolution

! Mersenne Twister RNG (grnd())
include 'mt.f90'
