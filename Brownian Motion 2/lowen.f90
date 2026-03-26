!===============================================================
! 
! Goal: Generate u(t) from AOUP SDE (Eq. 5) and verify Eqs. 6–8
!===============================================================
program AOUP_u_statistics
    implicit none
    
    integer,  parameter :: Nsteps = 2000000      ! total time steps
    real*8,   parameter :: dt     = 0.001d0      ! time step
    real*8,   parameter :: taup   = 0.50d0       ! persistence time τ_p
    real*8,   parameter :: v0     = 1.00d0       ! stationary std (self-propulsion speed scale)
    real*8,   parameter :: gamma  = 0.8d0        ! friction coeff (only used to build M(t))
    integer,  parameter :: Kfrac  = 20           ! max lag = Nsteps/Kfrac (keep memory & speed ok)
    integer,  parameter :: seed   = 1002         ! RNG seed

    
    integer :: Kmax
    real*8  :: a, b, D_theory

    
    real*8, allocatable :: u(:)
    real*8, allocatable :: C_est(:), C_theory(:)
    real*8, allocatable :: M_est(:), M_theory(:)

    
    integer :: n, k
    real*8  :: um, var, s, tlag, sumC, D_est
    real*8  :: tmp

    
    call sgrnd(seed)     
    a = dexp(-dt/taup)   
    b = v0*dsqrt(1.0d0 - a*a)   
    D_theory = v0*v0*taup

     allocate(u(Nsteps))
    Kmax = Nsteps / Kfrac
    allocate(C_est(0:Kmax), C_theory(0:Kmax))
    allocate(M_est(0:Kmax), M_theory(0:Kmax))

    ! ---------- Generate u(t) ----------
    u(1) = 0.0d0     ! start in mean; stationarity builds quickly
    do n = 2, Nsteps
        call gaussian(s)         
        u(n) = a*u(n-1) + b*s    ! exact OU update (preserves var = v0^2)
    end do

    ! ---------- Mean & variance ----------
    um = 0.0d0
    do n = 1, Nsteps
        um = um + u(n)
    end do
    um = um / dble(Nsteps)

    var = 0.0d0
    do n = 1, Nsteps
        tmp = u(n) - um
        var = var + tmp*tmp
    end do
    var = var / dble(Nsteps)

    ! ---------- Autocorrelation C[l] ----------
    do k = 0, Kmax
        C_est(k) = 0.0d0
        do n = 1, Nsteps - k
            C_est(k) = C_est(k) + (u(n)-um)*(u(n+k)-um)
        end do
        C_est(k) = C_est(k) / dble(Nsteps - k)
        C_theory(k) = v0*v0 * a**k           ! Eq. (6) prediction
    end do

    ! ---------- Kernel M[l] and white-noise D check ----------
    do k = 0, Kmax
        M_est(k)    = gamma*gamma * C_est(k)        ! Eq. (7)
        M_theory(k) = gamma*gamma * C_theory(k)
    end do

    sumC = 0.0d0
    do k = 0, Kmax
        if (k == 0) then
            sumC = sumC + 0.5d0*C_est(k)   ! trapezoid: half weight at ends
        else if (k == Kmax) then
            sumC = sumC + 0.5d0*C_est(k)
        else
            sumC = sumC + C_est(k)
        end if
    end do
    D_est = dt * sumC            ! should be ~ v0^2 * taup

    ! ---------- Write outputs ----------
    open(unit=11, file="u_stats.txt")
    write(11,'(A,1X,ES20.12)') "# mean_est", um
    write(11,'(A,1X,ES20.12)') "# var_est",  var
    write(11,'(A,1X,ES20.12)') "# v0^2_theory", v0*v0
    write(11,'(A,1X,ES20.12)') "# D_est (∫C dt)", D_est
    write(11,'(A,1X,ES20.12)') "# D_theory (v0^2*taup)", D_theory
    close(11)

    open(unit=12, file="u_autocorr.txt")
    write(12,*) "# tau   C_est   C_theory"
    do k = 0, Kmax
        tlag = dble(k)*dt
        write(12,'(3(1X,ES20.12))') tlag, C_est(k), C_theory(k)
    end do
    close(12)

    open(unit=13, file="M_kernel.txt")
    write(13,*) "# tau   M_est   M_theory"
    do k = 0, Kmax
        tlag = dble(k)*dt
        write(13,'(3(1X,ES20.12))') tlag, M_est(k), M_theory(k)
    end do
    close(13)

    print *,"Done. Files: u_stats.txt, u_autocorr.txt, M_kernel.txt"
    print *,"Check: var≈v0^2, C(τ)≈v0^2*exp(-τ/τp), D≈v0^2*taup."

    deallocate(u, C_est, C_theory, M_est, M_theory)
end program AOUP_u_statistics

! ------------------------------------------------------------
! RNG & Gaussian (your style)
include 'mt.f90'

! Box–Muller transform to get N(0,1) using your grnd()
subroutine gaussian(s)
    implicit none
    real*8 :: x1, x2, w, s, grnd
    w = 2.0d0
    do while (w > 1.0d0 .or. w == 0.0d0)
        x1 = 1.0d0 - 2.0d0*grnd()
        x2 = 1.0d0 - 2.0d0*grnd()
        w = x1*x1 + x2*x2
    end do
    s = dsqrt(-2.0d0 * dlog(w)/w) * x1
end subroutine gaussian
