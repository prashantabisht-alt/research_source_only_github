program Brownian_motion
  implicit none
  integer :: i, j, timestep, Tmax, seed
  real*8 :: gaussianvariable, eta, position, diff, dt
  character(len=100) :: filename

  ! Initialize random seed
  seed = 1002
  call sgrnd(seed)

  ! Parameters
  dt = 0.001d0
  diff = 1.0d0
  Tmax = 100000

  ! Loop over trajectories
  do j = 1, 1000
    print*, "Trajectory:", j

    position = 0.0d0
    write(filename, "('Brownian_motion_dt0.001_traj', I4.4, '.txt')") j
    !open(unit=42, file=filename)

    open(unit=42, file="Brownian_motion_dt0.001.txt")

    do timestep = 1, Tmax
      call gaussian(gaussianvariable)
      eta = gaussianvariable * dsqrt(2.0d0 * diff * dt)
      position = position + eta
      write(42,*) timestep * dt, position
    end do

    close(42)
  end do

end program Brownian_motion

include 'mt.f90'

!=============================================================================================
!This subroutine generates a unit gaussian random variable using the polar Box-Muller transform
!with distribution 1/sqrt(2 pi) exp (-x^2/2)
  
  subroutine gaussian(s)
  real*8 x1,x2,w,s,grnd
  
  w = 2.0d0
  do while(w.gt.1.0d0)  
    x1 = 1.0d0 - 2.0d0*grnd()
    x2 = 1.0d0 - 2.0d0*grnd()
    w = x1 * x1 + x2 * x2
  enddo
  
  s = dsqrt(-2.0d0*(dlog(w)/w))*x1
  
  end
!========================================
!=============================================================================================
! Write normalized PDF to file
subroutine write_pdf(hist, nbins, xmin, xmax, ntraj)
  implicit none
  integer, intent(in) :: nbins, ntraj
  real*8, intent(in) :: hist(nbins), xmin, xmax
  real*8 :: bin_width, norm, x
  integer :: i
  open(unit=43, file="pdf.txt")

  bin_width = (xmax - xmin)/nbins
  norm = sum(hist) * bin_width

  do i = 1, nbins
    x = xmin + (i - 0.5d0) * bin_width
    write(43,*) x, hist(i) / norm
  end do

  close(43)
end subroutine write_pdf
!=============================================================================================
