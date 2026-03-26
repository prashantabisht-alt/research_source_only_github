program Brownian_motion
implicit none

integer i,j,k,seed,timestep,Tmax
real*8 gaussianvariable,eta,grnd,s,position,diff,dt


!******************
seed = 1002
call sgrnd(seed)	
!******************

!Parameters
!******************
dt = 0.001d0
diff = 1.0d0
Tmax =100000
!******************


!****PRINT SIMPLE RANDOM NUMBERS********
!do i =1,100
!	print*, grnd()
!enddo
!***************************************


!open(unit=42,file="random_numbers.txt")
!****PRINT GAUSSIAN RANDOM NUMBERS******
!do i = 1,10000
!		call gaussian(s)
!		!s = grnd()
!		print*, s
!		write(42,*) i,s
!	enddo
!***************************************
!close(42)


do j = 1,1000

 print*,j

position = 0.0d0

open(unit=42,file="Brownian_motion_dt0.001.txt")

do timestep = 1,Tmax

	call gaussian(gaussianvariable)

	eta = gaussianvariable*dsqrt(2.0d0*diff/dt)	

	position =  position + eta*dt

	write(42,*) timestep*dt, position
	
enddo
close(42)

enddo

end

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
!===============================================================================================