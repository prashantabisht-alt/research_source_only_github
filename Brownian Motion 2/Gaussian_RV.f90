	program GaussianRV
	implicit none
	!=================================================================
	
	real*8 grnd,s
	integer i

	do i = 1,100
		call gaussian(s)
		print*, 1/s**2
	enddo

	
	end
	include 'mt.f90'	
	!========================================================================================================	
	!This subroutine generates a unit gaussian random variable using the polar Box-Muller transform
	
	subroutine gaussian(s)
	real*8 x1,x2,w,s,grnd
	
	w = 2.0d0
	do while(w.gt.1.0d0)  
		x1 = 1.0d0 - 2.0d0*grnd()
		x2 = 1.0d0 - 2.0d0*grnd()
		w = x1 * x1 + x2 * x2
	enddo
	
	s = sqrt(-2.0d0*(log(w)/w))*x1
	
	end
	!========================================================================================================