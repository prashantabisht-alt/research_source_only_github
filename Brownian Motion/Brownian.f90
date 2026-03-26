!----------------------------------------------------	
	module Particles
	
	implicit none
	
	integer N,T,Nbins	

	!max number of particles
	parameter(N=100000)
	!number of time steps
    parameter(T = 100)
	parameter(Nbins = 100000)
	
	real*8 ParticlePos(N)
         
	end module Particles
!----------------------------------------------------
!============================================================================================
	program BranchingWalks
	
	use Particles
	
	implicit none

	integer seed,i,j,k,m,timestep,totalsteps,rang,Nparticles,
     c	Nparticlestemp,bin_no,prog,time_bin,Tmax,diff_bin1,diff_bin2
	real*8 dx,dt,diff,branch,annih,gaussianvariable,eta,variable,
     c	grnd,rand,cumprob1,cumprob2,cumprob3,totalbin,sm,sm1,sm2,
     c	span,span_upto_t,span_at_t,probspan1(0:Nbins),probspan2(0:Nbins),
     c	Ntot1,Ntot2,max_pos,min_pos
     
	seed = 1000
	call sgrnd(seed)	
	
	!Parameters
	!****************************************
	dt = 0.1d0
	diff = 1.0d0
	branch = 1.0d0
	annih = 1.0d0
	Tmax = T
	!****************************************

	
	!do k = 1,8
	  !rand = k
	  !print*, k,2.0d0**k
	   
	  !do m = 1,T
	    !total_bin(m,k) = 0.0d0
	  !enddo
	  do i = 0,Nbins
	    probspan1(i) = 0.0d0
	    probspan2(i) = 0.0d0
	  enddo
	  Ntot1 = 0.0d0
	  Ntot2 = 0.0d0
	!enddo

	
	do prog = 1,1000
	print*, prog
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	do totalsteps = 1,1000
	!print*, totalsteps
	!print*, totalsteps
	!##########################################################################	
	
	!Time Step zero
	!*********PSI*******************
	NParticles = 1
	ParticlePos(1) = 0.0d0
	span_upto_t = 0.0d0
	max_pos = 0.0d0
	min_pos = 0.0d0
	!******************************

	do timestep = 1,Tmax
	!Time step 1
	!**************************************************************************
	
	Nparticlestemp = Nparticles

	do i = 1,Nparticles
	  !Updating Particles
	  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	  cumprob1 = annih*dt
	  cumprob2 = annih*dt + branch*dt
	  cumprob3 = 1.0d0
	
	  rand = grnd()	
	  !print*,rand
	
	  if (rand.lt.cumprob1) then
	    !ANNIHILATION
	    Nparticlestemp = Nparticlestemp - 1
	    do j = i+1,Nparticlestemp+1
	       ParticlePos(j-1) = ParticlePos(j)
	    enddo
	  else if (rand.lt.cumprob2) then
	    !BRANCHING
	       Nparticlestemp = Nparticlestemp + 1
	       ParticlePos(Nparticlestemp) = ParticlePos(i)
	  else
	    !DIFFUSION
	    !print*,"diffusion"
	    call gaussian(gaussianvariable)
	    eta = gaussianvariable*dsqrt(2.0d0*diff/dt)	
	    ParticlePos(i) =  ParticlePos(i) + eta*dt
	  endif	  
	  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		  	

	enddo

	Nparticles = Nparticlestemp

	if (Nparticles.eq.0) then
	   goto 10
	endif	
	if (Nparticles.gt.90000) then
	   goto 10
	endif
	
	
	 
	!Calculating the maxima and minima
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!print*, Nparticles
	if (Nparticles.gt.0) then

	  call particlesort(Nparticles)
	  
	  if (Particlepos(Nparticles).gt.max_pos) then
		max_pos = Particlepos(Nparticles)
	  endif
	
	  if (Particlepos(1).lt.min_pos) then
		min_pos = Particlepos(1)
	  endif

	 	
	endif
	
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	



	enddo

	!**************************************************************************
	!end of timestep loop

	
	
	!Measuring Probability
	!****************************************
	 span_upto_t = max_pos - min_pos
	 	
	  diff_bin1 = floor(10.0d0*span_upto_t)
	

	  if (diff_bin1.lt.Nbins) then

	    !total_bin(timestep,k,j) = total_bin(timestep) + 1.0d0

            probspan1(diff_bin1) = 
     c	    probspan1(diff_bin1) + 1.0d0

	    Ntot1 = Ntot1 + 1.0d0
 	  endif

	!****************************************
	  if (Nparticles.gt.0) then

	  call particlesort(Nparticles)
	  
	   span_at_t = particlepos(N) - particlepos(1)
	   diff_bin2 = floor(10.0d0*span_at_t)

	  if (diff_bin2.lt.Nbins) then

	    !total_bin(timestep,k,j) = total_bin(timestep) + 1.0d0

            probspan2(diff_bin2) = 
     c	    probspan2(diff_bin2) + 1.0d0

	    Ntot2 = Ntot2 + 1.0d0
 	  endif
	 	
	endif
	  
	!****************************************



	!##########################################################################
	!return of totalsteps loop 
10	enddo


	!!!!!!!!!!!!!!WRITE INTERMEDIATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!##########################################################################
	open(unit = 44, file = "probspan_upto_t10_b1.0_dt0.1")	
	open(unit = 45, file = "probspan_at_10_b1.0_dt0.1")
	
	!Ntot = 1.0d0*prog*totalsteps
	
	do i =0,Nbins
	  write(44,*) 0.1d0*i,10.0d0*probspan1(i)/(Ntot1), Ntot1
	enddo
	
	!do k = 2,10
	do i =0,Nbins
	  write(45,*) 0.1d0*i,10.0d0*probspan2(i)/(Ntot2), Ntot2
	enddo
	
	!enddo

100 	format (f10.4)
101 	format (e14.6)	
	
	close(44)
	close(45)
	
	!##########################################################################



	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	!return of progress loop
	enddo
	
	

	end

	
	include 'mt.f'
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
!===============================================================================================
	subroutine particlesort(Nparticles)
	use particles
		
	integer i,j,Nparticles
	real*16 temp

	do i = 1,Nparticles
	  do j = i+1, Nparticles
	
		if (ParticlePos(i).gt.ParticlePos(j)) then
		
		temp = ParticlePos(i)
		ParticlePos(i) = ParticlePos(j)
		ParticlePos(j) = temp

		endif	
		
	  enddo
	enddo
	

	end
!========================================================================================================
	
