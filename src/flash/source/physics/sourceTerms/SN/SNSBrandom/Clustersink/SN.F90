!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!! - CB April-May 2012 based on MKRJ original SN feedback code
!! added ridiculous amounts of comments 
!! kept most of MKRJ's comments
!! cleaned code up a bit but still needs polish
!! explosion times calculation more clear
!! merged Heat_block into Heat
!! rewrote routine so heating is always applied
!! exp_heat switch to reduce timestep was removed and replaced with more adaptive method
!! related to above, it should not be needed to know there will be a SN in the next timestep to control it
!! instead of counting number of zones and averaging by number, volumes are calculated as the zones might differ in size
!! restructured code so that radiation cooling of SN energy will be applied one timestep after SN explosion
!! ASSUMES CUBIC ZONES INSIDE BLOCK
!! TODO better 2d/3d switch
!!! Every process calculates the position of every global random explosion, this requires
!!! all random number generators to run with the same results, which is ok but not very safe
!!! a better approach would be to let every local domain compute its local random SN field
!!! only the SB would need global coordination, although this could also be rectified by pre-
!!! calculating SB locally and just working through them.
!!!
!!! - JCIM March-April 2015
!!!	  Added the sink SN implementation.
!!!	  Sink particles are local to each processor but SN computation is global, so all processors need to
!!!   know that a sink SN is blowing up and should know the position of the explosion. The routine is now
!!!   generalized for multiple SNs happening in many sinks in the same or in different processors
!!!   at the same timestep.
!!!
!! TODO FINAL CLEANUP
!! TODO he_ should be sn_
!! TODO change snap to grid to be refinement dependent
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***
!#define DEBUG_SN

#define Juan_Debug

subroutine SN (blockCount,blockList,dt,time)

! "only", I think I could just use the whole thing, worth contemplating, also very high priority (not actually)
	use SN_data!, Only: he_smallpres, he_smalldens,	&
!		 he_tsn1, he_tsn2,  he_nsndt,		&
!		 he_hstar1, he_hstar2,			&
!		 he_seed, he_seedsize, he_meshMe,	&
!		 he_imin, he_imax, he_jmin,		&
!		 he_jmax, he_kmin, he_kmax,		&
!		 he_r_init, he_exp_energy, 		&
!		 he_Mejc,				&
!		 he_nSN, he_stratifySN,  he_exp_flag,	&
!		 he_SNminstep, he_newDt, he_useSN, he_useSNrandom, &
!		 he_r_exp_max, he_SNmapToGrid
  use mtmod

	use SB_data, Only: useSB, sb_nSN, sbPOS, sb_useSNsink, sn_outflowFrac, &
					   sb_life, sn_sinkBulkMotion
					
	use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
				   									Grid_releaseBlkPtr, Grid_getDeltas

	use Timers_interface, ONLY : Timers_start, Timers_stop
	use Eos_interface, ONLY : Eos_wrapped, Eos_getAbarZbar
	! for the sinks 
	use Particles_sinkData, ONLY : particles_global, localnp, useSinkparticles, particles_local
	! just for output
	use Driver_data, ONLY : dr_nStep, dr_globalMe
	use tree, ONLY : lrefine, bnd_box, lrefine_max
	use Grid_data, ONLY : gr_delta

	implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Flash_mpi.h"

! arguments
	integer,intent(IN)			:: blockCount
	integer,dimension(blockCount),intent(IN):: blockList
	real,intent(IN) 			:: dt,time

	! block data
	integer					:: blockID, thisBlock
	real, pointer, dimension(:,:,:,:)	:: solnData
	real, allocatable, dimension(:)		:: xCoord, yCoord, zCoord
	real, allocatable,dimension(:)		:: dx, dy, dz
	integer					:: xSizeCoord, ySizeCoord, zSizeCoord
	logical					:: getGuardCells = .true.
	integer, dimension(2,MDIM)		:: blkLimits, blkLimitsGC, pointLimit
	real, dimension(MDIM)			:: del

	! communication and mass calculation
	! TODO make nms a read in parameter
	!integer, parameter	:: nms = 50! 'nms' is the number of shells used to compute r_exp where M~60 M_sun
	real, DIMENSION(he_nms)	:: vol_sum_loc, vol_sum_tot
	real, DIMENSION(he_nms)	:: mass_sum_tot, mass_sum_loc
	integer, DIMENSION(he_nms)	:: nzones_sum_loc,nzones_sum_tot
	real, DIMENSION(he_nms)	:: rho_sum
	real			:: rho_avg, vx, vy, vz, rate
	real			:: blockBounds(LOW:HIGH,1:MDIM)

	! two error integers? DECADENCE!
	integer ::  error!, ierr

	! iterators
	! TODO see if they can be reused
	integer :: imass, i, j, k, l, m, im_exp, isum

	! SN module scratch variables, only local
	integer :: nsn1, nsn2, nsn2dt, nsn1dt
	real    :: t0heat, tstar1, tstar2, theat1, theat2
	real    :: tnext1, tnext2
	integer :: num_exp, nSBw, nSN, nSNsink

	integer, DIMENSION(he_nsndt)	:: SNtype, SinkTag
	real, DIMENSION(he_nsndt)	:: x0h, y0h, z0h, r0h
	real, DIMENSION(he_nsndt)	:: vx0h, vy0h, vz0h

	! variables for the SN explosion
	real	:: r_exp, v_exp, outshf, r_shell, mass
	real	:: xx, yy, zz, xb, yb, zb
	real	:: x0heat, y0heat, z0heat

	! heating and cooling variables
	real :: tranheat, argm, sdot
	real :: tmp, rho, ei, ek, qheat

	! SB specific stuff
	real :: SBdt
	integer :: sbIndex, nSBr, SNtmpType, SNtag, zoneindex

    ! count local number of affected zones, current offset for unique particle tags
    integer :: NlocalZones, tagcounter, ncount, ntracer

	! spacing
    real, dimension(NDIM)	:: gridspace
	
	! Necessary for the MC tracer particles.
	logical :: injectLeftovers 

	! SNsink, cloud mass and stellar mass returned to the ISM
	real :: sn_star_ejecta, sn_sink_ejecta 
	
	! SinkSN. temporary flag and arrays to communicate the explosions in each proc.
	real :: cluster_age
	logical :: exp_flag_temp
	integer :: sink_exp_here, sink_exp_total, nsnSinkdt, num_exp_old
	integer, allocatable,dimension(:) :: SinkSN_Tag_here, SinkSN_Tag_final
	integer, dimension(he_meshNumProcs) :: SinkSN_Tag_total, SinkSN_Tag_disp
	

! 	JCIM - April 2015.
!	! This is how I plan to get the sink particles tags that are hosting SN in this timestep.
!
!	logical :: mpi_test
!	integer :: num_here, num_total
!	
!	integer, allocatable,dimension(:) :: array_here, array_final
!	integer, dimension(he_meshNumProcs) :: array_total, array_disp
!
!	! Communicate the Tags of the sink particles that are hosting a SN in this timestep.
!	mpi_test = .TRUE.
!	num_here = 0
!	
!	if (he_meshMe .eq. 1) num_here = 3
!	if (he_meshMe .eq. 3) num_here = 1
!	if (he_meshMe .eq. 4) num_here = 2
!	if (he_meshMe .eq. 5) num_here = 1
!
!
!	! Just writing some mock tags in the different procs.
!	if(mpi_test) then
!		allocate(array_here(num_here))
!		if (num_here .gt. 0 ) then
!			array_here(1) = -1
!			do j=1,num_here
!				array_here(j) = 10*he_meshMe + j
!			enddo
!			PRINT*, ""
!			PRINT *,"Proc", he_meshMe, 'Has', num_here,'numbers here and the following local array =', array_here
!		endif
!	endif
!
!	! Get the total number of SN per processor and organize it in an array.
!	call MPI_Gather (num_here, 1, MPI_INTEGER, array_total, 1, MPI_INTEGER, &
!				& MASTER_PE, MPI_Comm_World, error)
!			
!	if( he_meshMe .eq. MASTER_PE ) then	
!		PRINT *,"Number of SN per procs is", array_total
!	endif
!	
!	! Get the total number of SN in all procs and communicate this to everyone.
!	call MPI_Reduce(num_here, num_total, 1, MPI_INTEGER, MPI_Sum, MASTER_PE, MPI_Comm_World, error)
!	call MPI_Bcast(num_total, 1, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
!	
!	! Compute the array of displacements for each proc.
!	! For example: 6 procs and 7 SN.
!	! 			   proc  | 0 | 1 | 2 | 3 | 4 | 5 |
!	! 			   ------|---|---|---|---|---|---|
!	! 			   nSN   | 0 | 3 | 0 | 1 | 2 | 1 |
!	! 			   disp  | 0 | 0 | 3 | 3 | 4 | 6 |
!	if( he_meshMe .eq. MASTER_PE ) then	
!		! Now sort the displacements
!		array_disp(1) = 0
!		do i = 2, he_meshNumProcs
!			array_disp(i) = array_disp(i-1) + array_total(i-1)
!		enddo
!		
!		PRINT *,"The array of displacements is", array_disp
!	endif
!	
!	!! Tell each processor the displacement each has to know to write the sink tags array.
!	!call MPI_Scatter(array_disp, 1, MPI_INTEGER, my_disp, 1, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
!			
!	! Allocate the array for the Sink tags of the total number of SN explosions in this timestep.
!	allocate(array_final(num_total))
!	array_final(:) = 0
!
!	! Organize the array with the SN tags in the master processor.
!	call MPI_Gatherv (array_here, num_here, MPI_INTEGER, array_final, array_total, array_disp, MPI_INTEGER, &
!				& MASTER_PE, MPI_Comm_World, error)
!				
!	! Communicate the array to all processors.
!	call MPI_Bcast(array_final, num_total, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
!
!	if (he_meshMe .eq. MASTER_PE) then
!		PRINT *, "The final array of SNs is", array_final
!	endif
	
	

if (.not. he_useSN) return

	!  reset timestep limiter
	he_newDt   = 1e99
	! start the timer ticking
	call Timers_start("SN")

! gravitational acceleration parameters, not used, should come from heat_data
! not used but in comments, should all go to heat_data
! real, save :: a_parm1, a_parm2, a_parm3
! real, save :: p_ambient, rho_ambient

	!  reset explosion flag
	he_exp_flag = .FALSE.

  !====================================================
  ! SN explosion times for given rates
  !====================================================

	!  reset explosion and stellar wind numbers, just to be sure
	num_exp = 0
	! /////////////////
	! // change for SB
	! ////////////////
	! for SB execution
	nSBr 	= 0
	sbIndex	= 0
	! for output
	nSBw 	= 0
	nSN 	= 0
	nSNsink	= 0 

	! for SinkSN.
	sink_exp_here  = 0
	sink_exp_total = 0
	num_exp_old    = 0

  	! call SB
  	if (useSB) then
	!	print*,'calling SB'
		SBdt = 1e99
		call SB(he_exp_flag, num_exp, dt, time, SBdt, nSBw)
		nSBr   = nSBw
		he_nSN = he_nSN + nSBr
  	endif

  	if(he_useSNrandom) then
		!  for general purpose, look back one time step, extend time interval to check backwards for overlap
		t0heat = time - dt
    	
		!  for field supernovae of Type I
		nsn1   = time/he_tsn1 ! total number of type I SN since start of simulation
		tstar1 = nsn1*he_tsn1 ! current explosion time for type I SN
    	
		!  for field supernovae of Type II
		nsn2   = time/he_tsn2 ! total number of type II SN since start of simulation
		tstar2 = nsn2*he_tsn2 ! current explosion time for type II SN
		
		!  nsn2dt = t0heat/he_tsn2   ! 'nsn2dt' is the # of isolated Type II expl'ns in a single 'dt'
		!  TODO check if integer conversion is working
		nsn2dt = dt/he_tsn2
		nsn1dt = dt/he_tsn1

  		!====================================================
  		! check if explosion occurs in this time step
  		!====================================================

		!  MKRJ - use the current seeds to generate random numbers in the next timestep
		!		  call random_seed(GET=he_seed(1:he_seedsize))
    	
		!  next explosion times, assumes SB rate is < SNI and SNII rate
		tnext1 = tstar1+he_tsn1
		tnext2 = tstar2+he_tsn2
    	
		!  calculate new step to land before explosion add some slop
		he_newDt = 0.5*(tnext1 - 0.9*he_SNminstep - time)
    	
		!  explosion is inside minimum timestep already
		!  if(he_newDt .le. 0d0) then
		if(he_newDt .le. he_SNminstep) then
			he_newDt = he_SNminstep
		endif
    	
		tnext1 = he_newDt
		he_newDt = 0.5*(tnext2 - 0.9*he_SNminstep - time)
    	
		!	if(he_newDt .le. 0d0) then
		if(he_newDt .le. he_SNminstep) then
			he_newDt = he_SNminstep
		endif 

		if(tnext1 .lt. he_newDt) then
			he_newDt = tnext1
		endif

!====================================================
! random SN field
!====================================================
	!  first, check if it's time for a field supernova of Type I
	if (t0heat < tstar1) then
		he_exp_flag = .TRUE.
		do i = num_exp+1, num_exp + 1  + nsn1dt 

			if(he_radialSN) then
!				call random_number(x0h(i))
!				call random_number(y0h(i))

        		x0h(i) = grnd()
				y0h(i) = grnd()

				r0h(i) = sqrt(x0h(i)**2d0 + y0h(i)**2d0)
				x0h(i) = x0h(i)/r0h(i)
				y0h(i) = y0h(i)/r0h(i)
				r0h(i) = he_erstar1*log(r0h(i)) ! use inverse transform, could also be 1-r0h, just to exclude 0
				x0h(i) = x0h(i)*r0h(i)
				y0h(i) = y0h(i)*r0h(i)
			else
				
				!  determine x location of random SN
				!  call random_number(x0h(i))
        		x0h(i) = grnd()
				x0h(i) = x0h(i)*(he_imax-he_imin)
				x0h(i) = x0h(i) + he_imin
				
				!  determine y location of random SN
				!  call random_number(y0h(i))
				y0h(i) = grnd()
				y0h(i) = y0h(i)*(he_jmax-he_jmin)
				y0h(i) = y0h(i) + he_jmin
			endif

!			call random_number(z0h(i))
			z0h(i) = grnd()
!			rndCalls = rndCalls + 3

			!  MKRJ - exponentially decling supernova rate in the vertical direction
			!         the scale height, hstar1, is an input parameter
			if(he_stratifySN) then
				z0h(i) = z0h(i)*2.0 - 1.0
				z0h(i) = ( abs(z0h(i))/z0h(i) ) * log(abs(z0h(i)))
				z0h(i) = z0h(i)*he_hstar1
			else
				z0h(i) = z0h(i)*(he_kmax-he_kmin)
				z0h(i) = z0h(i) + he_kmin
			endif

			PRINT *," JCIM - This is line 399. I'm setting SN position at 0, 0, 0."
			PRINT *,"        I need to remove this after testing."
			
			x0h(i) = 0
			y0h(i) = 0
			z0h(i) = 0

			SNType(i) = 1
			he_nSN    = he_nSN + 1
			nSN 	  = nSN + 1
		enddo

		num_exp = num_exp + nsn1dt + 1
	endif

	!  second, check if it's time for a field supernova of Type II
	if (t0heat < tstar2) then
		he_exp_flag = .TRUE.
		do i = num_exp+1, num_exp + nsn2dt +1  ! can have multiple Type II isolated expl'ns

			if(he_radialSN) then
!			  call random_number(x0h(i))
!			  call random_number(y0h(i))
				x0h(i) = grnd()
				y0h(i) = grnd()

			  	r0h(i) = sqrt(x0h(i)**2d0 + y0h(i)**2d0)
			  	x0h(i) = x0h(i)/r0h(i)
			  	y0h(i) = y0h(i)/r0h(i)
			  	r0h(i) = he_erstar2*log(r0h(i)) ! use inverse transform, could also be 1-r0h, just to exclude 0
			  	x0h(i) = x0h(i)*r0h(i)
			  	y0h(i) = y0h(i)*r0h(i)
			else
!			  call random_number(x0h(i))
				x0h(i) = grnd()
			  	x0h(i) = x0h(i)*(he_imax-he_imin)
			  	x0h(i) = x0h(i) + he_imin
!			  call random_number(y0h(i))
				y0h(i) = grnd()
			  	y0h(i) = y0h(i)*(he_jmax-he_jmin)
			  	y0h(i) = y0h(i) + he_jmin
			endif

!			call random_number(z0h(i))
!			rndCalls = rndCalls + 3
			z0h(i) = grnd()
			!  MKRJ - exponentially decling supernova rate in the vertical direction
			!        the scale height, hstar2, is an input parameter
			if(he_stratifySN) then
				z0h(i) = z0h(i)*2.0 - 1.0
				z0h(i) = ( abs(z0h(i))/z0h(i) ) * log(abs(z0h(i)))
				z0h(i) = z0h(i)*he_hstar2
			else
				z0h(i) = z0h(i)*(he_kmax-he_kmin)
				z0h(i) = z0h(i) + he_kmin
			endif

			SNType(i) = 2
			he_nSN    = he_nSN + 1
			nSN = nSN + 1
		enddo

		num_exp = num_exp + nsn2dt +1 
	endif
  endif !SNrandom switch


!/////////////////////////////////////////////////////////
!///////////  Correlated Sink SN explosions //////////////
!/////////////////////////////////////////////////////////

	! check up-to-date global sink particle array for SN explosion and add it to SN to execute
	
	num_exp_old = num_exp
	sink_exp_here = 0
 	
	if(SINK_PART_TYPE .gt. 0 .and. sb_useSNsink .and. useSinkparticles) then	

#ifdef DEBUG_SN
	print*,'Checking sinks for SN:'
#endif

    do i = 1, localnp
		
		! JCI -  Particles local should be used always instead of global.
		cluster_age = time - particles_global(CREATION_TIME_PART_PROP, i) 

		nsn1   = cluster_age / (particles_global(TSN_PART_PROP, i))
		nsn1dt = dt/particles_global(TSN_PART_PROP,i)

#ifdef Juan_Debug			
		!PRINT *, "JC **** cluster age, time to next SN: ", cluster_age , (nsn1+1.0)*particles_global(TSN_PART_PROP,i) - cluster_age 
		!PRINT *, "JC    **** nsn1, nsn1dt :", nsn1, nsn1dt
#endif
			
			
#ifdef DEBUG_SN
			print*,'time, creation time',time,particles_global(CREATION_TIME_PART_PROP, i)
			print*,'k-blammo?',nsn1
#endif


! JC - Check if sink should be turned off
      	if(nsn1 .gt. particles_global(NSN_PART_PROP, i) ) then 
        	print*, '    Turning off Sink. No more SNs left inside.'
        
#ifdef DEBUG_SN
	print*,'Turning off sink',i 
#endif
			cycle
		endif
				
        tstar1 = nsn1*particles_global(TSN_PART_PROP,i) + particles_global(TCRT_PART_PROP, i) 
        tnext1 = tstar1 + particles_global(TSN_PART_PROP, i)
		
! calculate new step to land before explosision add some slop
        !he_newDt = 0.5*(tnext1 - 0.9*he_SNminstep - time)
 		tnext1 = 0.5*(tnext1 - 0.9*he_SNminstep - time)

! JCIM - What's this line ?? where is he_newDt defined ? why is it commented above ?

! compare to SB creation derived timestep
        if(tnext1 .le. he_newDt) then
          he_newDt = tnext1
        endif

        if(tnext1 .le. he_SNminstep) then
          he_newDt = he_SNminstep
        endif

#ifdef Juan_Debug			
		!PRINT*, "t0heat and tstar1", t0heat, tstar1
#endif

! reuse tstar1
			if(t0heat < tstar1) then

#ifdef Juan_Debug			
				PRINT *,'----------------------------- lets add the sink SN explosion ----------------------------  proc', he_meshMe, 'says'
				PRINT *,'Im Sink particle tagged as =', int(particles_global(TAG_PART_PROP,i))
#endif

				he_exp_flag = .TRUE.
				
				nsnSinkdt = dt / particles_global(TSN_PART_PROP, i)
					
!				do j = num_exp+1, num_exp + 1  + nsn1dt
				! JCI - can there be more than one SN explosion per sink per timestep ?
				do j = num_exp+1, num_exp + 1  + nsnSinkdt
					
					! JCIM - num_exp+1 - num_exp + 1 + nsnSinkdt 
					!		 nsnSinkdt -> number of SN sinks in this time !!	
					
					x0h(j) = particles_global(POSX_PART_PROP, i)
					y0h(j) = particles_global(POSY_PART_PROP, i)
					z0h(j) = particles_global(POSZ_PART_PROP, i)
					
#ifdef Juan_Debug			
	print *, "Sink particle position:"
	print *, x0h(j), y0h(j), z0h(j)
#endif
					
					!x0temp(j)= particles_global(POSX_PART_PROP, i)
					!y0temp(j)= particles_global(POSY_PART_PROP, i)
					!z0temp(j)= particles_global(POSZ_PART_PROP, i)
                    
				   	if(sn_sinkBulkMotion) then
				   		vx0h(j) = particles_global(VELX_PART_PROP, i)
				   		vy0h(j) = particles_global(VELY_PART_PROP, i)
				   		vz0h(j) = particles_global(VELZ_PART_PROP, i)
				   	endif
				enddo
				
				SinkTag(i) =  int(particles_global(TAG_PART_PROP,i))
			   	!he_nSN    = he_nSN + 1
			   	!nSN = nSN + 1
			   	!num_exp = num_exp + nsn1dt + 1
			
				! JCI - New variables individual for each proc. to be communicated globally later.
				sink_exp_here = sink_exp_here + 1
			endif
    	enddo
	endif

	allocate(SinkSN_Tag_here(sink_exp_here))
	
	if (sink_exp_here .gt. 0) then
		do i = 1, sink_exp_here
			SinkSN_Tag_here(i) = SinkTag(i)
		enddo
	endif
	
	! JCI - Only the processor containing the Sink particle knows that it is time to explode a SN.
	!		 This needs to be communicated to the master proc to then be broadcasted to all procs.
	call MPI_Reduce (he_exp_flag, exp_flag_temp, 1, MPI_LOGICAL, MPI_Lor, &
				& MASTER_PE, MPI_Comm_World, error)
	
	he_exp_flag = exp_flag_temp
	
	! JCI - If one, or more procs are hosting a Sink SN explosion, all procs are informed about this.
	call MPI_Bcast(he_exp_flag, 1, MPI_LOGICAL, MASTER_PE, MPI_Comm_World, error)

	! Put a conditional here, like if he_exp_flag is true then look for stuff.

	! add all the individual explosions from sinks in the different procs and tell it to everyone.
	call MPI_Reduce(sink_exp_here, sink_exp_total, 1, MPI_INTEGER, MPI_SUM, MASTER_PE, MPI_Comm_World, error)
	call MPI_Bcast (sink_exp_total, 1, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)

	! JCI - add all the sink SN explosions to the SN counters (global, local and explosion loop.)
	he_nSN  = he_nSN  + sink_exp_total
   	nSN     = nSN     + sink_exp_total
   	num_exp = num_exp + sink_exp_total

	! JCI - All procs should know that the new SNs are of type SinkSN.
	if (sink_exp_total .gt. 0) then
		do j = num_exp_old+1, num_exp + 1
			SNType(j) = 4
		enddo
	endif

	!if( he_meshMe .eq. MASTER_PE ) PRINT *,"Number of SN in all procs is", sink_exp_total

	! Get the total number of SN per processor and organize it in an array.
	call MPI_Gather (sink_exp_here, 1, MPI_INTEGER, SinkSN_Tag_total, 1, MPI_INTEGER, &
					 & MASTER_PE, MPI_Comm_World, error)

	!if(he_meshMe .eq. MASTER_PE) PRINT *,"The total number of SN per proc is:", SinkSN_Tag_total

	! Compute the array of displacements for each proc.
	! For example: 6 procs and 7 SN.
	! 			   proc  | 0 | 1 | 2 | 3 | 4 | 5 |
	! 			   ------|---|---|---|---|---|---|
	! 			   nSN   | 0 | 3 | 0 | 1 | 2 | 1 |
	! 			   disp  | 0 | 0 | 3 | 3 | 4 | 6 |
	if( he_meshMe .eq. MASTER_PE ) then	
		
		! Now sort the displacements
		SinkSN_Tag_disp(1) = 0
		do i = 2, he_meshNumProcs
			SinkSN_Tag_disp(i) = SinkSN_Tag_disp(i-1) + SinkSN_Tag_total(i-1)
		enddo
		
		!PRINT *,"The array of displacements is", SinkSN_Tag_disp
	endif

	! Allocate the array for the Sink tags of the total number of SN explosions in this timestep.
	allocate(SinkSN_Tag_final(sink_exp_total))
	SinkSN_Tag_final(:) = 0
		
	! Organize the array with the SN tags in the master processor.
	call MPI_Gatherv (SinkSN_Tag_here, sink_exp_here, MPI_INTEGER, SinkSN_Tag_final, SinkSN_Tag_total, SinkSN_Tag_disp, MPI_INTEGER, &
				& MASTER_PE, MPI_Comm_World, error)
				
	! Communicate the array to all processors.
	call MPI_Bcast(SinkSN_Tag_final, sink_exp_total, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
	
	deallocate(SinkSN_Tag_here)

! To position the sink SN exposions,
! All processors need to know: 	he_exp_flag - check.
!							  	num_exp		- check.
!								SinkTag		- testing. I think it's working now. But I dont need it !!
!								SNType		- check.


!/////////////////////////////////////////////////////////
!///////////  Super Bubble (SB) explosions ///////////////
!/////////////////////////////////////////////////////////

!  MKRJ - Since multiple explosions are possible in a give timestep, 
!        loop over all explosions within the current timestep (handle only one explosion at a time)
!        in heat_block. All explosions are equal in terms of how & how much energy is added
! update timestep with SB
	if (useSB) then
		if(SBdt .lt. he_newDt) then
			he_newDt = SBdt
		endif
	endif
!====================================================
! execute SN, this involves MPI communication
! this should go into its own subroutine, this is just for the SN
!====================================================

! TODO switch for SN
	if (he_exp_flag) then 

	! Only one proc knows about the SN explosion.
	!PRINT *,"Hello from proc", he_meshMe

! MKRJ - outshf*r_exp gives the initial outer radius of dense shell 
!        outshf is set so that r_init*(outshf**(nms-1)) = r_exp_max (now 50 pc)
	!outshf = 1.0679e0
	outshf = (he_r_exp_max / he_r_init)**(1.0/(he_nms-1.0))

! move make space for added tracer particles
! TODO this could lead to trouble if there are a lot of simultaneous explosions
! move this to explosion loop if that's the case
!  loop over all SN
		do l = 1, num_exp

		  	NlocalZones = 0
			!  reset mass calculation
			mass_sum_loc = 0.e0  ! sum of mass within initial explosion sphere on this proc
			mass_sum_tot = 0.e0  ! sum of local masses from all procs
			vol_sum_loc  = 0.e0  ! sum of zone volumes in initial explosion sphere on this proc
			vol_sum_tot  = 0.e0  ! sum of local zone volumes in initial explosion sphere from all procs
			nzones_sum_loc = 0
! get position of heat source
! first work through SBs

! /////////////////
! // change for SB
! ////////////////
! prepare send buffer
		if(nSBr .gt. 0) then
			do i = sbIndex+1, sb_nSN
			    if(abs(sbPOS(1,i)) .gt. 1d-15) then
			      	x0heat = sbPOS(1,i) 
			      	y0heat = sbPOS(2,i)
			      	z0heat = sbPOS(3,i)
                  	
			      	SNtag  = sbPOS(4,i)
! clean sbPOS	  	  	
			      	sbPOS(1,i) = 0d0
			      	sbPOS(2,i) = 0d0
			      	sbPOS(3,i) = 0d0
			      	sbPOS(4,i) = 0d0
			      	sbIndex = i
			      	nSBr = nSBr - 1
					SNType(l) = 3
			    endif
			enddo
		else
			x0heat = x0h(l)
			y0heat = y0h(l)
			z0heat = z0h(l)
			
		endif

		! to ensure a nice bubble, this maps to the nearest highest refinement position good for uniform boxes, not so great if SN goes off in lower refined region
		if(he_SNmapToGrid) then

			! use block information from sink particle
		  	if((SNType(l) .eq. 4) .or. SNType(l) .eq. 3) then
				
				! guaranteed to be local
			   	if(SNType(l) .eq. 4) then
					
					!print *, "Im proc", he_meshMe, "And the sink tag I'm looking at is:", SinkTag(l)
					
					
					! JCIM - Ok, I'm doing something wrong here. do I look for the blockID with the sink tag ?
					!		 or with some other index ??
					! 		 Well I communicated the sink tag everywhere. Now I need to find the position of
					!        that sink so all procs can synchronize the explosion center (x,y,z)0heat.
										

					!blockID = particles_global(BLK_PART_PROP,SinkSN_Tag_final(l))

					!x0heat = particles_global(POSX_PART_PROP, SinkSN_Tag_final(l))
					!y0heat = particles_global(POSY_PART_PROP, SinkSN_Tag_final(l))
					!z0heat = particles_global(POSZ_PART_PROP, SinkSN_Tag_final(l))

					!blockID = particles_global(BLK_PART_PROP,SinkSN_Tag_final(l))
					
					
					! Loopong over l is not the correct answer because every processor has a different realization
					! of the Sink particle array, starting with the local particles.
					
					x0heat = x0h(l)
					y0heat = y0h(l)
					z0heat = z0h(l)
					
					if (he_meshMe .eq. MASTER_PE) then
						print *, "======= JC: I'm positioning SN center at (0,0,0) ........................... "
						print *, "            This is for testing purposes only. Get the proper Sink positions"
					endif
					x0heat = 0
					y0heat = 0
					z0heat = 0
					
					! How do I get the proper particle position ? If I want them mapped ?
					
					call gr_findBlock(blockList,blockCount,(/x0heat, y0heat, z0heat/),blockID)
										
					! assume it goes off in highest refinement region
					! This is always true if we are always refining on Sink Particles.
					gridspace = gr_delta(1:MDIM,lrefine_max)
					
					! JCI -  If I have the block information I could use this one.
                    !gridspace = gr_delta(1:MDIM,lrefine(blockID)) 

					zoneindex = floor((x0heat) / gridspace(1))
					x0heat = gridspace(1)*(zoneindex + 0.5d0)
                    
					zoneindex = floor((y0heat) / gridspace(2))
					y0heat = gridspace(2)*(zoneindex + 0.5d0)
                    
					zoneindex = floor((z0heat) / gridspace(3))
					z0heat = gridspace(3)*(zoneindex + 0.5d0)
					

#ifdef Juan_Debug
		if (he_meshMe .eq. MASTER_PE) then
				PRINT *, "JC: adding the SN explosion in the sink particle    **************** CATAPLUMMMMMMMM **************** "		
		endif	
#endif
					
					! can I retrieve the position information here ?
		    	else 
					call gr_findBlock(blockList,blockCount,(/x0heat, y0heat, z0heat/),blockID)

					gridspace = gr_delta(1:MDIM,lrefine(blockID))
					
					zoneindex = floor((x0heat) / gridspace(1))
					x0heat = gridspace(1)*(zoneindex + 0.5d0)
                	
					zoneindex = floor((y0heat) / gridspace(2))
					y0heat = gridspace(2)*(zoneindex + 0.5d0)
                	
			  		zoneindex = floor((z0heat) / gridspace(3))
			  		z0heat = gridspace(3)*(zoneindex + 0.5d0)
			
				endif

!          		  blockBounds = bnd_box(:,:,blockID)
!				  zoneindex = floor((x0heat-blockBounds(LOW ,IAXIS)) / gridspace(1))
!				  x0heat = blockBounds(LOW ,IAXIS) + gridspace(1)*(zoneindex + 0.5d0)

!				  zoneindex = floor((y0heat-blockBounds(LOW ,JAXIS)) / gridspace(2))
!				  y0heat = blockBounds(LOW ,JAXIS) + gridspace(2)*(zoneindex + 0.5d0)

!				  zoneindex = floor((z0heat-blockBounds(LOW ,KAXIS)) / gridspace(3))
!				  z0heat = blockBounds(LOW ,KAXIS) + gridspace(3)*(zoneindex + 0.5d0)

!				  zoneindex = floor((x0heat) / gridspace(1))
!				  x0heat = gridspace(1)*(zoneindex + 0.5d0)

!				  zoneindex = floor((y0heat) / gridspace(2))
!				  y0heat = gridspace(2)*(zoneindex + 0.5d0)

!				  zoneindex = floor((z0heat) / gridspace(3))
!				  z0heat = gridspace(3)*(zoneindex + 0.5d0)

			else
! dont have block information so assume highest refinement, should work if refinement is forced to be maximal at explosion
! point
!					call gr_findBlock(blockList,blockCount,(/x0heat, y0heat, z0heat/),blockID)
!				  gridspace = gr_delta(1:MDIM,lrefine(blockID))

! assume it goes off in highest refinement region
				gridspace = gr_delta(1:MDIM,lrefine_max)
                
				zoneindex = floor((x0heat) / gridspace(1))
				x0heat = gridspace(1)*(zoneindex + 0.5d0)
                
				zoneindex = floor((y0heat) / gridspace(2))
				y0heat = gridspace(2)*(zoneindex + 0.5d0)
                
				zoneindex = floor((z0heat) / gridspace(3))
				z0heat = gridspace(3)*(zoneindex + 0.5d0)
			endif
		endif
				
		print*, 'SN of type',SNtype(l),' at position.'
        print*, x0heat,y0heat,z0heat

!  this loop just gathers the mass information
!  TODO make this much faster by using the tree and only looking up neighbours
			do thisBlock = 1, blockCount
				blockID = blockList(thisBlock)

				call Grid_getBlkPtr(blockID,solnData)
				call Grid_getDeltas(blockID, del) !grid spacing dx dz dy
				call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

				xSizeCoord = blkLimitsGC(HIGH,IAXIS)
				ySizeCoord = blkLimitsGC(HIGH,JAXIS)
				zSizeCoord = blkLimitsGC(HIGH,KAXIS)

		! allocate space for dimensions
				allocate(xCoord(xSizeCoord))
				allocate(yCoord(ySizeCoord))
				allocate(zCoord(zSizeCoord))

				allocate(dx(xSizeCoord))
				allocate(dy(ySizeCoord))
				allocate(dz(zSizeCoord))

	  	 ! actually only one dx would be needed, -SizeCoord might not be right size, for cubic zones irrelephant
				dx(:) = del(IAXIS)
				dy(:) = del(JAXIS)
				dz(:) = del(KAXIS)

				call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
				call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
				call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

				rho_sum = 0.0
				
				! loop over all zones in block
				! for 2d k is just 1
				do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
					zz = zCoord(k)
					do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
						yy = yCoord(j)
						do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
							xx = xCoord(i)

							rho = solnData(DENS_VAR,i,j,k)

!  TODO make periodic switch, all periodic right now, does not matter if there never is an explosion close to a boundary
!  think z exponential with big box
!  MKRJ 11/22/04 - periodic boundary condition used in x & y, so take that into account

!  TODO add mass from a sink to innermost shell	
!  TODO: add mass contribution of the star exploding as SN in the innermost shell
!        e.g. Sink_mass / (NSN - nsn1) * sn_outflowfrac.

! Or should NSN per sink change ?? It would be easier, I think.
! 	this is mandatory to explode SN from Sink particles !!! 

							argm   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							argm   = argm + min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							argm   = argm + min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )

! is it anywhere inside the maximum radius
						    if (argm*argm <= he_maxradius*he_maxradius) then
								NlocalZones = NlocalZones + 1
              			    endif

!  compute the mass inside outshf*r_exp before the explosion
							do imass = 1, he_nms
							!  check which shell to add the zone
								if (sqrt(argm) <= (outshf**(imass-1))*he_r_init) then
									rho_sum(imass) = rho_sum(imass) + rho ! add density
									vol_sum_loc(imass) = vol_sum_loc(imass) + dx(i)*dy(j)*dz(k)
									nzones_sum_loc(imass) = nzones_sum_loc(imass) + 1
								endif
							enddo

						enddo !x
					enddo !y
				enddo !z

!  compute the mass inside explosion radius
!  TODO this does not check if enough mass was gathered
				do imass = 1, he_nms
					mass_sum_loc(imass) = mass_sum_loc(imass) + rho_sum(imass)*vol_sum_loc(imass)
				enddo
									
!  clean up memory 
        call Grid_releaseBlkPtr(blockID,solnData)

				deallocate(xCoord)
				deallocate(yCoord)
				deallocate(zCoord)

				deallocate(dx)
				deallocate(dy)
				deallocate(dz)
			enddo ! block loop

  !====================================================
  ! communicate found mass shell pieces to everyone
  ! also for the SN subroutine
  !====================================================

! TODO make faster by only talking to neighbouring cpus, or using FLASH4 communication routines 

!  MKRJ 7/3/04
!      Just for mass redistribution within initial explosion sphere, 
!      heat.F90 currently makes 4 MPI calls / timestep (2 MPI_Reduce's & 2 MPI_Bcast's) !!PER SN CALL!! -CB May 2012

!  add mass_s (nc_int) of all processors to get mass_sum_tot (nc_int_sum_tot)
!  TODO this should really go into a seperate communication routine
!  TODO use a better communcation routine than an reduce call, as most SN information just has to be send to neighbours

		call MPI_Reduce (mass_sum_loc, mass_sum_tot, he_nms, MPI_Double_Precision, MPI_Sum, &
				& MASTER_PE, MPI_Comm_World, error)

		call MPI_Reduce (vol_sum_loc, vol_sum_tot, he_nms, MPI_Double_Precision, MPI_Sum, &
				& MASTER_PE, MPI_Comm_World, error)

!  MKRJ - broadcast mass_sum from Master_PE to all processors

		call MPI_Bcast(mass_sum_tot, he_nms, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_Comm_World, error)

!  MKRJ - broadcast nc_int_sum as well to all processors!
!         call MPI_Bcast(nc_int_sum, nms, MPI_Integer, MASTER_PE, MPI_Comm_World, error)
!  CB - broadcast volumes to all procs kinda to much

		call MPI_Bcast(vol_sum_tot, he_nms, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_Comm_World, error)
!         nc_int_loc = nc_int_sum


! JC 2015
!!!!!!! Need To check this lines
!!!!!!! Actual injection of the SN.

		!!  CB - find out how many tracer particles will be injected
		if(he_useSNTracer) then
			call MPI_Reduce(nzones_sum_loc, nzones_sum_tot, he_nms, MPI_Integer, MPI_Sum, &
					& MASTER_PE, MPI_Comm_World, error)
			injectLeftovers = .false.
			!! broadcast to all processors
			call MPI_Bcast(nzones_sum_tot, he_nms, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
		else
			ntracer = 0
		endif
        

		!!  CB - broadcast volumes to all procs kinda to much
		call MPI_Bcast(vol_sum_tot, he_nms, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_Comm_World, error)

		! JCIM - Add sink mass to the shells. 
		!        could have a more sophisticated model based on age of sink
		if(SNType(l) .eq. 4) then
			
			! This needs to be checked closely.
			
			
			! JC 2015
			! TODO: Link the subgrid model to know the mass of the exploding star to be added here.
			
			! This is terrible !!!!! The mass of the star should be included when the mass per shell is being computed.
			! And the mininmum explosion radius should be given by the highest level of refinement.
			! SN explosions should happend always at the highest resolution. mmmmm not really because in zoom-in simulations that 
			! would definitely break everything.
			
			! Ok something to think about later...
			
			! I shouldn't have left this for later !! now that I have no time !! thanks past me !.

			! JCI: Ask if my sink particle is a single star or a cluster
			if (particles_global(NSN_PART_PROP,l) .le. 1.0) then
				! Sink has only one SN star inside.
				sn_sink_ejecta = sn_outflowFrac*particles_global(MASS_PART_PROP,l) 			
				mass_sum_tot(1) = mass_sum_tot(1) + sn_sink_ejecta
			else
				! JC: Simplest subgrid model ever. Each star returns to the ISM some
				!     an equal mass fraction to the ISM. => Mejecta = Cluster_Mass / Nstars
				! sn_star_ejecta  = 8.0 * 1.9884e33 
				sn_sink_ejecta = sn_outflowFrac*particles_global(MASS_PART_PROP,l) / &
								 & particles_global(NSN_PART_PROP,l)
				mass_sum_tot(1) = mass_sum_tot(1) + sn_sink_ejecta		
			endif
			
!			mass_sum_tot(:) = mass_sum_tot(:) + sn_outflowFrac*particles_global(MASS_PART_PROP,SinkTag(l))
! remove mass from sink, should be done on all cpus so no prob there.

#ifdef DEBUG_SN
			!print*,'sink mass before', particles_global(MASS_PART_PROP,SinkTag(l))
#endif

#ifdef Juan_Debug
	print *, 'sink mass before', particles_global(MASS_PART_PROP,l)
	!Print *, " *^*^*^*^ Need to check this computation line 1055 SN.F90"
#endif

!			particles_global(MASS_PART_PROP,SinkTag(l)) = particles_global(MASS_PART_PROP,SinkTag(l)) - sn_outflowFrac*particles_global(MASS_PART_PROP,SinkTag(l))
			
			! JCI - This needs to be done.
			if (he_meshMe .eq. MASTER_PE) then
				print *," ---- Subtracting mass from sink ----"
			endif 
			! Subtract the Mass of the SN explosion from the Sink. Only one proc needs to do this.
			!if (he_meshMe .eq. MASTER_PE) then
				!particles_global(MASS_PART_PROP,SinkSN_Tag_final(l)) = particles_global(MASS_PART_PROP,SinkSN_Tag_final(l)) - &
				!	& sn_outflowFrac * particles_global(MASS_PART_PROP,SinkSN_Tag_final(l))
			!endif
			! Subtract one SN from the NSN in the particle

#ifdef Juan_Debug
	print*,'sink mass after', particles_global(MASS_PART_PROP,l)
	!print*,'mass shells',     mass_sum_tot
#endif

#ifdef DEBUG_SN
			print*,'sink mass after', particles_global(MASS_PART_PROP,l)
#endif


#ifdef DEBUG_SN
			print*,'mass shells', mass_sum_tot
#endif

		endif
		

		! MKRJ - choose explosion radius
		im_exp = 0

sloop: do isum = 1, he_nms
			if (mass_sum_tot(isum) > he_mejc) then
				  im_exp = isum
				exit sloop
			endif
    enddo sloop

		if (im_exp == 0) then 
			im_exp = he_nms
      		if(he_meshMe == MASTER_PE) then 
				print*,'did not find enough mass:',mass_sum_tot(he_nms)/he_mejc
      		endif
		endif

		rho_avg = mass_sum_tot(im_exp)/(vol_sum_tot(im_exp))
		
		if(he_meshMe == MASTER_PE) then 
			print*,'##average density:',rho_avg
     	endif

		if(rho_avg .ne. rho_avg) then
			print*, 'NaN?', rho_avg, he_meshMe
	    	call Driver_abortFlash("[SN.F90] average density is NaN")
		endif
		
		r_exp    = (outshf**(im_exp-1))*he_r_init	! [cm]
    	mass     = rho_avg*vol_sum_tot(im_exp)		! [g]

		v_exp    = vol_sum_tot(im_exp)				! [cm^3]
		qheat    = dt*v_exp							! [cm^3/s]
		qheat    = he_exp_energy/qheat				! [erg*s/cm^3]

! JC this part is completely wrong, is old and not well maintained and should be updated
! calculate total number of new tracer
		if(he_useSNTracer) then ! if the tracer particles are following the velocity, inject them everywhere
			if(he_meshMe == MASTER_PE) then
				PRINT *, "WARNING: Dont trust this SN tracers, they are part of an old legacy code that "
				PRINT *, "		   needs to be updated."
			endif
  		if(he_veltracer) then 
			! calculate number of tracers per zone
			ncount = he_TracerPerSN/nzones_sum_tot(im_exp)
			!ncount = NlocalZones*he_tracerPerZone**3
      	else 
			ncount = he_TracerPerSN/nzones_sum_tot(im_exp)
			!ncount = NlocalZones*he_tracerPerZone
		endif
			
		! update tag offset, this does a couple of MPI calls >(
      	call pt_findTagOffset(ncount,tagcounter)
    endif

!print*,mass_sum_tot(im_exp),im_exp,he_meshMe
!print*,vol_sum_tot(im_exp),im_exp, he_meshMe
!print*,x0heat,y0heat,z0heat,r_exp
!print*,he_imax,he_imin,he_jmax,he_jmin, he_kmax,he_kmin

  !====================================================
  ! insert SN energy
  ! here should be an interface for the actual SN feedback implementation
  ! at the moment only thermal based
  !====================================================


#ifdef DEBUG_SN
			print*,'SN position',x0heat,y0heat,z0heat
			print*,'Average density in SN region',rho_avg
#endif

!  TODO record the relevant blockIDs beforehand and just loop over those
			do thisBlock = 1, blockCount
				blockID = blockList(thisBlock)
				call Grid_getBlkPtr(blockID,solnData)
				call Grid_getDeltas(blockID, del) !grid spacing dx dz dy
				call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
				xSizeCoord = blkLimitsGC(HIGH,IAXIS)
				ySizeCoord = blkLimitsGC(HIGH,JAXIS)
				zSizeCoord = blkLimitsGC(HIGH,KAXIS)

			! allocate space for dimensions
				allocate(xCoord(xSizeCoord))
				allocate(yCoord(ySizeCoord))
				allocate(zCoord(zSizeCoord))

				allocate(dx(xSizeCoord))
				allocate(dy(ySizeCoord))
				allocate(dz(zSizeCoord))

	  	 ! actually only one dx would be needed, -SizeCoord might not be right size, for cubic zones irrelephant
				dx(:) = del(IAXIS)
				dy(:) = del(JAXIS)
				dz(:) = del(KAXIS)
	
				call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
				call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
				call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

		 ! loop over all zones in block
		 ! for 2d k is just 1
				do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
					zz = zCoord(k)
					! get zone-local coordinate system
					zb = bnd_box(LOW,KAXIS,blockID) + dz(k)*(k-blkLimits(LOW,KAXIS))
					do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
						yy = yCoord(j)
						! get zone-local coordinate system
						yb = bnd_box(LOW,JAXIS,blockID) + dy(j)*(j-blkLimits(LOW,JAXIS))
						do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
							xx = xCoord(i)
							! get zone-local coordinate system
							xb = bnd_box(LOW,IAXIS,blockID) + dx(i)*(i-blkLimits(LOW,IAXIS))


!  MKRJ 11/22/04 - periodic boundary condition used in x & y, so take that into account
!  radial distance to SN = heat source

							argm   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							!argmx   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							!if (ndim >= 2) then
							argm   = argm + min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							!argmy   = min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							!endif
							!if (ndim == 3) then
							argm   = argm + min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )
							!argmz   = min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )
							!endif
							r_shell = sqrt(argm)
							
							if (r_shell <= r_exp) then
								!  get zone variables
								tmp = solnData(TEMP_VAR,i,j,k)
								rho = solnData(DENS_VAR,i,j,k)
								ei  = solnData(EINT_VAR,i,j,k)
								rate= solnData(PHHE_VAR,i,j,k)
								
!								vx  = solnData(VELX_VAR,i,j,k)
!								vy  = solnData(VELY_VAR,i,j,k)
!								vz  = solnData(VELZ_VAR,i,j,k)
!								tranheat = 0.0

#ifdef DEBUG_SN
!			print*,'shell radius',r_shell,r_exp,xx,yy,zz
#endif
								!print*,'x,y,z',argmx,argmy,argmz,he_meshMe
								!print*,'x,y,z',xx,yy,zz
								!v_exp    = 1.333*PI*(r_exp**3)

!								qheat    = dt*v_exp
!								qheat    = he_exp_energy/qheat

!  overwrites local values with SN
!								tranheat = qheat
								!tmp      = tmp*rho/rho_avg
								!ei       = ei*rho/rho_avg
								!rho      = rho_avg
								
#ifdef DEBUG_SN
			print*,'density input',rho,xx,yy,zz
#endif

								if(he_useSNTracer) then
							! check if zone is on most outer shell
!								if(r_shell - r_exp .lt. dx(i)) then
							! insert tracer particles 
									if(nSBr .gt. 0) then
								    	SNtmpType = 3
								  	else
								    	SNtmpType = SNType(l)
								    	!random SN have 0 tag
								    	SNtag     = he_nSN-(nSN-1)
								  	endif

								  	gridspace = gr_delta(1:MDIM,lrefine(blockID))

									if(he_veltracer) then 
								    call sn_createSNtracer(xb,yb,zb,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
								       									 tmp,rho,time,SNtag,SNtmpType,blockID, gridspace,tagcounter)
									else
								    call sn_createSNtracer(xx,yy,zz,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
								       									 tmp,rho,time,SNtag,SNtmpType,blockID, gridspace,tagcounter)
									endif
								endif


!  safey check
								!if (tmp > 1.e9) then
								!	tmp = 1.e9
								!	write(*,*) 'Limit temperature to a BILLION DEGREES'
								!endif

								!  radiation cooling and all the other terms were already applied, radiation cooling of SN energy will be done next timestep
								!sdot = tranheat/rho
								!  apply SN heating 
								! MKRJ - 8/26/04 at the timestep of a SN explosion, use explicit method
								!        (since implicit method seems to produce rather low peak temps - for a yet unknown reason)
								!  switch to explicit scheme for high thermal energy input
								!ei = ei + dt*sdot

								! set bulk motion
								if(sn_sinkBulkMotion) then
									solnData(VELX_VAR,i,j,k) = vx0h(l)
									solnData(VELY_VAR,i,j,k) = vy0h(l)
									solnData(VELZ_VAR,i,j,k) = vz0h(l)
								endif

								ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
								&           solnData(VELY_VAR,i,j,k)**2 + &
								&           solnData(VELZ_VAR,i,j,k)**2)


								! add SN heating to heating rate, sdot is in [erg/(s cm^3)]
								solnData(PHHE_VAR,i,j,k)  = rate + qheat*rho/rho_avg

								!solnData(TEMP_VAR,i,j,k)  = tmp
								!solnData(DENS_VAR,i,j,k)  = rho
								!solnData(EINT_VAR,i,j,k)  = ei
								!solnData(ENER_VAR,i,j,k)  = ei + ek
								
								pointLimit(:,IAXIS) = i
					 			pointLimit(:,JAXIS) = j
					 			pointLimit(:,KAXIS) = k

! check if radiation transport is included by looking for ionisation rate field
#ifdef PHIO_VAR
				! neutral and ionised hydrogen fraction
				! set all gas to ionized, this changes temperature effectively through eos call to 1/2
				solnData(IHA_SPEC,i,j,k) = 0d0
				solnData(IHP_SPEC,i,j,k) = 1d0
                call Eos_wrapped(MODE_DENS_EI,pointLimit,blockID)
#endif								
								
					
							endif
						enddo !z
					enddo !y
				enddo !x

				!  crank changed state variables through EOS, check it this is the right EOSMODE, again only changed blockIDs should be called
				!call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)


#ifdef DEBUG_SN
! check variables after eos call
		do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
					zz = zCoord(k)
					zb = bnd_box(LOW,KAXIS,blockID) + dz(k)*(k-blkLimits(LOW,KAXIS))
					do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
						yy = yCoord(j)
						yb = bnd_box(LOW,JAXIS,blockID) + dy(j)*(j-blkLimits(LOW,JAXIS))
						do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
							xx = xCoord(i)
							xb = bnd_box(LOW,IAXIS,blockID) + dx(i)*(i-blkLimits(LOW,IAXIS))


!  MKRJ 11/22/04 - periodic boundary condition used in x & y, so take that into account
!  radial distance to SN = heat source
							argm   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							!argmx   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							!if (ndim >= 2) then
							argm   = argm + min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							!argmy   = min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							!endif
							!if (ndim == 3) then
							argm   = argm + min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )
							!argmz   = min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )
							!endif
							r_shell = sqrt(argm)
							
							if (r_shell <= r_exp) then
								print*,solnData(DENS_VAR,i,j,k)
							endif

						enddo !z
					enddo !y
				enddo !x
#endif

!  clean up memory 
				call Grid_releaseBlkPtr(blockID,solnData)

				deallocate(xCoord)
				deallocate(yCoord)
				deallocate(zCoord)

				deallocate(dx)
				deallocate(dy)
				deallocate(dz)

		!		if(sb_useSNTracer) then
		!		  deallocate(dx)
		!		endif

			enddo ! block loop

!  write basic SN data to 
!  skip SB output
			if(nSBw .gt. 0) then
				call WriteSNFeedback(he_nSN-nSBw+1, 3, dr_nstep, time, x0heat, y0heat, z0heat, r_exp, mass)
				nSBw = nSBw - 1
				cycle
			else
				nSN = nSN - 1
! could have multiple calls fix
				call WriteSNFeedback(he_nSN-nSN, SNType(l), dr_nstep, time, x0heat, y0heat, z0heat, r_exp, mass)
			endif
		enddo ! SN loop
	endif !he_exp_flag

  !====================================================
  ! finalize
  !====================================================

!	TODO redo the whole tracer particle approach

! MKRJ - use the current seeds to generate random numbers in the next timestep
! TODO RESTARTING DOES NOT WORK AS RANDOM SEED IS RESET FIX THAT !!
!	call random_seed(GET=he_seed(1:he_seedsize))

!	TODO only call when SN are used, SN switch here
	call MPI_Barrier (MPI_Comm_World, error)
	call Timers_stop("SN")
	return

contains

  !====================================================
  ! output routine for supernovae, writes to SNfeedback.dat
  !====================================================

  subroutine WriteSNFeedback(nSN,type,ndt,time,x,y,z,radius,mass)

    use SN_data, ONLY :  he_meshMe
!  use Driver_data, ONLY : dr_globalMe

    implicit none

#include "constants.h"
#include "Flash.h"

    real, intent(IN)   :: time, x,y,z,radius,mass
    integer, intent(IN):: nSN,type,ndt ! 1 is I 2 is II and 3 is SN in SB

    integer, parameter :: funit_evol = 15
    character(len=80)  :: outfile = "SNfeedback.dat"
    integer            :: i
    logical, save      :: firstCall = .TRUE.

!    if (firstCall) then
!
!      if (he_meshMe .eq. MASTER_PE) then
!        open(funit_evol, file=trim(outfile), position='APPEND')
!        write(funit_evol,'(9(1X,A16))') '[00]n_SN', '[01]type', '[02]n_timestep', '[03]time', '[04]posx', &
!                                      '[05]posy', '[06]posz', '[07]radius'
!        close(funit_evol)
!      endif
!
!      firstCall = .false.
!
!    endif

!   as all processors calculate the SN in lockstep, only the Master process has to output data.
    if (he_meshMe .NE. MASTER_PE) return

    open(funit_evol, file=trim(outfile), position='APPEND')
!do i = 1, localnpf

    write(funit_evol,'(3(1X,I16),6(1X,ES16.9))') &
			 nSN,		&
	     type,  &
	     ndt,   &
	     time,  &
	     x,     &
	     y,     &
	     z,     &
	     radius,&
       mass

    close(funit_evol)

    return
  end subroutine WriteSNFeedback

end subroutine SN
