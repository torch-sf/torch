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
!!! TODO every process calculates the position of every global random explosion, this requires
!!! all random number generators to run with the same results, which is ok but not very safe
!!! a better approach would be to let every local domain compute its local random SN field
!!! only the SB would need global coordination, although this could also be rectified by pre-
!!! calculating SB locally and just working through them.
!! TODO FINAL CLEANUP
!! TODO he_ should be sn_
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

!!===============================
! leftovers are not injected, search lefty
!!===============================

!#DFINE VERBOSE
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
!		 he_r_exp_max, he_SNmapToGrid, he_erstar1, he_erstar2, &
!		 he_radialSN, he_useSNTracer
  use mtmodSN

! for output at right time
  use IO_data, ONLY : io_justCheckpointed, io_checkpointFileNumber

	use SB_data, Only: useSB, sb_nSN, sbPOS
	use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
	    Grid_releaseBlkPtr, Grid_getDeltas

	use Timers_interface, ONLY : Timers_start, Timers_stop
	use Eos_interface, ONLY : Eos_wrapped, Eos_getAbarZbar

	! just for output
	use Driver_data, ONLY : dr_nStep, dr_simTime
	use tree, ONLY : lrefine, bnd_box, lrefine_max
	use Grid_data,		ONLY : gr_delta

! for tracer offset calculation
  use Particles_data, ONLY : pt_startTagNumber

! just for output to file

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
	integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
	real, dimension(MDIM)			 :: del

	! communication and mass calculation
!	integer, parameter	  :: nms = 50! 'nms' is the number of shells used to compute r_exp where M~60 M_sun
	real, DIMENSION(he_nms)	:: vol_sum_loc, vol_sum_tot
	real, DIMENSION(he_nms)	:: mass_sum_tot, mass_sum_loc
	integer, DIMENSION(he_nms)	:: nzones_sum_loc,nzones_sum_tot

	real			:: rho_avg
	real			:: blockBounds(LOW:HIGH,1:MDIM)
	! two error integers? DECADENCE!
	integer ::  error!, ierr

	! iterators
	integer :: imass, i, j, k, l, m, im_exp, isum, im_MC

	! SN module scratch variables, only local
	integer :: nsn1, nsn2, nsn2dt, nsn1dt
	real    :: t0heat, tstar1, tstar2, theat1, theat2
	real    :: tnext1, tnext2
	integer :: num_exp, nSBw, nSN

	integer, DIMENSION(he_nsndt)	:: SNtype
	real, DIMENSION(he_nsndt)	:: x0h, y0h, z0h, r0h

	! variables for the SN explosion
	real	:: r_exp, v_exp, outshf, r_shell, mass
	real	:: xx, yy, zz, xb, yb, zb
	real	:: x0heat, y0heat, z0heat

	! heating and cooling variables
	real :: tranheat, argm, sdot
	real :: tmp, rho, ei, ek, qheat, convf

	! SB specific stuff
	real :: SBdt
	integer :: sbIndex, nSBr, SNtmpType, SNtag, zoneindex
  
  ! count local number of affected zones, current offset for unique particle tags
  integer :: tagcounter, ncount, ntracer, leftovers
	! outermost radius for MC tracer particle injection
	real	  :: MCradius, checkRadius
	logical :: injectLeftovers 

	! spacing
  real, dimension(NDIM)	:: gridspace

  character (len=4) ::convert
	

  if (.not. he_useSN) return

!	do i = 1, 1000*3+
!		SBdt = call grndSN()
!		print*,i,SBdt
!	enddo
!  call mtsavef('RNG_state_SN_0011','u')


! write RNG state if there was a checkpoint
  if(io_justCheckpointed) then
! save random number state, independent of flash4 random number stream
! use only main thread for this, as all other should have the same state
		if(he_meshMe == MASTER_PE) then 
      write(convert, '(i4.4)') io_checkpointFileNumber - 1
      call mtsavef('RNG_state_SN_' // trim(convert),'u')
		endif
	endif

!  reset timestep limiter
	he_newDt   = 1e99

	! start the timer ticking
	call Timers_start("SN")

! gravitational acceleration parameters, not used, should come from heat_data
! not used but in comments, should all go to heat_data

! reset explosion flag
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
	nSBr 		= 0
	sbIndex	= 0
! for output
	nSBw 		= 0
	nSN 		= 0 

! call SB
  if (useSB) then
		SBdt = 1e99
! checks for new SB creation and updates number of SN from SB
! by changing num_exp and nSBw, also calculates timestep from SB SNs
		call SB(he_exp_flag, num_exp, dt, time, SBdt, nSBw)
! nSBr and nSBw defined for output loop and iteration over execution
! clumsy but works
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
!	call random_seed(GET=he_seed(1:he_seedsize))

!  next explosion times, assumes SB rate is < SNI and SNII rate
		tnext1 = tstar1+he_tsn1
		tnext2 = tstar2+he_tsn2

!  calculate new step to land before explosion add some slop
		he_newDt = 0.5*(tnext1 - 0.9*he_SNminstep - time)

!  explosion is inside minimum timestep already
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
					x0h(i) = grndSN()
					y0h(i) = grndSN()

			  	r0h(i) = sqrt(x0h(i)**2d0 + y0h(i)**2d0)
			  	x0h(i) = x0h(i)/r0h(i)
			  	y0h(i) = y0h(i)/r0h(i)
			  	r0h(i) = he_erstar1*log(r0h(i)) ! use inverse transform, could also be 1-r0h, just to exclude 0
			  	x0h(i) = x0h(i)*r0h(i)
			  	y0h(i) = y0h(i)*r0h(i)
				else
!  determine x location of random SN
        	x0h(i) = grndSN()
			  	x0h(i) = x0h(i)*(he_imax-he_imin)
			  	x0h(i) = x0h(i) + he_imin
!  determine y location of random SN
					y0h(i) = grndSN()
			  	y0h(i) = y0h(i)*(he_jmax-he_jmin)
			  	y0h(i) = y0h(i) + he_jmin
				endif

				z0h(i) = grndSN()

!  MKRJ - exponentially decling supernova rate in the vertical direction
!        the scale height, hstar1, is an input parameter
				if(he_stratifySN) then
					z0h(i) = z0h(i)*2.0 - 1.0 !-1 to 1 
					z0h(i) = ( abs(z0h(i))/z0h(i) ) * log(abs(z0h(i))) !use inverse transform
					z0h(i) = z0h(i)*he_hstar1 !log(z)*h
				else
					z0h(i) = z0h(i)*(he_kmax-he_kmin)
					z0h(i) = z0h(i) + he_kmin
				endif

				SNType(i) = 1
				he_nSN    = he_nSN + 1
				nSN	  = nSN + 1
			enddo

			num_exp = num_exp + nsn1dt + 1
		endif

	!  second, check if it's time for a field supernova of Type II
		if (t0heat < tstar2) then
			he_exp_flag = .TRUE.
			do i = num_exp+1, num_exp + nsn2dt +1  ! can have multiple Type II isolated expl'ns

				if(he_radialSN) then
					x0h(i) = grndSN()
					y0h(i) = grndSN()
				  r0h(i) = sqrt(x0h(i)**2d0 + y0h(i)**2d0)
				  x0h(i) = x0h(i)/r0h(i)
				  y0h(i) = y0h(i)/r0h(i)
				  r0h(i) = he_erstar2*log(r0h(i)) ! use inverse transform, could also be 1-r0h, just to exclude 0
				  x0h(i) = x0h(i)*r0h(i)
				  y0h(i) = y0h(i)*r0h(i)
				else
					x0h(i) = grndSN()
				  x0h(i) = x0h(i)*(he_imax-he_imin)
				  x0h(i) = x0h(i) + he_imin
					y0h(i) = grndSN()
	 			  y0h(i) = y0h(i)*(he_jmax-he_jmin)
				  y0h(i) = y0h(i) + he_jmin
				endif

				z0h(i) = grndSN()
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

! MKRJ - outshf*r_exp gives the initial outer radius of dense shell 
!        outshf is set so that r_init*(outshf**(nms-1)) = r_exp_max (now 50 pc)
! CB - changed to actual calculation
		outshf = (he_r_exp_max / he_r_init)**(1.0/(he_nms-1.0))

!  loop over all SN
		do l = 1, num_exp

!  reset mass calculation
			mass_sum_loc = 0.e0  ! sum of mass within initial explosion sphere on this proc
			mass_sum_tot = 0.e0  ! sum of local masses from all procs
			vol_sum_loc  = 0.e0  ! sum of zone volumes in initial explosion sphere on this proc
			vol_sum_tot  = 0.e0  ! sum of local zone volumes in initial explosion sphere from all procs
			nzones_sum_loc = 0
! first work through SBs

! /////////////////
! // change for SB
! /////////////////
! TODO add type and tag fields for tracers
			if(nSBr .gt. 0) then
				do i = sbIndex+1, sb_nSN
! stop if it's not empty, so ugly, would be nicer with logical array, but would require more MPI
					if(abs(sbPOS(1,i)) .gt. 1d-14) then
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

				if(SNType(l) .eq. 3) then
! should be mapped to grid in SB.F90

! guaranteed to be local
!					call gr_findBlock(blockList,blockCount,(/x0heat, y0heat, z0heat/),blockID)
!					gridspace = gr_delta(1:MDIM,lrefine(blockID))
!					gridspace = gr_delta(1:MDIM,lrefine_max)

! blockbounds not needed
!					zoneindex = floor((x0heat) / gridspace(1))
!					x0heat = gridspace(1)*(zoneindex + 0.5d0)

!				  zoneindex = floor((y0heat) / gridspace(2))
!				  y0heat = gridspace(2)*(zoneindex + 0.5d0)

!				  zoneindex = floor((z0heat) / gridspace(3))
!				  z0heat = gridspace(3)*(zoneindex + 0.5d0)
				else
! get grid spacing, assumes cubic box and 8^3 zones per block 
!			  gridspace = (he_imax-he_imin)/(2**(lrefine_max+2))

! this does not work for MPI runs
!				call gr_findBlock(blockList,blockCount,(/x0heat,y0heat,z0heat/),blockID)

! assume it goes off in highest refinement region, as it is not local (would need mpi to find out)
			  	gridspace = gr_delta(1:MDIM,lrefine_max)

! find nearest zone position
!			  x0heat = DNINT(x0heat/gridspace)*gridspace
!			  y0heat = DNINT(y0heat/gridspace)*gridspace
!			  z0heat = DNINT(z0heat/gridspace)*gridspace

					zoneindex = floor((x0heat) / gridspace(1))
					x0heat = gridspace(1)*(zoneindex + 0.5d0)
					zoneindex = floor((y0heat) / gridspace(2))
					y0heat = gridspace(2)*(zoneindex + 0.5d0)
					zoneindex = floor((z0heat) / gridspace(3))
					z0heat = gridspace(3)*(zoneindex + 0.5d0)
       endif
			endif

			if(he_meshMe == MASTER_PE) then 
				print*, 'SN of type',SNtype(l),' at position'
				print*, x0heat,y0heat,z0heat
			endif

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

				! loop over all zones in block
				! for 2d k is just 1
				do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
					zz = zCoord(k)
					do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
						yy = yCoord(j)
						do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
							xx = xCoord(i)

							rho = solnData(DENS_VAR,i,j,k)

!  TODO make periodic switch, all periodic right now
!  MKRJ 11/22/04 - periodic boundary condition used in x & y, so take that into account
!  radial distance to SN = heat source

							argm   = min( (xx-x0heat)**2, (xx-(x0heat+he_imax-he_imin))**2, (xx-(x0heat-he_imax+he_imin))**2 )
							argm   = argm + min( (yy-y0heat)**2, (yy-(y0heat+he_jmax-he_jmin))**2, (yy-(y0heat-he_jmax+he_jmin))**2 )
							argm   = argm + min( (zz-z0heat)**2, (zz-(z0heat+he_kmax-he_kmin))**2, (zz-(z0heat-he_kmax+he_kmin))**2 )

! is it anywhere inside the maximum radius
!  compute the mass inside outshf*r_exp before the explosion
							do imass = 1, he_nms
							!  check which shells to add the zone density to
							!  it is added to all shells containing the zone
								if (sqrt(argm) <= (outshf**(imass-1))*he_r_init) then
									mass_sum_loc(imass) = mass_sum_loc(imass) + rho*dx(i)*dy(j)*dz(k) ! add density
									vol_sum_loc(imass) = vol_sum_loc(imass) + dx(i)*dy(j)*dz(k)
									nzones_sum_loc(imass) = nzones_sum_loc(imass) + 1
								endif
							enddo
						enddo !x
					enddo !y
				enddo !z

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

!  TODO use a better communcation routine than an reduce call, as most SN information just has to be send to neighbours
			call MPI_Reduce (mass_sum_loc, mass_sum_tot, he_nms, MPI_Double_Precision, MPI_Sum, &
					& MASTER_PE, MPI_Comm_World, error)

			call MPI_Reduce (vol_sum_loc, vol_sum_tot, he_nms, MPI_Double_Precision, MPI_Sum, &
					& MASTER_PE, MPI_Comm_World, error)

!  MKRJ - broadcast mass_sum from Master_PE to all processors
			call MPI_Bcast(mass_sum_tot, he_nms, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_Comm_World, error)
	
!  CB - find out how many tracer particles will be injected
			if(he_useSNTracer) then

				call MPI_Reduce(nzones_sum_loc, nzones_sum_tot, he_nms, MPI_Integer, MPI_Sum, &
						& MASTER_PE, MPI_Comm_World, error)
				injectLeftovers = .false.
! broadcast to all processors
				call MPI_Bcast(nzones_sum_tot, he_nms, MPI_INTEGER, MASTER_PE, MPI_Comm_World, error)
			else
				ntracer = 0
			endif

!  CB - broadcast volumes to all procs kinda to much
			call MPI_Bcast(vol_sum_tot, he_nms, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_Comm_World, error)
	
! MKRJ - choose explosion radius
			im_exp = 0

sloop: do isum = 1, he_nms
				if (mass_sum_tot(isum) > he_mejc) then
					im_exp = isum
			  	exit sloop
				endif
       enddo sloop

! maybe reroll to a different location?
			if (im_exp == 0) then
				im_exp = he_nms
      	if(he_meshMe == MASTER_PE) then 
				  print*,'did not find enough mass:',mass_sum_tot(he_nms)/he_mejc
      	endif
			endif

			rho_avg = mass_sum_tot(im_exp)/(vol_sum_tot(im_exp)) 

    	if(he_meshMe == MASTER_PE) then 
			  print*,'average density:',rho_avg
     	endif

			r_exp   = (outshf**(im_exp-1))*he_r_init
    	mass    = rho_avg*vol_sum_tot(im_exp)
			qheat   = dt*vol_sum_tot(im_exp)
			qheat   = he_exp_energy/qheat

! calculate total number of new tracer
			if(he_useSNTracer) then

! if the tracer particles are following the velocity, inject them everywhere

				if(he_veltracer) then
! calculate number of tracers per zone
					ncount = he_TracerPerSN/nzones_sum_tot(im_exp)
				else
! inject tracer particles in -1 shells
					! find shell 1 pc away from edge
					! loop over all shells, backwards, stop at first inside shell that fullfills distance to SN region edge

						if(r_exp .gt. sn_MCTracerMaxRad) then
							checkRadius = sn_MCTracerMaxRad
							if(he_meshMe == MASTER_PE) then
								print*,'injection and expl. radius, shell index'
								print*, MCradius, r_exp, im_MC
							endif
						else
							checkRadius = r_exp
						endif

					  im_MC = he_nms
					  do im_MC = he_nms, 1, -1
					    MCradius = (outshf**(im_MC-1))*he_r_init
					    if(MCradius .le. (checkRadius - sn_MCTracerShellDis)  ) then
								if(he_meshMe == MASTER_PE) then
									print*,'injection and expl. radius, shell index'
									print*, MCradius, r_exp, im_MC
								endif
					     exit 
					    endif
					  enddo
					
					if(nzones_sum_tot(im_MC) .ge. 1) then
						! find shell 1 pc away from edge
						ncount = he_TracerPerSN/nzones_sum_tot(im_MC)
						leftovers = he_TracerPerSN- ncount*nzones_sum_tot(im_MC)

						if(he_meshMe == MASTER_PE) then
							print*,'number of particles per zone, leftovers, zones in shell, total zones'
							print*,ncount,leftovers,nzones_sum_tot(im_MC),nzones_sum_tot(im_exp)
						endif
					else
						ncount = he_TracerPerSN/nzones_sum_tot(im_exp)
						MCradius = (outshf**(im_exp-1))*he_r_init
            leftovers = he_TracerPerSN- ncount*nzones_sum_tot(im_exp)
						im_MC = im_exp
						if(he_meshMe == MASTER_PE) then
							print*,'no zones far enough from sn region, using the whole bubble'
							print*,'number of particles per zone, leftovers'
							print*,ncount,leftovers,nzones_sum_tot(im_exp)
						endif
					endif

				endif

	    	if(he_meshMe == MASTER_PE) then
					if(he_veltracer) then
! lefty 
 						ntracer = nzones_sum_tot(im_exp)*ncount!+leftovers
					else
						if(nzones_sum_tot(im_MC) .ge. 1) then
! lefty 
							ntracer = nzones_sum_tot(im_MC)*ncount!+leftovers
						else
						!	ntracer = nzones_sum_tot(im_exp)*ncount+leftovers
							print*,'no zones found in SN.F90 tracer particle injection'
						endif
					endif
  	  	endif

! update tag offset, this does a couple of MPI calls >(

				if(he_veltracer) then

				  he_oldStartTag = pt_startTagNumber
! calculate number of tracers per zone
					call pt_findTagOffset(ncount*nzones_sum_loc(im_exp),tagcounter)
					injectLeftovers = .false.

					if((he_TracerPerSN - (tagcounter-he_oldStartTag) - ncount*nzones_sum_loc(im_exp)) .eq. leftovers) then
! lefty
!						injectLeftovers = .true.
						injectLeftovers = .false.
					endif
					he_oldStartTag = pt_startTagNumber + leftovers
					pt_startTagNumber  = pt_startTagNumber + leftovers
				else

					if(nzones_sum_loc(im_MC) .ge. 1) then
   				  he_oldStartTag = pt_startTagNumber
						call pt_findTagOffset(ncount*nzones_sum_loc(im_MC),tagcounter)
						injectLeftovers = .false.

            if((he_TracerPerSN - (tagcounter-he_oldStartTag) - ncount*nzones_sum_loc(im_MC)) .eq. leftovers) then
! lefty
!						injectLeftovers = .true.
							injectLeftovers = .false.
						endif
					else
						call pt_findTagOffset(0,tagcounter)
					endif
! lefty
! adjust tagnumber to account for leftovers (leftovers known on each process)
!					he_oldStartTag = pt_startTagNumber + leftovers
!					pt_startTagNumber  = pt_startTagNumber + leftovers


				endif
! add them leftovers
	endif

  !====================================================
  ! insert SN energy
  ! here should be an interface for the actual SN feedback implementation
  ! at the moment only thermal based
  !====================================================

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
!  get zone variables
							  tmp = solnData(TEMP_VAR,i,j,k)
							  rho = solnData(DENS_VAR,i,j,k)
							  ei  = solnData(EINT_VAR,i,j,k)

! conversion only working for values consistent with eos
								convf = tmp/ei
                tranheat = 0.0

!  overwrites local values with SN
								tranheat = qheat
								tmp      = tmp*rho/rho_avg
								ei       = ei*rho/rho_avg
								rho      = rho_avg

!								vx  = solnData(VELX_VAR,i,j,k)
!								vy  = solnData(VELY_VAR,i,j,k)
!								vz  = solnData(VELZ_VAR,i,j,k)
#ifdef DEBUG_SN
			print*,'density input',rho,xx,yy,zz
#endif

								if(he_useSNTracer) then
							! check if zone is on most outer shell
							! insert tracer particles 


								  if(nSBr .gt. 0) then
								    SNtmpType = 3

								  else
								    SNtmpType = SNType(l)
								    !random SN start with 0 tag
								    SNtag     = he_nSN-(nSN-1)
								  endif

								  gridspace = gr_delta(1:MDIM,lrefine(blockID))

! here the passed temperature value is wrong, as there was no EOS call
! but these particles are just passive so no biggie.
 									if(he_veltracer) then 
										if(injectLeftovers) then
											call sn_createSNtracer(xb,yb,zb,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
														tmp,rho,time,SNtag,SNtmpType,blockID, gridspace, tagcounter, ncount+leftovers)
											injectleftovers = .false.
										else
											call sn_createSNtracer(xb,yb,zb,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
															tmp,rho,time,SNtag,SNtmpType,blockID, gridspace, tagcounter, ncount)
										endif
									else
! extra check if in outer most shell or not
										if (r_shell <= MCradius) then
											if(injectLeftovers) then
               
												call sn_createSNtracer(xx,yy,zz,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
 												tmp,rho,time,SNtag,SNtmpType,blockID, gridspace, tagcounter, ncount+leftovers)
                        print*,'inject',ncount+leftovers,leftovers
												injectLeftovers = .false.
											else
												call sn_createSNtracer(xx,yy,zz,solnData(VELX_VAR,i,j,k),solnData(VELY_VAR,i,j,k),solnData(VELZ_VAR,i,j,k),&
 												tmp,rho,time,SNtag,SNtmpType,blockID, gridspace, tagcounter, ncount)
											endif
										endif
									endif ! tracer
								endif ! shell

							!  radiation cooling and all the other terms were already applied, radiation cooling of SN energy will be done next timestep
								sdot = tranheat/rho
							!  apply SN heating 
								ei = ei + dt*sdot

!  safety check
		            if(ei*convf .ge. sn_max_temp ) then
! limit temperature to 1e9 K
! print's should go somewhere else...

								open(sn_funit_evol, file=trim(sn_outfile), position='APPEND')

								write(sn_funit_evol,'(2(1X,I16),6(1X,ES16.9),5(1X,I16))') &
									dr_nStep, SNtag, dr_simTime,(ei*convf), sn_max_temp,(ei*convf)/sn_max_temp, &
									dt*sdot*dy(j)*dx(i)*dz(k)*rho/he_exp_Energy,mass_sum_tot(he_nms)/he_mejc,i,j,k,blockID,he_meshMe
								close(sn_funit_evol)

#ifdef VERBOSE
			            print*, '\\temp. too high in SN, changing internal energy input'
			            print*, '\\factor of ',(ei*convf)/sn_max_temp,'too high'
			            print*, '\\lost energy',dt*sdot*dy(j)*dx(i)*dz(k)*rho/he_exp_Energy ,'%'
			            print*, '\\core ID, block ID, zone ID xyz',he_meshMe,blockID,i,j,k
#endif
			            ei = sn_max_temp/(ei*convf)*ei
		            endif

								ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
								&           solnData(VELY_VAR,i,j,k)**2 + &
								&           solnData(VELZ_VAR,i,j,k)**2)

								solnData(TEMP_VAR,i,j,k)  = tmp
								solnData(DENS_VAR,i,j,k)  = rho
								solnData(EINT_VAR,i,j,k)  = ei
								solnData(ENER_VAR,i,j,k)  = ei + ek
							endif
						enddo !z
					enddo !y
				enddo !x

!  crank changed state variables through EOS, check it this is the right EOSMODE, again only changed blockIDs should be called
				call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
!  clean up memory 
				call Grid_releaseBlkPtr(blockID,solnData)

				deallocate(xCoord)
				deallocate(yCoord)
				deallocate(zCoord)

				deallocate(dx)
				deallocate(dy)
				deallocate(dz)

			enddo ! block loop

! SB output, always worked on first
			if(nSBw .gt. 0) then
				call WriteSNFeedback(he_nSN-nSBw+1, 3, dr_nstep, time, x0heat, y0heat, z0heat, r_exp, mass, ntracer)
				nSBw = nSBw - 1
! no idea why there is a cycle here anymore (write more documentation)
! should not be needed 
				cycle
			else
				nSN = nSN - 1
! could have multiple calls fix
				call WriteSNFeedback(he_nSN-nSN, SNType(l), dr_nstep, time, x0heat, y0heat, z0heat, r_exp, mass, ntracer)
			endif

		enddo ! SN loop

	endif !he_exp_flag

  !====================================================
  ! finalize
  !====================================================

! MKRJ - use the current seeds to generate random numbers in the next timestep
!	call random_seed(GET=he_seed(1:he_seedsize))

!	TODO only call when SN are used, SN switch here
	call MPI_Barrier (MPI_Comm_World, error)
	call Timers_stop("SN")
	return

contains

  !====================================================
  ! output routine for supernovae, writes to SNfeedback.dat
  !====================================================

  subroutine WriteSNFeedback(nSN,type,ndt,time,x,y,z,radius,mass,ntracer)

    use SN_data, ONLY :  he_meshMe
!  use Driver_data, ONLY : dr_globalMe

    implicit none

#include "constants.h"
#include "Flash.h"

    real, intent(IN)   :: time, x,y,z,radius,mass
    integer, intent(IN):: nSN,type,ndt,ntracer ! 1 is I 2 is II and 3 is SN in SB

    integer, parameter :: funit_evol = 15
    character(len=80)  :: outfile = "SNfeedback.dat"
    integer            :: i

!   as all processors calculate the SN in lockstep, only the Master process has to output data.
    if (he_meshMe .NE. MASTER_PE) return

    open(funit_evol, file=trim(outfile), position='APPEND')
!do i = 1, localnpf

    write(funit_evol,'(4(1X,I16),6(1X,ES16.9))') &
			 nSN,		  &
	     type,    &
	     ndt,     &
			 ntracer, &
	     time,    &
	     x,       &
	     y,       &
	     z,       &
	     radius,  &
       mass

    close(funit_evol)

    return
  end subroutine WriteSNFeedback

end subroutine SN
