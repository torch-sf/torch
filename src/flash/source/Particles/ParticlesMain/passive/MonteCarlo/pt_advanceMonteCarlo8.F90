!!****if* source/Particles/ParticlesMain/passive/Euler/pt_advanceEuler_passive
!!
!! NAME
!!
!!  pt_advanceEuler_passive
!!
!! SYNOPSIS
!!
!!  call pt_advanceMonteCarlo_passive(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!! (1) grab mass fluxes in different directions.
!! (2) loop over local tracer particle and call random number
!! (3) compare each normalized mass flux to the random number
!! (4) if there is a tracer particle transfer, move particle
!! (5) check if neighbour is different refinement, if so take measures
!! (6) finish loop and do MPI step or write integer stepping so the Flash4 routines take care of it
!! possible issues: randomize which face is chosen first
!! 								 if higher refinement distribution of tracers is not clear
!!  
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind   -- index for type into pt_typeInfo data structure
!!  
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!***

!! TODO make MPI safe, i.e. only use local information | done
!! TODO clean up
!! TODO properly test high dynamic range simulations with lots of de/refinement going on
!! TODO Trim down for only outgoing fluxes (throw away neg. velocity branches)
!! TODO make interface for general MC function based on hydro variables
!! TODO TODO at the moment uses variables from SN driving, write not SN independent implementation 
!! TODO check signs in outflows

!===============================================================================
!#define DEBUG_MC
subroutine pt_advanceMonteCarlo (dtOld, dtNew, particles, p_count, ind)
    
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
       								pt_posAttrib, pt_velNumAttrib,pt_velAttrib, &
											pt_xmin, pt_xmax, pt_ymin, pt_ymax, pt_zmin, pt_zmax, &
											therm_NumAttrib, therm_Attrib

  use Grid_interface, ONLY : Grid_mapMeshToParticles, &
														 Grid_getPointData, &
														 Grid_putPointData, Grid_getBlkIndexLimits, &
														 Grid_getFluxData, Grid_releaseBlkPtr, &
														 Grid_getBlkPtr

	use Grid_data, ONLY : gr_delta,gr_xflx, gr_yflx, gr_zflx
  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes

	use tree, ONLY : lrefine, bnd_box, coord!, lnblocks
	use pt_advanceMonteCarlo_data
  use mtmod

  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  	:: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles
  integer :: i, nstep, j, numNegh, newBID
	integer :: k, blockType, blockproc
  real, pointer, dimension(:,:,:,:) :: solnData
  integer                    :: part_props = NPART_PROPS
	integer                    :: blockID, currID
	integer                    :: face, dir
  integer, dimension(MDIM)   :: facedir
  integer, dimension(2*MDIM) :: order
	real 											 :: dens, mflx, temp!, stepsize,
	real 											 :: x, y, z
	real 											 :: xpos, ypos, zpos
	integer                    :: xID, yID, zID
	real, dimension(NDIM)      :: vel
	integer, dimension(NDIM)   :: xyzID
	logical                    :: move
  real, dimension(NDIM)	     :: gridspace
  integer, dimension(MDIM,4) :: negh
  integer, dimension(MDIM,4) :: neghCornerID
	real											 :: mtot, Pm, Prnd, area 
	integer 									 :: p1, p2, lref, rnd
	real 											 :: lowBx,lowBy,lowBz, d1, d2, ds
	integer 									 :: child, swp
  integer 									 :: iSize, jSize, kSize
  integer 									 :: xIoffset,yIoffset,zIoffset
	integer 									 :: direction

! should be set from parameters, or looked up properly 
	logical 									 :: xperiodic = .false., yperiodic  = .false., zperiodic = .false.

!	integer,dimension(lnblocks):: blockList
  integer, dimension(MDIM)   :: dataSize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimitsGC, blkLimits

	integer, dimension(2*MDIM) :: tmpA !indexpool
	integer, dimension(2*MDIM) :: tmpB !reduced indexpool

  real, pointer, dimension(:,:,:,:) :: blockData

!  real, allocatable,dimension(:,:,:) :: totXFlx
  real, allocatable,dimension(:,:,:) :: totXFlx
  real, allocatable,dimension(:,:,:) :: totYFlx
  real, allocatable,dimension(:,:,:) :: totZFlx
  real, allocatable,dimension(:,:,:) :: leftMass

!!------------------------------------------------------------------------------  
  integer :: mapType
!!------------------------------------------------------------------------------

! TODO re-calculate new zone centered position inside the newblock to not accrue any num. errors
! use indices and block position and gridspacing for that (probably not needed but juuuust to 
! make sure?? | done


! access density flux and face size
! +- dt * velxyz * density * facesize  [s*cm/s*g/cm^3*cm^2] = [g_xzy]
! OR use flux data from hydro solver (hope it's up to date)

! get hydrovariables and calculate fluxes, just three cardinal directions
! velocity, density and face area needed

  ! Don't do anything if runtime parameter isn't set
  if (.not. useParticles) return
  
  	!call gr_ensureValidNeighborInfo(10)
		mapType=pt_typeInfo(PART_MAPMETHOD,ind)
		currID = -1 !particles(BLK_PART_PROP,1)
		mtot = -1d0
     
  ! Update the particle positions.
		do i = 1, p_count

! get blockid
			blockID   = particles(BLK_PART_PROP,i)

! new block
			if(blockID .ne. currID ) then

! this is real ugly, but there's no way to know how many particles are in the same block,
! or memory allocation would have to be done for each particle in the same block
! could also be done once if size is known, which usually does not change but eh
				if(i .gt. 1) then
	  			deallocate(totXFlx)
	  			deallocate(totYFlx)
				  deallocate(totZFlx)
				  deallocate(leftMass)
	        call Grid_releaseBlkPtr(blockID,blockData,CENTER)
				endif

				gridspace = gr_delta(1:MDIM,lrefine(blockID))
				call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

! Get block pointer
				call Grid_getBlkPtr(blockID,blockData,CENTER)

				iSize	=	blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1 !number of zones in x
				jSize	=	blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1 !number of zones in y
				kSize	=	blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1 !number of zones in z

				xIoffset = blkLimits(LOW,IAXIS)-1
				yIoffset = blkLimits(LOW,JAXIS)-1
				zIoffset = blkLimits(LOW,KAXIS)-1				

! space
! actually just need density flux but Grid_getFluxData returns all fluxes
				allocate(totXFlx(iSize,jSize,kSize))
				allocate(totYFlx(iSize,jSize,kSize))
				allocate(totZFlx(iSize,jSize,kSize))

				totXFlx = 0.0
				totYFlx = 0.0
				totZFlx = 0.0

! should all be zero
				allocate(leftMass(iSize,jSize,kSize))

! hope that gridspace is right
				leftMass = blockData( DENS_VAR,:,:,:)*gridspace(1)*gridspace(2)*gridspace(3)
				totXFlx  = blockData(XMFLX_VAR,:,:,:)
				totYFlx  = blockData(YMFLX_VAR,:,:,:)
				totZFlx  = blockData(ZMFLX_VAR,:,:,:)

!				dataSize(1) = iSize
!				dataSize(2) = jSize
!				dataSize(3) = kSize
! get flux data once for all particles in block instead of several times for each particle

! get X fluxes
!				call Grid_getFluxData(blockID,IAXIS,totXFlx,dataSize)
!				print*,blockData(XMFLX_VAR,7,7,5:10)*gridspace(2)*gridspace(3)*dtnew
!				print*,totXFlx(DENS_FLUX,7,7,5:10)*gridspace(2)*gridspace(3)*dtnew

!! get Y fluxes
!				call Grid_getFluxData(blockID,JAXIS,totYFlx,dataSize)
! get Z fluxes
!				call Grid_getFluxData(blockID,KAXIS,totZFlx,dataSize)
				currID = blockID
				mtot = -1d0
			endif

#ifdef DEBUG_MC
		print*,'At particle', particles(TAG_PART_PROP,i)
#endif

			newBid = blockID

			gridspace = gr_delta(1:MDIM,lrefine(blockID))

			lowbx = bnd_box(LOW,IAXIS,blockID)
			lowby = bnd_box(LOW,JAXIS,blockID)
			lowbz = bnd_box(LOW,KAXIS,blockID)

! block bounds
			xID = particles(VELX_PART_PROP,i)
			yID = particles(VELY_PART_PROP,i)
			zID = particles(VELZ_PART_PROP,i)

! check integer zone position from refinement level stored in mass property
! if there is discrepancy adjust position and integer values
! this is caused by changing refinement
! could also be done in AMR routines but would be more complicated
			lref = int(particles(LREF_PART_PROP,i))

			if( lref .ne. lrefine(blockID) ) then

#ifdef DEBUG_MC
			print*,'different ref'
			print*,lref,particles(LREF_PART_PROP,i)
			print*,lrefine(blockID),x,y,z

			print*,'refine',particles(TAG_PART_PROP,i)
			print*,'ids',xid,yid,zid
#endif

				if( lref .lt. lrefine(blockID) ) then
! throw dice to see which child block gets the tracer particle
					child = floor(grnd()*8) + 1
! TODO implement tracer particle splitting
	  		  gridspace = gridspace/2d0

! this maps index to new indices, always lower one
! from 1->1,2;  2->3,4; 3->5,6; 4->7,8 for le 4
! from 5->1,2;  6->3,4; 7->5,6; 8->7,8 for ge 5
			  if     (xID .le. 4 .and. yID .le. 4 .and. zID .le. 4 ) then
! x 1..4 to 1...8, y 1..4 to 1...8, z 1..4 to 1...8
				  xID = xID*2-1
					yID = yID*2-1
					zID = zID*2-1

#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .ge. 5 .and. yID .le. 4 .and. zID .le. 4 ) then
! x 5..8 to 1...8, y 1..4 to 1...8, z 1..4 to 1...8
! map x 	
				  xID = (xID-4)*2-1
				  yID = yID*2-1
				  zID = zID*2-1

#ifdef DEBUG_MC
				  print*,'mapped ids',xID,yID,zID
#endif
		    else if(xID .ge. 5 .and. yID .ge. 5 .and. zID .le. 4 ) then
! x 5..8 to 1...8, y 5..8 to 1...8, z 1..4 to 1...8
					xID = (xID-4)*2-1
					yID = (yID-4)*2-1
					zID = zID*2-1
#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .le. 4 .and. yID .ge. 5 .and. zID .le. 4 ) then
! x 1..4 to 1...8, y 5..8 to 1...8, z 1..4 to 1...8
					xID = xID*2-1
					yID = (yID-4)*2-1
					zID = zID*2-1

#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .le. 4 .and. yID .le. 4 .and. zID .ge. 5 ) then
! x 1..4 to 1...8, y 1..4 to 1...8, z 5..8 to 1...8
					xID = xID*2-1
					yID = yID*2-1
					zID = (zID-4)*2-1

#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .ge. 5 .and. yID .le. 4 .and. zID .ge. 5 ) then
! x 5..8 to 1...8, y 1..4 to 1...8, z 5..8 to 1...8
					xID = (xID-4)*2-1
					yID = yID*2-1
					zID = (zID-4)*2-1							
#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .ge. 5 .and. yID .ge. 5 .and. zID .ge. 5 ) then
! x 5..8 to 1...8, y 5..8 to 1...8, z 5..8 to 1...8
					xID = (xID-4)*2-1
					yID = (yID-4)*2-1
					zID = (zID-4)*2-1
#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				else if(xID .le. 4 .and. yID .ge. 5 .and. zID .ge. 5 ) then
! x 1..4 to 1...8, y 5..8 to 1...8, z 5..8 to 1...8
					xID = xID*2-1
					yID = (yID-4)*2-1
					zID = (zID-4)*2-1
#ifdef DEBUG_MC
					print*,'mapped ids',xID,yID,zID
#endif
				endif

#ifdef DEBUG_MC
				print*,'child',child
#endif

! this change physical position so the particle will be at a center of a refined zone
	      select case (child)
					case(1)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

					case(2)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

						xID = xID + 1
					case(3)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

						yID = yID + 1
					case(4)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

						xID = xID + 1
						yID = yID + 1
					case(5)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

						zID = zID + 1
					case(6)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

						xID = xID + 1
						zID = zID + 1
					case(7)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

						yID = yID + 1
						zID = zID + 1
					case(8)
						particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
						particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
						particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

						xID = xID + 1
						yID = yID + 1
						zID = zID + 1
				end select			

#ifdef DEBUG_MC
				print*,(particles(POSX_PART_PROP,i)-bnd_box(LOW,IAXIS,blockID))/(gridspace(1)*2d0)
				print*,(particles(POSY_PART_PROP,i)-bnd_box(LOW,JAXIS,blockID))/(gridspace(2)*2d0)
				print*,(particles(POSZ_PART_PROP,i)-bnd_box(LOW,KAXIS,blockID))/(gridspace(3)*2d0)

				print*,particles(POSX_PART_PROP,i),gridspace(1)*2d0*(xID-0.5)+lowbx
				print*,particles(POSY_PART_PROP,i),gridspace(2)*2d0*(yID-0.5)+lowby
				print*,particles(POSZ_PART_PROP,i),gridspace(3)*2d0*(zID-0.5)+lowbz
#endif
	
				particles(VELX_PART_PROP,i) = xID
				particles(VELY_PART_PROP,i) = yID
				particles(VELZ_PART_PROP,i) = zID
				particles(LREF_PART_PROP,i) = lrefine(blockID)

! reset gridspace
		    gridspace = gr_delta(1:MDIM,lrefine(blockID))

! lower ref block
			else if( lref .gt. lrefine(blockID) ) then

#ifdef DEBUG_MC
				print*,'oldids',xid,yid,zid
#endif

#ifdef DEBUG_MC
				print*,'derefine',particles(TAG_PART_PROP,i)
#endif

! recalculate indices and use them to adjust position

				x = particles(POSX_PART_PROP,i) - bnd_box(LOW,IAXIS,blockID)
				y = particles(POSY_PART_PROP,i) - bnd_box(LOW,JAXIS,blockID)
				z = particles(POSZ_PART_PROP,i) - bnd_box(LOW,KAXIS,blockID)

! this is already coarser gridspace
				xid = floor(x/gridspace(1))+1
				yid = floor(y/gridspace(2))+1
				zid = floor(z/gridspace(3))+1

#ifdef DEBUG_MC
				print*,x/gridspace(1)
				print*,y/gridspace(2)
				print*,z/gridspace(3)
				print*,'newids',xid,yid,zid
#endif

				particles(VELX_PART_PROP,i) = xID
				particles(VELY_PART_PROP,i) = yID
				particles(VELZ_PART_PROP,i) = zID

				particles(POSX_PART_PROP,i) = bnd_box(LOW,IAXIS,blockID) + (xid-0.5)*gridspace(1)
				particles(POSY_PART_PROP,i) = bnd_box(LOW,JAXIS,blockID) + (yid-0.5)*gridspace(2)
				particles(POSZ_PART_PROP,i) = bnd_box(LOW,KAXIS,blockID) + (zid-0.5)*gridspace(3)

				particles(LREF_PART_PROP,i) = lrefine(blockID)
!		    gridspace = gr_delta(1:MDIM,lrefine(blockID))
			endif
		endif

		x = particles(POSX_PART_PROP,i)
		y = particles(POSY_PART_PROP,i)
		z = particles(POSZ_PART_PROP,i)

		xpos = x
		ypos = y
		zpos = z

! IDs should be up to date
		xyzID(1) = xID
 		xyzID(2) = yID
	  xyzID(3) = zID

!		call Grid_getPointData(blockID, CENTER, DENS_VAR, INTERIOR, xyzID, dens)
!		call Grid_getPointData(blockID, CENTER, TEMP_VAR, INTERIOR, xyzID, temp)
		call Grid_getPointData(blockID, CENTER, VELX_VAR, INTERIOR, xyzID, vel(1))
		call Grid_getPointData(blockID, CENTER, VELY_VAR, INTERIOR, xyzID, vel(2))
		call Grid_getPointData(blockID, CENTER, VELZ_VAR, INTERIOR, xyzID, vel(3))

#ifdef DEBUG_MC
		if(abs(particles(POSX_PART_PROP,i)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,newBID))) .gt. 1e5 .or. & 
		   abs(particles(POSY_PART_PROP,i)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,newBID))) .gt. 1e5 .or. &
		   abs(particles(POSZ_PART_PROP,i)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,newBID))) .gt. 1e5) then

		  print*,'position wrong before'
		  print*,particles(POSX_PART_PROP,i)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,newBID))
		  print*,particles(POSY_PART_PROP,i)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,newBID))
		  print*,particles(POSZ_PART_PROP,i)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,newBID))						
			print*,particles(TAG_PART_PROP,i)

		  print*,(particles(POSX_PART_PROP,i)-bnd_box(LOW,IAXIS,newBID))/gridspace(1),xid
		  print*,(particles(POSY_PART_PROP,i)-bnd_box(LOW,JAXIS,newBID))/gridspace(2),yid
		  print*,(particles(POSZ_PART_PROP,i)-bnd_box(LOW,KAXIS,newBID))/gridspace(3),zid
			print*,gridspace
		  stop
		endif
#endif

!=========================================================================================
!=========================================================================================
! after check refinement


! use centered density
! first tracer in block
!    if(mtot .lt. 0) then
!  		mtot = dens*gridspace(1)*gridspace(2)*gridspace(3)
!			leftMass = dens*gridspace(1)*gridspace(2)*gridspace(3)
!    endif

! update hydro variables independent if particle moves or not.
!		particles(NSN_PART_PROP,i)  = temp
!		particles(MASS_PART_PROP,i) = dens

! generate a randomized order for checking which face to check first
! scale to 1 to 6
! not sure if order has to be randomized but can't hurt
		order = 0
		do j=0, 4
			if(j .eq. 0) then
				tmpA = (/1,2,3,4,5,6/)
			endif

! position in face array
			rnd 				= floor(grnd()*(6-j)) + 1
			order(j+1)	= tmpA(rnd)
! move value to the edge and zero out
			swp 		  =	tmpA(6-j)
			tmpA(rnd) = swp  !cshift(tmpA, shift=int(rnd))
			tmpA(6-j)= 0

!		tmpB(1:5-j) = tmpA(1:6-j-1)
		enddo
		order(6) = tmpA(1)

		mtot = leftMass(xID+xIoffset,yID+yIoffset,zID+zIoffset)

		do j = 1, 2*NDIM
			dir  = order(j)

			select case (dir)
! neg x 
			  case (1)

! dens flux is g/cm^2*cm^2*s = g 
					area = gridspace(2)*gridspace(3)
					mflx = totXFlx(xID+xIoffset,yID+yIoffset,zID+zIoffset)*area*dtNew

! determine which face it moves out of for later switches
					direction = -1
					dir = 1

! check sign for outflow
					if(mflx .gt. 0d0) then
						cycle
					endif
! pos x 
			  case (2)
! dens flux is g/cm^2*cm^2*s = g 
					area = gridspace(2)*gridspace(3)
					mflx = totXFlx(xID+xIoffset+1,yID+yIoffset,zID+zIoffset)*area*dtNew
					direction = 1
					dir = 1

! check sign for outflow
					if(mflx .lt. 0d0) then
						cycle
					endif
! neg y 
			  case (3)
					area = gridspace(1)*gridspace(3)
					mflx = totYFlx(xID+xIoffset,yID+yIoffset,zID+zIoffset)*area*dtNew  ! units [g/cm^3*cm/s*s*cm^2] = [g]
! determine which face it moves out of for later switches
					direction = -1
					dir = 2

! check sign for outflow
					if(mflx .gt. 0d0) then
						cycle
					endif
! pos y
			  case (4)
					area = gridspace(1)*gridspace(2)
					mflx = totYFlx(xID+xIoffset,yID+yIoffset+1,zID+zIoffset)*area*dtNew  ! units [g/cm^3*cm/s*s*cm^2] = [g]

! determine which face it moves out of for later switches
					direction = 1
					dir = 2

! check sign for outflow
					if(mflx .lt. 0d0) then
						cycle
					endif
! neg z
			  case (5)

					area = gridspace(1)*gridspace(2)
					mflx = totZFlx(xID+xIoffset,yID+yIoffset,zID+zIoffset)*area*dtNew  ! units [g/cm^3*cm/s*s*cm^2] = [g]          

! determine which face it moves out of for later switches
					direction = -1
					dir = 3

! check sign for outflow
					if(mflx .gt. 0d0) then
						cycle
					endif

! pos z
			  case (6)
					area = gridspace(1)*gridspace(2)
					mflx = totZFlx(xID+xIoffset,yID+yIoffset,zID+zIoffset+1)*area*dtNew  ! units [g/cm^3*cm/s*s*cm^2] = [g]
! determine which face it moves out of for later switches
					direction = 1
					dir = 3

! check sign for outflow
					if(mflx .lt. 0d0) then
						cycle
					endif
			end select

!			mflx = dens*vel(dir)*dtNew*area  ! units [g/cm^3*cm/s*s*cm^2] = [g]

! switch face as there's inflow
! other face should have outflow 
! in and outflow from same zone should not be possible
			Pm   = abs(mflx/mtot)
!			print*, Pm, mflx,mtot

! update mass in zone, next transports more likely
! didn't ride the first transfer?
! try again

! outflow
			mtot = mtot - abs(mflx)
!			leftMass(xID+xIoffset,yID+yIoffset,zID+zIoffset) = leftMass(xID+xIoffset,yID+yIoffset,zID+zIoffset) - abs(mflx)

! throw some dice
			Prnd = grnd()

! Prnd < Pm, so Pm = 0.1 then only 10% chance of particle moving
			if(Prnd .le. Pm) then

! find zone in which the particle resides
! use velocity info for that, so blk + vel gives full location information
  			xyzID(1) = particles(VELX_PART_PROP,i)
	  		xyzID(2) = particles(VELY_PART_PROP,i)
			  xyzID(3) = particles(VELZ_PART_PROP,i)

! if edge of block check for neighbours, could just be done with eq 1 instead of le
! this also selects the face
!//////////////// neg. direction
				if( direction .lt. 0 .and. xyzID(dir) .le. 1) then

! not totally sure of face direction vectors, 3 should point to negative direction
				  select case (dir)
					  case (1)
						  facedir = (/1,2,2/)
					  case (2)
						  facedir = (/2,1,2/)
					  case (3)
						  facedir = (/2,2,1/)
				  end select

					numnegh = 0

!	find neighbour information and move
  	      call gr_ptFindNegh(blockID,facedir,negh,neghCornerID,numNegh)

! higher refined region
		      if(numNegh .gt. 1) then

#ifdef DEBUG_MC
		        if(numNegh .ne. 4) then
						  print*,'not 4 neigbours'
						else
						  print*,'high ref',	particles(TAG_PART_PROP,i)
						endif
#endif
						particles(LREF_PART_PROP,i) = lrefine(blockID)+1

! have to decide between 4 new zones
!	throw some dice
						face = floor(grnd()*4) + 1
			      newBid = -1
			      select case (dir)
! neg. x
					    case (1)
!// choose new block based on zone index
! y low:  1..4; y high: 2..8
! z low:  1..4; z high: 2..8
! 4 combinations
! z
! |----|----|
! | 3  |  4	|
! |----|----| --> j selects right face
! | 1  |  2 |
! |----|----|y

! also weak tests, could be out of lower or upper bounds
							  if     (yID .le. 4 .and. zID .le. 4) then ! 1 
					        newBID    = negh(1,1)! 
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(yID .ge. 5 .and. zID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(yID .le. 4 .and. zID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  else if(yID .ge. 5 .and. zID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								x = x - gridspace(1)/2d0

! make gridspacing smaller, independent of tree information
								gridspace = gridspace/2d0

! update lower bound
								lowbx = lowbx - 8*gridspace(1)

								x = x - gridspace(1)/2d0					
								xpos = x

!								if(xperiodic) then
!								  if(xpos .lt. pt_xmin ) then
!#ifdef DEBUG_MC
!									  print*, 'periodic high ref boundary x'
!#endif
!										xpos = xpos + (pt_xmax-pt_xmin)
!										x = xpos
!									endif
!								endif
! move to local coordinate system
! tree independent
								x = x - lowbx
!								xID = floor(x/gridspace(1))+1
                xID = 8

!// choose new zone 
			      		select case (face)
! 3 dimensions, 1 is fixed 
! the other two give 4 options:
! low low, max low, max max, low max 
					        case (1)
! y and z axis
								    y = y - gridspace(2)/2d0
								    z = z - gridspace(3)/2d0
										ypos = y
										zpos = z
! move to local coordinate system
! tree independent
										y = y - lowby
										z = z - lowbz

! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x-',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
										print*,'case --'
#endif
					        case (2)
! y and z axis
									  y = y + gridspace(2)/2d0
									  z = z - gridspace(3)/2d0
										ypos = y
										zpos = z
! tree independent
										y = y - lowby
										z = z - lowbz
						
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x-',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
										print*,'case +-'
#endif
					        case (3)
! y and z axis
								    y = y + gridspace(2)/2d0
								    z = z + gridspace(3)/2d0
										ypos = y
										zpos = z
! tree independent
										y = y - lowby
										z = z - lowbz
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x-',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
										print*,'case --'
#endif

					        case (4)
! y and z axis
									  y = y - gridspace(2)/2d0
									  z = z + gridspace(3)/2d0
										ypos = y
										zpos = z
! tree independent
										y = y - lowby
										z = z - lowbz										
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x-',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
										print*,'case -+'
#endif
								end select
!//

! /3,4,7,8/, /xmin,zmin,xmax,zmin,xmin,zmin,xmax,zmax/
! neg. y
					    case (2)
!// choose new block based on zone index
! x low:  1..4; x high: 2..8
! z low:  1..4; z high: 2..8
! also weak tests, could be out of lower or upper bounds
							  if     (xID .le. 4 .and. zID .le. 4) then ! 1
					        newBID    = negh(1,1)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(xID .ge. 5 .and. zID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(xID .le. 4 .and. zID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  else if(xID .ge. 5 .and. zID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								y = y - gridspace(2)/2d0
! new gridspacing of higher refined grid
								gridspace = gridspace/2d0

! update lower bound
								lowby = lowby - 8*gridspace(2)

								y = y - gridspace(2)/2d0
   							ypos = y

!								if(yperiodic) then
!								  if(ypos .lt. pt_ymin ) then
!#ifdef DEBUG_MC
!									  print*, 'periodic high ref boundary y'
!#endif
!										ypos = ypos + (pt_ymax-pt_ymin)
!										y = ypos
!									endif
!								endif

! move to local coordinate system
								y = y - lowby
!								yID = floor(y/gridspace(2))+1
								yID = 8

! new gridspacing
! 3 dimensions, 1 is fixed 
! the other two give 4 options 
			      		select case (face)
					        case (1)
! x and z axis
										x = x - gridspace(1)/2d0
										z = z - gridspace(3)/2d0
										xpos = x
										zpos = z

										x = x - lowbx
										z = z - lowbz
! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(3),'->',zID,'z'
										print*,'case --'
#endif

					        case (2)
! y and z axis
										x = x + gridspace(1)/2d0
										z = z - gridspace(3)/2d0
										xpos = x
										zpos = z

										x = x - lowbx
										z = z - lowbz		

! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(3),'->',zID,'z'
										print*,'case +-'
#endif
					        case (3)
! y and z axis
										x = x + gridspace(1)/2d0
										z = z + gridspace(3)/2d0
										xpos = x
										zpos = z
										x = x - lowbx
										z = z - lowbz
! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1

#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(3),'->',zID,'z'
										print*,'case ++'
#endif

					        case (4)
! x and z axis
										x = x - gridspace(1)/2d0
										z = z + gridspace(3)/2d0
										xpos = x
										zpos = z
										x = x - lowbx
										z = z - lowbz

! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(3),'->',zID,'z'
										print*,'case -+'
#endif
								end select
! /5,6,7,8/, /xmin,ymin,xmax,ymin,xmin,zmax,xmax,zmax/
! neg. z
					    case (3)
!// choose new block based on zone index
! x low:  1..4; x high: 2..8
! y low:  1..4; y high: 2..8
! also weak tests, could be out of lower or upper bounds
							  if     (xID .le. 4 .and. yID .le. 4) then ! 1 
					        newBID    = negh(1,1)
									lowbz = bnd_box(LOW,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
							  else if(xID .ge. 5 .and. yID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowbz = bnd_box(LOW,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowby = bnd_box(LOW,JAXIS,blockID)
							  else if(xID .le. 4 .and. yID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowbz = bnd_box(LOW,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
							  else if(xID .ge. 5 .and. yID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowbz = bnd_box(LOW,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								z = z - gridspace(3)/2d0

! new gridspacing of higher refined grid
								gridspace = gridspace/2d0

! update lower bound
								lowbz = lowbz - 8*gridspace(3)

								z = z - gridspace(3)/2d0
								zpos = z

!								if(zperiodic) then
!								  if(zpos .lt. pt_zmin ) then
!
!#ifdef DEBUG_MC
!									  print*, 'periodic high ref boundary z'
!#endif
!										zpos = zpos + (pt_zmax-pt_zmin)
!										z = zpos
!									endif
!								endif

! move to local coordinate system
								z = z - lowbz
!								zID = floor(z/gridspace(3))+1
! fix to be 8
								zID = 8

! new gridspacing
! 3 dimensions, 1 is fixed 
! the other two give 4 options 
			      	  select case (face)
					        case (1)
! x and y axis
										x = x - gridspace(1)/2d0
										y = y - gridspace(2)/2d0
										xpos = x
										ypos = y
										x = x - lowbx
										y = y - lowby

! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1

#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
										print*,'case --'
#endif

					        case (2)
! y and z axis
										x = x + gridspace(1)/2d0
										y = y - gridspace(2)/2d0
										xpos = x
										ypos = y
										x = x - lowbx
										y = y - lowby

! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
										print*,'case +-'
#endif

					        case (3)
! y and z axis
										x = x + gridspace(1)/2d0
										y = y + gridspace(2)/2d0
										xpos = x
										ypos = y
										x = x - lowbx
										y = y - lowby
! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1

#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
										print*,'case ++'
#endif

					        case (4)
! x and z axis
										x = x - gridspace(1)/2d0
										y = y + gridspace(2)/2d0
										xpos = x
										ypos = y
										x = x - lowbx
										y = y - lowby

! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(z/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z-',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
										print*,'case -+'
#endif

						  end select
				    end select
				  endif

! same or lower res
				  if(numnegh .eq. 1) then
				    newBID = negh(1,1)
! move to face
			      select case (dir)
! neg. x
					    case (1)
					      x = x - gridspace(1)/2d0
								ds = bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID)
! neg. y
					    case (2)
							  y = y - gridspace(2)/2d0
								ds = bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID)
! neg. z
					    case (3)
							  z = z - gridspace(3)/2d0
								ds = bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID)
						end select

! not sure this works with MPI
						if( negh(3,1) .lt. 0 )then
#ifdef DEBUG_MC
							print*,'low res neg'
#endif
	 						particles(LREF_PART_PROP,i) = lrefine(blockID)-1
! find new indices for the other two dimensions
			        select case (dir)
! neg. x
					      case (1)
! divide block position by the block size
! if odd then left block in parent, if even then right block in parent
! there is probably an easier by using the hiblockIDert curve
									ds = 2d0*(bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID))
!									lowbx = bnd_box(LOW,IAXIS,blockID) - 2d0*d1

									d1 = bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID)
									d2 = bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID)

									p1 = mod(floor((coord(2,blockID)-pt_ymin)/d1)+1,2)
									p2 = mod(floor((coord(3,blockID)-pt_zmin)/d2)+1,2)

									if(p1 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block y'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y

										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										yid = ceiling(yid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block y'
#endif

										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										lowby = bnd_box(LOW,JAXIS,blockID) - d1

										yid = ceiling(yid/2d0)+4
									endif

									if(p2 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										zid = ceiling(zid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										lowbz = bnd_box(LOW,KAXIS,blockID) - d2
										zid = ceiling(zid/2d0)+4
									endif
! neg. y
					      case (2)
									ds = 2d0*(bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID))
!									lowby = bnd_box(LOW,JAXIS,blockID) - 2d0*d1

									d1 = bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID)
									d2 = bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID)

									p1 = mod(floor((coord(1,blockID)-pt_xmin)/d1)+1,2)
									p2 = mod(floor((coord(3,blockID)-pt_zmin)/d2)+1,2)

									if(p1 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block x'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										xid = ceiling(xid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block x'
#endif

										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										lowbx = bnd_box(LOW,IAXIS,blockID) - d1
										xid = ceiling(xid/2d0)+4
									endif

									if(p2 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										zid = ceiling(zid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										lowbz = bnd_box(LOW,KAXIS,blockID) - d2
										zid = ceiling(zid/2d0)+4
									endif

! neg. z
					      case (3)
									ds = 2d0*(bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID))

									d2 = bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID)
									d1 = bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID)

									p2 = mod(floor((coord(1,blockID)-pt_xmin)/d2)+1,2)
									p1 = mod(floor((coord(2,blockID)-pt_ymin)/d1)+1,2)

									if(p1 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block y'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										yid = ceiling(yid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block y'
#endif

										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										lowby = bnd_box(LOW,JAXIS,blockID) - d2
										yid = ceiling(yid/2d0)+4
									endif

									if(p2 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block x'
#endif
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										xid = ceiling(xid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block x'
#endif
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										lowbx = bnd_box(LOW,IAXIS,blockID) - d1
										xid = ceiling(xid/2d0)+4
									endif

						  end select	
! new bigger gridspacing												
							gridspace = gridspace*2d0
						endif
! move face
			      select case (dir)
! neg. x
					    case (1)

							  x   = x - gridspace(1)/2d0
								xpos = x
								lowbx = bnd_box(LOW,IAXIS,blockID) - ds

								if(xperiodic) then
								  if(xpos .lt. pt_xmin ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary x'
#endif
										xpos = xpos + (pt_xmax-pt_xmin)
										x = xpos
									endif
								endif
! don't need this actually
! move to local coordinate system
								xID = 8				
#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(1),'->',xID,'x-',i,pt_xmax,pt_xmin
#endif

								if(abs(xyzID(1) - xID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif

! neg. y
					    case (2)
							  y    = y - gridspace(2)/2d0
								ypos = y
								lowby = bnd_box(LOW,JAXIS,blockID) - ds

								if(yperiodic) then
								  if(ypos .lt. pt_ymin ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary y'
#endif
										ypos = ypos + (pt_ymax-pt_ymin)
										y = ypos

									endif
								endif

								yID = 8
#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(2),'->',yID,'y-',i,pt_ymax,pt_ymin
#endif
								if(abs(xyzID(2) - yID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif

! neg. z
					    case (3)
							  z    = z - gridspace(3)/2d0
								zpos = z
								lowbz = bnd_box(LOW,KAXIS,blockID) - ds

								if(zperiodic) then
								  if(zpos .lt. pt_zmin ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary z'
#endif
										zpos = zpos + (pt_zmax-pt_zmin)
										z = zpos
									endif
								endif
								zID = 8
#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(3),'->',zID,'z-',i,pt_zmax,pt_zmin
#endif

								if(abs(xyzID(3) - zID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif
									
						end select

			    else if(numnegh .eq. 0) then
#ifdef DEBUG_MC
		        print*,'no neighbour'
#endif			
				    newBid = -20
			    endif
!////////////////
!//////////////// pos. direction
! this also selects the face
				else if( direction .gt. 0 .and. xyzID(dir) .ge. 8) then

#ifdef DEBUG_MC
		    if(numNegh .ne. 4) then
				  print*,'pos boundary reached'
				endif
#endif

				  select case (dir)
					  case (1)
						  facedir = (/3,2,2/)
					  case (2)
						  facedir = (/2,3,2/)
					  case (3)
						  facedir = (/2,2,3/)
				  end select

					numnegh = 0

!	find neighbour information and move
  	      call gr_ptFindNegh(blockID,facedir,negh,neghCornerID,numNegh)

! higher refined region
		      if(numNegh .gt. 1) then
#ifdef DEBUG_MC
		        if(numNegh .ne. 4) then
						  print*,'not 4 neigbours'
						else
						  print*,'high ref',	particles(TAG_PART_PROP,i)
						endif

#endif
	 					particles(LREF_PART_PROP,i) = lrefine(blockID)+1
! have to decide between 4 new zones
!	throw some dice
						face = floor(grnd()*4) + 1

			      newBid = -1
			      select case (dir)
! pos. x
					    case (1)
!// choose new block based on zone index
! y low:  1..4; y high: 2..8
! z low:  1..4; z high: 2..8
! 4 combinations
! z
! |----|----|
! | 3  |  4	|
! |----|----| --> j selects right face
! | 1  |  2 |
! |----|----|y

! also weak tests, could be out of lower or upper bounds
							  if     (yID .le. 4 .and. zID .le. 4) then ! 1 
					        newBID    = negh(1,1)! 
									lowbx = bnd_box(HIGH,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(yID .ge. 5 .and. zID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowbx = bnd_box(HIGH,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(yID .le. 4 .and. zID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowbx = bnd_box(HIGH,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  else if(yID .ge. 5 .and. zID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowbx = bnd_box(HIGH,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								x = x + gridspace(1)/2d0
! new gridspacing of higher refined grid
								gridspace = gridspace/2d0
								x    = x + gridspace(1)/2d0
								xpos = x

								if(xperiodic) then
								  if(xpos .gt. pt_xmax ) then
#ifdef DEBUG_MC
									  print*, 'periodic high ref boundary x'
#endif
										xpos = xpos - (pt_xmax-pt_xmin)
										x = xpos
									endif
								endif

! move to local coordinate system
								x = x- lowbx
								xID = floor(x/gridspace(1))+1

!// choose new zone 
			      		select case (face)
! 3 dimensions, 1 is fixed 
! the other two give 4 options:
! low low, max low, max max, low max
					        case (1)
! y and z axis
								    y = y - gridspace(2)/2d0
								    z = z - gridspace(3)/2d0
										ypos = y
										zpos = z
! move to local coordinate system
										y = y - lowby
										z = z - lowbz
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x+',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
#endif

					        case (2)
! y and z axis
									  y = y + gridspace(2)/2d0
									  z = z - gridspace(3)/2d0										
										ypos = y
										zpos = z

! move to local coordinate system
										y = y - lowby
										z = z - lowbz

! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x+',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
#endif
					        case (3)
! y and z axis
								    y = y + gridspace(2)/2d0
								    z = z + gridspace(3)/2d0				
										ypos = y
										zpos = z

! move to local coordinate system
										y = y - lowby
										z = z - lowbz
						
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x+',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
#endif

					        case (4)
! y and z axis
									  y = y - gridspace(2)/2d0
									  z = z + gridspace(3)/2d0										
										ypos = y
										zpos = z

! move to local coordinate system
										y = y - lowby
										z = z - lowbz
! find new indices 
								    yID = floor(y/gridspace(2))+1
								    zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(1),'->',xID,'x+',i
										print*,xyzID(2),'->',yID,'y'
										print*,xyzID(3),'->',zID,'z'
#endif

								end select
!//

! /3,4,7,8/, /xmin,zmin,xmax,zmin,xmin,zmin,xmax,zmax/
! pos. y
					    case (2)
!// choose new block based on zone index
! x low:  1..4; x high: 2..8
! z low:  1..4; z high: 2..8
! also weak tests, could be out of lower or upper bounds
							  if     (xID .le. 4 .and. zID .le. 4) then ! 1
					        newBID    = negh(1,1)
									lowby = bnd_box(HIGH,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(xID .ge. 5 .and. zID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowby = bnd_box(HIGH,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID)
							  else if(xID .le. 4 .and. zID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowby = bnd_box(HIGH,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  else if(xID .ge. 5 .and. zID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowby = bnd_box(HIGH,JAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowbz = bnd_box(LOW,KAXIS,blockID) + (bnd_box(HIGH,KAXIS,blockID)- bnd_box(LOW,KAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								y = y + gridspace(2)/2d0
								print*,y/gridspace(2),gridspace(2)
! new gridspacing of higher refined grid
								gridspace = gridspace/2d0
								y = y + gridspace(2)/2d0
								ypos = y

								if(yperiodic) then
								  if(ypos .gt. pt_ymax ) then
#ifdef DEBUG_MC
									  print*, 'periodic high ref boundary y'
#endif
										ypos = ypos - (pt_ymax-pt_ymin)
										y = ypos
									endif
								endif

								y = y -lowby
#ifdef DEBUG_MC
								print*,y/gridspace(2),gridspace(2)
#endif
								yID = floor(y/gridspace(2))+1
! new gridspacing
! 3 dimensions, 1 is fixed 
! the other two give 4 options 
			      		select case (face)
					        case (1)
! x and z axis
										x = x - gridspace(1)/2d0
										z = z - gridspace(3)/2d0
										xpos = x
										zpos = z

! move to local coordinate system
										x = x - lowbx
										z = z - lowbz
! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y+',i,y/gridspace(2)
										print*,xyzID(1),'->',xID,'x',x/gridspace(1)
										print*,xyzID(3),'->',zID,'z',z/gridspace(2)
#endif

					        case (2)
! y and z axis
										x = x + gridspace(1)/2d0
										z = z - gridspace(3)/2d0
										xpos = x
										zpos = z

! move to local coordinate system
										x = x - lowbx
										z = z - lowbz

! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y+',i,y/gridspace(2)
										print*,xyzID(1),'->',xID,'x',x/gridspace(1)
										print*,xyzID(3),'->',zID,'z',z/gridspace(2)
#endif

					        case (3)
! y and z axis
										x = x + gridspace(1)/2d0
										z = z + gridspace(3)/2d0
										xpos = x
										zpos = z
! move to local coordinate system
										x = x - lowbx
										z = z - lowbz

! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y+',i,y/gridspace(2)
										print*,xyzID(1),'->',xID,'x',x/gridspace(1)
										print*,xyzID(3),'->',zID,'z',z/gridspace(2)
#endif

					        case (4)
! x and z axis
										x = x - gridspace(1)/2d0
										z = z + gridspace(3)/2d0
										xpos = x
										zpos = z
! move to local coordinate system
										x = x - lowbx
										z = z - lowbz
! find new indices 
										xID = floor(x/gridspace(1))+1
										zID = floor(z/gridspace(3))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(2),'->',yID,'y+',i,y/gridspace(2)
										print*,xyzID(1),'->',xID,'x',x/gridspace(1)
										print*,xyzID(3),'->',zID,'z',z/gridspace(2)
#endif

								end select
! /5,6,7,8/, /xmin,ymin,xmax,ymin,xmin,zmax,xmax,zmax/
! pos. z
					    case (3)
!// choose new block based on zone index
! x low:  1..4; x high: 2..8
! z low:  1..4; z high: 2..8
! also weak tests, could be out of lower or upper bounds
							  if     (xID .le. 4 .and. yID .le. 4) then ! 1 
					        newBID    = negh(1,1)
									lowbz = bnd_box(HIGH,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID)
							  else if(xID .ge. 5 .and. yID .le. 4) then ! 2 
					        newBID    = negh(1,2)
									lowbz = bnd_box(HIGH,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowby = bnd_box(LOW,JAXIS,blockID)
							  else if(xID .le. 4 .and. yID .ge. 5) then ! 3
					        newBID    = negh(1,3)
									lowbz = bnd_box(HIGH,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID)
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
							  else if(xID .ge. 5 .and. yID .ge. 5) then ! 4
					        newBID    = negh(1,4)
									lowbz = bnd_box(HIGH,KAXIS,blockID)
									lowbx = bnd_box(LOW,IAXIS,blockID) + (bnd_box(HIGH,IAXIS,blockID)- bnd_box(LOW,IAXIS,blockID))/2d0
									lowby = bnd_box(LOW,JAXIS,blockID) + (bnd_box(HIGH,JAXIS,blockID)- bnd_box(LOW,JAXIS,blockID))/2d0
							  endif
!//
! move to boundary
								z = z + gridspace(3)/2d0
! new gridspacing of higher refined grid
								gridspace = gridspace/2d0
								z = z + gridspace(3)/2d0
								zpos = z

								if(zperiodic) then
								  if(zpos .gt. pt_zmax ) then
#ifdef DEBUG_MC
									  print*, 'periodic high ref boundary z'
#endif
										zpos = zpos - (pt_zmax-pt_zmin)
										z = zpos
									endif
								endif

! move to local coordinate system
								z = z - lowbz
								zID = floor(z/gridspace(3))+1
! new gridspacing
! 3 dimensions, 1 is fixed 
! the other two give 4 options 
			      		select case (face)
					        case (1)
! x and y axis
										x = x - gridspace(1)/2d0
										y = y - gridspace(2)/2d0
										xpos = x
										ypos = y
! move to local coordinate system
										x = x - lowbx
										y = y - lowby
! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z+',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
#endif

					        case (2)
! y and z axis
										x = x + gridspace(1)/2d0
										y = y - gridspace(2)/2d0
										xpos = x
										ypos = y
! move to local coordinate system
										x = x - lowbx
										y = y - lowby

! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z+',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
#endif

					        case (3)
! y and z axis
										x = x + gridspace(1)/2d0
										y = y + gridspace(2)/2d0
										xpos = x
										ypos = y
! move to local coordinate system
										x = x - lowbx
										y = y - lowby
! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z+',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
#endif
					        case (4)
! x and z axis
										x = x - gridspace(1)/2d0
										y = y + gridspace(2)/2d0
										xpos = x
										ypos = y
! move to local coordinate system
										x = x - lowbx
										y = y - lowby

! find new indices 
										xID = floor(x/gridspace(1))+1
										yID = floor(y/gridspace(2))+1
#ifdef DEBUG_MC
										print*,particles(TAG_PART_PROP,i)
										print*,xyzID(3),'->',zID,'z+',i
										print*,xyzID(1),'->',xID,'x'
										print*,xyzID(2),'->',yID,'y'
#endif

								end select
				    end select
				  endif

! same or lower res
				  if(numnegh .eq. 1) then
				    newBID = negh(1,1)
! move face
			      select case (dir)
! pos. x
					    case (1)
					      x = x + gridspace(1)/2d0
! pos. y
					    case (2)
							  y = y + gridspace(2)/2d0
! pos. z
					    case (3)
							  z = z + gridspace(3)/2d0
						end select

! TODO change to different check 
						if( negh(3,1) .lt. 0 )then

#ifdef DEBUG_MC
							print*,'low res pos'
#endif
	 						particles(LREF_PART_PROP,i) = lrefine(blockID)-1

			        select case (dir)
! neg. x
					      case (1)
! divide block position by the block size
! if odd then left block in parent, if even then right block in parent
! there is probably an easier by using the hiblockIDert curve
									d1 = bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID)
									d2 = bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID)

									p1 = mod(floor((coord(2,blockID)-pt_ymin)/d1)+1,2)
									p2 = mod(floor((coord(3,blockID)-pt_zmin)/d2)+1,2)
									
									lowbx = bnd_box(HIGH,IAXIS,blockID)

									if(p1 .eq. 1) then

#ifdef DEBUG_MC
										print*,'odd block y'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										yid = ceiling(yid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block y'
#endif

										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										lowby = bnd_box(LOW,JAXIS,blockID) - d1

										yid = ceiling(yid/2d0)+4
									endif


									if(p2 .eq. 1) then

#ifdef DEBUG_MC
										print*,'odd block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										zid = ceiling(zid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										lowbz = bnd_box(LOW,KAXIS,blockID) - d2

										zid = ceiling(zid/2d0)+4
									endif
! neg. y
					      case (2)
									d1 = bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID)
									d2 = bnd_box(HIGH,KAXIS,blockID)-bnd_box(LOW,KAXIS,blockID)

									p1 = mod(floor((coord(1,blockID)-pt_xmin)/d1)+1,2)
									p2 = mod(floor((coord(3,blockID)-pt_zmin)/d2)+1,2)

									lowby = bnd_box(HIGH,JAXIS,blockID)

									if(p1 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block x'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										xid = ceiling(xid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block x'
#endif

										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										lowbx = bnd_box(LOW,IAXIS,blockID) - d1
										xid = ceiling(xid/2d0)+4
									endif

									if(p2 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										zid = ceiling(zid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block z'
#endif
										if(mod(zid,2) .eq. 0) then 
										  z = z - gridspace(3)/2d0
											zpos = z
										endif

										if(mod(zid,2) .eq. 1) then 
										  z = z + gridspace(3)/2d0
											zpos = z
										endif

										lowbz = bnd_box(LOW,KAXIS,blockID) - d2

										zid = ceiling(zid/2d0)+4
									endif
! neg. z
					      case (3)
									d2 = bnd_box(HIGH,IAXIS,blockID)-bnd_box(LOW,IAXIS,blockID)
									d1 = bnd_box(HIGH,JAXIS,blockID)-bnd_box(LOW,JAXIS,blockID)

									p2 = mod(floor((coord(1,blockID)-pt_xmin)/d2)+1,2)
									p1 = mod(floor((coord(2,blockID)-pt_ymin)/d1)+1,2)

									lowbz = bnd_box(HIGH,KAXIS,blockID)

									if(p1 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block y'
#endif
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										yid = ceiling(yid/2d0)
									endif

									if(p1 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block y'
#endif

										if(mod(yid,2) .eq. 0) then 
										  y = y - gridspace(2)/2d0
											ypos = y
										endif

										if(mod(yid,2) .eq. 1) then 
										  y = y + gridspace(2)/2d0
											ypos = y
										endif

										lowby = bnd_box(LOW,JAXIS,blockID) - d1
										yid = ceiling(yid/2d0)+4
									endif

									if(p2 .eq. 1) then
#ifdef DEBUG_MC
										print*,'odd block x'
#endif
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then 
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										xid = ceiling(xid/2d0)
									endif

									if(p2 .eq. 0) then
#ifdef DEBUG_MC
										print*,'even block x'
#endif
										if(mod(xid,2) .eq. 0) then 
										  x = x - gridspace(1)/2d0
											xpos = x
										endif

										if(mod(xid,2) .eq. 1) then
										  x = x + gridspace(1)/2d0
											xpos = x
										endif

										lowbx = bnd_box(LOW,IAXIS,blockID) - d2
										xid = ceiling(xid/2d0)+4
									endif

						  end select	

! new bigger gridspacing												
							gridspace = gridspace*2d0
					  endif
! move face
			      select case (dir)
! pos. x
					    case (1)
						    x   = x + gridspace(1)/2d0
								xpos = x

								lowbx = bnd_box(HIGH,IAXIS,blockID)

								if(xperiodic) then
								  if(xpos .gt. pt_xmax ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary x'
#endif
										xpos = xpos- (pt_xmax-pt_xmin)
										x = xpos
									endif
								endif

								xID = 1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)				
								print*,xyzID(1),'->',xID,'x+',i,pt_xmax,pt_xmin
#endif
								if(abs(xyzID(1) - xID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif

! pos. y
					    case (2)
							  y = y + gridspace(2)/2d0
								ypos = y

								lowby = bnd_box(HIGH,JAXIS,blockID)

								if(yperiodic) then
								  if(ypos .gt. pt_ymax ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary y'
#endif
										ypos = ypos - (pt_ymax-pt_ymin)
										y = ypos
									endif
								endif
								yID = 1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(2),'->',yID,'y+',i,pt_ymax,pt_ymin
#endif
								if(abs(xyzID(2) - yID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif

! pos. z
					    case (3)
							  z = z + gridspace(3)/2d0
								zpos = z

								lowbz = bnd_box(HIGH,KAXIS,blockID)

								if(zperiodic) then
								  if(zpos .gt. pt_zmax ) then
#ifdef DEBUG_MC
									  print*, 'periodic boundary z'
#endif
										zpos = zpos - (pt_zmax-pt_zmin)
										z = zpos
									endif
								endif

								zID = 1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(3),'->',zID,'z+',i,pt_zmax,pt_zmin
#endif
								if(abs(xyzID(3) - zID) .ne. 7 ) then 
									print*,'broky broke',i
									stop
								endif

						end select
			    else if(numnegh .eq. 0) then
				    newBid = -20
			    endif			  
			  else
!////////////////
!////////////////
! calculate new physical position of MC tracer inside block
#ifdef DEBUG_MC
						  print*,'in block'
#endif							

					if(direction .gt. 0 ) then
            select case (dir)
!pos. x
					    case (1)
					      x = x + gridspace(1)
								xpos = x
								x = x - bnd_box(LOW,IAXIS,blockID)

						    xID = floor(x/gridspace(1))+1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(1),'->',xID,'x+',newbid,blockid,i
#endif

								if(abs(xyzID(1)-xid) .ne. 1) then
								  print*,'broky broke',i
									stop
								endif
								ypos = y
								zpos = z
!pos. y
					    case (2)
					      y = y + gridspace(2)
								ypos = y
								y = y - bnd_box(LOW,JAXIS,blockID)
						    yID = floor(y/gridspace(2))+1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(2),'->',yID,'y+',newbid,blockid,i
#endif

								if(abs(xyzID(2)-yid) .ne. 1) then
								  print*,'broky broke',i									
									stop
								endif

								zpos = z
								xpos = x

!pos. z
					    case (3)
	              z = z + gridspace(3)
								zpos = z
								z = z - bnd_box(LOW,KAXIS,blockID)
					      zID = floor(z/gridspace(3))+1

#ifdef DEBUG_MC
								print*,particles(TAG_PART_PROP,i)
								print*,xyzID(3),'->',zID,'z+',newbid,blockid,i
#endif
								if(abs(xyzID(3)-zid) .ne. 1) then
								  print*,'broky broke',i
									stop
								endif

								xpos = x
								ypos = y

						end select
					endif

					if(direction .lt. 0 ) then
            select case (dir)
!neg. x
					  case (1)
					    x   = x - gridspace(1)
							xpos = x
							x   = x - bnd_box(LOW,IAXIS,blockID)
						  xID = floor(x/gridspace(1))+1
#ifdef DEBUG_MC
							print*,particles(TAG_PART_PROP,i)
							print*,xyzID(1),'->',xID,'x-',newbid,blockid,i
#endif

							if(abs(xyzID(1)-xid) .ne. 1) then
							  print*,'broky broke',i
								stop
							endif

							ypos = y
							zpos = z
!neg. y
					  case (2)
					    y   = y - gridspace(2)

							ypos = y
							y   = y - bnd_box(LOW,JAXIS,blockID)
						  yID = floor(y/gridspace(2))+1

#ifdef DEBUG_MC
							print*,particles(TAG_PART_PROP,i)
							print*,xyzID(2),'->',yID,'y-',newbid,blockid,i
#endif
							if(abs(xyzID(2)-yid) .ne. 1) then
								print*,'broky broke',i
								stop
							endif

							xpos = x
							zpos = z
!neg. z
					  case (3)
	            z    = z - gridspace(3)
							zpos = z
							z    = z - bnd_box(LOW,KAXIS,blockID)
					    zID  = floor(z/gridspace(3))+1

#ifdef DEBUG_MC
							print*,particles(TAG_PART_PROP,i)
							print*,xyzID(3),'->',zID,'z-',newbid,blockid,i
#endif
							if(abs(xyzID(3)-zid) .ne. 1) then
							  print*,'broky broke',i
								stop
							endif
							xpos = x
							ypos = y
					  end select
					endif
			  endif
! hopefully the particle move routine moves it to its new location if block is on new domain
! change particle properties
! new IDs 
			   particles(POSX_PART_PROP,i) = xpos
				 particles(POSY_PART_PROP,i) = ypos
				 particles(POSZ_PART_PROP,i) = zpos

! find zone in which the particle resides
! use velocity info for that, so blk + vel gives full location information
  			 particles(VELX_PART_PROP,i) = xID
	  		 particles(VELY_PART_PROP,i) = yID
			   particles(VELZ_PART_PROP,i) = zID

! correct numerical errors
				 particles(POSX_PART_PROP,i) = gridspace(1)*(xid-0.5)+lowbx
				 particles(POSY_PART_PROP,i) = gridspace(2)*(yid-0.5)+lowby
				 particles(POSZ_PART_PROP,i) = gridspace(3)*(zid-0.5)+lowbz

! update carried information
!				 particles(NSN_PART_PROP,i)  = temp
!				 particles(MASS_PART_PROP,i) = dens
! save volume

!				 gridspace =  gr_delta(1:MDIM,lrefine(newBID))

! test only works for 1 thread simulation
#ifdef DEBUG_MC
 				 if(abs(particles(POSX_PART_PROP,i)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,newBID))) .gt. 1e6 .or. & 
 				   abs(particles(POSY_PART_PROP,i)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,newBID))) .gt. 1e6 .or. &
				   abs(particles(POSZ_PART_PROP,i)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,newBID))) .gt. 1e7) then

				   print*,'position wrong'
				   print*,particles(POSX_PART_PROP,i)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,newBID))
				   print*,particles(POSY_PART_PROP,i)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,newBID))
				   print*,particles(POSZ_PART_PROP,i)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,newBID))						
					 print*,particles(TAG_PART_PROP,i)
           print*, xid,yid,zid
				 endif
#endif

#ifdef DEBUG_MC
				if(xID .lt. 1 .or. yID .lt. 1 .or. zID .lt. 1 .or. xID .gt. 8 .or. yID .gt. 8 .or. zID .gt. 8 ) then
					print*,'id out of range',xID,yID,zID
					stop
				endif
#endif

! move to next particle
				exit 
			endif
		enddo ! face loop
  enddo ! particle loop


! deallocate last block
  if(p_count .gt. 0) then
		deallocate(totXFlx)
		deallocate(totYFlx)
 		deallocate(totZFlx)
		deallocate(leftMass)
	  call Grid_releaseBlkPtr(blockID,blockData,CENTER)
	endif

#ifdef DEBUG_SN
	print*,'done MC-ing'
#endif

  ! update thermodynamic variables, could have different mapType 
  ! temperature and density, mapped to tracer particles which overload SB particle properties
  ! attribute mapping is done in pt_setDataStructures
	! NSN  = temp
  ! MASS = dens
  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,therm_NumAttrib,therm_Attrib,mapType)
  
  return
!!------------------------------------------------------------------------------

end subroutine pt_advanceMonteCarlo
