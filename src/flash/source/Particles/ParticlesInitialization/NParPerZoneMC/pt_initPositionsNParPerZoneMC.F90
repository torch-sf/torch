!!****if* source/Particles/ParticlesInitialization/Lattice/pt_initPositionsLattice
!!
!! NAME
!!    pt_initPositionsLattice
!!
!! SYNOPSIS
!!
!!    call pt_initPositionsLattice(integer(in)  :: blockID,
!!                                 logical(OUT) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations. 
!! 
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!! PARAMETERS
!!
!!***
subroutine pt_initPositionsNParPerZoneMC (blockID, success)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
														 pt_typeInfo, pt_meshMe, pt_posAttrib, pt_velNumAttrib,pt_velAttrib, &
														 therm_NumAttrib, therm_Attrib
       
  use Particles_data, ONLY : pt_parPerSide,pt_tracerPosX,pt_tracerPosY,pt_tracerPosZ
	use tree, only : bnd_box, lrefine

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_mapMeshToParticles, Grid_getBlkPtr, & 
														 Grid_releaseBlkPtr, Grid_getFluxData
	use Grid_data, ONLY : gr_delta

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  integer :: i, j, k, l, p, startID,x,y,z
	integer :: xSizeCoord, ySizeCoord, zSizeCoord
  integer :: part_props=NPART_PROPS

  real, allocatable, dimension(:,:,:) :: tracerPosX
  real, allocatable, dimension(:,:,:) :: tracerPosY
  real, allocatable, dimension(:,:,:) :: tracerPosZ
	real, pointer, dimension(:,:,:,:)	  :: solnData

	real :: xx, yy, zz
	real :: xpos, ypos, zpos
	real :: rho, temp
  integer :: xID, yID, zID, tag


	real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
	real, allocatable,dimension(:)  :: dx, dy, dz

	integer, dimension(2,MDIM)			:: blkLimits, blkLimitsGC
	real, dimension(MDIM)			 			:: del
	logical													:: getGuardCells = .true.

	real :: blockBounds(LOW:HIGH,1:MDIM)
  real, dimension(NDIM)	:: gridspace

!----------------------------------------------------------------------

  p = pt_numLocal
! debug
	blockBounds = bnd_box(:,:,blockID)

! assume it goes off in highest refinement region
  gridspace = gr_delta(1:MDIM,lrefine(blockID))

	call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
	call Grid_getBlkPtr(blockID,solnData)

	xSizeCoord = blkLimitsGC(HIGH,IAXIS)
  ySizeCoord = blkLimitsGC(HIGH,JAXIS)
	zSizeCoord = blkLimitsGC(HIGH,KAXIS)

	! allocate space for dimensions
	allocate(xCoord(xSizeCoord))
	allocate(yCoord(ySizeCoord))
	allocate(zCoord(zSizeCoord))

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

				startID = pt_numLocal+1
				do l = 0, pt_parPerSide-1

					p = startID + l

          if (p > pt_maxPerProc) then
             print *,' '
             print *,'PARAMETER pt_maxPerProc is set to ',pt_maxPerProc
             print *,'  To avoid this crash, redimension bigger in your flash.par'
             call Driver_abortFlash &
                   ("pt_initPositionsNParPerZone:  Exceeded max # of particles/processor!")
          endif

					rho = solnData(DENS_VAR,i,j,k)
          temp = solnData(TEMP_VAR,i,j,k)

!   same block as source initially
          particles(BLK_PART_PROP,p)  = blockID
          particles(PROC_PART_PROP,p) = real(pt_meshMe)

! unique tag set afterwards
          particles(TAG_PART_PROP,p)  = 0d0

#ifdef TYPE_PART_PROP
          particles(TYPE_PART_PROP,p) = PASSIVE_PART_TYPE
#endif

! multiply position array by zonesize and move to right zone
          particles(POSX_PART_PROP,p) = xx
          particles(POSY_PART_PROP,p) = yy
          particles(POSZ_PART_PROP,p) = zz

! refinement level tracer sits at
				  particles(LREF_PART_PROP,p) = lrefine(blockID)

! move to local coordinate system
				  xpos = xx - bnd_box(LOW,IAXIS,blockID)! + gridspace(1)*0.5d0
				  ypos = yy - bnd_box(LOW,JAXIS,blockID)! + gridspace(2)*0.5d0
				  zpos = zz - bnd_box(LOW,KAXIS,blockID)! + gridspace(3)*0.5d0

! get the IDs from zone position
   				xID = floor(xpos/gridspace(1))+1
		  		yID = floor(ypos/gridspace(2))+1
			  	zID = floor(zpos/gridspace(3))+1

! set from zone indices
          particles(VELX_PART_PROP,p) = xID
          particles(VELY_PART_PROP,p) = yID
          particles(VELZ_PART_PROP,p) = zID

! tag 
          particles(TSN_PART_PROP,p)  = 0d0

! explosion time 
          particles(TCRT_PART_PROP,p) = 0d0

! free for whatever, temperature 
          particles(NSN_PART_PROP,p)  = temp

! nee done more field for density, could be mass if sinks are in etc. check with some preprocessor code maybe
          particles(MASS_PART_PROP,p) = rho

!        	if(abs(particles(POSX_PART_PROP,p)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,blockID))) .gt. 1e6 .or. & 
!		        abs(particles(POSY_PART_PROP,p)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,blockID))) .gt. 1e6 .or. &
!		        abs(particles(POSZ_PART_PROP,p)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,blockID))) .gt. 1e6) then

!	      	  print*,'position wrong in init'
! 		      print*,particles(POSX_PART_PROP,p)-(gridspace(1)*(xid-0.5)+bnd_box(LOW,IAXIS,blockID))
!		        print*,particles(POSY_PART_PROP,p)-(gridspace(2)*(yid-0.5)+bnd_box(LOW,JAXIS,blockID))
!      		  print*,particles(POSZ_PART_PROP,p)-(gridspace(3)*(zid-0.5)+bnd_box(LOW,KAXIS,blockID))						
!			      print*,particles(TAG_PART_PROP,p)

!       		print*,(particles(POSX_PART_PROP,p)-bnd_box(LOW,IAXIS,blockID))/gridspace(1),xid
!		        print*,(particles(POSY_PART_PROP,p)-bnd_box(LOW,JAXIS,blockID))/gridspace(2),yid
!		        print*,(particles(POSZ_PART_PROP,p)-bnd_box(LOW,KAXIS,blockID))/gridspace(3),zid
!			      print*,gridspace
!           print*,zCoord(2)-zCoord(1)
!           print*,zCoord(4)-zCoord(3)
!           print*,zCoord(6)-zCoord(5)
!		        stop
!		      endif

				enddo

				pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) + pt_parPerSide

				! is just last occupied slot
				pt_numLocal = pt_numLocal + pt_parPerSide

			enddo
		enddo
	enddo
    
  ! Set the particle database local number of particles.
  success=.true.  

  deallocate(xCoord)
	deallocate(yCoord)
	deallocate(zCoord)

!  clean up memory 
  call Grid_releaseBlkPtr(blockID,solnData)

	return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsNParPerZoneMC
