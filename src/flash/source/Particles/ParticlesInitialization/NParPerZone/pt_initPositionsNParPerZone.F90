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
subroutine pt_initPositionsNParPerZone (blockID, success)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
														 pt_typeInfo, pt_meshMe, pt_posAttrib, pt_velNumAttrib,pt_velAttrib, &
														 therm_NumAttrib, therm_Attrib
       
  use Particles_data, ONLY : pt_parPerSide,pt_tracerPosX,pt_tracerPosY,pt_tracerPosZ
	use tree, only : bnd_box

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_mapMeshToParticles, Grid_getDeltas

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  integer :: tracerDimN, i, j, k, l, p, startID,x,y,z
	integer :: xSizeCoord, ySizeCoord, zSizeCoord
  integer :: part_props=NPART_PROPS

  real, allocatable, dimension(:,:,:) :: tracerPosX
  real, allocatable, dimension(:,:,:) :: tracerPosY
  real, allocatable, dimension(:,:,:) :: tracerPosZ

	real :: space, xx, yy, zz

	real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
	real, allocatable,dimension(:)  :: dx, dy, dz

	integer, dimension(2,MDIM)			:: blkLimits, blkLimitsGC
	real, dimension(MDIM)			 			:: del
	logical													:: getGuardCells = .true.

	real :: blockBounds(LOW:HIGH,1:MDIM)


!----------------------------------------------------------------------

  p = pt_numLocal
	tracerDimN = pt_parPerSide**3 

! debug
	blockBounds = bnd_box(:,:,blockID)

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

	dx(:) = del(IAXIS)
	dy(:) = del(JAXIS)
	dz(:) = del(KAXIS)

	call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
	call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
	call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)
	! loop over all zones in block
	! for 2d k is just 1
	do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
		zz = blockBounds(LOW,KAXIS) + dz(k)*(k-blkLimits(LOW,KAXIS))
		do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
		  yy = blockBounds(LOW,JAXIS) + dy(j)*(j-blkLimits(LOW,JAXIS))
			do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
			  xx = blockBounds(LOW,IAXIS) + dx(i)*(i-blkLimits(LOW,IAXIS))

				startID = pt_numLocal+1
				do l = 0, tracerDimN-1

					p = startID + l

          if (p > pt_maxPerProc) then
             print *,' '
             print *,'PARAMETER pt_maxPerProc is set to ',pt_maxPerProc
             print *,'  To avoid this crash, redimension bigger in your flash.par'
             call Driver_abortFlash &
                   ("pt_initPositionsNParPerZone:  Exceeded max # of particles/processor!")
          endif

! same block as source initially
					particles(BLK_PART_PROP,p)  = real(blockID)
					particles(PROC_PART_PROP,p) = real(pt_meshMe)
! SB or sink tag, otherwise 0
					particles(TAG_PART_PROP,p)  = l
					particles(TYPE_PART_PROP,p) = PASSIVE_PART_TYPE

! quick hack to get from linear array to 3d array indices
					x = mod(l, pt_parPerSide) + 1
					y = mod((l - mod(l,pt_parPerSide))/pt_parPerSide, pt_parPerSide) + 1
					z = mod((l - mod(l, pt_parPerSide*pt_parPerSide))/(pt_parPerSide*pt_parPerSide),&
								 pt_parPerSide*pt_parPerSide) + 1

! multiply position array by zonesize and move to right zone
					particles(POSX_PART_PROP,p) = pt_tracerPosX(x, y, z)*dx(i) + xx
					particles(POSY_PART_PROP,p) = pt_tracerPosY(x, y, z)*dy(j) + yy
					particles(POSZ_PART_PROP,p) = pt_tracerPosZ(x, y, z)*dz(k) + zz

! test if particle is inside
!					if(particles(POSX_PART_PROP,p) - blockBounds(LOW,1) .lt. 0d0 .or. &
!						 blockBounds(HIGH,1) - particles(POSX_PART_PROP,p) .lt. 0d0) then
!
!						print*,blockBounds(LOW:HIGH,1)
!						print*,particles(POSX_PART_PROP,p)
!						print*,pt_tracerPosX(x, y, z),pt_tracerPosY(x, y, z),pt_tracerPosZ(x, y, z)
!						print*,dx(i),dy(j),dz(k)
!						print*,xx,yy,zz
!						stop
!						print*,'not inside block!'
!					endif

!					if(particles(POSY_PART_PROP,p) - blockBounds(LOW,2) .lt. 0d0 .or. &
!						 blockBounds(HIGH,2) - particles(POSY_PART_PROP,p) .lt. 0d0) then
!
!						print*,blockBounds(LOW:HIGH,1)
!						print*,particles(POSY_PART_PROP,p)
!						print*,'not inside block!'
!					stop
!					endif

!					if(particles(POSZ_PART_PROP,p) - blockBounds(LOW,3) .lt. 0d0 .or. &
!						 blockBounds(HIGH,3) - particles(POSZ_PART_PROP,p) .lt. 0d0) then
!
!						print*,blockBounds(LOW:HIGH,1)
!						print*,particles(POSZ_PART_PROP,p)
!
!						print*,'not inside block!'
!						stop
!					endif

! set from zone 
					particles(VELX_PART_PROP,p) = 0d0
					particles(VELY_PART_PROP,p) = 0d0
					particles(VELZ_PART_PROP,p) = 0d0
! type of origin SN
					particles(TSN_PART_PROP,p)  = 0d0
! explosion time 
					particles(TCRT_PART_PROP,p) = 0d0
! free for whatever, temperature 
					particles(NSN_PART_PROP,p)  = 0d0
! need one more field for density, could be mass if sinks are in etc. check with some preprocessor code maybe
					particles(MASS_PART_PROP,p) = 0d0
				enddo

				pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) + tracerDimN

				! is just last occupied slot
				pt_numLocal = pt_numLocal + tracerDimN

			enddo
		enddo
	enddo

  ! map velocities to the particles
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, pt_numLocal,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,pt_typeInfo(PART_MAPMETHOD,PASSIVE_PART_TYPE))

  ! update thermodynamic variables, could have different mapType 
  ! temperature and density, mapped to tracer particles which overload SB particle properties

  call Grid_mapMeshToParticles(particles, &
       part_props, BLK_PART_PROP,pt_numLocal, &
	     pt_posAttrib,therm_NumAttrib,therm_Attrib,pt_typeInfo(PART_MAPMETHOD,PASSIVE_PART_TYPE))
    
  ! Set the particle database local number of particles.
  success=.true.  

  deallocate(xCoord)
	deallocate(yCoord)
	deallocate(zCoord)

	deallocate(dx)
	deallocate(dy)
	deallocate(dz)

	return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsNParPerZone
