!!****if* source/Particles/ParticlesInitialization/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!   call Particles_initPositions( logical(inout) :: partPosInitialized,
!!                                 logical(out)   :: updateRefine)
!!
!!
!! DESCRIPTION
!!
!!    Initialize particle locations. This routine calls pt_initPositions
!!    which a Particles unit's local routine to initialize the positions
!!    on leaf blocks. The routine also creates tags for all the particles
!!    This routine will initialize based on Lattice or with Density 
!!    distribution depending upon which of the two is selected. 
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized.
!!            On entry, a value of TRUE is taken to mean that position
!!            initialization has already been completed previously,
!!            so the routine returns immediately (leaving partPosInitialized
!!            TRUE).
!!            If particles are disabled (as per runtime parameter useParticles),
!!            this implementation also returns immediately with
!!            partPosInitialized set to TRUE.
!!            Otherwise, partPosInitialized will be to TRUE if all particles
!!            have been placed in the domain successfully. A return value
!!            of FALSE may indicate that only some particles have been placed
!!            in the domain, perhaps because of space limitations in the
!!            particles array in some MPI tasks; partially initialized
!!            particles data of this kind may still be useful during FLASH
!!            initialization, in particular for the the purpose of providing
!!            refinement criteria for the initial Grid construction if the
!!            runtime parameter refine_on_particle_count is TRUE, but a
!!            fully successful Particles_initPositions invocation is still
!!            required before the simulation is allowed to proceed with
!!            its main evolution loop.
!!
!!  updateRefine : is set to TRUE if the routine wishes to indicate that during
!!                 the initial iterative construction of an AMR Grid (see
!!                 Grid_initDomain and gr_expandDomain), the initialization of
!!                 particle positions need not be repeated for each iteration if
!!                 all particles have already been placed in the domain in a
!!                 previous iteration.  Under that condition, subsequent calls
!!                 to Particles_initPositions from the Grid construction loop
!!                 will have partPosInitialized=.TRUE.  so will return
!!                 immediately, and Particles_updateRefinement will be called in
!!                 each iteration when the number of blocks or the block
!!                 distribution may have changed, to make sure that the retained
!!                 particles get moved to the correct block (hence the name of
!!                 the dummy argument).
!!
!!                 This implementation always returns FALSE.
!!                 Alternative implementations may wish to return TRUE
!!                 instead if initialization is very expensive.
!!                 
!!
!! NOTES
!!
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates particles data that is private to the Particles unit.
!!
!!  May modify an internal flag (pt_posInitialized) that keeps track of
!!  whether initialization of particle positions is complete.
!!
!! SEE ALSO
!!
!!  Driver_initFlash
!!  Grid_initDomain
!!  Particles_initData
!!***

#define VORAMR !Must be commented out to turn off VorAMR - SCL
!!#define DEBUG_PARTICLES

subroutine Particles_initPositions (partPosInitialized,updateRefine)


  use Grid_interface, ONLY : Grid_getListOfBlocks
  use Driver_interface, ONLY : Driver_abortFlash
  use pt_interface, ONLY : pt_initPositions,pt_createTag,&
       pt_initVoronoiPositions
  use Particles_data, ONLY : pt_posInitialized,pt_numLocal,useParticles,&
       pt_typeInfo, particles, pt_meshNumProcs, pt_meshMe

  use pt_interface, ONLY :  pt_initLocal

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine

  integer       :: i, j, k, b
  integer       :: p
  integer       :: numLocalThisType, numNewLocalThisType, numLocalPreviousTypes
  integer       :: numPreviousLocal
  integer       :: blockID
  logical       :: IsInBlock
  real          :: xpos, ypos, zpos, bxl, byl, bzl, bxu, byu, bzu
  real          :: xvel, yvel, zvel

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  real, dimension(2,MDIM):: boundBox
  integer :: blkCount
  integer,dimension(MAXBLOCKS) :: blkList
!----------------------------------------------------------------------
#ifdef VORAMR
  if (pt_meshMe == MASTER_PE) then
     print*,"Calling VorAMR routine."
  end if
  ! VorAMR subroutine, VORAMR compile flag is set in pt_initVoronoiPositions.F90
  call pt_initVoronoiPositions(partPosInitialized,updateRefine)
#endif
#ifndef VORAMR
  if (pt_meshMe == MASTER_PE) then
     print*,"Skipping VorAMR."
  end if
#endif
  
  if(.not.useParticles) then
     partPosInitialized = .true.
  end if
  if(partPosInitialized) return

  !CD: We need to move the call to pt_initLocal to this level. 
  !Otherwise, if it is in pt_initPositions, we get a deadlock
  !when the number of blocks are not the same on all processors.
  call pt_initLocal()

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  ! Distribute initial positions

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
#ifdef DEBUG_PARTICLES
  if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
     print*,pt_meshMe,': am initializing particles on number of blocks=',blkCount
  end if
#endif

  updateRefine = .FALSE.

  partPosInitialized=.true.

  numLocalPreviousTypes = 0
  if(.not.updateRefine) then
     pt_numLocal = 0
  else
     pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
  end if
  do i = 1,NPART_TYPES
     if(.not.updateRefine) pt_typeInfo(PART_LOCAL,i) = 0
     numLocalThisType = pt_typeInfo(PART_LOCAL,i)
     
     b=0
     numNewLocalThisType = 0
     numPreviousLocal = pt_numLocal
     do while((b<blkCount).and.partPosInitialized)
        b=b+1
        blockID = blkList(b)
!!        print*,'pt_initPositions',blockID, pt_typeInfo(PART_INITMETHOD,i)
        select case(pt_typeInfo(PART_INITMETHOD,i))
        case(LATTICE)
           call pt_initPositionsLattice(blockID,partPosInitialized)
        case(WITH_DENSITY, CELLMASS, REJECTION)
           call pt_initPositionsWithDensity(blockID,partPosInitialized)
        case(NPARPERZONE)
           call pt_initPositionsNParPerZone(blockID,partPosInitialized)
        case(NPARPERZONEMC)
           call pt_initPositionsNParPerZoneMC(blockID,partPosInitialized)
        case(CUSTOM)
           call pt_initPositions(blockID,partPosInitialized)
        case default
           call Driver_abortFlash("Particles_initPosition: no valid initialization method")
        end select
        numNewLocalThisType = pt_numLocal - numPreviousLocal
        pt_typeInfo(PART_LOCAL,i) = numNewLocalThisType + numLocalThisType
     enddo
#ifdef TYPE_PART_PROP
     particles(TYPE_PART_PROP, &
          pt_numLocal-numNewLocalThisType+1:pt_numLocal) = pt_typeInfo(PART_TYPE,i) 
#endif
     numLocalThisType=pt_typeInfo(PART_LOCAL,i)
     numLocalPreviousTypes = numLocalPreviousTypes + numLocalThisType
  end do
  pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
!!  print*,'Particles_initPositions: pt_numLocal now is',pt_numLocal

  pt_posInitialized = partPosInitialized

#ifdef DEBUG_PARTICLES
  if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
     print*,'Particles_initPositions on processor', pt_meshMe, 'done, pt_numLocal=',pt_numLocal
  end if
#endif

  call pt_createTag()
  return
  
end subroutine Particles_initPositions
