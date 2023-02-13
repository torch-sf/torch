!!****if* source/Particles/ParticlesMain/Particles_updateRefinement
!!
!! NAME
!!
!!  Particles_updateRefinement
!!
!! SYNOPSIS
!!
!!  Particles_updateRefinement(real(inout) :: oldLocalNumBlocks)
!!
!! DESCRIPTION
!!   This routine provides a hook into the particle data structure
!!   for the Grid. It is called during Grid_updateRefinement processing
!!   by the Grid. The routine passes the control right back to
!!   Grid, with Particles-specific data structures in the argument
!!   list, so that Grid can operate on them.
!!
!! ARGUMENTS
!!
!!    oldLocalNumBlocks :   number of blocks on a processor before 
!!                          refinement. 
!!
!! PARAMETERS
!!  
!!
!!***
#include "constants.h"
#include "Particles.h"
#include "Flash.h"

subroutine Particles_updateRefinement(oldLocalNumBlocks)

  use Particles_data, ONLY : particles, pt_numLocal, pt_maxPerProc,useParticles, &
       pt_posInitialized, pt_logLevel, pt_meshMe,&
       pt_indexList, pt_indexCount, pt_keepLostParticles, pt_numLost
  use Grid_interface,ONLY : Grid_moveParticles, Grid_sortParticles
  use Logfile_interface,ONLY : Logfile_stamp
  use Particles_interface, only: Particles_sinkMoveParticles
  use pt_interface, ONLY: pt_updateTypeDS
  
  implicit none 
  integer,intent(INOUT) :: oldLocalNumBlocks
  logical, parameter :: regrid = .true.

! JCI - added this line present in a previous release.
  integer       :: pfor,pbak, lostNow,i 
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  if(.not.useParticles)return
  if(.not.pt_posInitialized) then
     if (pt_logLevel > PT_LOGLEVEL_WARN_USE) then
        if (pt_meshMe==MASTER_PE) then
           print*,'WARNING: Particles_updateRefinement was called while particles positions are not yet initialized!'
        end if
        call Logfile_stamp( &
             'WARNING: Called while particles positions are not yet initialized!','Particles_updateRefinement')
     end if
     return
  end if
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc, pt_numLocal, &
       pt_indexList, pt_indexCount,&
       regrid)
  call Particles_sinkMoveParticles(regrid)
    
! added to allow particle advancement outside Particle advancement call, i.e. radiation transport
#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
  
  if(pt_keepLostParticles) then
     pfor=pt_numLocal
     do while(particles(BLK_PART_PROP,pt_numLocal)==LOST)
        pt_numLocal=pt_numLocal-1
     end do
     lostNow=pfor-pt_numLocal
     pt_numLost=pt_numLost+lostNow
     pbak=pt_maxPerProc-pt_numLost
     if(pbak<pt_numLocal)call Driver_abortFlash("no more space for lost particles")
     do i = 1,lostNow
        particles(:,pbak+i)=particles(:,pt_numLocal+i)
     end do
  end if

  ! Now update the pt_typeInfo data structure
  call pt_updateTypeDS(particlesPerBlk)

  return
end subroutine Particles_updateRefinement
