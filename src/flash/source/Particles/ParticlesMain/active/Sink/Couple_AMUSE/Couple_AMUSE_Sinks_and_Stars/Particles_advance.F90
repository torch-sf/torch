!!****if* source/Particles/ParticlesMain/Particles_advance
!!
!! NAME
!!
!!  Particles_advance
!!
!! SYNOPSIS
!!
!!  Particles_advance(real(in) :: dtOld,
!!                    real(in) :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  Calls passive and active versions
!!  
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles.
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!***

!#define DEBUG_VPARTICLES

!===============================================================================

subroutine Particles_advance (dtOld,dtNew)
  
#ifdef DEBUG_PARTICLES
  use Driver_data, ONLY: dr_globalMe
#endif
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, useParticles, & 
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_meshMe, pt_typeInfo,&
       pt_indexList, pt_indexCount, pt_keepLostParticles, pt_numLost, &
       pt_reduceGcellFills

  use pt_interface, ONLY: pt_advancePassive, pt_advanceActive,  pt_updateTypeDS
  use Grid_interface, ONLY : Grid_moveParticles, Grid_fillGuardCells, &
                             Grid_mapMeshToParticles, Grid_sortParticles
  use Particles_interface, ONLY: Particles_sinkMoveParticles, Particles_sinkSortParticles, &
                                 Particles_sinkAdvanceParticles, &
                                 Particles_sinkCreateAccrete, Particles_wind
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
#include "GridParticles.h"
  
  real, INTENT(in)  :: dtOld, dtNew

  integer       :: i,nstep,kk
  !integer       :: totalPassive, totalActive
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  !real          :: jumpx,jumpy,jumpz
  integer       ::p_begin,p_end, p_count
  !logical       :: allTypesNotDone
  logical,parameter :: regrid=.false.
  logical,save      :: gcMaskLogged = .FALSE.
  integer       :: pfor,pbak, lostNow
  real          :: dm, starMass
!!------------------------------------------------------------------------------
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return

#ifdef DEBUG_VPARTICLES
  do kk=1,pt_numLocal
     write(*,*)'local particle ',kk,'In proc ',int(particles(PROC_PART_PROP,kk)),'on blk=',int(particles(BLK_PART_PROP,kk))
     write(*,*)'Xpos= ',particles(POSX_PART_PROP,kk),'Ypos=',particles(POSY_PART_PROP,kk)
  enddo
#endif 

  ! Prepare guardcell data needed for particle interpolation.
  !
  ! Experimentation with passive particles (with the old way of advancing particles)
  ! has shown that at least 2 layers of guardcells need to be filled
  ! with updated data for vel[xyz] and density, in order to get the
  ! same results as for a full guardcell fill, when using native grid interpolation. - KW
  ! With "monotonic" interpolation, even more layers are needed. - KW

#ifndef PRNT_PART_PROP
  if (pt_reduceGcellFills) then
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,unitReadsMeshDataOnly=.true.)
  else
     call Grid_fillGuardCells( CENTER, ALLDIR,&
          maskSize=pt_gcMaskSizeForAdvance,mask=pt_gcMaskForAdvance,&
          doLogMask=.NOT.gcMaskLogged)
  end if
#endif

  ! Sort particles there so we only have to move the minimum of them.
  !  Sort by type, and then within each sort by block

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
  

  !! Now do actual movement, advance particles in time
  !! Here we assume that all passive particles have identical
  !! integration. The active particles may also chose to have
  !! identical integration, but they also have to option of picking
  !! different integration methods for different types.

!  do i = 1, NPART_TYPES
!     p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
!     p_count=pt_typeInfo(PART_LOCAL,i)
!     p_end=p_begin+p_count-1
!     select case(pt_typeInfo(PART_ADVMETHOD,i))
!     case(RUNGEKUTTASTATIC) 
!        call pt_advanceRKstatic(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(RUNGEKUTTATHERMO) 
!        call pt_advanceRKthermo(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(RUNGEKUTTA) 
!        call pt_advanceRK(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(MIDPOINT)
!        call pt_advanceMidpoint(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(ESTI)
!        call pt_advanceEsti(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(EULER_TRA)
!        call pt_advanceEuler_passive(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(EULER_MAS)
!        call pt_advanceEuler_active(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(LEAPFROG)
!        call pt_advanceLeapfrog(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(LEAPFROG_COSMO)
!        call pt_advanceLeapfrog_cosmo(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(CHARGED)
!        call pt_advanceCharged(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count)
!     case(DPD)
!        call pt_advanceDPD(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!     case(CUSTOM)
!        call pt_advanceCustom(dtOld,dtNew,particles(:,p_begin:p_end),&
!             p_count, i)
!	 case(MONTECARLO)
!	    call pt_advanceMonteCarlo(dtOld,dtNew,particles(:,p_begin:p_end),&
!	         p_count, i)
!     end select

!  end do

  ! sink particle routines

#ifdef DEBUG_PARTICLES
  if (dr_globalMe .eq. MASTER_PE) print *, 'Particles_advance: entering sink routines.'
#endif
  ! advance sink particles based on velocity and acceleration
  call Particles_sinkAdvanceParticles(dtNew)

#ifdef AMUSE_STARS
  call stochastic(dtNew)
#else
  ! check for creation and accretion
  call Particles_sinkCreateAccrete(dtNew)
#endif

  ! move sink particles to the right blocks/procs
  call Particles_sinkMoveParticles(regrid)

  ! sort sink particles
  call Particles_sinkSortParticles()

#ifdef DEBUG_PARTICLES
  if (dr_globalMe .eq. MASTER_PE) print *, 'Particles_advance: done with sink routines.'
#endif


#ifdef DEBUG_PARTICLES
  print*,' ready to move Particles'
#endif


#ifdef DEBUG_VPARTICLES
  do kk=1,pt_numLocal
     write(*,*)'local particle ',kk,'In proc ',int(particles(PROC_PART_PROP,kk)),'on blk=',int(particles(BLK_PART_PROP,kk))
     write(*,*)'Xpos= ',particles(POSX_PART_PROP,kk),'Ypos=',particles(POSY_PART_PROP,kk)
  enddo
#endif 
  
  ! Put the particles in the appropriate blocks if they've moved off
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       pt_indexList, pt_indexCount, regrid) 


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

  ! inject stellar winds
  !call Particles_sinkWind(dtNew)
  call Particles_wind(dtNew)

#ifdef DEBUG_PARTICLES
  print*,' back from Grid_moveParticles'
#endif
  ! If predictive routines are used, they will need to sort and prepare for the
  !  next time step.  Since sorting is so expensive, we suffer code duplication
  !  and do it in the pt_preparePassive routines.
  ! Many algorithms use the stub routines.

  do i=1,NPART_TYPES
     if (pt_typeInfo(PART_ADVMETHOD,i)==ESTI)then
        p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
        p_count=pt_typeInfo(PART_LOCAL,i)
        p_end=p_begin+p_count-1
        call pt_prepareEsti(dtOld,dtNew,particles(:,p_begin:p_end),p_count,i)
     end if
  end do

  gcMaskLogged = .TRUE.
  
! just does star particle creation, this means they are created at the end 
! of the advancement, they are integrated with one of the standard advancement
! schemes, e.g. Leapfrog in the case block above in this routine
!  call Particles_starsAdvance(dtNew)

  return

  !!-----------------------------------------------------------------------
end subroutine Particles_advance


