!!****if* source/Particles/ParticlesMain/passive/RungeKutta/pt_advanceRK
!!
!! NAME
!!
!!  pt_advanceRK
!!
!! SYNOPSIS
!!
!!  call pt_advanceRK(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version implements the improved Euler method, also called Heun's method,
!!  for integration in time. The improved Euler methos is one of the family of
!!  2-stage, second-order Runge-Kutta methods.  (It can probably also be framed
!!  as a simple Predictor-Corrector method, with a first-order predictor and
!!  corrector.)
!!
!!  In detail:
!!      x*(t+dtNew) = x(t)  + dtNew *        v(x(t),t)
!!      x(t+dtNew)  = x(t)  + dtNew * 1/2* [ v(x(t),t) + v(x*(t+dtNew),t+dtNew) ]
!!      v(t+dtNew)  =  v(x(t+dtNew),t+dtNew)
!!  where x* is ephemeral (a predicted position that then is corrected after
!!  another evaluation there).
!!
!!  Implementation detail:
!!  This can be rewritten to save on memory for intermediate result storage:
!!      x*(t+dtNew) = x (t)        + dtNew *          v(x(t),t)
!!      x(t+dtNew)  = x*(t+dtNew)  + dtNew * 1/2* [ - v(x(t),t) + v(x*(t+dtNew),t+dtNew) ]
!!      v(t+dtNew)  =  v(x(t+dtNew),t+dtNew)
!!  where x* can be stored in the same location as the previous and final x.
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind -- index into pt_typeInfo
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
!!
!!***

!===============================================================================

subroutine pt_advanceRKthermo (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: pt_numLocal, pt_maxPerProc, &
       useParticles, pt_typeInfo, &
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_meshMe, &
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib, therm_NumAttrib, therm_Attrib

  use Grid_interface, ONLY : Grid_mapMeshToParticles
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

  integer       :: i,particleTypes
  
  integer,dimension(MAXBLOCKS, 1) :: perBlk
  real          :: jumpx,jumpy,jumpz
  real,allocatable :: origVel(:,:)
  integer :: part_props=NPART_PROPS

  integer :: mapType 

!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return

  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

  ! Update the particle positions to temporary ("predicted") values
  do i = 1, p_count
 
     jumpx = dtNew * particles(VELX_PART_PROP,i)
     particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + jumpx

     if(NDIM >1) then
        jumpy = dtNew * particles(VELY_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + jumpy
        if(NDIM >2) then
           jumpz = dtNew * particles(VELZ_PART_PROP,i)
           particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + jumpz
        end if
     end if
  enddo

  ! Now save the original velocity values
  allocate(origVel(p_count,MDIM))
  origVel(:,1) = particles(VELX_PART_PROP,1:p_count)
  origVel(:,2) = particles(VELY_PART_PROP,1:p_count)
  origVel(:,3) = particles(VELZ_PART_PROP,1:p_count)

  ! Map the updated gas velocity field at the temporary positions to
  ! obtain a second estimate of velocities;

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP, p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! Adjust particle positions, using the second point velocities
  do i = 1, p_count
     particles(POSX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELX_PART_PROP,i) - origVel(i,1))
     if(NDIM>1)&
          particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELY_PART_PROP,i) - origVel(i,2))
     if(NDIM>2)&
          particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
          dtNew * 0.5*(particles(VELZ_PART_PROP,i) - origVel(i,3))
  enddo

  ! done with this temporary storage
  deallocate(origVel)

  ! Map the updated gas velocity field onto the current particle positions to
  ! obtain the updated particle velocities - for the next integration step
  ! as well as for particle plot files etc.

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  ! update thermodynamic variables, could have different mapType 
  ! temperature and density, mapped to tracer particles which overload SB particle properties

  call Grid_mapMeshToParticles(particles,&
       part_props, BLK_PART_PROP,p_count,&
       pt_posAttrib,therm_NumAttrib,therm_Attrib,mapType)

  return
!!------------------------------------------------------------------------------
  
end subroutine pt_advanceRKthermo


