!!****if* source/Particles/ParticlesMain/passive/RungeKutta/pt_advanceRK
!!
!! NAME
!!
!!  pt_advanceRK
!!
!! SYNOPSIS
!!
!!  call pt_advanceRKstatic(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count,
!!                         integer(in):: ind)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  CB: does not call Grid_mapMeshToParticles
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
!! EPHEMERAL AND FLEETING
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

subroutine pt_advanceRKstatic (dtOld,dtNew,particles,p_count, ind)
    
  implicit none

#include "Flash.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

!!------------------------------------------------------------------------------
  
end subroutine pt_advanceRKstatic


