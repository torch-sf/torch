!!****if* source/Particles/ParticlesMain/active/Sink/Couple_AMUSE/Particles_moveAndSort
!!
!! NAME
!!
!!  Particles_moveAndSort
!!  a modification of Particles_advance that just sorts and updates stuff for
!!  sinks
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

!! All I want to do is update particle proc locations and sort them on the processors.
!! Actual integration is done by AMUSE. - JW

!===============================================================================

subroutine Particles_moveAndSort(regrid)

#include "Flash.h"

  implicit none

  logical, intent(in) :: regrid

  return

end subroutine Particles_moveAndSort


