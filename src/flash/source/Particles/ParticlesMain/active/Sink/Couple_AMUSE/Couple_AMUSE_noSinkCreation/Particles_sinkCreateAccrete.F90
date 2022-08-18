!!****f* source/Particles/Particles_sinkCreateAccrete
!!
!! NAME
!!
!!  Particles_sinkCreateAccrete
!!
!! SYNOPSIS
!!
!!  call Particles_sinkCreateAccrete(real, intent(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Handles creation of sink particles and accretion of gas onto sink particles.
!!  It also calls routines for sink particle merging and for dumping sink
!!  particle properties every time step (pt_sinkDumpParticles).
!!
!! ARGUMENTS
!!
!!   dt - the current simulation time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2014
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   renamed and cleaned by Christoph Federrath, 2013
!!
!!***

!! I'm a stub, I don't want to allow sink creation or accretion in this
!! module. If you want that, you need to use sink/Couple_AMUSE instead. - JW

subroutine Particles_sinkCreateAccrete(dt)
  implicit none
  real, intent(IN) :: dt
end subroutine Particles_sinkCreateAccrete
