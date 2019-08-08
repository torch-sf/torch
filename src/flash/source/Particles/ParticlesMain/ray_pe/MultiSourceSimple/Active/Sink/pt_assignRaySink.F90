!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014

!! Modified by J. Wall to read in from sink particles, 2015.
!! Note, handling of the Nion and Eion will be done in a separate
!! code that calculates stellar evolution. - JW

!! Description:
!!   loops over all source particles and generates initial rays
!!   
!!
!! Input: 
!!    dt: current simulation timestep
subroutine pt_assignRaySink(i, sourcedata)
! Include Flash.h first so you can access the #defines in there. - JW
#include "Flash.h"

#ifdef SINK_PART_TYPE
  use Particles_sinkdata, ONLY : particles_local, localnp
  use pt_sinkInterface, ONLY : pt_sinkGatherGlobal
#endif
  use Particles_data,     ONLY : pt_typeInfo
  use Particles_interface, ONLY : Particles_sinkMoveParticles, &
      Particles_sinkSortParticles
  implicit none

#include "constants.h"
#include "Particles.h"

! index into particle array 
  integer, intent(inout)    :: i
  integer ::  sourceStr, sourceNum

! too big but not important 
  real, dimension(NPART_PROPS), intent(out) :: sourcedata
! Only compile with this stuff if sinks are present. - JW  

#ifdef SINK_PART_TYPE
#ifndef ACTIVE_PART_TYPE
! This isn't necessary, we already know how many particles are on the local
! processor. - JW
!  sourceStr = pt_typeInfo(PART_TYPE_BEGIN, SINK_PART_TYPE) + i
!  sourceNum = pt_typeInfo(PART_TYPE_BEGIN, SINK_PART_TYPE) + pt_typeInfo(PART_LOCAL, SINK_PART_TYPE) - 1

! i trails by one
!  if(i .lt. sourceNum) then
! Indeed, the way its set up i starts at 0 and goes to localnp-1, so
! we want to bump it by one *before* we copy info here... - JW
  if (i .lt. localnp) then
    i = i + 1
    sourcedata(BLK_PART_PROP)   = particles_local(BLK_PART_PROP,  i)
    sourcedata(TAG_PART_PROP)   = particles_local(TAG_PART_PROP,  i)
    sourcedata(PROC_PART_PROP)  = particles_local(PROC_PART_PROP, i)

    sourcedata(NION_PART_PROP)  = particles_local(NION_PART_PROP, i)
    sourcedata(EION_PART_PROP)  = particles_local(EION_PART_PROP, i)

    sourcedata(SIGH_PART_PROP)  = particles_local(SIGH_PART_PROP, i)

    sourcedata(POSX_PART_PROP)  = particles_local(POSX_PART_PROP, i)
    sourcedata(POSY_PART_PROP)  = particles_local(POSY_PART_PROP, i)
    sourcedata(POSZ_PART_PROP)  = particles_local(POSZ_PART_PROP, i)
    
#ifdef PE_HEAT
! # of PE photons.
    sourcedata(NPEP_PART_PROP) = particles_local(NPEP_PART_PROP, i)
! average energy of PE photons.
    sourcedata(EPEP_PART_PROP) = particles_local(EPEP_PART_PROP, i)
! cross section of dust to PE photons.
    sourcedata(SPEP_PART_PROP) = particles_local(SPEP_PART_PROP, i)
#endif
! But not here - JW
! next source
!    i = i + 1
  else
! no source on local processor
    i = -1
  endif
#endif
#endif
return
end subroutine pt_assignRaySink
