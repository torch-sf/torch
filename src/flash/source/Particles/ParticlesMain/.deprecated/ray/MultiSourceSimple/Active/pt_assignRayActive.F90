!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   loops over all source particles and generates initial rays
!!   
!!
!! Input: 
!!    dt: current simulation timestep
subroutine pt_assignRayActive(i, sourcedata)

  use Particles_data,ONLY : particles, pt_typeInfo
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

! index into particle array 
  integer, intent(inout)    :: i
  integer ::  sourceStr, sourceNum

! too big but not important 
  real, dimension(NPART_PROPS), intent(out) :: sourcedata
! Only compile with this stuff if active particles are present. -JW 
#ifdef ACTIVE_PART_TYPE
  sourceStr = pt_typeInfo(PART_TYPE_BEGIN, ACTIVE_PART_TYPE) + i
  sourceNum = pt_typeInfo(PART_TYPE_BEGIN, ACTIVE_PART_TYPE) + pt_typeInfo(PART_LOCAL, ACTIVE_PART_TYPE) - 1

! i trails by one
  if(i .lt. sourceNum) then
    sourcedata(BLK_PART_PROP)   = particles(BLK_PART_PROP,    sourceStr)
    sourcedata(TAG_PART_PROP)   = particles(TAG_PART_PROP,    sourceStr)
    sourcedata(PROC_PART_PROP)  = particles(PROC_PART_PROP,   sourceStr)

    sourcedata(NION_PART_PROP)  = particles(NION_PART_PROP,   sourceStr)
    sourcedata(EION_PART_PROP)  = particles(EION_PART_PROP,   sourceStr)

    sourcedata(SIGH_PART_PROP)  = particles(SIGH_PART_PROP,   sourceStr)

    sourcedata(POSX_PART_PROP)  = particles(POSX_PART_PROP,   sourceStr)
    sourcedata(POSY_PART_PROP)  = particles(POSY_PART_PROP,   sourceStr)
    sourcedata(POSZ_PART_PROP)  = particles(POSZ_PART_PROP,   sourceStr)
! next source
    i = i + 1
  else
! no source on local processor
    i = -1
  endif
#endif
end subroutine pt_assignRayActive
