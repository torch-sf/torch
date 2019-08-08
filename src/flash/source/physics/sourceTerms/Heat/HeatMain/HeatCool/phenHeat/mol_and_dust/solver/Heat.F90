!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!  Apply the stat+gauss source term operator to a block of
!!  zones. The energy generation rate is used to update the
!!  internal energy in the zone. The phonomenological heating
!!  rate is described as a 3-D Gauss function.
!!
!!  After we call stat+gauss, call the eos to update the
!!  pressure and temperature based on the phenomenological
!!  heating.
!!  
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine Heat (blockCount,blockList,dt,time)

#include "Flash.h"

#ifdef IHP_SPEC
use rt_data, only : rt_heatInRad
use RadTrans_data, only : rt_useRadTrans
#endif
use Heat_interface, ONLY : RadHeat

!
!==============================================================================
!
  implicit none
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  
#ifdef IHP_SPEC
  if ((.not. rt_heatInRad) .or. (.not. rt_useRadTrans)) then
    call RadHeat(blockCount,blockList,dt,time)
  end if
#else
  call RadHeat(blockCount,blockList,dt,time)
#endif

  return
end subroutine Heat


