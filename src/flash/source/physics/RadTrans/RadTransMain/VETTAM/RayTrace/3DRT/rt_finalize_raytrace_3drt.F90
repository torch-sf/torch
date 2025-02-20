!!****if* source/physics/RadTrans/RadTransMain/HybridChar3DRTFLD/Radtrans_finalize
!!
!! NAME
!!
!!  Radtrans_finalize
!!
!! SYNOPSIS
!!
!!  call Radtrans_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the Radtrans unit.
!!
!! ARGUMENTS
!!
!!***
subroutine rt_finalize_raytrace_3drt
  use RadTrans_RayTrace_3DRT, ONLY : rad_deallocate
  implicit none
#include "Flash.h"
#include "constants.h"
  call rad_deallocate()

  return
end subroutine rt_finalize_raytrace_3drt
