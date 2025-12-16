!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/raytrace_finalize
!!
!! NAME
!!
!!  raytrace_finalize
!!
!! SYNOPSIS
!!
!!  call raytrace_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the rt unit.
!!
!! ARGUMENTS
!!
!!***
subroutine raytrace_finalize ()

  implicit none
#include "Flash.h"
#include "constants.h"
#ifdef RAYTRACE_3DRT
  call rt_finalize_raytrace_3drt
#endif

  return
end subroutine raytrace_finalize
