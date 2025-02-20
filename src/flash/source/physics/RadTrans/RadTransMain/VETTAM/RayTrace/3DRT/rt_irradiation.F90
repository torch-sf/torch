!!****if* source/physics/RadTrans/RadTransMain/HybridChar3DRTFLD/RayTrace/3DRT/rt_irradiation
!!
!! NAME
!!
!!  rt_irradiation
!!
!! SYNOPSIS
!!
!!  call rt_irradiation ()
!!
!! DESCRIPTION
!!
!!  Set the irradiation in the specified direction.
!!  This can be overwritten by a user implementation of the routine
!!  to achieve different irradiation from different directions.
!!
!! ARGUMENTS
!!
!!***
subroutine rt_irradiation(direction,intensity)
#include "Flash.h"
#include "constants.h"

  use raytrace_data
#ifdef ERAD_VAR
  use RadTrans_data, ONLY: current_band
#endif
  implicit none

  REAL, DIMENSION(NDIM), INTENT(IN) :: direction
  REAL, INTENT(OUT) :: intensity

  intensity = irradiation
  !Set vacuum BCs for UV bands
  if(current_band .ne. 'IR') &
    & intensity = 0.0
  return
end subroutine rt_irradiation
