!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/raytrace_init.F90
!!
!!  NAME 
!!
!!  raytrace_init
!!
!!  SYNOPSIS
!!
!!  call raytrace_init()
!!
!!  DESCRIPTION 
!!    Initialize local data for the ray-tracer
!!
!!***
subroutine raytrace_init()

  use RadTrans_data, ONLY: rt_meshMe, rt_acrossme
  use raytrace_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  call RuntimeParameters_get('useRayTrace', rt_useRayTrace)
  call RuntimeParameters_get('rt_nPhi', rt_nPhi)
  call RuntimeParameters_get('rt_nTheta', rt_nTheta)

     !
  call RuntimeParameters_get('rt_dirX', rt_dirX)
  call RuntimeParameters_get('rt_dirY', rt_dirY)
  call RuntimeParameters_get('rt_dirZ', rt_dirZ)
     !
     !
  call RuntimeParameters_get('rt_ALI', rt_ALI)
  call RuntimeParameters_get('rt_epsilon', rt_epsilon)
  IF(rt_epsilon.eq.1..AND.rt_ALI.eq.1) THEN
    ! Accelerated lambda iteration (ALI) with epsilon=1 automatically
    ! reduces to normal lambda iteration. Therefore we can deactive
    ! ALI to save some computations.
    rt_ALI=0
  END IF

  call RuntimeParameters_get('rt_maxNrOfBoundIter', rt_maxNrOfBoundIter)
  call RuntimeParameters_get("rt_nrOfAngleGroups", rt_nrOfAngleGroups)

  call RuntimeParameters_get("rt_irradiation", irradiation) 

#ifdef RT_HEALPIX
  call RuntimeParameters_get("rt_healpix_nSide", rt_healpix_nSide)
#endif

  call RuntimeParameters_get("rt_healpix_randomize", rt_healpix_randomize)
     !
  
#ifdef RAYTRACE_3DRT
  call rt_init_raytrace_3drt
#endif

end subroutine raytrace_init
