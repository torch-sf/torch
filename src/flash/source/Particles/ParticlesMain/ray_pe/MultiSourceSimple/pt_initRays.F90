!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   allocates memory for ray tracing
!!   MPI and local buffers
!!
!! Input:
!!
!!
!! TODO different sizes for MPI and local buffers

subroutine pt_initRays()

!  use Particles_data, only : pt_meshMe
  use Particles_rayData, only : ph_maxNRays, &
      ph_numIntProp, ph_numRealProp, raysIntProp, raysRealProp

  use pt_rayAsyncComm, only : ph_initComm
  implicit none

#include "Flash.h"
#include "constants.h"

! local buffer	
!  allocate(raysIntProp (ph_numIntProp,  ph_maxNRays))
!  allocate(raysRealProp(ph_numRealProp, ph_maxNRays))

! MPI buffer 
!  allocate(rayDestBuf  (ph_transProp, ph_maxNRays))
!  allocate(raySourceBuf(ph_transProp, ph_maxNRays))

!  allocate(rayDestBufTarget(1:ph_maxNRays))

! mark as empty
  raysIntProp  = -1
  raysRealProp = -1
!	rayDestBuf   = -1 
!	raySourceBuf = -1 
!	rayDestBufTarget = -1

  call ph_initComm()

! update source information 
!#ifdef WINDSOURCE
!  call pt_sourceUpdate()
!#endif

  return
end subroutine pt_initRays
