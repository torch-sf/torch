!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   deallocates memory for ray tracing
!!   and finalizes communication

subroutine pt_finalizeRays()

  use Particles_rayData, only : raysIntProp, raysRealProp, ph_localRays

  use Driver_interface, ONLY : Driver_abortFlash

  use pt_rayAsyncComm, only : ph_finalizeComm
  implicit none

#include "Flash.h"
#include "constants.h"

  if(ph_localRays .gt. 0) then
    print*,'rays left',ph_localRays
    call Driver_abortFlash ("raytracing:  missed rays!")
  endif

! local buffer	
!  deallocate(raysIntProp)
!  deallocate(raysRealProp)

! MPI buffer 
!  deallocate(rayDestBuf)
!  deallocate(raySourceBuf)
!  deallocate(rayDestBufTarget)

  call ph_finalizeComm()

 return
end subroutine pt_finalizeRays
