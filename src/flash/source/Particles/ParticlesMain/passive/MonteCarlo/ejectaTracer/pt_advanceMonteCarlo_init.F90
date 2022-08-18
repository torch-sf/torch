!===============================================================================

subroutine pt_advanceMonteCarlo_init()
!===============================================================================

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get  
	use pt_advanceMonteCarlo_data

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  !print*,'There are HEALPix particles! Hooray!'
  call RuntimeParameters_get("ej_TempUnfreeze",  ej_TempUnfreeze)
  call RuntimeParameters_get("ej_TimeUnfreeze",  ej_TimeUnfreeze)
  call RuntimeParameters_get("ej_shockUnfreeze", ej_shockUnfreeze)

!-------------------------------------------------------------------------------

	return
end subroutine pt_advanceMonteCarlo_init
