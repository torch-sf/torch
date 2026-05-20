!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! subroutine Particles_sinkWind()
!!
!! Calls the wind injection subroutine over all the sink
!! particles. 
!!
!! Joshua Wall 5/26/2016
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Particles_sinkWind(dt)

use Particles_interface, ONLY: Particles_wind
use Particles_sinkData, ONLY: localnp, particles_local

use Particles_sort

#include "Flash.h"
#include "constants.h"

real, intent(in) :: dt

real    :: dmdt
real    :: xloc, yloc, zloc
integer :: i



! I'm a stub.

end subroutine Particles_sinkWind
