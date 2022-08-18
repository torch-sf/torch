!!****h* source/Particles/Particles_interface
!!
!! This is the header file for the particles module that defines its
!! public interfaces.
!!***

Module SC_interface

  implicit none

#include "Flash.h"
#include "constants.h"

	interface
	  subroutine SC_init (restart)
	    logical, intent(IN) :: restart
	  end subroutine SC_init
	end interface

	interface
	  subroutine SC_clusterSample (Mass, seed)
	    real, intent(IN) :: Mass, seed
	  end subroutine SC_clusterSample
	end interface



end Module SC_interface
