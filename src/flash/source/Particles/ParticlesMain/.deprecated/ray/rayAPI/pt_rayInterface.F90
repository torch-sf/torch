!!****ih* source/Particles/localAPI/pt_interface
!!
!! Used to give interfaces to face and point source routines
!! the stub is compiled in depending on the the type that is not used
!! It defines the interfaces with subunit scope.
!! 
!!***
Module pt_rayInterface
#include "constants.h"
#include "Flash.h"

  interface
     subroutine pt_generateRaysFace(dt)
			 real,intent(in)       :: dt
     end subroutine pt_generateRaysFace
  end interface

  interface
     subroutine pt_generateRaysPoint(dt)
			 real,intent(in)       :: dt
     end subroutine pt_generateRaysPoint
  end interface

  interface
     subroutine pt_advanceRaysPoint(dt)
			 real,intent(in)       :: dt
     end subroutine pt_advanceRaysPoint
  end interface

  interface
     subroutine pt_advanceRaysFace(dt)
			 real,intent(in)       :: dt
     end subroutine pt_advanceRaysFace
  end interface

! only needed for async comm
  interface
    subroutine pt_startComm()
    end subroutine pt_startComm
  end interface

end Module pt_rayInterface
