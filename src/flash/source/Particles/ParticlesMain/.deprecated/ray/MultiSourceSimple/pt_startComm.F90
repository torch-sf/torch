
subroutine pt_startComm()

  use pt_rayAsyncComm, only : ph_startComm
  implicit none

#include "Flash.h"
#include "constants.h"

! starts the counter and allocates bits
! well thats not named well
    call ph_startComm()
  return
end subroutine pt_startComm
