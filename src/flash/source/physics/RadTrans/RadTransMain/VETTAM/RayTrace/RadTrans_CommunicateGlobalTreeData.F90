!*******************************************************************************

! Routine:      RadTrans_CommunicateGlobalTreeData

! Description:  Makes sure all processors have a copy of the structure of the
!               AMR mesh hierarchy.


module ModuleRadTransCommunicateGlobalTreeData

contains

  subroutine RadTrans_CommunicateGlobalTreeData

!===============================================================================

#include "Flash.h"
#include "constants.h"
    
  use RadTrans_TreeCommData 
  use tree, ONLY: lrefine, bnd_box
  
  implicit none

#include "Flash_mpi.h"

  integer                            :: ierr

!===============================================================================


  call mpi_allgather (lrefine, MAX_LOCBLOCKS, FLASH_INTEGER, &
                      g_lrefine, MAX_LOCBLOCKS, FLASH_INTEGER, &
                      MPI_COMM_WORLD, ierr)

! Gather the block boundig boxes.

  call mpi_allgather (bnd_box, 2*MAX_NDIM*MAX_LOCBLOCKS, FLASH_REAL, &
                      g_bndbox, 2*MAX_NDIM*MAX_LOCBLOCKS, FLASH_REAL, &
                      MPI_COMM_WORLD, ierr)
  return

!===============================================================================

end subroutine RadTrans_CommunicateGlobalTreeData

end module ModuleRadTransCommunicateGlobalTreeData
