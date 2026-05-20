!!
!! NAME
!!  
!!
!! SYNOPSIS
!!  
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***
subroutine RadTrans_update_radiation_field(nblk, blklst)
!
!==============================================================================
  !
#include "Flash.h"
#include "constants.h"
  
  use tree, ONLY                      : nodetype
  use Driver_data, ONLY               : dr_globalMe
  use Timers_interface, ONLY          : Timers_start, Timers_stop
  use RadTrans_hybridCharModule, ONLY : dMeanMax, io_conv
!

  implicit none
  
  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)

#include "Flash_mpi.h"

      integer block_no, lb, ierr

      
!==============================================================================
!

      call Timers_start("lambda step")
!
      dMeanMax = 0.0
!
      do block_no = 1, nblk
         lb = blklst(block_no)
         if (nodetype(lb) .eq. LEAF) then
            call rt_mean_intensity_block (lb)
         end if
      enddo
!
      call Timers_stop("lambda step")
!
!     Exchange the global maximun 
!
      call MPI_Allreduce(MPI_IN_PLACE,dMeanMax,1,       &
                         FLASH_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
!
!      call MPI_Barrier (MPI_Comm_World, ierr)
!
      if(dr_globalMe.eq.MASTER_PE) write(io_conv,*) dMeanMax
!
!==============================================================================
!
      return
    end subroutine RadTrans_update_radiation_field
