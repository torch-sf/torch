!!****if* source/IO/IOMain/IO_readCheckpoint
!!
!! NAME
!!
!!  IO_readCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_readCheckpoint()
!!
!!
!!
!! DESCRIPTION
!!
!!  IO_readCheckpoint is a generic subroutine that retrieves
!!  the unklabels and then calls io_readData
!!  which is specific to pnetcdf, hdf5 and
!!  the necessary grid, UG, paramesh etc.
!!  io_readData reads a checkpoint file and reinitializes
!!  grid and scalar values to resume the run
!!  
!!
!!
!! ARGUMENTS
!!  
!!
!!
!!***

subroutine IO_readCheckpoint()

  use IO_data, ONLY : io_checkpointFileNumber, io_unklabels, io_chkptFileID
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use IO_interface, ONLY : IO_readParticles,IO_readUserArray
  implicit none
#include "Flash.h"
#include "constants.h"

  integer :: i
  
  call Timers_start("IO_readCheckpoint")

  !get the string names for the unknown variables
  do i= UNK_VARS_BEGIN,UNK_VARS_END
     call Simulation_mapIntToStr(i, io_unklabels(i),MAPBLOCK_UNK)
  end do

  call io_readData()

  call io_rescaleCellBoxes()

 ! call Grid_fixupCheckpointData
 
  call IO_readParticles()

  call IO_readUserArray()

  ! close the file
  call io_closeFile( io_chkptFileID)

  !increment the checkpoint number
  io_checkpointFileNumber = io_checkpointFileNumber + 1

  call Timers_stop("IO_readCheckpoint")

  return
end subroutine IO_readCheckpoint

