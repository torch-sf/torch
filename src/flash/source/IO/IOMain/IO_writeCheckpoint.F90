!!****if* source/IO/IOMain/IO_writeCheckpoint
!!
!! NAME
!!
!!  IO_writeCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  call IO_writeCheckpoint() 
!!
!!
!!
!! DESCRIPTION
!!
!!  This is a generic call to write the important simulation data to a
!!  checkpoint file.  A checkpoint file writes a few different types of
!!  data to a file, first the physical data like, temperature, pressure, density
!!  etc. which are stored in each cell on the grid.  
!!  Second, in order to recreate the simulation from a checkpoint file a
!!  number of other single quantities are needed as well.  We call these
!!  scalar values which include simTime, dt, nstep, globalNumBlocks etc.
!!  We also store descriptive strings that describe the simulation run.
!!
!!  The same IO_writeCheckpoint routine is called regardless of the type of
!!  file being written, (such as hdf5 parallel, hdf5 serial or pnetcdf)
!!  IO_writeCheckpoint prepares the Grid_ioData (like getting the
!!  globalNumBlocks) and collects the scalars wanting to be checkpointed
!!  from each unit. 
!!  IO_writeCheckpoint then calls four methods, io_initFile, io_writeData, 
!!  IO_writeParticles and io_closeFile.  Each of these routines _is_ specific
!!  to the type of io library used and have their own implementation.  
!!  In addition, io_writeData has its own
!!  implementation for io library and type of grid (UG, Paramesh, or other)
!!
!!
!!  In FLASH IO_writeCheckpoint is called from IO_output (or IO_outputInitial or
!!  IO_outputFinal) IO_output checks to see if enough wall clock time,
!!  simTim, or nsteps has passed to checkpoint.
!!
!!  We have put IO_writeCheckpoint in the API because a user may want to write
!!  a checkpoint at another time or for another reason without having to go through
!!  IO_output.  For most flash users IO_writeCheckpoint will only ever be
!!  called through IO_output and possibly IO_outputFinal and IO_outputInitial.
!!
!! ARGUMENTS
!! 
!!
!! NOTES
!!
!!  For those familiar with FLASH2, breaking up the checkpoint routine into
!!  these four different methods is a change.  Because FLASH3 now supports
!!  different grid packages and we are committed to supporting both
!!  hdf5 and parallel netCDF having each grid and io library writing its
!!  own checkpoint file proved to be a lot of code duplication.  We believe
!!  that while dividing up the checkpoint routines created more files it 
!!  will in the end be easier to maintain.
!!
!!***


subroutine IO_writeCheckpoint()

  use IO_data, ONLY : io_checkpointFileNumber, io_unklabels, &
       io_doublePrecision, io_globalMe,&
       io_outputSplitNum, io_chkptFileID, io_alwaysComputeUserVars
  use Grid_interface, ONLY : Grid_restrictAllLevels, &
    Grid_computeUserVars
  use Logfile_interface, ONLY : Logfile_stamp
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use IO_interface, ONLY : IO_writeParticles,IO_writeUserArray,IO_updateScalars
  use mtmod  ! Torch add-on
  
  implicit none
  
#include "constants.h"
#include "Flash.h"  
#include "io_flash.h"
  include "Flash_mpi.h"

  
  integer :: i 

  character (len=MAX_STRING_LENGTH) :: filename
  character (len=MAX_STRING_LENGTH) :: rngFilename  ! Torch add-on
  
 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  logical :: restrictNeeded

  !double precision for all checkpoint files
  io_doublePrecision = .true.

  !collect data from all the units for outputting
  call IO_updateScalars()


  call io_prepareSimInfo()
  

  call Timers_start("writeCheckpoint")

  call io_restrictBeforeWrite( restrictNeeded)
  if (restrictNeeded .eqv. .true.) then
     call Grid_restrictAllLevels()
  end if
  if(io_alwaysComputeUserVars) call Grid_computeUserVars()

  !---------------------------------------------------------------
  ! open the file -- checkpoints are never 'forced'
  !---------------------------------------------------------------
  IO_TIMERS_START("create file")
  call io_initFile( io_checkpointFileNumber, io_chkptFileID, filename, "_chk_", .false.)
  IO_TIMERS_STOP("create file")


  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[IO_writeCheckpoint] open")
  end if

  
  do i= UNK_VARS_BEGIN,UNK_VARS_END
     call Simulation_mapIntToStr(i, io_unklabels(i),MAPBLOCK_UNK)
  end do
  
  call io_writeData( io_chkptFileID)

  !if particles are  not included in this simulation this
  !function will be empty

  call IO_writeParticles( .true.)

  call IO_writeUserArray()

  !----------------------------------------------------------------------
  ! close the file
  !----------------------------------------------------------------------
  IO_TIMERS_START("close file")
  call io_closeFile( io_chkptFileID)
  IO_TIMERS_STOP("close file")

  !----------------------------------------------------------------------
  ! Torch add-on: save random number generator state
  !----------------------------------------------------------------------
  if (io_globalMe == MASTER_PE) then
     call io_getOutputName( io_checkpointFileNumber, "txt", "_RNG_state_", rngFilename, .false.)
     call mtsavef(rngFilename,'f')  ! formatted file is plaintext, so machine portable
  end if

  !increment the checkpoint number unless it is a dump checkpoint file
  !DUMP_IOFILE_NUM typically 9999 or some large number, located in constants.h 
  if(io_checkpointFileNumber /= DUMP_IOFILE_NUM) then
     io_checkpointFileNumber = io_checkpointFileNumber + 1
  end if

  call Timers_stop("writeCheckpoint")


  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[IO_writeCheckpoint] close")
     print *, '*** Wrote checkpoint file to ', trim(filename),  ' ****' 
     write (strBuff(1,1), "(A)") "file"
     write (strBuff(1,2), "(A)") trim(filename)
  end if

  return
end subroutine IO_writeCheckpoint




