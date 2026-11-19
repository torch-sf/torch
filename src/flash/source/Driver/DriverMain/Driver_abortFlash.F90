!!****if* source/Driver/DriverMain/Driver_abortFlash
!!
!! NAME
!!
!!  Driver_abortFlash
!!
!! SYNOPSIS
!!
!!  Driver_abortFlash(character(len=*)(IN) :: errorMessage)
!!
!! DESCRIPTION
!!
!!  Write an error message to the logfile and abort FLASH.
!!  Attempts to shut down all processes (using MPI_Abort()).
!!  If you wish to call Driver_abortFlash from a 'c' routine
!!  use the API function Driver_abortFlashC
!!
!!
!! ARGUMENTS
!!
!!  errorMessage -    A string to write to the logfile (presumably 
!!                    indicating what went wrong).
!!
!!***

subroutine Driver_abortFlash (errorMessage, withOutput)
 
  use iso_fortran_env, ONLY : error_unit
  use Driver_data, ONLY : dr_globalMe,dr_globalComm, dr_eachProcWritesOwnAbortLog, &
       dr_abortPause
  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp, Logfile_open, Logfile_close, &
       Logfile_getDateTimeStr
  use IO_interface, ONLY : IO_writeCheckpoint
  implicit none

  include "Flash_mpi.h"
#include "constants.h"          
#include "Flash.h"

  logical, intent(in), optional :: withOutput
  character(len=*), intent(in) :: errorMessage
  character(len=MAX_STRING_LENGTH*2) :: errorMessagePlusAbort
  character(len=28)                  :: dateTimeStr
  logical, parameter :: forceStampMessage = .false. !disabled alternative
  logical :: writeOutput

  ! To be passed to MPI_Abort
  integer, parameter::      errcode = 1
  integer :: logUnit
  integer :: ierr

  ! Use to generate a message giving the processor #.
  character(len=MAX_STRING_LENGTH):: peMessage = &
            & "[DRIVER_ABORT] Driver_abort() called by PE "

  
  writeOutput = .false.
  if (present(withOutput)) writeOutput = withOutput

  if (dr_globalMe .eq. MASTER_PE) then
     print *, " "
     if (dr_eachProcWritesOwnAbortLog) then
        print *, "Driver_abort called. See log files for details."
     else
        print *, "Driver_abort called. See log file for details."
     end if
     print *, "Error message is ", errorMessage  

     if (dr_abortPause .LE. 0) then
        print *, "Calling MPI_Abort() for immediate shutdown!"
     else
101     format(" Calling MPI_Abort() for shutdown in", I4," seconds!")
        write(*,101) dr_abortPause
     end if
     print  *, ''
  else
     print *, "DRIVER_ABORT: ", errorMessage  
  end if

  write (peMessage(len_trim(peMessage)+1:), *) dr_globalMe

  ! make a complete string out of error message
  errorMessagePlusAbort = "abort_message " // errorMessage

  ! Write first line to the regular log file
  if (forceStampMessage) then
     call Logfile_stampMessage( peMessage, forceStampMessage)
  else
     call Logfile_stamp( peMessage(16:), "[DRIVER_ABORT]")
  end if

  if (dr_eachProcWritesOwnAbortLog) then
     call Logfile_getDateTimeStr(dateTimeStr)
     ! Write two lines to per-PE log file
     call Logfile_open(logUnit,logUnitLocal=.TRUE.)
100  format(" ",A," ",A)
     write (logUnit, 100) dateTimeStr, trim(peMessage)
     write (logUnit, *) trim(errorMessagePlusAbort)
     call Logfile_close(logUnitLocal=.TRUE.)
  end if

  ! Write second line to the regular log file
  if (forceStampMessage) then
     call Logfile_stampMessage( errorMessagePlusAbort, forceStampMessage)
  else
     call Logfile_stamp( errorMessage, "abort_message")
  end if

  call Logfile_close()

  ! Print error message
  write(error_unit, '(A,I0,2A)') &
       "[DRIVER_FATAL] PE ", dr_globalMe, ": ", trim(errorMessage)
  flush(error_unit)

#ifdef __GFORTRAN__
  ! Print traceback (only GNU Fortran)
  call backtrace()
#endif

  if (writeOutput) then
     if (dr_globalMe .eq. MASTER_PE) then
        print *, "Writing checkpoint and plot files"
     end if
     call IO_writeCheckpoint()
     call IO_writePlotfile(.true.)  ! mark forced
  endif

  if (dr_abortPause > 0) then
     call dr_sleep(dr_abortPause)
  end if
  
  ! Terminate all MPI processes
  call MPI_Abort(dr_globalComm, errcode, ierr)

  stop 1
  
  return
end subroutine Driver_abortFlash
