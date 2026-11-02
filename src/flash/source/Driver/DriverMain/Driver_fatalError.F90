
!!****if* source/Driver/DriverMain/Driver_fatalError
!!
!! NAME
!!
!!  Driver_fatalError
!!
!! SYNOPSIS
!!
!!  Driver_fatalError(character(len=*)(IN) :: errorMessage)
!!
!! DESCRIPTION
!!
!!  Writes a minimal error, a traceback and exits the code.
!!
!!
!! ARGUMENTS
!!
!!  errorMessage -    A string to write to the error log file.
!!
!!***


subroutine Driver_fatalError(errorMessage)

  use Driver_data, ONLY : dr_globalMe, dr_globalComm
  use iso_fortran_env, ONLY : error_unit
  implicit none

  include "Flash_mpi.h"

  character(len=*), intent(in) :: errorMessage

  integer, parameter :: errcode = 1
  integer :: ierr

  ! Print error message
  write(error_unit, '(A,I0,2A)') &
       "[DRIVER_FATAL] PE ", dr_globalMe, ": ", trim(errorMessage)
  flush(error_unit)

#ifdef __GFORTRAN__
  ! Print traceback (only GNU Fortran)
  call backtrace()
#endif

  ! Terminate all MPI processes
  call MPI_Abort(dr_globalComm, errcode, ierr)

  stop 1

end subroutine Driver_fatalError
