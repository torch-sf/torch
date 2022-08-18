!!****f* source/physics/sourceTerms/SN/SN_finalize
!!
!! NAME
!!
!!  SN_finalize
!!
!!
!! SYNOPSIS
!!
!!  call SN_finalize
!!
!!
!! DESCRIPTION
!!
!!  Cleans up for the SN module.
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!***
subroutine SN_finalize ()

  use mtmodSN, ONLY : grndSN, mtsavef

  use SN_data, ONLY : he_meshMe, he_outputDir, he_useSN
  use IO_data, ONLY : io_checkpointFileNumber, io_summaryOutputOnly

  implicit none

#include "constants.h"

  character (len=MAX_STRING_LENGTH) :: convert

  if (.not. he_useSN) return

  ! Logic adapted from IO_outputFinal().  Always dump RNG state, regardless of
  ! whether checkpoint was dumped by IO_output or IO_outputFinal in last
  ! Driver_evolveFlash loop.
  if (.not. io_summaryOutputOnly) then

    ! save random number state, independent of flash4 random number stream
    ! use only main thread for this, as all other should have the same state
    if(he_meshMe == MASTER_PE) then
      write(convert, '(i4.4)') io_checkpointFileNumber - 1
      call mtsavef(trim(he_outputDir) // 'RNG_state_SN_' // trim(convert),'u')
    endif

  endif

  return

end subroutine SN_finalize
