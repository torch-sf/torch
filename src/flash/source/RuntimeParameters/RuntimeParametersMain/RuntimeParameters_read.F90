!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_read
!!
!!
!! NAME
!!  RuntimeParameters_read
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_read( character(in), (len=MAX_STRING_LEN) :: parmfile)
!!
!!
!! DESCRIPTION
!!
!!
!!         Parses the parameter file parmfile, which contains
!!         job-dependent parameter definitions.  Syntax of
!!         parameter file lines is:
!!
!!                       # comment               rest of line is a comment
!!                       variable = value        set variable to value
!!                       strvar   = "word"     set string variable strvar
!!                                                to word
!!
!!         Some error checking is performed:  if a variable is
!!         unrecognized, it is ignored.  Syntax errors are signaled
!!         along with the offending lines.  However, type mismatch
!!         errors force termination.  Also, no checking is done to
!!         determine whether all of the variables have received some
!!         value.  In case of repeated definitions of the same
!!         variable, the last definition overrides the others.
!!
!! ARGUMENTS
!!
!!   parmfile :       the name of the parameter file to read
!! 
!! NOTES
!!
!!   This routine is called during FLASH initialization to read the flash.par file.
!!   In general it would not be used by users.  Instead, use the restart = .true.
!!   capability within a flash.par to restart a run from a checkpoint file,
!!   using a different flash.par configuration.
!!
!!
!!
!!***



subroutine RuntimeParameters_read (parmfile)

  use RuntimeParameters_data, ONLY : parameter
  use Driver_interface, ONLY : Driver_fatalError
  use RuntimeParameters_interface, ONLY : RuntimeParameters_set, &
    RuntimeParameters_mapStrToInt

  use nameValueLL_data, ONLY : name_real, name_int, name_log, name_str, &
       nameValueLL_getType


  implicit none

#include "constants.h"        
#include "Flash.h"

  character(len=MAX_STRING_LENGTH), intent(in)    :: parmfile
  character(len=MAX_STRING_LENGTH+16)       :: inline  ! some extra space for 'varname =' + long string
  character(len=MAX_STRING_LENGTH)          :: name, value

  logical                                   :: quoted, inQuote, tempLog, debug
  integer                                   :: ioStatus, first(MAX_STRING_LENGTH)
  integer                                   :: last(MAX_STRING_LENGTH)
  integer                                   :: nwords, tempInt, lineLen
  integer                                   :: nameType, i
  real                                      :: tempReal
  character, parameter                      :: quote = '"', space = ' ', & 
       &                                             hash = '#', & 
       &                                             equals = '='
  
  character(len=MAX_STRING_LENGTH)          :: abortMessage
  
  
  !  Debugging option allows us to check parser.
  
  debug = .false.
  
  !  Open parameter file.

    
  open (1, file=parmfile, status='old', iostat=ioStatus)


  
  if (ioStatus /= 0) then
     write (*,*) 'read :  unable to open parameter file', & 
          &              ' "', parmfile(1:max(len_trim(parmfile),1)), '"'
     abortMessage = "Error: unable to open parameter file " // parmfile(1:max(len_trim(parmfile),1))
     call Driver_fatalError(abortMessage)
     stop
  endif
  
  do                ! While not EOF, read lines.
     
     read (1, '(A96)', iostat=ioStatus) inline
     if (ioStatus > 0) then
        write (*,*) 'read_parameters:  I/O error ', ioStatus, & 
             &                'while reading parameter file'
        abortMessage = "Error: I/O error while reading parameter file"
        call Driver_fatalError(abortMessage)
        stop
     endif
     if (ioStatus < 0) exit
     
     !                       Ignore empty, blank, or comment (#) lines.
     
     i = index( inline, hash )
     if (i == 0) then
        lineLen = len_trim(inline)
     else
        lineLen = len_trim(inline(:i-1))
     endif
     if (lineLen == 0) cycle
     
     !                       Convert nonquoted tabs to spaces.
     
     inQuote = .false.
     do i = 1, lineLen
        if (inline(i:i) == quote) inQuote = .not. inQuote
        if ((inline(i:i) == TAB_CHAR) .and. (.not. inQuote))  & 
               &          inline(i:i) = space
     enddo
     if (len_trim(inline(1:lineLen)) == 0) cycle
       
       !             Find first and last character positions of all words.
       
     call findWords (inline(1:lineLen), first, last, & 
            &                     MAX_STRING_LENGTH, nwords)
       
     if (debug) then
        write (*,*) 'line:'
        do i = 1, nwords
           write (*,*) '  ', first(i), last(i), & 
                  &                    '<', inline(first(i):last(i)), '>'
        enddo
     endif
       
       !                       Catch syntax errors; accepted syntax:
       !                         variable = value
       !                       Value can be multiple words if double-quoted.
       
       ! fewer than three words
     if (nwords < 3) then
        call nameSyntaxError (inline)
        cycle
     endif
       ! second word not an =
     if (inline(first(2):last(2)) /= equals) then
        call nameSyntaxError (inline)
        cycle
     endif
       ! unmatched quotes
     if ( (inline(first(3):first(3)) == quote) .neqv. & 
            &         (inline(last(nwords):last(nwords)) == quote) ) then
        call nameSyntaxError (inline)
        cycle
     endif
       ! no quotes if more than three words
     quoted = (inline(first(3):first(3)) == quote)
     if ( (nwords > 3) .and. (.not. quoted) ) then
        call nameSyntaxError (inline)
        cycle
     endif
       
       !                       Syntax OK.  Isolate the parameter name.
       
     name = inline(first(1):last(1))
       
       !                       Isolate the value.
       
     if (quoted) then
        value = inline(first(3)+1:last(nwords)-1)
     else
        value = inline(first(3):last(3))
     endif
       
       ! Find out whether the name specified really
       ! exists, and if so, what its type is.  If the
       ! name is valid, set its value.
       
     call nameValueLL_getType(parameter, name, nameType)
              
#ifndef STRICT_PARAMS
#define STRICT_PARAMS 0
#endif
     select case (nameType)
          
       case (name_real)
          read (value, *, iostat=ioStatus) tempReal
          if (ioStatus == 0) then
             call RuntimeParameters_set(name, tempReal)
          else
             call nameSyntaxError (inline)
          endif
          
       case (name_int)
          read (value, *, iostat=ioStatus) tempInt
          if (ioStatus == 0) then
             call RuntimeParameters_set (name, tempInt)
          else
             call RuntimeParameters_mapStrToInt(value, tempInt)
             if (tempInt == NONEXISTENT) then
                write (*,*) 'RuntimeParameters_read: "', &
                     trim(name), '" found with unrecognized value ','"' &
                     // trim(value), '",',' aborting...'
                abortMessage = 'RuntimeParameters_read: "' &
                     // name(1:len_trim(name)) // '" with unrecognized value.'
                call Driver_fatalError(abortMessage)
             endif
             call RuntimeParameters_set (name, tempInt)
          endif

          
       case (name_str)
          if (index(value, quote) == 0) then
             call RuntimeParameters_set (name, value)
          else
             call nameSyntaxError (inline)
          endif
          
       case (name_log)
          read (value, *, iostat=ioStatus) tempLog
          if (ioStatus == 0) then
             call RuntimeParameters_set (name, tempLog)
          else
             call nameSyntaxError (inline)
          endif
          
       case default

          if (STRICT_PARAMS.ne.0) then
            abortMessage = 'RuntimeParameters_read:  encountered unknown parameter "' &
                                      // name(1:len_trim(name)) // '" aborting...'
            call Driver_fatalError(abortMessage)
          else
            write (*,*) 'RuntimeParameters_read:  ignoring unknown parameter ' & 
                 &                    ,'"', name(1:len_trim(name)), '"...'
          endif
          

          call rp_storeIgnoredParams(name)


     end select
       
       ! Done with line.
       
  enddo
    
    ! Done with file.
    
  close (1)
    


  return
end subroutine RuntimeParameters_read
