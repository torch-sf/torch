!!****f* source/physics/sourceTerms/Heat/Heat_init
!!
!! NAME
!!  
!!  SN_init
!!
!!
!! SYNOPSIS
!! 
!!  call SN_init()
!!
!!  
!! DESCRIPTION
!!
!!		Perform various initializations (apart from the problem-dependent ones)
!!	  for the heat module.
!!	  Ported from FLASH2.5, 4.4.2012, the code was part of the heat.f90 module,   
!!	  specifically the first_call block
!!	  precalculates SN explosion positions and times 
!!	  also adorned with generous amount of comments and clarifications
!!
!!	  ARGUMENTS
!!	  no arguments from me
!!  
!!
!!***
subroutine SB_init()

	use SB_data
  use mtmod

	use RuntimeParameters_interface, ONLY : RuntimeParameters_get
	use Driver_data, ONLY : dr_dt, dr_restart, dr_simTime
	use SN_data, ONLY: he_meshMe

	implicit none

! some local variables
    logical :: restart
    integer :: m,nside, readstat
    real, dimension(3) :: dir
! for writing headers
    integer, parameter :: funit_evol = 15
    character(len=80)  :: outfile = "SBfeedback.dat"
    character(len=80)  :: outfile2 = "SBcreate.dat"
    integer*8 :: i
! for restart read in, only first column is needed
    real 	:: c0,c1,c2,c3,c4,c5
    real 	:: c6,c7,c8

#include "constants.h"
#include "Flash.h"

! call get_parm_from_context(global_parm_context,"cpnumber", cpnumber)         

! SN rates and numbers
	call RuntimeParameters_get('tsb', sb_tsb)
	call RuntimeParameters_get('nsnmax', sb_nsnmax)
	call RuntimeParameters_get('nsnmin', sb_nsnmin)
	call RuntimeParameters_get('hstarb', sb_hstarb)
	call RuntimeParameters_get('sblife', sb_life)
	call RuntimeParameters_get('sbMaxV', sb_MaxV)
	call RuntimeParameters_get('useSB', useSB)
	call RuntimeParameters_get('useSBrandom', sb_useSBrandom)
	call RuntimeParameters_get('sbTrackV', sb_trackV)
	call RuntimeParameters_get('SBmax', sb_SBmax)
	!call RuntimeParameters_get('accTracer', sb_accTracer)
	!call RuntimeParameters_get('useSNsink', sb_useSNsink)
	!call RuntimeParameters_get('cloudMassPerStar', sn_cloudMassPerStar)
	!call RuntimeParameters_get('outflowFrac', sn_outflowFrac)
	!call RuntimeParameters_get('tracerAccFrac', sn_tracerAccFrac)
	!call RuntimeParameters_get('SFdelay', sn_SFdelay)
	!call RuntimeParameters_get('sinkBulkMotion', sn_sinkBulkMotion)

! first index is 0, minimum value 
	sb_edge(1) = 0.e0

! defaults are nsnmin 4, nsnmax 40
! this actually calculates the cummulative distribution function, hence it actually is ~ n^-1 as one power is integrated out
! probability function is max/(max-min)*(min/(min+n)^2) ~ n^-2
  do m = 2, sb_nsnmax-sb_nsnmin+1
	  sb_edge(m) = float(sb_nsnmax-sb_nsnmin)*(float(sb_nsnmin+m)-1.5) ! -1.5 moves to the center of the second probability bin
	  sb_edge(m) = float(sb_nsnmax)*(float(m)-1.5)/sb_edge(m)
  enddo

! last index is 1, maximum value
	sb_edge(sb_nsnmax-sb_nsnmin+2) = 1.e0

	if(.not. dr_restart) then
		sb_nSN   =  0

		sbPOS(:,:) = 0d0
		sbPOSbuff(:,:) = 0d0

		! write headers
		if(he_meshMe .eq. MASTER_PE) then
			open(funit_evol, file=trim(outfile2), position='APPEND')
			write(funit_evol,'(10(1X,A16))') '[00]SBid', '[01]totSN', '[02]time', '[03]SNinvT','[04]posx', &
				                              '[05]posy', '[06]posz', '[07]velx','[08]vely', &
				                              '[09]velz'
			close(funit_evol)
		endif

		if (he_meshMe .eq. MASTER_PE) then
			open(funit_evol, file=trim(outfile), position='APPEND')
			write(funit_evol,'(9(1X,A16))') '[00]SBid', '[01]nSN', '[02]time', '[03]posx', &
				                            '[04]posy', '[05]posz', '[06]velx','[07]vely', &
				                            '[08]velz'
			close(funit_evol)
		endif
	else
	  ! read in from ascii file 
		open(funit_evol, file=trim(outfile2), status="old", action="read")
! do one empty read to skip header
		read(funit_evol,*,IOSTAT=readstat)

    if(readstat .lt. 0 ) then 
			print*, 'End of File reached'
			sb_nSN = 0
    else
! read through all SN
			do
				read(funit_evol,*,IOSTAT=readstat) c0,c1,c2,c3,c4,c5, &
														 c6,c7,c8

				if(readstat .lt. 0 ) then
					sb_nSN = 0
					print*, 'End of File reached'
					exit
				endif

				if(c3 .gt. dr_simTime) then
					sb_nSN   = c0-1
					exit
		    endif
    	enddo
		endif
		close(funit_evol)
	endif
! this is used for the MPI buffer
! the max allowed number per core is 1000
	sbPOS(:,:) = 0d0
	sbPOSbuff(:,:) = 0d0

	return

end subroutine SB_init
