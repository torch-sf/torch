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
subroutine SN_init()

! only few parameters are needed 
! there really should be an except statement
! ah one more: day-umm
	use SN_data
  use mtmodSN

	use RuntimeParameters_interface, ONLY : RuntimeParameters_get
	use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
	use Driver_data, ONLY : dr_dt, dr_restart, dr_simTime
  use IO_data, ONLY : io_checkpointFileNumber 

	implicit none

! some local variables
	logical :: restart
	integer :: tracerDimN,i,j,k,m,readstat
	real    :: space
! for writing headers
  integer, parameter :: funit_evol = 15
  character(len=80)  :: outfile = "SNfeedback.dat"
  real 	:: c0,c1,c2,c3,c4,c5
  real 	:: c6,c7,c8,c9
  character(len=80)  :: convert

#include "constants.h"
#include "Flash.h"

! SN rates and numbers
	call RuntimeParameters_get('tsn1', he_tsn1)
	call RuntimeParameters_get('tsn2', he_tsn2)
	call RuntimeParameters_get('hstar1', he_hstar1)
	call RuntimeParameters_get('hstar2', he_hstar2)
	call RuntimeParameters_get('nsndt', he_nsndt)
	call RuntimeParameters_get('SNminstep', he_SNminstep)
	call RuntimeParameters_get('erstar1', he_erstar1)
	call RuntimeParameters_get('erstar2', he_erstar2)

! comp. domain
	call RuntimeParameters_get('xmin',he_imin)
	call RuntimeParameters_get('xmax',he_imax)
	call RuntimeParameters_get('ymin',he_jmin)
	call RuntimeParameters_get('ymax',he_jmax)
	call RuntimeParameters_get('zmin',he_kmin)
	call RuntimeParameters_get('zmax',he_kmax)

! SN 
	call RuntimeParameters_get('r_init',he_r_init)
	call RuntimeParameters_get('r_exp_max',he_r_exp_max)
	call RuntimeParameters_get('mejc',he_mejc)
	call RuntimeParameters_get('exp_energy',he_exp_energy)
	call RuntimeParameters_get('SNmapToGrid',he_SNmapToGrid)
!	call RuntimeParameters_get('SNmapToGridCorner',he_SNmapToGridCorner)
	call RuntimeParameters_get('nms',he_nms)
	call RuntimeParameters_get('sn_max_temp',sn_max_temp)

! switches 
	call RuntimeParameters_get('useSN', he_useSN)
	call RuntimeParameters_get('useSNrandom', he_useSNrandom)
	call RuntimeParameters_get('stratifySN', he_stratifySN)
	call RuntimeParameters_get('injMass', he_injMass)
	call RuntimeParameters_get('injVol',  he_injVol)

! tracers
	call RuntimeParameters_get('useSNTracer',  he_useSNTracer)
	call RuntimeParameters_get('TracerPerSN', he_tracerPerSN)
	call RuntimeParameters_get('veltracer',  he_veltracer)
	call RuntimeParameters_get('MCTracerShellDis', sn_MCTracerShellDis)

! Errrybody should know these
	call Driver_getMype(MESH_COMM,he_meshMe)
	call Driver_getNumProcs(MESH_COMM,he_meshNumProcs)

  if(he_useSNTracer) then
    he_oldStartTag = 0 
    if (he_meshMe .eq. MASTER_PE) then
      print*,'mapping SB particle properties to ejecta particle properties' 
      print*,'mass prop-> dens'
      print*,'nsn prop -> temp'
      print*,'tsn prop -> SN number'
      print*,'sbid prop -> frozen? flag'
    endif
  endif

!====================================================
! this precalculates possible locations and times
! for SB explosions over the simulation time 
! also assigns number of SN per SB  
!====================================================
!!	he_seedsize = 2
!	allocate (he_seed(he_seedsize))
!	call random_seed(SIZE = he_seedsize)
!	he_seed(1:he_seedsize) = (/123456789, 987654321/)
!	call random_seed(PUT = he_seed(1:he_seedsize))

	if(he_injMass .and. he_injVol) then
    if (he_meshMe .eq. MASTER_PE) then
        print*,'Mass AND volume injection criterion not possible, choosing mass'
				he_injVol = .false.
    endif
	endif

! when restarting just call the random number he_nSN times, to get to where we should be
! it's not the most elegant method but eh.
  if(dr_restart) then
		open(funit_evol, file=trim(outfile),status="old",action="read")
		read(funit_evol,*,IOSTAT=readstat) 
    do
			read(funit_evol,*,IOSTAT=readstat) c0,c1,c2,c3,c4,c5, &
															 c6,c7,c8,c9

			if(readstat .lt. 0) then
				exit
			endif

      if(c4 .gt. dr_simTime) then
				he_nSN   = c0-1
				exit
      endif
				
    enddo

! if checkpoint file time equals last SN time
		if(he_nSN .eq. 0 .and. c0 .gt. 0) then
			he_nSN = c0
		endif 

    if (he_meshMe .eq. MASTER_PE) then
        print*,'Restarting from SN number', he_nSN
    endif

		he_newDt = 1e99
		close(funit_evol)

    write(convert, '(i4.4)') io_checkpointFileNumber - 1
    call mtgetf('RNG_state_SN_' // trim(convert),'u')
  else
	  he_nSN = 0
	  he_newDt = 1e99
! write header

    if (he_meshMe .eq. MASTER_PE) then
      open(funit_evol, file=trim(outfile), position='APPEND')
      write(funit_evol,'(11(1X,A16))') '[00]n_SN', '[01]type', '[02]n_timestep', '[03]n_tracer','[04]time', '[05]posx', &
                                    '[06]posy', '[07]posz', '[08]radius', '[09]mass'
      close(funit_evol)
    endif
  endif  

! nms, number of shells
  he_maxradius = ((he_r_exp_max / he_r_init)**( 1.0/ (he_nms-1.0)))**(he_nms-1.0)*he_r_init

! for the SB
	call SB_init()

  write(convert, '(i4.4)') he_meshMe
  sn_outfile = trim(sn_outfile) //trim(convert)// trim('.log')


	if ( .not. dr_restart) then
 
		open(sn_funit_evol, file=trim(sn_outfile), position='APPEND')
		write(sn_funit_evol,'(13(1X,A16))') '[00]step','[01]SNtag' ,'[02]time', '[03]currT', '[04]maxT', '[05]facT', &
                     '[06]lostE','[07]massFrac' ,'[08]i', '[09]j', '[10]k','[11]blockID', '[12]proc'
		close(sn_funit_evol)
	endif 



! create 3d position array for tracers in (0..1)^3 box
!  if(he_useSNTracer .and. he_veltracer) then
!    allocate(he_tracerPosX(1:he_tracerPerZone,1:he_tracerPerZone,1:he_tracerPerZone))
!    allocate(he_tracerPosY(1:he_tracerPerZone,1:he_tracerPerZone,1:he_tracerPerZone))
!    allocate(he_tracerPosZ(1:he_tracerPerZone,1:he_tracerPerZone,1:he_tracerPerZone))
! fill array, 4 because we want some distance from the boundaries
!		space = 1.0/(he_tracerPerZone)
!    do i = 1, he_tracerPerZone
!      do j = 1, he_tracerPerZone
!        do k = 1, he_tracerPerZone
!		      he_tracerPosX(i,j,k) = (i-0.5)*space
!		      he_tracerPosY(i,j,k) = (j-0.5)*space
!		      he_tracerPosZ(i,j,k) = (k-0.5)*space
!        enddo
!      enddo
!    enddo
!  endif

	return
end subroutine SN_init
