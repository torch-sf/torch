!!****f* source/physics/sourceTerms/SN/Simple/SN_init
!!
!! NAME
!!  
!!  SN_init
!!
!! SYNOPSIS
!! 
!!  call SN_init()
!!  
!! DESCRIPTION
!!
!!    Perform various initializations (apart from the problem-dependent ones)
!!    for the SN module.
!!    Ported from FLASH2.5, 4.4.2012, the code was part of the heat.f90 module,   
!!    specifically the first_call block
!!    precalculates SN explosion positions and times 
!!
!! ARGUMENTS
!!
!!    none
!!
!!***
subroutine SN_init()

  use SN_data
  use mtmodSN, ONLY : grndSN, mtgetf

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_abortFlash
  use Driver_data, ONLY : dr_restart, dr_simTime
  use IO_data, ONLY : io_checkpointFileNumber 

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: i
  real    :: rateadj1, rateadj2
  ! for SN log read/write
  real    :: c0,c1,c2,c3,c4,c5
  real    :: c6,c7,c8,c9
  integer :: pos, readstat
  logical :: exst
  character(len=MAX_STRING_LENGTH) :: convert

  ! SN rates and numbers
  call RuntimeParameters_get('tsn1',      sn_tsn1)
  call RuntimeParameters_get('tsn2',      sn_tsn2)
  call RuntimeParameters_get('sn_tstop',  sn_tstop)
  call RuntimeParameters_get('hstar1',    sn_hstar1)
  call RuntimeParameters_get('hstar2',    sn_hstar2)
  call RuntimeParameters_get('nsndt',     sn_nsndt)
  call RuntimeParameters_get('SNminstep', sn_SNminstep)

  ! comp. domain
  call RuntimeParameters_get('xmin',sn_imin)
  call RuntimeParameters_get('xmax',sn_imax)
  call RuntimeParameters_get('ymin',sn_jmin)
  call RuntimeParameters_get('ymax',sn_jmax)
  call RuntimeParameters_get('zmin',sn_kmin)
  call RuntimeParameters_get('zmax',sn_kmax)

  ! SN 
  call RuntimeParameters_get('r_init',     sn_r_init)
  call RuntimeParameters_get('r_exp_max',  sn_r_exp_max)
  call RuntimeParameters_get('mejc',       sn_mejc)
  call RuntimeParameters_get('exp_energy', sn_exp_energy)
  call RuntimeParameters_get('SNmapToGrid',sn_SNmapToGrid)
  call RuntimeParameters_get('nms',        sn_nms)
  call RuntimeParameters_get('sn_max_temp',sn_max_temp)

  ! switches 
  call RuntimeParameters_get('useSN',       sn_useSN)
  call RuntimeParameters_get('stratifySN',  sn_stratifySN)
  call RuntimeParameters_get('sn_kinetic',  sn_kinetic)

  ! more controls for field SN distribution
  call RuntimeParameters_get('sn_fieldMode',sn_fieldMode)
  call RuntimeParameters_get('sn_single_x', sn_single_x)
  call RuntimeParameters_get('sn_single_y', sn_single_y)
  call RuntimeParameters_get('sn_single_z', sn_single_z)

  ! I/O
  call RuntimeParameters_get('output_directory', sn_outputDir)  ! IO/IOMain

  call Driver_getMype(MESH_COMM, sn_meshMe)

  !====================================================
  ! adjust SN rates for z stratification
  !====================================================
  ! strategy: include only SNe in z-interval when drawing SNe, discard all SNe
  ! outside of z-interval, and re-draw until SNe falls within domain.
  ! -AT, 2019 May 03
  if (sn_stratifySN) then
    if (sn_kmin .lt. 0 .and. sn_kmax .gt. 0) then  ! nominal case,
      rateadj1 = 0.5*(1-exp(-1*abs(sn_kmin)/sn_hstar1)) + 0.5*(1-exp(-1*abs(sn_kmax)/sn_hstar1))
      rateadj2 = 0.5*(1-exp(-1*abs(sn_kmin)/sn_hstar2)) + 0.5*(1-exp(-1*abs(sn_kmax)/sn_hstar2))
    else if (sn_kmin .ge. 0) then
      rateadj1 = 0.5*(exp(-1*abs(sn_kmin)/sn_hstar1) - exp(-1*abs(sn_kmax)/sn_hstar1))
      rateadj2 = 0.5*(exp(-1*abs(sn_kmin)/sn_hstar2) - exp(-1*abs(sn_kmax)/sn_hstar2))
      if(sn_meshMe == MASTER_PE) then
        ! not tested as of 2019 May 07 - AT
        print*,'[SN_init] WARNING: SN rate adj untested for zmin >= 0'
      endif
    else if (sn_kmax .le. 0) then
      rateadj1 = 0.5*(exp(-1*abs(sn_kmax)/sn_hstar1) - exp(-1*abs(sn_kmin)/sn_hstar1))
      rateadj2 = 0.5*(exp(-1*abs(sn_kmax)/sn_hstar2) - exp(-1*abs(sn_kmin)/sn_hstar2))
      if(sn_meshMe == MASTER_PE) then
        ! not tested as of 2019 May 07 - AT
        print*,'[SN_init] WARNING: SN rate adj untested for zmax <= 0'
      endif
    end if
    if(sn_meshMe == MASTER_PE) then 
      print*,'[SN_init] adjusting stratified SN rates for finite z-domain'
      print*,'[SN_init] decreasing Ia rate by factor', rateadj1
      print*,'[SN_init] decreasing CC rate by factor', rateadj2
    endif
    sn_tsn1 = sn_tsn1 / rateadj1
    sn_tsn2 = sn_tsn2 / rateadj2
  endif

  !====================================================
  ! set directory for RNG saves
  !====================================================
  ! code adapted from io_getOutputName.F90 - AT, 2019 May 07
  ! if sn_outputDir isn't empty, and it doesn't have a
  ! directory separator at the end, add the directory seperator.
  pos = index(sn_outputDir, ' ')  ! strings are whitespace padded to MAX_STRING_LENGTH
  if (pos.gt.1) then
    if (sn_outputDir(pos-1:pos-1).ne."/") then
      sn_outputDir(pos:pos) = "/"
    end if
  end if
  if (sn_outputDir == "./") then
    sn_outputDir = ""
  end if

  !====================================================
  ! SN numbering
  !====================================================

  sn_nSN = 0
  sn_newDt = 1e99

  if(dr_restart) then

    inquire(file=trim(sn_outfile), exist=exst)
    if (.not. exst) then
      call Driver_abortFlash('[SN_init] Error, cannot find: ' // trim(sn_outfile))
    endif

    ! find last SN number prior to restart checkpoint time
    open(sn_funit, file=trim(sn_outfile),status="old",action="read")
    read(sn_funit,*,IOSTAT=readstat)  ! skip header row
    do
      read(sn_funit,*,IOSTAT=readstat) c0,c1,c2,c3,c4,c5,c6,c7,c8,c9
      if (readstat .lt. 0) then
        exit
      endif

      if (c4 .le. dr_simTime) then
        sn_nSN = c0
      else
        exit
      endif

    enddo
    close(sn_funit)

    if (sn_meshMe .eq. MASTER_PE) then
      print*,'[SN_init] Restarting from SN number', sn_nSN
    endif

    ! TODO (enhancement): perform sanity checks based on flash.par SN rate & SN
    ! number, because SN detonation times are pretty regularly spaced
    ! -AT 2019 May 11

    ! Read in random seed
    write(convert, '(i4.4)') io_checkpointFileNumber - 1
    call mtgetf(trim(sn_outputDir) // 'RNG_state_SN_' // trim(convert),'u')

    ! Call the RNG if user wants a new state
    if(sn_callRNG .gt. 0) then
      if (sn_meshMe .eq. MASTER_PE) then
        print*,'Calling RNG at restart', sn_callRNG, 'times'
      endif
      do i = 0, sn_callRNG
        c0 = grndSN()
      enddo
    endif

  else

    if (sn_meshMe .eq. MASTER_PE) then

      inquire(file=trim(sn_outfile), exist=exst)
      if (exst) then
        print*,'[SN_init] WARNING, clobbering existing SN data: ', trim(sn_outfile)
      endif

      open(sn_funit, file=trim(sn_outfile))
      write(sn_funit, '(10(1X,A16))') '[00]n_SN', '[01]type', &
        '[02]n_timestep', '[03]dud', '[04]time', '[05]posx', '[06]posy', &
        '[07]posz', '[08]radius', '[09]mass'
      close(sn_funit)

    endif

  endif  

  return
end subroutine SN_init
