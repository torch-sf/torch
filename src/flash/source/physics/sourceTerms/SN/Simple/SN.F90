!!****f* source/physics/sourceTerms/SN/Simple/SN
!!
!! NAME
!!  
!!  SN
!!
!! SYNOPSIS
!! 
!!  call SN (integer(IN) :: blockCount,
!!           integer(IN) :: blockList(blockCount),
!!           real(IN)    :: dt,
!!           real(IN)    :: time)
!!  
!! DESCRIPTION
!!
!! original work by MKRJ from 2004
!!
!! - CB April-May 2012 based on MKRJ original SN feedback code
!! added ridiculous amounts of comments 
!! kept most of MKRJ's comments
!! cleaned code up a bit but still needs polish
!! explosion times calculation more clear
!! merged Heat_block into Heat
!! rewrote routine so heating is always applied
!! exp_heat switch to reduce timestep was removed and replaced with more adaptive method
!! related to above, it should not be needed to know there will be a SN in the next timestep to control it
!! instead of counting number of zones and averaging by number, volumes are calculated as the zones might differ in size
!! restructured code so that radiation cooling of SN energy will be applied one timestep after SN explosion
!! ASSUMES CUBIC ZONES INSIDE BLOCK
!! TODO better 2d/3d switch
!!! TODO every process calculates the position of every global random explosion, this requires
!!! all random number generators to run with the same results, which is ok but not very safe
!!! a better approach would be to let every local domain compute its local random SN field
!!! only the SB would need global coordination, although this could also be rectified by pre-
!!! calculating SB locally and just working through them.
!! TODO FINAL CLEANUP
!!
!! AT 2019 May 02 cleanup -- strip to bare-bones.
!! remove tracers, SB, all the fancy stuff.
!! just random thermal SN driving
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine SN (blockCount,blockList,dt,time)

  use SN_data
  use mtmodSN, ONLY : mtsavef

  use IO_data, ONLY : io_justCheckpointed, io_checkpointFileNumber

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use GridInject_interface, ONLY : GridInject_thermalSN, GridInject_kineticSN

  ! just for output
  use Driver_data, ONLY : dr_nStep
  use tree, ONLY : lrefine_max
  use Grid_data,    ONLY : gr_delta
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Flash_mpi.h"

  integer, intent(IN)                         :: blockCount
  integer, dimension(blockCount), intent(IN)  :: blockList
  real, intent(IN)                            :: dt, time

  ! communication and mass calculation
  integer ::  ierr

  ! iterators
  integer :: i, j, k, l, m

  ! SN module scratch variables, only local
  integer :: nsn1, nsn2, nsn2dt, nsn1dt, num_exp
  real    :: t0heat, theat1, theat2, tnext1, tnext2, tstar1, tstar2

  integer, dimension(sn_nsndt)  :: SNtype
  real, dimension(sn_nsndt)     :: x0h, y0h, z0h, r0h

  ! variables for the SN explosion
  integer :: dud
  real  :: r_exp, mass
  real  :: xx, yy, zz
  real  :: x0heat, y0heat, z0heat
  ! energy deposition
  real :: rho_avg

  ! spacing
  integer :: zoneindex
  real, dimension(NDIM) :: gridspace

  character (len=MAX_STRING_LENGTH) :: convert

  if (.not. sn_useSN) return

  if(io_justCheckpointed) then
    ! save random number state, independent of flash4 random number stream
    ! use only main thread for this, as all other should have the same state
    if(sn_meshMe == MASTER_PE) then 
      write(convert, '(i4.4)') io_checkpointFileNumber - 1
      call mtsavef(trim(sn_outputDir) // 'RNG_state_SN_' // trim(convert),'u')
    endif
  endif

  call Timers_start("SN")

  !====================================================
  ! SN explosion times for given rates
  !====================================================

  t0heat = time - dt ! extend time interval back one step to check for overlap

  nsn1   = time/sn_tsn1 ! (int) expected total number of type I SN since start of simulation
  tstar1 = nsn1*sn_tsn1 ! current explosion time for type I SN
  tnext1 = tstar1 + sn_tsn1  ! next explosion time

  nsn2   = time/sn_tsn2 ! (int) expected total number of type II SN since start of simulation
  tstar2 = nsn2*sn_tsn2 ! current explosion time for type II SN
  tnext2 = tstar2 + sn_tsn2  ! next explosion time

  nsn2dt = dt/sn_tsn2  ! (int) 'nsn2dt' is the # of isolated Type II expl'ns in a single 'dt'
  nsn1dt = dt/sn_tsn1  ! expect to be zero for both, very often...

  if (t0heat > sn_tstop) then
    sn_newDt = 1e99  ! no more SN timestep adjustment needed
    call Timers_stop("SN")
    return
  endif

  !====================================================
  ! calculate new step to land before explosion, with some slop
  ! for SN_computeDt()
  !====================================================

  sn_newDt = min( 0.5*(tnext1 - 0.9*sn_SNminstep - time), &
                  0.5*(tnext2 - 0.9*sn_SNminstep - time) )
  sn_newDt = max(sn_newDt, sn_SNminstep)

  !====================================================
  ! any explosions in this time step?
  !====================================================

  num_exp = 0

  ! time for a field Type I supernova?
  if (t0heat < tstar1) then
    do i = num_exp + 1, num_exp + 1  + nsn1dt 

      SNType(i) = 1
      sn_nSN    = sn_nSN + 1  ! must update before sn_nextFieldPosition(...)
      num_exp   = num_exp + 1

      call sn_nextFieldPosition(SNType(i), xx, yy, zz)
      x0h(i) = xx
      y0h(i) = yy
      z0h(i) = zz
    enddo
  endif

  ! time for a field Type II supernova?
  if (t0heat < tstar2) then
    do i = num_exp+1, num_exp + nsn2dt +1

      SNType(i) = 2
      sn_nSN    = sn_nSN + 1  ! must update before sn_nextFieldPosition(...)
      num_exp   = num_exp + 1

      call sn_nextFieldPosition(SNType(i), xx, yy, zz)
      x0h(i) = xx
      y0h(i) = yy
      z0h(i) = zz
    enddo
  endif

  !====================================================
  ! execute SN, this involves MPI communication
  !====================================================

  !  loop over all SN
  do l = 1, num_exp

    x0heat = x0h(l)
    y0heat = y0h(l)
    z0heat = z0h(l)

    ! to ensure a nice bubble, this maps to the nearest highest refinement position good for uniform boxes, 
    ! not so great if SN goes off in lower refined region
    if(sn_SNmapToGrid) then 

      ! get grid spacing, assumes cubic box and 8^3 zones per block 
      !gridspace = (sn_imax-sn_imin)/(2**(lrefine_max+2))
      ! this does not work for MPI runs
      !call gr_findBlock(blockList,blockCount,(/x0heat,y0heat,z0heat/),blockID)

      ! assume it goes off in highest refinement region, as it is not local (would need mpi to find out)
      gridspace = gr_delta(1:MDIM,lrefine_max)

      zoneindex = floor((x0heat) / gridspace(1))
      x0heat = gridspace(1)*(zoneindex + 0.5d0)
      zoneindex = floor((y0heat) / gridspace(2))
      y0heat = gridspace(2)*(zoneindex + 0.5d0)
      zoneindex = floor((z0heat) / gridspace(3))
      z0heat = gridspace(3)*(zoneindex + 0.5d0)
    endif

    if (sn_kinetic) then
      call GridInject_kineticSN(x0heat, y0heat, z0heat, sn_exp_energy, &
                                sn_mejc, .false.)
      r_exp = 0 ! dummy placeholders
      mass = 0 ! dummy placeholders
      rho_avg = 0 ! dummy placeholders
    else
      call GridInject_thermalSN(x0heat, y0heat, z0heat, sn_exp_energy, &
              sn_mejc, sn_r_init, sn_r_exp_max, sn_nms, r_exp, mass, rho_avg)
    endif

    if(sn_meshMe == MASTER_PE) then 
      print*, '##SN of type',SNtype(l),' at position'
      print*, x0heat,y0heat,z0heat
      print*,'##average density:',rho_avg
    endif

    if(rho_avg .ne. rho_avg) then
      print*, 'NaN?', rho_avg, sn_meshMe
      call Driver_abortFlash("[SN.F90] average density is NaN")
    endif

    dud = 0
    if (mass .lt. sn_mejc) then
      dud = 1
      if(sn_meshMe == MASTER_PE) then
        print*, '##did not find enough mass:', mass/sn_mejc
      endif
    end if

    ! TODO do we need to keep this warning or restore max_temp cap?
    !! this should only happen if we can't find enough mass.
    !! but might be affected by discretization noise, too.
    !if (ei*convf .ge. sn_max_temp) then
    !  print*, '##temp. too high in SN, changing internal energy input'
    !  print*, '##factor of ',(ei*convf)/sn_max_temp,'too high'
    !  print*, '##lost energy',qheat*dt*dy(j)*dx(i)*dz(k)/sn_exp_energy ,'%'
    !  print*, '##core ID, block ID, zone ID xyz',sn_meshMe,blockID,i,j,k
    !  ei = sn_max_temp/convf
    !end if

    ! could have multiple calls fix
    call WriteSNFeedback(sn_nSN-num_exp+l, SNType(l), dr_nstep, time, x0heat, y0heat, z0heat, r_exp, mass, dud)

  enddo ! SN loop

  !====================================================
  ! finalize
  !====================================================

  call MPI_Barrier (MPI_Comm_World, ierr)
  call Timers_stop("SN")
  return

contains

  !====================================================
  ! output routine for supernovae, writes to SNfeedback.dat
  !====================================================

  subroutine WriteSNFeedback(nSN,type,ndt,time,x,y,z,radius,mass,dud)

    use SN_data, ONLY :  sn_meshMe, sn_funit, sn_outfile
    !use Driver_data, ONLY : dr_globalMe

    implicit none

#include "constants.h"

    real, intent(IN)   :: time, x,y,z,radius,mass
    integer, intent(IN):: nSN,type,ndt,dud ! 1 is type I, 2 is type II
    integer            :: i

    ! as all processors calculate the SN in lockstep, only the Master process has to output data.
    if (sn_meshMe .ne. MASTER_PE) return

    open(sn_funit, file=trim(sn_outfile), position='APPEND')
    write(sn_funit,'(4(1X,I16),6(1X,ES16.9))') &
       nSN,     &
       type,    &
       ndt,     &
       dud,     &
       time,    &
       x,       &
       y,       &
       z,       &
       radius,  &
       mass
    close(sn_funit)

    return
  end subroutine WriteSNFeedback

end subroutine SN
