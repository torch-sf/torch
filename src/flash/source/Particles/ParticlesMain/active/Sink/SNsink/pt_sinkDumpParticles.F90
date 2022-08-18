!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkDumpParticles
!!
!! NAME
!!
!!  pt_sinkDumpParticles
!!
!! SYNOPSIS
!!
!!  call pt_sinkDumpParticles(real, intent(IN) :: simtime
!!                            character(len=MAX_STRING_LENGTH), intent(IN) :: outfilename)
!!
!! DESCRIPTION
!!
!!  Dumps a file called outfilename, containing the time evolution of the most 
!!  important sink particle properties (tag, positions, velocities, accelerations
!!  spin, mass, accretion rate, etc.).
!!
!! ARGUMENTS
!!
!!   simtime - the current simulation time
!!
!!   outfilename - the file name of the dump file
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2014
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!
!!	 CB added SNsink Module to FLASH 3/4
!!	 JC ported and modified the SNsink Module to FLASH4.2.2 2015
!!
!!***

subroutine pt_sinkDumpParticles(simtime, outfilename)

  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Particles_data, ONLY : pt_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  real, intent(IN) :: simtime
  character(len=*), intent(IN) :: outfilename

  integer, parameter :: nfilesmax = 2
  character(len=MAX_STRING_LENGTH), dimension(nfilesmax), save :: outfilename_stored
  integer, save :: n_known_files
  logical :: header_written

  integer, parameter :: funit_evol = 15
  integer            :: i
  integer, save      :: MyPE, MasterPE, ipax, ipay, ipaz
  logical, save      :: firstCall = .TRUE.

  ! JCI -three new properties added for the SNsink particles (NSN, TSN, TCRT)
  integer, parameter :: gather_nprops = 19
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: TAG_PART_PROP, MASS_PART_PROP, ACCR_RATE_PART_PROP, &
                  POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, &
                  VELX_PART_PROP, VELY_PART_PROP, VELZ_PART_PROP, &
                  ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP, &
                  X_ANG_PART_PROP, Y_ANG_PART_PROP, Z_ANG_PART_PROP, &
                  CREATION_TIME_PART_PROP, NSN_PART_PROP, TSN_PART_PROP, &
 				  TCRT_PART_PROP /)

  if (firstCall) then

     MyPE     = pt_globalMe
     MasterPE = MASTER_PE
     ipax     = ACCX_PART_PROP
     ipay     = ACCY_PART_PROP
     ipaz     = ACCZ_PART_PROP

     n_known_files = 0
     outfilename_stored(:) = ""

     firstCall = .false.

  endif

  ! exchange particle information across CPUs (this needs to be called by all processors)
  call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

  ! only the master processor dumps the data of all particles (global list perticlest)
  if (MyPE .NE. MasterPE) return

  ! open file for write
  open(funit_evol, file=trim(outfilename), position='APPEND')

  ! check whether we have to write a header for this file
  header_written = .false.
  do i = 1, nfilesmax
     if (trim(outfilename) .eq. trim(outfilename_stored(i))) then
        header_written = .true.
        exit ! the loop
     endif
  enddo

  ! write header if not yet written to the passed file
  if (.not. header_written) then

     ! print *, 'pt_sinkDumpParticles: writing header to >'//trim(outfilename)//'<'

     open(funit_evol, file=trim(outfilename), position='APPEND')
     write(funit_evol,'(20(1X,A16))') '[00]part_tag', '[01]time', '[02]posx', '[03]posy', '[04]posz', &
                                      '[05]velx', '[06]vely', '[07]velz', '[08]accelx', '[09]accely', &
                                      '[10]accelz', '[11]anglx', '[12]angly', '[13]anglz', '[14]mass', &
                                      '[15]mdot', '[16]ptime', '[17]nSN', '[18]tSN', '[19]Creationt'

     ! and mark the file name as known
     n_known_files = n_known_files + 1
     outfilename_stored(n_known_files) = outfilename

  endif ! write a header

  ! now write the actual data to the file
  do i = 1, localnpf

     write(funit_evol,'(1(1X,I16),21(1X,ES16.9))') &
!     write(funit_evol,'(1(1X,I16),18(1X,ES16.9))') &
          int(particles_global(iptag,i)), &
          simtime, &
          particles_global(ipx,i), &
          particles_global(ipy,i), &
          particles_global(ipz,i), &
          particles_global(ipvx,i), &
          particles_global(ipvy,i), &
          particles_global(ipvz,i), &
          particles_global(ipax,i), &
          particles_global(ipay,i), &
          particles_global(ipaz,i), &
          particles_global(iplx,i), &
          particles_global(iply,i), &
          particles_global(iplz,i), &
          particles_global(ipm,i), &
          particles_global(ipmdot,i), &
          particles_global(ipt,i), &
  		  particles_global(NSN_PART_PROP,i), &
          particles_global(TSN_PART_PROP,i), &
          particles_global(TCRT_PART_PROP,i)

  enddo

  close(funit_evol)

  return

end subroutine pt_sinkDumpParticles
