!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkMergingAfterCreation
!!
!! NAME
!!
!!  pt_sinkMergingAfterCreation
!!
!! SYNOPSIS
!!
!!  call pt_sinkMergingAfterCreation(real, intent(IN) :: delta_at_lrefmax)
!!
!! DESCRIPTION
!!
!!  If the FLASH setup exhibits symmetries (e.g., a spherical collapse or a disk, etc.),
!!  then the sink particle must be created exactly in the point of symmetry or on the
!!  symmetry plane or axis. This routine achieves that, by merging the multiple
!!  zero-mass sinks created in such a symmetric case into one particle. If such symmetries
!!  occur in a setup, then Particles_sinkCreateAccrete() allows for the creation of
!!  multiple zero-mass sinks that are merged into the symmetry center by this routine.
!!
!! ARGUMENTS
!!
!!   delta_at_lrefmax - cell size at the maximum level of AMR
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2014
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine pt_sinkMergingAfterCreation(delta_at_lrefmax)

  use Particles_sinkData
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal, pt_sinkCorrectForPeriodicBCs
  use Driver_data, ONLY : dr_globalMe
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "GridParticles.h"
#include "Particles.h"

  include "Flash_mpi.h"
  
  real, intent(IN)    :: delta_at_lrefmax

  integer, parameter  :: funit = 42
  logical, save       :: first_call = .true.
  integer, save       :: MyPE, MasterPE

  real                :: pdist, search_radius, merged_xpos, merged_ypos, merged_zpos, dx, dy, dz
  integer             :: i, j, np, ip1, ip2, ip2index, n_merged, ierr
  integer             :: n_delete, n_delete_red, n_keep, n_keep_red
  logical             :: merged_already, merged_with_itself, sink_merging_happened, sink_merging_happened_red

  integer, parameter  :: maxpart = 200
  integer, dimension(maxpart) :: pindex, plistcount, merged_list
  integer, dimension(maxpart) :: delete_id, keep_id, delete_id_red, keep_id_red
  real, dimension(maxpart) :: keep_posx, keep_posy, keep_posz, keep_posx_red, keep_posy_red, keep_posz_red
  integer, dimension(maxpart, maxpart) :: plist

  character(len=80), save :: grav_boundary_type

  logical, parameter  :: debug = .false.

  if (first_call) then

     MyPE = dr_globalMe
     MasterPE = MASTER_PE
     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
     first_call = .false.

  end if

  if (debug) print *, 'SinkMergingAfterCreation entering.'

  sink_merging_happened = .false.

  n_delete = 0
  n_keep   = 0
  delete_id (:) = 0
  keep_id   (:) = 0
  keep_posx (:) = 0.0
  keep_posy (:) = 0.0
  keep_posz (:) = 0.0

  ! only the master proc takes care of merging the particles
  ! (only work on the global particles list)
  if (MyPE == MasterPE) then

     ! find all particles that were just created (zero mass).
     np = 0
     do i = 1, localnpf
        if (particles_global(ipm,i) .le. 0) then
           !print*, "Reported mass is =", particles_global(ipm,i)
           ! only merge particles with zero mass - this guarantees they are fresh out of creation
           np = np + 1
           if (np .gt. maxpart) call Driver_abortFlash("pt_sinkMergingAfterCreation: Error, list too small, increase size")
           pindex(np) = i
        end if
     end do
     ! np is the total number of mass-less particles

     ! for each particle: compute the distance to all other particles (including itself)
     ! for each particle fill a list of particle indexes that are inside the search radius;
     ! the search radius is set to 1.1*sqrt(3)*cellsize
     search_radius = 1.1*sqrt(3.0)*delta_at_lrefmax

     do ip1 = 1, np
        ip2index = 0
        do ip2 = 1, np
           dx = particles_global(ipx,pindex(ip2)) - particles_global(ipx,pindex(ip1))
           dy = particles_global(ipy,pindex(ip2)) - particles_global(ipy,pindex(ip1))
           dz = particles_global(ipz,pindex(ip2)) - particles_global(ipz,pindex(ip1))
           if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)
           pdist = sqrt(dx**2 + dy**2 + dz**2)
           if (pdist .le. search_radius) then
              ip2index = ip2index + 1
              plist(ip1,ip2index) = ip2
              plistcount(ip1) = ip2index
           end if

        end do
     end do

! Note: This is all well and good for even numbered symmetric situations, but in
! a turbulent gas you can get formation of sinks that are very close and odd
! numbered like:
!   1    2    3
!   o    o    o 
! Where #2 has two neighbors that are within 1.1*dx while #1 and #3 only have
! one each. Then this check fails even though I'd like these to merge. So I'm
! basically saying who needs sanity? - JW

#ifdef need_sanity
     ! Sanity check...
     do ip1 = 1, np
        do ip2 = 1, plistcount(ip1)
           if (plistcount(ip1) .NE. plistcount(plist(ip1,ip2))) then
              call Driver_abortFlash('pt_sinkMergingAfterCreation: Error, newly formed sink particles &
    & have unexpected relative positions!')
           endif
        end do
     end do
#endif

     ! merge zero mass particles that are inside of the search radius of one another
     n_merged = 0
     merged_list(:) = 0
     do ip1 = 1, np
        merged_already = .false.
        merged_with_itself = .true.
        do i = 1, n_merged
           if (ip1 .eq. merged_list(i)) merged_already = .true.
        end do
        if (.not. merged_already) then
           merged_xpos = 0.0
           merged_ypos = 0.0
           merged_zpos = 0.0
           do ip2 = 1, plistcount(ip1)
              ! compute new position for the merged particle
              dx = particles_global(ipx,pindex(plist(ip1,ip2))) - particles_global(ipx,pindex(ip1))
              dy = particles_global(ipy,pindex(plist(ip1,ip2))) - particles_global(ipy,pindex(ip1))
              dz = particles_global(ipz,pindex(plist(ip1,ip2))) - particles_global(ipz,pindex(ip1))
              if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)
              merged_xpos = merged_xpos + particles_global(ipx,pindex(ip1))+dx
              merged_ypos = merged_ypos + particles_global(ipy,pindex(ip1))+dy
              merged_zpos = merged_zpos + particles_global(ipz,pindex(ip1))+dz
              n_merged = n_merged + 1
              merged_list(n_merged) = plist(ip1,ip2)
              ! delete merged particles; only keep ip1 (to which we merge)
              if (ip1 .ne. plist(ip1,ip2)) then
                 n_delete = n_delete+1
                 delete_id(n_delete) = int(particles_global(iptag,pindex(plist(ip1,ip2))))
                 write(*,'(A,I8,A,I8,A,3(3X,ES16.9))') 'SinkMergingAfterCreation: zero-mass sink merged TO tag: ', &
                    int(particles_global(iptag,pindex(ip1))), ' FROM tag: ', &
                    int(particles_global(iptag,pindex(plist(ip1,ip2)))), &
                    ' from pos: ', particles_global(ipx,pindex(plist(ip1,ip2))), &
                                   particles_global(ipy,pindex(plist(ip1,ip2))), &
                                   particles_global(ipz,pindex(plist(ip1,ip2)))
                 sink_merging_happened = .true.
                 merged_with_itself = .false.
              end if
           end do !ip2

           ! move ip1 to the new position (remember that we just merge zero mass particles in
           ! this incredibly complicated subroutine, so that only the position has to be changed
           ! to conserve all relevant physical quantities)

           if (.not. merged_with_itself) then

              n_keep = n_keep + 1
              keep_id(n_keep) = int(particles_global(iptag,pindex(ip1)))

              keep_posx(n_keep) = merged_xpos/real(plistcount(ip1))
              keep_posy(n_keep) = merged_ypos/real(plistcount(ip1))
              keep_posz(n_keep) = merged_zpos/real(plistcount(ip1))

              particles_global(ipx,pindex(ip1)) = merged_xpos/real(plistcount(ip1))
              particles_global(ipy,pindex(ip1)) = merged_ypos/real(plistcount(ip1))
              particles_global(ipz,pindex(ip1)) = merged_zpos/real(plistcount(ip1))

              write(*,'(A,I8,A,3(3X,ES16.9))') 'SinkMergingAfterCreation: new position of tag ', &
                int(particles_global(iptag,pindex(ip1))), ': ', particles_global(ipx,pindex(ip1)), &
                                                                particles_global(ipy,pindex(ip1)), &
                                                                particles_global(ipz,pindex(ip1))
              write(*,'(A)') &
                'SinkMergingAfterCreation: Its first accretion step moves this sink to the actual center of mass.'

           end if ! not merged_with_itself

        end if ! not merged_already 

     end do ! ip1

     if (np .ne. n_merged) call Driver_abortFlash("pt_sinkMergingAfterCreation: Error, something went wrong")

  endif ! MasterPE

  ! Spread the word that merging happened
  sink_merging_happened_red = .false.
  call MPI_ALLREDUCE(sink_merging_happened, sink_merging_happened_red, 1, FLASH_LOGICAL, & 
       FLASH_LOR, MPI_COMM_WORLD, ierr)
  sink_merging_happened = sink_merging_happened_red

  ! all processors do this if sink_merging_happened = .true.
  ! modify the local particle list (delete, keep)
  ! Note that this is only done if merging after creation is neccessary,
  ! which basically never happens or only a few times in a real application,
  ! so we don't worry about the massage size (maxpart) for the MPI_ALLREDUCE here.
  if (sink_merging_happened) then

     ! reduce the indices and ids to be deleted and kept
     n_delete_red  = 0
     n_keep_red    = 0
     delete_id_red (:) = 0
     keep_id_red   (:) = 0
     keep_posx_red (:) = 0.
     keep_posy_red (:) = 0.
     keep_posz_red (:) = 0.
     call MPI_ALLREDUCE(n_delete , n_delete_red,   1,       FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(n_keep   , n_keep_red,     1,       FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(delete_id, delete_id_red,  maxpart, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(keep_id  , keep_id_red,    maxpart, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(keep_posx, keep_posx_red,  maxpart, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(keep_posy, keep_posy_red,  maxpart, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(keep_posz, keep_posz_red,  maxpart, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     n_delete  = n_delete_red
     n_keep    = n_keep_red
     delete_id = delete_id_red
     keep_id   = keep_id_red
     keep_posx = keep_posx_red
     keep_posy = keep_posy_red
     keep_posz = keep_posz_red

     ! loop over particles to keep
     do j = 1, n_keep
        do i = 1, localnp
           if (int(particles_local(iptag,i)) .eq. keep_id(j)) then
              particles_local(ipx,i) = keep_posx(j)
              particles_local(ipy,i) = keep_posy(j)
              particles_local(ipz,i) = keep_posz(j)
           end if
        end do ! i
     end do ! j

     ! loop over the particles to delete
     do j = 1, n_delete
        do i = 1, localnp
           if (int(particles_local(iptag,i)) .eq. delete_id(j)) then
              ! I can't find anywhere else in Flash this this used. So
              ! comment it out. -JW
              !NumParticlesPerBlock(int(particles_local(ipblk,i))) = &
              !    & NumParticlesPerBlock(int(particles_local(ipblk,i))) - 1
              particles_local(:,i) = particles_local(:,localnp) ! copy the last in list here
              particles_local(ipblk,i) = NONEXISTENT ! delete the last in list
              n_empty = n_empty + 1 ! delete
              localnp = localnp - 1 ! reduce local number of particles by one
              exit ! found the one to delete, so go to next j
           end if
        end do ! i
                            
        !!! Additions for capturing new particles from Flash to
        !!! move to AMUSE. Note we have to do this after merging
        !!! of new particles, or we get too many sinks. - JW
        do i = 1, number_new_sinks
           if (int(new_sink_tags(i)) .eq. delete_id(j)) then
              new_sink_tags(i) = new_sink_tags(number_new_sinks)
              number_new_sinks = number_new_sinks - 1
              exit
           end if
        end do
     end do ! j

     call pt_sinkGatherGlobal()

  end if ! sink_merging_happened

  return

end subroutine pt_sinkMergingAfterCreation
