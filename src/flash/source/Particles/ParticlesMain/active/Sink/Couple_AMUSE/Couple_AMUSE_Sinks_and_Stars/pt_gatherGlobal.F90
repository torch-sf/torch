!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkGatherGlobal
!!
!! NAME
!!
!!  pt_gatherGlobal
!!
!! SYNOPSIS
!!
!!  call pt_gatherGlobal(integer, dimension(:), intent(in), optional :: propinds,
!!                           integer, intent(in), optional :: nprops)
!!
!! DESCRIPTION
!!
!!  Fills the global particle list with a copy of all sink particles on all
!!  processors. Note that the global sink particle list is not ordered in the same
!!  way on all processors. It always contains the local sink particles in the first
!!  1:localnp indexes after the call is finished. If optional arguments are given,
!!  only specified particle properties are communicated across all processors.
!!
!! ARGUMENTS
!!
!!  propinds (optional) - integer list of particle properties to communicate
!!
!!  nprops (optional) - number of particle properties to communicate (length of propinds)
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2014
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   modified to use on regular particles by Josh Wall 2017
!!
!!***

subroutine pt_gatherGlobal(propinds, nprops)

   ! Get particle information from all processors
   ! i.e., update allproc_particles from particles

   use Particles_Data, ONLY : particles, allproc_particles, pt_numLocal, pt_numGlobal
   use Particles_sinkData, ONLY : recv_buff, send_buff, MAX_MSGS
                            !& pt_sinkParticleProps, sink_maxSinks
   use Particles_data, ONLY : pt_globalMe, pt_globalNumProcs, pt_globalComm
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none

   include "Flash_mpi.h"

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

   integer, dimension(:), intent(in), optional :: propinds
   integer, intent(in), optional :: nprops

   integer       :: pcount_recv, pcount_send, pmsgcount_recv, pmsgcount_send
   integer       :: localnpt, jproc, i, k, ierr, np_offset, reqr
   integer       :: statr(MPI_STATUS_SIZE), nprops_to_communicate, prop_index
   logical       :: send_receive
   integer, save :: MyPE, NumPEs
   logical, save :: first_call = .TRUE.

   integer, dimension(NPART_PROPS), save :: prop_indices_all
   integer, dimension(NPART_PROPS) :: prop_indices

   if (first_call) then

     MyPE = pt_globalMe
     NumPEs = pt_globalNumProcs

     ! make a save index array with all particle property indices
     do k = 1, NPART_PROPS
       prop_indices_all(k) = k
     enddo

     first_call = .false.

   end if

   ! see if we only have to communicate specified particle properties
   if (present(propinds) .and. present(nprops)) then
     nprops_to_communicate = nprops
     prop_indices(1:nprops) = propinds(1:nprops)
   else ! communicate them all
     nprops_to_communicate = NPART_PROPS
     prop_indices(1:NPART_PROPS) = prop_indices_all(1:NPART_PROPS)
   endif

   ! the global list always contains the local particles in the first indices (1:localnp)
   if (pt_numLocal .gt. 0) then
     do k = 1, nprops_to_communicate
       prop_index = prop_indices(k)
       allproc_particles(prop_index,1:pt_numLocal) = particles(prop_index,1:pt_numLocal)
     end do
   endif

   ! Get the total number of particles on each processor and the total number overall
   pt_numGlobal = 0

   ! Get the total number of sink particles on all processors
   call MPI_ALLREDUCE(pt_numLocal, pt_numGlobal, 1, FLASH_INTEGER, MPI_SUM, pt_globalComm, ierr)
   
!   if (localnpf .gt. sink_maxSinks) then
!     call Driver_abortFlash("pt_sinkGatherGlobal: Number of all sink particles exceeds sink_maxSinks. Increase.")
!   end if

   np_offset = pt_numLocal

   !-----------------------------------------------------------------------------
   ! loop over all of the processors.  All the data is moved to local processor 
   ! using MPI sends and receives.
   !-----------------------------------------------------------------------------

   do jproc = 0, NumPEs-1

     if (jproc .ne. MyPE) then

        call MPI_IRECV(localnpt,1,FLASH_INTEGER,jproc,889,pt_globalComm,reqr,ierr)
        call MPI_SSEND(pt_numLocal ,1,FLASH_INTEGER,jproc,889,pt_globalComm,ierr)

        call MPI_WAIT(reqr,statr,ierr)

        ! if there are no particles on this processor
        ! and jproc, do not go through with send/receive
        if ((localnpt .eq. 0) .and. (pt_numLocal .eq. 0)) send_receive = .false.

        ! Let's send data; one particle property at a time.
        do k = 1, nprops_to_communicate

           prop_index = prop_indices(k)

           pcount_recv = 0
           pcount_send = 0
           send_receive = .true.

           do while (send_receive)

              ! Do not receive more than 12 (=MAX_MSGS) at a time
              pmsgcount_recv = min(localnpt - pcount_recv, MAX_MSGS)
              if (pmsgcount_recv .gt. 0) then

                 call MPI_IRECV(recv_buff, pmsgcount_recv, FLASH_REAL, &
                      jproc, 4711, pt_globalComm, reqr, ierr)

              end if

              pmsgcount_send = min(pt_numLocal - pcount_send, MAX_MSGS)
              if (pmsgcount_send .gt. 0) then
                 do i = 1, pmsgcount_send
                    send_buff(i) = particles(prop_index, pcount_send+i)
                 end do

                 call MPI_SSEND(send_buff, pmsgcount_send, FLASH_REAL, &
                      jproc, 4711, pt_globalComm, ierr)

                 pcount_send = pcount_send + pmsgcount_send
              end if

              if (pmsgcount_recv .gt. 0) then

                 call MPI_WAIT(reqr, statr, ierr)
                 do i = 1, pmsgcount_recv
                    allproc_particles(prop_index, np_offset+pcount_recv+i) = recv_buff(i)
                 end do
                 pcount_recv = pcount_recv + pmsgcount_recv
              end if

              if ((pcount_recv .ge. localnpt) .and. (pcount_send .ge. pt_numLocal)) &
                   send_receive = .false.
           end do   ! while send_receive
        end do  ! particle properties

        np_offset = np_offset + localnpt

     end if
   end do

end subroutine pt_gatherGlobal
