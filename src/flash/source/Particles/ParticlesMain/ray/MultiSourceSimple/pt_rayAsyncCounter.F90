!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014

!! Description:
!!   issues asynchronous all reduce statements to figure out 
!!   if rays are still active
!!   needs MPI-3, for asynchronous collective calls
!!   

#define NBC_OK 0
#define NBC_CONTINUE 3
#define NBC_REQUEST_NULL -1
#define FLASH_REQUEST_NULL NBC_REQUEST_NULL

module pt_rayAsyncCounter 
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getComm

  implicit none

#include "constants.h"
#include "Flash.h"
  include 'Flash_mpi.h'

! local to counter
  integer, dimension(3), save :: ph_localCountMsg
  integer, dimension(3), save :: ph_sharedCountMsg
  integer, save :: ph_counterRequest

! sends
  integer, save :: ph_localSends
! recvs
  integer, save :: ph_localRecvs

! dublicate in pt_rayAsyncComm
  integer, save :: ph_meshComm
  integer, save :: ph_rank
  integer, save :: ph_size

  logical, save :: ph_isCounterDone      = .true.
  logical, save :: ph_isCounterTargetMet = .false.
  logical, save :: ph_testIfTargetMet = .false.

contains

!============================================
!==== main counter routiness
!============================================
  subroutine ph_initCounter()
    implicit none
  end subroutine ph_initCounter
!############################################
  subroutine ph_startCounter()!counterTarget)
! direct access is not nice from a code structure point of view but meh
    use Particles_rayData, only : ph_localRays
    implicit none
    integer :: ierr

! get comm info
    call Driver_getComm(  MESH_COMM, ph_meshComm)

! get rank
    call MPI_Comm_rank(ph_meshComm, ph_rank, ierr)

! get size
    call MPI_Comm_size(ph_meshComm, ph_size, ierr)

    ph_counterRequest = FLASH_REQUEST_NULL

! post first iallreduce
    if (ph_size > 1) then

!No outstanding MPI_Iallreduce so we post an MPI_Iallreduce.
      ph_sharedCountMsg = 0
! check send buffers and local rays for rays, this should be robust
! if receives append directly to local ray arrrray
      if(ph_localRays .gt. 0) then
        ph_localCountMsg(1) = 1 
      else
        ph_localCountMsg(1) = 0
      endif
      ph_localSends = 0
      ph_localRecvs = 0

      ph_localCountMsg(2) = ph_localSends
      ph_localCountMsg(3) = ph_localRecvs

#ifdef FLASH_LIBNBC
         call NBC_Iallreduce(ph_localCountMsg, ph_sharedCountMsg, &
              3, MPI_INT, MPI_SUM, ph_meshComm, ph_counterRequest, ierr)

         if (ierr /= NBC_OK) call Driver_abortFlash("NBC_Iallreduce error")
#else
         call MPI_Iallreduce(ph_localCountMsg, ph_sharedCountMsg, &
              3, FLASH_INTEGER, MPI_SUM, ph_meshComm, ph_counterRequest, ierr)
#endif
    endif

    ph_isCounterTargetMet = .false.
    ph_isCounterDone = .false.
  end subroutine ph_startCounter
!############################################
 subroutine ph_progressCounter(isTargetMet, commInProgress)
    use Particles_rayData, only : ph_localRays
    implicit none
    logical, intent(OUT) :: isTargetMet
    logical, intent(IN)  :: commInProgress
    integer :: ierr
    logical :: testIfTargetMet

    if (.not. ph_isCounterDone) then
       ph_isCounterTargetMet = .false.

       if (ph_size > 1) then
         testIfTargetMet = .false.
         if (ph_counterRequest .ne. FLASH_REQUEST_NULL .and. .not. ph_testIfTargetMet) then

#ifdef FLASH_LIBNBC
           call NBC_Test(ph_counterRequest, ierr)

           if (ierr /= NBC_OK .and. ierr /= NBC_CONTINUE) call Driver_abortFlash("NBC_Test error")           

           ph_testIfTargetMet = (ph_counterRequest == FLASH_REQUEST_NULL)
#else
           call MPI_Test(ph_counterRequest, ph_testIfTargetMet, MPI_STATUS_IGNORE,ierr)
#endif

         endif
       else
          ph_localCountMsg(:)  = ph_localRays
          ph_sharedCountMsg(:) = ph_localRays
          ph_testIfTargetMet = .true.
       end if

! did the allreduces complete?
       if (ph_testIfTargetMet) then
! reset
         ph_testIfTargetMet = .false.

         ph_isCounterTargetMet = (  (ph_sharedCountMsg(1) .eq. 0) .and. (ph_sharedCountMsg(2) .eq. ph_sharedCountMsg(3) ))

         if (ph_isCounterTargetMet) then
! returns asynchronous
           call ph_stopCounter()
         else
! immediatly post new receive
! this is done to not allow any gap in captured evolution of the ray number state
! 's a bit spammy
           ph_sharedCountMsg = 0

           if(ph_localRays .gt. 0 .or. commInProgress) then
             ph_localCountMsg(1) = 1 
           else
             ph_localCountMsg(1) = 0
           endif

           ph_localCountMsg(2) = ph_localSends
           ph_localCountMsg(3) = ph_localRecvs

           if (ph_size > 1) then
! LOCAL
#ifdef FLASH_LIBNBC
             call NBC_Iallreduce(ph_localCountMsg, ph_sharedCountMsg, &
                3, MPI_INT, MPI_SUM, ph_meshComm, ph_counterRequest, ierr)
             if (ierr /= NBC_OK) call Driver_abortFlash("NBC_Iallreduce error")
#else
             call MPI_Iallreduce(ph_localCountMsg, ph_sharedCountMsg, &
                  3, FLASH_INTEGER, MPI_SUM, ph_meshComm, ph_counterRequest, ierr)
#endif
           endif
        endif
      endif
    end if

    isTargetMet = ph_isCounterTargetMet
  end subroutine ph_progressCounter
!############################################
  subroutine ph_stopCounter()
    use Particles_rayData, only : ph_localRays
    implicit none
    integer :: i, ierr
    ph_isCounterDone = .true.
  end subroutine ph_stopCounter
!############################################
  subroutine ph_finalizeCounter()
    implicit none
  end subroutine ph_finalizeCounter
!############################################
end module pt_rayAsyncCounter
