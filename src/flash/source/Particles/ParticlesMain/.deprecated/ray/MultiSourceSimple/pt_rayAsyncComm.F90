!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014
!!
!!
!! Description:
!!   routines for asynchronous communication 
!!   3 stages:
!!   init: posts speculative recvs, inits buffers and MPI info
!!   loop: checks for rays locally, in send buffers
!!         posts sends and re-issues speculative recvs, saves rays locally after receive
!!   finish: once counter is fulfilled
!!           posts 0-byte sends to fulfill open receives
!!           deallocate all buffers
!!

! based on energy deposition asynchronous communication, and UTPipeline

! changes compared to energy deposition async comm.
! recv. rays are directly saved to ray array, no intermediate buffer
! this is done as we do not need a sorting step to prepare rays for 

!! TODO change lookups into processor list to something else than a loop
!!      for large number of procs this becomes costly

! total asso comm, ey

#include "constants.h"
#include "Flash.h"
!#define DEBUG
!#define PROC 0
module pt_rayAsyncComm
  use Driver_interface, ONLY : Driver_abortFlash
  use Particles_rayData, ONLY : ph_transProp

  implicit none

  include 'Flash_mpi.h'
! general comm variables
  integer, save :: ph_meshComm
  integer, save :: ph_numNeigh

! not used in async comm. module
  integer, save :: ph_meshMe
  integer, save :: ph_rank
  integer, save :: ph_size
  integer, save :: ph_globalNumProcs, ph_globalComm, ph_globalMe, ph_meshNumProcs

  integer, allocatable, save, dimension(:)   :: ph_neighProcList

! hardcoded for the moment
!	integer,parameter :: ph_transProp = 6

! so ph_raysToBundle*neighbours is total size
  integer,save :: ph_raysToBundle, ph_CommCheckInterval

! taggin da streets
  integer, parameter :: ph_tag = 1235
! 2D
  integer, allocatable, save, dimension(:,:) :: ph_recvStatus

  integer, allocatable, save, dimension(:)   :: ph_recvRequest
  integer, allocatable, save, dimension(:)   :: ph_recvIndex
  integer, allocatable, save, dimension(:)   :: ph_recvCount
! 2D 
  integer, allocatable, save, dimension(:,:) :: ph_sendStatus

  integer, allocatable, save, dimension(:)   :: ph_sendRequest
  integer, allocatable, save, dimension(:)   :: ph_sendIndex
  integer, allocatable, save, dimension(:)   :: ph_sendCount

  real, allocatable, save, dimension(:,:,:) :: ph_sendBuf
  real, allocatable, save, dimension(:,:,:) :: ph_recvBuf

  logical, save :: ph_isSendCommDone
  logical, save :: ph_isRecvCommDone
  logical, save :: ph_isCommDone = .false.

  integer, allocatable, save, dimension(:) :: ph_sendState

  integer, parameter :: OPEN_STATE = -3000
  integer, parameter :: PROMISE_TO_CLOSE_STATE = -4000
  integer, parameter :: WAITING_TO_CLOSE_STATE = -5000
  integer, parameter :: CLOSE_STATE = -6000

! needed as speculative recvs are posted to cancel those 
! and finish comm when all rays are done
  integer, save      :: ph_commState
  integer, parameter :: COMM_ACTIVE = -100
  integer, parameter :: COMM_CLOSING = -200
  integer, parameter :: COMM_INACTIVE = -300

  integer, save      :: ph_counterState
  integer, parameter :: COUNTER_ACTIVE = -400
  integer, parameter :: COUNTER_INACTIVE = -500

!#define DEBUG
contains 
!============================================
!==== main routiness
!============================================
  subroutine ph_progressComm(force, NoMoRays, doCounter)
    use pt_rayAsyncCounter, only : ph_progressCounter    
    use Particles_rayData, only : ph_localRays
    implicit none
    logical, intent(IN)  :: force, doCounter
    logical, intent(OUT) :: NoMoRays
    logical :: raysInComm
    logical :: isCommDone, isTargetMet
    integer :: i, ierr, opencount, tmp

! checks local rays and rays in sendBuffer
    if(doCounter) then

      if(ph_numNeigh .gt. 0) then
        raysInComm = .false.

        if( any(ph_sendCount(:) .gt. 0) ) then
          raysInComm = .true.
        endif

        call ph_progressCounter(isTargetMet, raysInComm)
      else
        call ph_progressCounter(isTargetMet, .false.)
      endif

! signal that counter is done and we should wrap up comm
      if(isTargetMet) then 

        ph_counterState = COUNTER_INACTIVE
      endif
    else
      isTargetMet = .false.
    endif

! update counter 
    select case (ph_commState)
      case (COMM_ACTIVE)

! these routines check for communication status and update status 
! according the communication state of each thread
        call ph_progressRecvComm()
        call ph_progressSendComm()

! send what we have without waiting for full buffer, buffers are checked in ph_sendRay
         if(force) then
! check if there is anything to send
           call ph_sendFullestChannel()
         endif

! check counter
        if(COUNTER_INACTIVE == ph_counterState) then


#ifdef DEBUG
if(ph_meshMe .eq. PROC) then
          if(ph_localRays .gt. 0) then
            print*,'rays left :(',ph_localRays,ph_sendCount,ph_recvCount,ph_meshMe
            call Driver_abortFlash("local rays left, abort")
          endif
endif
#endif

! check if communication step is done
! clean up requests
          call ph_commDone(isCommDone)

          if (isCommDone) ph_commState = COMM_CLOSING
        endif

      case (COMM_CLOSING)
! reuse logical
! check if global communication is done
        call ph_finishComm(isCommDone)

        if (isCommDone) then
          ph_commState = COMM_INACTIVE
          ph_isCommDone = .true.
        end if

      case (COMM_INACTIVE)
        continue

      end select

! this makes sure that finish comm code is called, canceling opened comm requests
! just using met counter is not enough
     NoMoRays = ( COUNTER_INACTIVE == ph_counterState .and. &
                     COMM_INACTIVE == ph_commState )
    
   return
  end subroutine ph_progressComm
!#############################################
  subroutine ph_initComm()
    use pt_rayAsyncCounter, only : ph_initCounter
    use Grid_interface, ONLY : Grid_getNeighProcList

    implicit none
    integer :: i
    integer, pointer :: tmpNeighProcList (:)
! neighbour processer list, should not change over simulation run time
! boolean flags if calling processor should be included
    call Grid_getNeighProcList (.false., tmpNeighProcList, ph_numNeigh)

    if (ph_numNeigh > 0) then
! associated checks pointer pair
      if (.not. associated (tmpNeighProcList)) then
        call Driver_abortFlash ("[pt_initRays] Neigh proc list broken?")
      end if

! allocate right sized array
      allocate (ph_neighProcList (1:ph_numNeigh))
      ph_neighProcList (:) = tmpNeighProcList (:)
    else
! no neighbours
      allocate (ph_neighProcList(1))
      ph_neighProcList (1) = -1
    end if

! destroy pointers and free associated memory
    if (associated (tmpNeighProcList)) deallocate (tmpNeighProcList)
    nullify (tmpNeighProcList)

    if (ph_numNeigh > 0) then
! not sure this is needed, buffer for all received rays
! might just drop in directly into local array

! communication status
      allocate(ph_sendStatus(MPI_STATUS_SIZE, ph_numNeigh))
      allocate(ph_recvStatus(MPI_STATUS_SIZE, ph_numNeigh))

! storage for actual communication, ray properties x number of rays to send/recv at once x 
! number of neighbours -> 3D array 
      allocate(ph_sendBuf(ph_transProp, ph_raysToBundle, ph_numNeigh))
      allocate(ph_recvBuf(ph_transProp, ph_raysToBundle, ph_numNeigh))

! communication status arrays to each neighbour
! send 
      allocate(ph_sendRequest(ph_numNeigh))
      allocate(ph_sendIndex(  ph_numNeigh))
      allocate(ph_sendCount(  ph_numNeigh))
      allocate(ph_sendState(  ph_numNeigh))

! receive
      allocate(ph_recvRequest(ph_numNeigh))
      allocate(ph_recvIndex(  ph_numNeigh))
      allocate(ph_recvCount(  ph_numNeigh))

! open to anything
      ph_sendRequest = MPI_REQUEST_NULL
      ph_recvRequest = MPI_REQUEST_NULL
    end if

    if (ph_numNeigh > 0) then
      ph_sendCount(:) = 0
      ph_recvCount(:) = 0

      ph_isCommDone     = .false.
      ph_isSendCommDone = .false.
      ph_isRecvCommDone = .false.

      do i = 1, ph_numNeigh
        ph_sendState(i) = OPEN_STATE
! send recv message, waiting to be filled in raytracing,
! every thread is listening

        call ph_postRecvMsg(i)
      end do
! don't do comm
    else
      ph_isCommDone     = .true.
      ph_isSendCommDone = .true.
      ph_isRecvCommDone = .true.
    end if

! turn on communication
    ph_counterState = COUNTER_ACTIVE
    ph_commState    = COMM_ACTIVE

! copies relevant comm info to counter module
! cause otherwise circular dependencies
    call ph_initCounter()

  end subroutine ph_initComm
!############################################
  subroutine ph_startComm()
    use pt_rayAsyncCounter, only : ph_startCounter
    implicit none
! get the count counting, determines when no rays are left
! one ray ah ah ah, two rays ah ah ah... (he'll be busy for a while)
    call ph_startCounter()
  end subroutine ph_startComm
!############################################
  subroutine ph_finishComm(isCommDone)
    implicit none
    integer :: i, ierr
    logical,intent(out) :: isCommDone

    if(ph_isSendCommDone .and. ph_isRecvCommDone) then

      if (ph_numNeigh > 0) then
! cancel speculative receives
         do i = 1, ph_numNeigh
            if (ph_recvRequest(i) /= MPI_REQUEST_NULL) then
               call MPI_Cancel(ph_recvRequest(i), ierr)
               call CheckMPIErrorCode(ierr)
            end if
         end do

         call MPI_Waitall(ph_numNeigh, ph_recvRequest, &
              ph_recvStatus, ierr)

         ph_isRecvCommDone = .true.

         do i = 1, ph_numNeigh
            if (ph_sendRequest(i) /= MPI_REQUEST_NULL) then
               call MPI_Cancel(ph_sendRequest(i), ierr)
               call CheckMPIErrorCode(ierr)
            end if
         end do

         call MPI_Waitall(ph_numNeigh, ph_sendRequest, &
              ph_sendStatus, ierr)
         call CheckMPIErrorCode(ierr)

         ph_sendState(:) = CLOSE_STATE
         ph_isSendCommDone = .true.
      end if
    else
      call ph_progressRecvComm()
      call ph_progressSendComm()
    endif

    isCommDone = .true.
    return
  end subroutine ph_finishComm
!############################################
  subroutine ph_finalizeComm()
    use pt_rayAsyncCounter, only : ph_finalizeCounter
    implicit none
    integer :: ierr, tmp, i
    logical :: raysInComm

    if (ph_numNeigh > 0) then

! deallocate stuff
! communication status
      deallocate(ph_sendStatus)
      deallocate(ph_recvStatus)

! storage for actual communication, ray properties x number of rays to send/recv at once x 
! number of neighbours -> 3D array 
      deallocate(ph_sendBuf)
      deallocate(ph_recvBuf)

! communication status arrays to each neighbour
! send
      deallocate(ph_sendRequest)
      deallocate(ph_sendIndex)
      deallocate(ph_sendCount)
      deallocate(ph_sendState)

! receive
      deallocate(ph_recvRequest)
      deallocate(ph_recvIndex)
      deallocate(ph_recvCount)

    endif

    deallocate(ph_neighproclist)

    
    call ph_finalizeCounter()

! return synchronized
    call MPI_Barrier(ph_meshComm, ierr)
  end subroutine ph_finalizeComm
!############################################
  !Caller should probably add the following: if (isHandled) item = NONEXISTENT
  subroutine ph_sendRay(ray, procID, isHandled)
    use pt_rayAsyncCounter, only : ph_localSends
    implicit none
    real, dimension(:), intent(IN) :: ray
    integer, intent(IN) :: procID

    logical, intent(OUT) :: isHandled
    integer :: channel, ptr, i

    !It may be necessary to change the utpipe_procList data structure
    !to make the lookup faster.
    channel = -1
! find target
    do i = 1, ph_numNeigh
       if (ph_neighProcList(i) == procID) then
          channel = i
          exit
       end if
    end do

    if (channel == -1) call Driver_abortFlash("Msg channel not found")

    !If there is a pending send in our desired channel we test all
    !send channels.  Request values are reset to MPI_REQUEST_NULL when
    !sends complete.
    if (ph_sendRequest(channel) /= MPI_REQUEST_NULL) then
       call ph_progressSendComm()
    end if

    !We can safetly add items to the send buffer if there is no pending send.
    if ( ph_sendState(channel) == OPEN_STATE .and. &
         ph_sendRequest(channel) == MPI_REQUEST_NULL ) then

       ptr = ph_sendCount(channel) + 1

       if (ptr > ph_raysToBundle) call Driver_abortFlash("Counting error")

       ph_sendBuf(:,ptr,channel) = ray(:)
       ph_sendCount(channel) = ptr !Array is needed in utpipe_postSendMsg
          
       if (ph_sendCount(channel) .eq. ph_raysToBundle) then

! enough for send
! add to counter
          ph_localSends = ph_localSends + ph_raysToBundle
          call ph_postSendMsg(channel)
! reset?

       end if

       isHandled = .true.
    else
       isHandled = .false.
    end if
  end subroutine ph_sendRay
!============================================
!==== helper routiness
!============================================
  subroutine ph_sendFullestChannel
    use pt_rayAsyncCounter, only : ph_localSends
    implicit none
    integer :: fullestChannel, bufSize, i

    if (ph_numNeigh > 0) then
       if (any(ph_sendCount(:) > 0)) then

          fullestChannel = -1
          bufSize = -1

! find fullest channel
          do i = 1, ph_numNeigh

             !Test for data that is not currently being sent
            if ( ph_sendState(i) == OPEN_STATE .and. &
                 ph_sendRequest(i) == MPI_REQUEST_NULL .and. &
                 ph_sendCount(i) > 0 ) then

! find highest amount of rays to send
              if (ph_sendCount(i) > bufSize) then
                fullestChannel = i
                bufSize = ph_sendCount(i)                 
              end if
            end if
          end do

          ! found anything?          
          if (fullestChannel >= 1 .and. fullestChannel <= ph_numNeigh) then
! add to counter
             ph_localSends = ph_localSends + bufsize
             call ph_postSendMsg(fullestChannel)
          end if

       end if
    end if
  end subroutine ph_sendFullestChannel
!############################################
 subroutine ph_progressSendComm()
    use Particles_rayData, only : ph_localRays
    implicit none
    integer :: opencount, ind, ierr, i

! ph_isSendCommDone flags if other communication process is 
    if (ph_numNeigh > 0 .and. .not. ph_isSendCommDone) then
       call ph_progressClosePromise()

! checks all neighbours for comm completion/free recv request
       call MPI_Testsome(ph_numNeigh, ph_sendRequest, &
            opencount, ph_sendIndex, ph_sendStatus, ierr)

       call CheckMPIErrorCode(ierr)
       ! loop over free recv requests
       do i = 1, opencount
          ind = ph_sendIndex(i)

          ph_sendCount(ind) = 0

! change status handle to closed
          if (ph_sendState(ind) == WAITING_TO_CLOSE_STATE) then
            ph_sendState(ind) = CLOSE_STATE
          end if
       end do

       !Check for completion
       ph_isSendCommDone = all(ph_sendState == CLOSE_STATE)

! all flags show closed but rays are still there to be send? -> something went wrong
       if (ph_isSendCommDone .and. &
            any(ph_sendRequest /= MPI_REQUEST_NULL .or. &
            ph_sendCount /= 0)) then
          call Driver_abortFlash('Bad shutdown')
       end if
    end if
  end subroutine ph_progressSendComm
!############################################
  subroutine ph_progressRecvComm()
    use Particles_rayData!, only : raysRealProp, ph_localRays, ph_maxNRays, raysIntProp, & 
                         !         ieion, ih2col, ihnum, inion, irad, ihcol, iblk, ihlev, & 
                         !         ithcol, ithnum, itinfo, ith2col, itrad, itnion, inih2, &
                         !         raysIntProp,iblk, ph_periodic, ph_BCcase, itposa, &
                         !         itposb, itposc, iposx, iposy, iposz, itid, isid, itblk, itnih2



    use pt_rayAsyncCounter, only : ph_localRecvs
    implicit none

    integer :: opencount, ind, ierr, i, j, msgLen, procID
    integer :: localCount

    localCount = ph_localRays

    if (ph_numNeigh > 0 .and. .not. ph_isRecvCommDone) then
       !Test all receive channels for new messages.  Save the corresponding
       !items and then post a new receive.
       opencount = 0
       call MPI_Testsome(ph_numNeigh, ph_recvRequest, opencount, &
            ph_recvIndex, ph_recvStatus, ierr)
       call CheckMPIErrorCode(ierr)

       do i = 1, opencount 

          ind = ph_recvIndex(i)
          procID = ph_recvStatus(MPI_SOURCE,i)

          if (procID /= ph_neighprocList(ind)) then
             call Driver_abortFlash("ProcID mismatch")
          end if

          call MPI_Get_count(ph_recvStatus(:,i), FLASH_REAL, msgLen, ierr)
          call CheckMPIErrorCode(ierr)
! calculate number of received rays
          ph_recvCount(ind) = msgLen / ph_transProp
! add to counter
          ph_localRecvs = ph_localRecvs + ph_recvCount(ind)

! save new rays to local ray array
          if (msgLen > 0) then
            do j = 1, ph_recvCount(ind)

              localCount=localCount+1
              if (localCount > ph_maxNRays) then
                call Driver_abortFlash("[pt_rayAsyncComm]: Insufficient space "//&
                   "in rays array: increase ph_maxNRays, German: nix platz" )
              end if

! sanity check if slot is actually free 
							if(raysRealProp(inion,localcount) .gt. 0) then
								call Driver_abortFlash &
											("recv:  occupied particle slot!")
							endif

! copy to local structure
! write only transported properties
              raysRealProp(inion,localCount)  = ph_recvBuf(itnion, j,ind)
              raysRealProp(ihnum,localCount)  = ph_recvBuf(ithnum, j,ind)
              raysRealProp(irad,localCount)   = ph_recvBuf(itrad,  j,ind)

! will be overwritten later, info always > 0, have to be careful here
              raysIntProp(isid,localCount)    = ph_recvBuf(itid,   j,ind)
              raysIntProp(iblk,localCount)    = ph_recvBuf(itblk,  j,ind)
              raysIntProp(ihlev,localCount)   = ph_recvBuf(itinfo, j,ind)
! add additional properties from send

              raysRealProp(ivelx,localCount)  = ph_recvBuf(itvelx,  j,ind) !10
              raysRealProp(ively,localCount)  = ph_recvBuf(itvely,  j,ind) !11
              raysRealProp(ivelz,localCount)  = ph_recvBuf(itvelz,  j,ind) !12

              raysRealProp(iposx,localCount)  = ph_recvBuf(itposx,  j,ind) !13
              raysRealProp(iposy,localCount)  = ph_recvBuf(itposy,  j,ind) !14
              raysRealProp(iposz,localCount)  = ph_recvBuf(itposz,  j,ind) !15

              raysRealProp(isigh,localCount)  = ph_recvBuf(itsigh,  j,ind) !18
              raysRealProp(ieion,localCount)  = ph_recvBuf(iteion,  j,ind) !24
              raysRealProp(istpd,localCount)  = ph_recvBuf(itstpd,  j,ind) !25

! assign local processor ID
              raysIntProp(iproc,localCount)   = ph_meshMe                  !26
            enddo

! update number of total rays 
            ph_localRays = localCount 
! changing ray number directly is bad

            ph_recvCount(ind) = 0
! send a new recv request
            call ph_postRecvMsg(ind)
          endif
       end do

! update total number of rays
       ph_localRays = localCount

! Check for completion
       ph_isRecvCommDone = all(ph_recvRequest == MPI_REQUEST_NULL .and. &
       ph_recvCount == 0)
    end if
  end subroutine ph_progressRecvComm
!############################################
  !Fulfill the close promise by sending a zero-byte notification message
  subroutine ph_progressClosePromise
    implicit none
    integer :: i

    do i = 1, ph_numNeigh
! after send status is set to promise_to_close_state, after counter hits target
       if ( ph_sendState(i) == PROMISE_TO_CLOSE_STATE .and. &
            ph_sendRequest(i) == MPI_REQUEST_NULL ) then

          ph_sendCount(i) = 0 !for a zero-byte message

          call ph_postSendMsg(i)
          ph_sendState(i) = WAITING_TO_CLOSE_STATE
       end if
    end do
  end subroutine ph_progressClosePromise
!############################################
  subroutine ph_commDone(isClosing)
    implicit none
    logical, intent(OUT) :: isClosing
    integer :: i

! this is a local operation, tell the local comm
! to finish up
    if (ph_numNeigh > 0) then
      do i = 1, ph_numNeigh
! if channel open and empty, mark as closed
        if (ph_sendState(i) == OPEN_STATE) then
            ph_sendState(i) = PROMISE_TO_CLOSE_STATE
        end if
      end do
! not done? get it done
      call ph_progressSendComm()
    end if

    isClosing = .true.
  end subroutine ph_commDone
!============================================
!==== MPI wrappers
!============================================
  subroutine ph_postSendMsg(ind)
    implicit none
! target
    integer, intent(IN) :: ind
    integer :: procID, msgSize, ierr

    msgSize = ph_sendCount(ind)

    if (msgSize >= 0) then
       procID = ph_neighprocList(ind)

!int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest,
!    int tag, MPI_Comm comm, MPI_Request *request)

!INPUT
!buf
!    Initial address of send buffer (choice). 
!count
!    Number of elements in send buffer (integer). 
!datatype
!    Datatype of each send buffer element (handle). 
!dest
!    Rank of destination (integer). 
!tag
!    Message tag (integer). 
!comm
!    Communicator (handle). 

!OUTPUT
!request
!    Communication request (handle). 
!IERROR
!    Fortran only: Error status (integer). 


! if msgSize is 0 then 0 byte is send
! only up to ph_transProp*ph_raysToBundle can be send
       call MPI_Isend(ph_sendBuf(1,1,ind), ph_transProp*msgSize, &
            FLASH_REAL, procID, ph_tag, ph_meshComm, &
            ph_sendRequest(ind), ierr)
       call CheckMPIErrorCode(ierr)
    end if

! wait for receipt
    ph_isSendCommDone = .false.
  end subroutine ph_postSendMsg
!############################################
  subroutine ph_postRecvMsg(ind)
    implicit none

! index is id of processor to receive rays from
    integer, intent(IN) :: ind
    integer :: procID, ierr

    procID = ph_neighProcList(ind)

!int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
!              int tag, MPI_Comm comm, MPI_Request *request)
! INTPUT
!buf
!    initial address of receive buffer (choice)
!count
!    number of elements in receive buffer (integer)
!datatype
!    datatype of each receive buffer element (handle)
!source
!    rank of source (integer)
!tag
!    message tag (integer)
!comm
!    communicator (handle)
!
! OUTPUT
!request
!    communication request (handle)

! always reads full recv buffer ph_raysToBundle*pt_transProp
    call MPI_Irecv(ph_recvBuf(1,1,ind), &
         ph_transProp*ph_raysToBundle, FLASH_REAL, procID, &
         ph_tag, ph_meshComm, ph_recvRequest(ind), ierr)
             call CheckMPIErrorCode(ierr)
         ph_isRecvCommDone = .false.
  end subroutine ph_postRecvMsg

  subroutine CheckMPIErrorCode(errorCode)
    implicit none
    integer, intent(IN) :: errorCode
    if (errorCode /= MPI_SUCCESS) call Driver_abortFlash('Error in MPI')
  end subroutine CheckMPIErrorCode

end module pt_rayAsyncComm
