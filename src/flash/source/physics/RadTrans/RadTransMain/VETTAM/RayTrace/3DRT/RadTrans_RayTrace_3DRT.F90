!===============================================================================
!
! Subroutine: RadTrans_raytrace_3DRT
! Path: source/physics/RadTrans/RadTransMain/VET/RayTrace/3DRT/RadTrans_raytrace_3DRT.F90
!
!===============================================================================
!
! Description: Ray-Trace step for sink particles and diffuse radiation, used with the VET cooling module. 
! Original Author: Manuel Jung (2018)
! Modified by Shyam Harimohan Menon to incorporate with VET module (2020-2021)
! Email: shyam.menon@anu.edu.au
!
!===============================================================================
!
! TO DO:
!
!===============================================================================
!
!===============================================================================
!
  !
  ! debug gives general geometrical information of the angle loops
  ! in the log-file. Each mpi-process will create it's own log-file.
  ! Verbose is more general and even spits out information on
  ! a cell-by-cell level, be careful with the size of log-files.
  ! If one wants the raytracer to be verbose, debug will automatically
  ! be true as well.
  !

#undef DEBUG_RT
#undef VERBOSE_RT

#ifdef VERBOSE_RT
#define DEBUG_RT
#endif

MODULE RadTrans_RayTrace_3DRT
#include "Flash.h"
#include "constants.h"

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs , dr_simGeneration

  use tree, only : lrefine, lnblocks, grid_changed

  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_fillGuardCells

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Logfile_interface, ONLY : Logfile_stamp

  use raytrace_data

  use rt_data_raytrace_3drt, ONLY : domainSizeX, domainSizeY, domainSizeZ, &
       nrOfAnglesPerGroup, nrOfAngles

  use RadTrans_hybridCharModule

  use ModuleRadTransCommunicateGlobalTreeData, ONLY : RadTrans_CommunicateGlobalTreeData

  use gr_interface,     ONLY : gr_findAllNeghID
  !
  implicit none
#include "Flash_mpi.h"

  ! -----------------------
  ! The heavy memory stuff:
  ! -----------------------
  ! These fields will be allocated according to the maximum
  ! number of leaf blocks of all mpi processes. This makes
  ! sure, that we use the available memory efficiently,
  ! even if MAXBLOCKS is much larger than the maximum number
  ! of leaf blocks.
  !

  integer, dimension(NDIM), parameter :: ib = NGUARD*KD
  integer, dimension(NDIM), parameter :: ie = NGUARD*KD+NB

  ! This is the list containing the id's of leaf blocks
  ! handled by this mpi process.
  !
  integer, dimension(MAXBLOCKS)     :: leafList

  integer              :: maxLevel, maxLevel_local
  integer              :: ierr

  logical             :: do_lambda  = .false. ! is true at first
  !                                             iteration step.
  !
  ! This is the number of leaf blocks this mpi process is
  ! handling locally.
  !
  integer             :: nrOfLeafBlocks
  !
  ! The maxmimum number of leaf blocks of all mpi processes
  !
  integer             :: maxNrOfLeafBlocks
  !
  real                               :: allocatedMemory
  
  ! Healpix randomization angles
  ! A random axis of rotation is chosen (2 random nos) + a random angle about that axis (1 no)
  real                            :: rand_angles(3) = 0.
  integer, save :: rand_counter = 0

  ! PUBLIC :: rad_allocate, rad_deallocate, rad_gatherAll, &
  !   rad_initData, RadTrans_raytrace_3DRT, allocGeneration, FirstStep

CONTAINS

SUBROUTINE rad_allocate(firstAllocation)
  implicit none
  logical, INTENT(IN) :: firstAllocation
  !
  ! Now allocate. This is now handled dynamically and depends on
  ! the maximum number of leaf blocks of all mpi processes.
  !

     allocatedMemory = &
       dr_globalNumProcs &
       + maxNrOfLeafBlocks * dr_globalNumProcs &
       + (NXB+1)*(NYB+1)*6 * nrOfAnglesPerGroup &
           * maxNrOfLeafBlocks * (dr_globalNumProcs + 1)

     allocatedMemory = allocatedMemory * 8.d0  ! double precision

     if(dr_globalMe.eq.MASTER_PE) then
        if(firstAllocation) then
           call writeMemoryBanner
        else
           write(io,*) ''
           write(io,*) 'REQUIRED MEMORY [MB]', allocatedMemory/1.d6
           write(io,*) ''
           write(io,*) 'Current Mesh settings'
           write(io,*) ''
           write(io,*) '  maxNrOfLeafBlocks    :', maxNrOfLeafBlocks
           write(io,*) '  nrOfFaceValues       :', (NXB+1)*(NYB+1)*6
           write(io,*) '  dr_globalNumProcs    :', dr_globalNumProcs
           write(io,*) '  nrOfAnglesPerGroup   :', nrOfAnglesPerGroup
           write(io,*) ''
        endif
     endif
     if(dr_globalMe.eq.MASTER_PE) then
        write(io,*) '--------------- TRYING TO ALLOCATE MEMORY... ----------------------'
     endif

     !    Get memory:
     !      - faceValueAll:       Values of local column density at all faces of
     !                            all blocks on all processors.
     call Timers_start("rad_allocate")
     allocate(nrOfLeafBlocksList(0:dr_globalNumProcs-1),&
              leafListAll(maxNrOfLeafBlocks,0:dr_globalNumProcs-1),&
              faceValueAll(ib(1):ie(1),ib(2):ie(2),3,2,maxNrOfLeafBlocks,nrOfAnglesPerGroup,0:dr_globalNumProcs-1),&
              stat=ierr)
     if(ierr.NE.0) CALL Driver_abortFlash("Could not allocate memory in Radtrans_Raytrace_3DRT.")
     call Timers_stop("rad_allocate")
     !
     if(dr_globalMe.eq.MASTER_PE) then
        write(io,*) '---------------        SUCCESS...             ----------------------'
        write(io,*) ''
     endif
     !
     !
END SUBROUTINE rad_allocate



SUBROUTINE rad_gatherAll()
  implicit none
  integer :: p, d, neghsize, dir, b, nrOfLeafBlocksSum, blockID, s
  integer, dimension(dr_globalNumProcs) :: displs,recvcounts
  integer, parameter, dimension(MDIM,2*NDIM) :: edge = &
     reshape (source = (/&
     LEFT_EDGE,  CENTER, CENTER, &
     RIGHT_EDGE, CENTER, CENTER, &
     CENTER, LEFT_EDGE,  CENTER, &
     CENTER, RIGHT_EDGE, CENTER, &
     CENTER, CENTER,  LEFT_EDGE, &
     CENTER, CENTER, RIGHT_EDGE  &
     /), shape = (/MDIM,2*NDIM/))
  !
  ! Create a list that contains the number for leaf blocks each processor owns
  !
#ifdef DEBUG_RT
  if(dr_globalMe.eq.MASTER_PE) write(*,*) 'gathering nrOfLeafBlocks...'
#endif
  call Timers_start("rad_MPI_ALLGATHER")
  call MPI_ALLGATHER(nrOfLeafBlocks, 1, MPI_INTEGER, &
                     nrOfLeafBlocksList, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD,ierr)
  call Timers_stop("rad_MPI_ALLGATHER")
#ifdef DEBUG_RT
  if(dr_globalMe.eq.MASTER_PE) write(*,*) 'gathering nrOfLeafBlocks... done'
#endif
  !
  ! Merge all leaf lists form different processors
  !
  leafListAll = 0
#ifdef VERBOSE_RT
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'gathering leaf-block lists...'
#endif
  call Timers_start("rad_MPI_ALLGATHER")
  call MPI_ALLGATHER(leafList, maxNrOfLeafBlocks, MPI_INTEGER, &
                     leafListAll, maxNrOfLeafBlocks, MPI_INTEGER, &
                     MPI_COMM_WORLD,ierr)
  call Timers_stop("rad_MPI_ALLGATHER")
  !
#ifdef VERBOSE_RT
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'gatherin leaf-block lists... done'
#endif
  !
  ! Get updated global hierarchy tree data for tree walk in create_cut_block_list
  !
#ifdef DEBUG_RT
  write(*,*) 'call RadTrans_CommunicateGlobalTreeData...'
#endif
  call Timers_start("rad_MPI_ALLGATHER_tree")
  call RadTrans_CommunicateGlobalTreeData
  call Timers_stop("rad_MPI_ALLGATHER_tree")
#ifdef DEBUG_RT
  write(*,*) 'call RadTrans_CommunicateGlobalTreeData... done'
#endif

  call Timers_start("rad_MPI_ALLGATHER")
  ! Find the global maximum refinement level.
  call MPI_ALLREDUCE(maxLevel_local, maxLevel, &
         1, MPI_INTEGER, &
         MPI_MAX, &
         MPI_COMM_WORLD, ierr)



  ! Given a Blkno and a Procno rblockListAll (reverse block list all)
  ! will return the global contigous number of that block
  ALLOCATE(rblockListAll(MAXBLOCKS,0:dr_globalNumProcs-1))
  rblockListAll(:,:) = -1
  DO b = 1, nrOfLeafBlocks
    blockID = leafList(b)
    rblockListAll(blockID, dr_globalMe) = 1 ! This is a leaf block
  END DO
  call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, &
                     rblockListAll, MAXBLOCKS, MPI_INTEGER, &
                     MPI_COMM_WORLD,ierr)
  ! All leaf blocks are marked with a "1", now make the cummulative Sum.
  ! In principle any other unique numbering would also be fine, it just has to
  ! be contingues and the same on all procs.
  s = 0
  DO p = 0, dr_globalNumProcs-1
    DO b = 1, MAXBLOCKS
      IF(rblockListAll(b,p).EQ.1) THEN
        rblockListAll(b,p) = s
!        print *,"b -> s: ", b, " -> ", s
        s = s + 1
      END IF
    END DO
  END DO

  ! SurrBlkSum holds information about block and proc number of neighbour leaf
  ! blocks for a leaf blocks of this proc.
  allocate (SurrBlkSum(nrOfLeafBlocks))

  ! Make sure surr_blks is updated.
  call gr_ensureValidNeighborInfo(10)

  DO b = 1, nrOfLeafBlocks
    blockID = leafList(b)
    call gr_findAllNeghID(blockID, SurrBlkSum(b))
  END DO

  nrOfLeafBlocksSum = SUM(nrOfLeafBlocksList)
  
  ! SurrBlkSumArray will hold information about block and proc number of
  ! neighbour leaf blocks for all leaf blocks on all procs - therefore global
  ! information. This is needed to follow the rays over blocks on the global grid.
  ! We only keep neighbours at faces and disregard the information of blocks
  ! which touch only edges or corners, since we do not need them.
  ALLOCATE(SurrBlkSumArray(BLKNO:PROCNO,2**(NDIM-1),2*NDIM,nrOfLeafBlocksSum), &
           SurrBlkSumNumNegh(2*NDIM,nrOfLeafBlocksSum))
  neghsize = (PROCNO-BLKNO+1)*(2**(NDIM-1))*(2*NDIM)

  ! This is the displacment of the local leaf block number.
  d = SUM(nrOfLeafBlocksList(0:dr_globalMe-1))

  ! Copy the local information from fortran custom type structures to arrays.
  DO b = 1, nrOfLeafBlocks
    DO dir = 1, 2*NDIM
      SurrBlkSumArray(BLKNO:PROCNO,1:2**(NDIM-1),dir,b+d) &
        = SurrBlkSum(b) &
            % regionInfo(edge(1,dir),edge(2,dir),edge(3,dir)) &
            % details(BLKNO:PROCNO,1:2**(NDIM-1))
      SurrBlkSumNumNegh(dir,b+d) &
        = SurrBlkSum(b) &
            % regionInfo(edge(1,dir),edge(2,dir),edge(3,dir)) &
            % numNegh
    END DO
  END DO

  ! Calculate displacments
  recvcounts(:) = nrOfLeafBlocksList(:)*(2*NDIM)
  d = 0
  DO p=1,dr_globalNumProcs
    displs(p) = d
    d = d + recvcounts(p)
  END DO

  ! Communicate neighbour informations from local blocks to all other procs.
  ! - number of neighbours.
  call MPI_ALLGATHERV( &
       MPI_IN_PLACE, &
       0, &
       MPI_DATATYPE_NULL, &
       SurrBlkSumNumNegh, &
       recvcounts, &
       displs, &
       FLASH_INTEGER, &
       MPI_COMM_WORLD,ierr)
  
  ! Calculate displacments
  recvcounts(:) = nrOfLeafBlocksList(:) * neghsize
  d = 0
  DO p=1,dr_globalNumProcs
    displs(p) = d
    d = d + recvcounts(p)
  END DO

  ! Communicate neighbour informations from local blocks to all other procs.
  ! - Blk/Proc numbers of neighbour blocks
  call MPI_ALLGATHERV( &
       MPI_IN_PLACE, &
       0, &
       MPI_DATATYPE_NULL, &
       SurrBlkSumArray, &
       recvcounts, &
       displs, &
       FLASH_INTEGER, &
       MPI_COMM_WORLD,ierr)

  call Timers_stop("rad_MPI_ALLGATHER")

end subroutine rad_gatherAll

SUBROUTINE rad_deallocate()
#include "Flash.h"
#include "constants.h"
  !use Driver_data, ONLY: dr_globalMe, dr_simGeneration
  !use Timers_interface, ONLY : Timers_start, Timers_stop
  !use RadTrans_hybridCharModule, ONLY : io
  implicit none
  !
  ! Deallocate all the heavy stuff
  !
  if(dr_globalMe.eq.MASTER_PE) then
     write(io,*) 'deallocating...'
  endif

  call Timers_start("rad_allocate")
  deallocate(nrOfLeafBlocksList,&
             leafListAll,&
             faceValueAll,&
             SurrBlkSumArray,&
             SurrBlkSumNumNegh,&
             SurrBlkSum,&
             rblockListAll,&
             stat=ierr)
  if(ierr.NE.0) CALL Driver_abortFlash("Could not deallocate memory in Radtrans_Raytrace_3DRT.")

  call Timers_stop("rad_allocate")

  if(dr_globalMe.eq.MASTER_PE) then
     write(io,*) 'deallocating... done'
  endif
END SUBROUTINE rad_deallocate

SUBROUTINE rad_initData(first)
  implicit none
  logical, intent(in) :: first
  INTEGER :: b, blk
  !
  ! Determine the number of leaf nodes on this processor and
  ! create their list.
  !
  leafList(:) = 0
#ifdef VERBOSE_RT
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'call create_leaf_list...'
#endif
  call Grid_getListOfBlocks(LEAF, leafList, nrOfLeafBlocks)
#ifdef VERBOSE_RT
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'call create_leaf_list... done'
#endif
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'current number of leaf blocks:',nrOfLeafBlocks
  !
  ! Find the maximum number of leaf blocks that a processor has, where
  ! nrOfLeafBlocks is the current number of leaf blocks on this processor.
  !
#ifdef VERBOSE_RT
  if(dr_globalMe.eq.MASTER_PE) write(io,*) 'call MPI_ALLREDUCE for nrOfLeafBlocks...'
#endif
  !
  call Timers_start("rad_MPI_ALLREDUCE")
  call MPI_ALLREDUCE(nrOfLeafBlocks, maxNrOfLeafBlocks, 1, MPI_INTEGER, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  call Timers_stop("rad_MPI_ALLREDUCE")

  call rad_allocate(first)

  maxLevel_local = 0
  !
  ! Get the current maximum refinement level
  ! for all local blocks.
  do b=1,nrOfLeafBlocks
     blk = leafList(b)
     maxLevel_local = MAX(maxLevel_local, lrefine(blk))
  enddo

END SUBROUTINE rad_initData

SUBROUTINE rad_resetData
  implicit none
  INTEGER :: b, blk

  call Timers_start("rad_resetData")

  do b=1,nrOfLeafBlocks
     !
     blk = leafList(b)

     call Grid_getBlkPtr(blk, solnData)
     !
     !  Save old radiation field
     !
#ifdef MEAO_VAR
     solnData(MEAO_VAR,:,:,:) = solnData(MEAN_VAR,:,:,:)
#endif
     !
     !  reset specific and mean intensity
     !
#ifdef MEAN_VAR
     solnData(MEAN_VAR,:,:,:) = 0.0
#endif
#ifdef RPXX_VAR
     solnData(RPXX_VAR:RPZZ_VAR,:,:,:) = 0.0
#endif

#ifdef LAMB_VAR
     if(do_lambda) then
        solnData(LAMB_VAR,:,:,:) = 0.0
     endif
#endif
     call Grid_releaseBlkPtr(blk, solnData)

  enddo ! b
  !
  call Timers_stop("rad_resetData")
  !
end subroutine rad_resetData


SUBROUTINE rad_resetsinkContainers
  implicit none
  INTEGER :: b, blk

  call Timers_start("rad_resetData")

  do b=1,nrOfLeafBlocks
     !
     blk = leafList(b)

     call Grid_getBlkPtr(blk, solnData)
     !
     !  Save old radiation field
     !
#ifdef FSPT_VAR
     solnData(FSPT_VAR,:,:,:) = 0.0
     solnData(STMX_VAR:STMZ_VAR,:,:,:) = 0.0
#endif

#ifdef LAMB_VAR
     if(do_lambda) then
        solnData(LAMB_VAR,:,:,:) = 0.0
     endif
#endif
     call Grid_releaseBlkPtr(blk, solnData)

  enddo ! b
  !
  call Timers_stop("rad_resetData")
  !
end subroutine rad_resetsinkContainers



subroutine writeMemoryBanner()
      !
      write(io,*) '-------------------------------------------------------------------'
      write(io,*) '------------- IMPORTANT NOTE ON MEMORY CONSUMPTION ----------------'
      write(io,*) '-------------------------------------------------------------------'
      write(io,*) ''
      write(io,*) 'To this point, the raytrace_3DRT routine of the hybrid_char_3DRT   '
      write(io,*) 'module has allocated the following amount of memory.'
      write(io,*) ''
      write(io,*) 'REQUIRED MEMORY [MB]', allocatedMemory/1.d6
      write(io,*) ''
      write(io,*) 'Actually, this is the memory, that one mpi process has allocated.  '
      write(io,*) 'If the code crashed, this is probably the point where to tweak your'
      write(io,*) 'runtime parameters. The most memory consuming fields are the face  '
      write(io,*) 'values for the optical depths and the specific intensities.        '
      write(io,*) 'The memory Mem (in single precision) that is allocated for these   '
      write(io,*) 'fields is computed by as follows'
      write(io,*) ''
      write(io,*) '            Mem = maxNrOfLeafBlocks * nrOfFaceValues * '
      write(io,*) '                  dr_globalNumProcs * nrOfAnglesPerGroup * 4 byte'
      write(io,*) ''
      write(io,*) 'with'
      write(io,*) ''
      write(io,*) '  maxNrOfLeafBlocks    : the maxmimum number of leaf blocks of     '
      write(io,*) '                         all mpi processes'
      write(io,*) '  nrOfFaceValues       : the number of face values on one block for each ray'
      write(io,*) '                         (in 3D: 3*(NXB+1)*(NYB+1)*2    '
      write(io,*) '                                 i.e. 3 faces of slice where ray exits, for 2 variables '
      write(io,*) '  dr_globalNumProcs               : the total number of mpi processes'
      write(io,*) '  nrOfAnglesPerGroup   : the number of angles per angle group'
      write(io,*) ''
      write(io,*) ''
      write(io,*) 'Your current settings are'
      write(io,*) ''
      write(io,*) '  maxNrOfLeafBlocks    :', maxNrOfLeafBlocks
      write(io,*) '  nrOfFaceValues       :', (NXB+1)*(NYB+1)*6
      write(io,*) '  dr_globalNumProcs    :', dr_globalNumProcs
      write(io,*) '  nrOfAnglesPerGroup   :', nrOfAnglesPerGroup
      write(io,*) ''
      write(io,*) 'If you have a lot of (leaf) blocks, that one mpi process has to    '
      write(io,*) 'handle, you need to reduce the number of angles per angle-group    '
      write(io,*) 'angle-group by either incresing the number of angle-groups or by   '
      write(io,*) 'using more mpi processes.                                          '
      write(io,*) ''
      !
endsubroutine writeMemoryBanner

! Run this before doing the first iteration
subroutine FirstStep()
  IMPLICIT NONE
  real :: u
  ! generate Healpix angle randomization
  IF(rt_healpix_randomize.GT.0) THEN
    if(dr_globalMe.eq.MASTER_PE) then
      IF(MOD(rand_counter,rt_healpix_randomize).EQ.0) THEN
        IF(rand_counter.eq.0) &
          CALL RANDOM_SEED()
        !3 random nos to be sampled

        !Generate theta of rotation axis using inverse sampling . 
        !Note: uniform theta does not sample the sphere correctly
        call RANDOM_NUMBER(u)
        rand_angles(1) = ACOS(1-2*u)

        !random phi of rotation axis
        call RANDOM_NUMBER(rand_angles(2))
        rand_angles(2) = rand_angles(2)*2*PI

        !+ an angle psi uniformly sampled representing the angle of rotation
        call RANDOM_NUMBER(rand_angles(3))
        rand_angles(3) = rand_angles(3)*2*PI

      END IF
      rand_counter = rand_counter + 1
    END IF
    call MPI_Bcast(rand_angles, 3, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  END IF
  IF(allocGeneration.eq.-1) THEN
    ! first run
    if(dr_globalMe.eq.MASTER_PE) &
      call Logfile_stamp("Allocate memory and initialize after start", "[3DRT] ")
    allocGeneration = dr_simGeneration
    call rad_initData(.true.)
    ! gather some arrays
    ! and take the global maximum of maxLevel
    call rad_gatherAll()

  ! Note : the grid_changed parameter remains non-zero for two steps after grid refinement
  ! so some unnecessary allocation occurs that is unavoidable as forcing grid_changed to 
  ! be zero breaks some solvers (e.g. the 8 wave solver)
  ELSE IF(grid_changed.NE.0) THEN
    ! refinement/derefinements have happened
    if(dr_globalMe.eq.MASTER_PE) &
      call Logfile_stamp("Reallocate memory and initialize after Re/Derefinment", "[3DRT] ")
    allocGeneration = dr_simGeneration
    call rad_deallocate() ! nothing might happen, if first run

    call rad_initData(.false.)
    ! gather some arrays
    ! and take the global maximum of maxLevel
    call rad_gatherAll()

  END IF
end subroutine FirstStep

END MODULE RadTrans_RayTrace_3DRT
!
!===============================================================================
!
