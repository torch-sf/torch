! For the arepo port
module flash_run

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

! Various interfaces and variables we modify directly in Flash. - JW

#define bisect

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

use Driver_interface, ONLY : Driver_initFlash, &
    Driver_evolveFlash, Driver_finalizeFlash, Driver_getComm, Driver_getNumProcs, &
    Driver_getSimTime, Driver_init, Driver_getMype, Driver_abortFlash

use Driver_data, ONLY : dr_nbegin, dr_nend, dr_dtInit, dr_tmax, &
    dr_globalMe, dr_globalNumProcs, dr_globalComm, dr_dt, dr_dtOld, &
    dr_dtAdvect, dr_restart, dr_shortenLastStepBeforeTMax, & ! This last one is important to keep code synced!
    dr_nstep, dr_nbegin, dr_simTime, dr_dtMin, dr_initialSimTime

use Grid_interface, ONLY : Grid_getBlkData, Grid_getPointData, &
    Grid_getBlkIDFromPos, Grid_getCellCoords, Grid_getBlkIndexLimits, &
    Grid_getListOfBlocks, Grid_getBlkCenterCoords, Grid_getBlkPhysicalSize, &
    Grid_getDeltas, Grid_fillGuardCells, Grid_getBlkCornerID, Grid_getBlkPtr, &
    Grid_getBlkBoundBox, Grid_releaseBlkPtr, Grid_sortParticles, &
    Grid_getSingleCellVol, Grid_getLocalNumBlks, Grid_getMaxRefinement, &
    Grid_notifySolnDataUpdate, Grid_mapMeshToParticles, Grid_updateRefinement


#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
    use Particles_interface, ONLY : Particles_sinkSyncWithParticles, &
    Particles_getGlobalNum, Particles_mapFromMesh, Particles_sinkMoveParticles, &
    Particles_longRangeForce, &
    Particles_sinkSortParticles, Particles_addNew
#endif

use Timers_interface, ONLY : Timers_start, Timers_stop

use IO_interface, ONLY : IO_output

use IO_data, ONLY: io_checkpointFileNumber, io_plotFileNumber, &
                   io_rollingCheckpoint

use ut_qsortInterface, ONLY : ut_qsort

!use Particles_interface, ONLY : Particles_wind


use ut_interpolationInterface, ONLY: ut_polint


use RuntimeParameters_interface, ONLY : RuntimeParameters_set, &
    RuntimeParameters_get

#ifdef SINK_PART_TYPE
    use pt_sinkInterface, ONLY : pt_sinkCreateParticle, &
    pt_sinkGatherGlobal

    use pt_sinkSort
#endif


#ifdef SINK_PART_TYPE
    use Particles_sinkdata
#endif

#if defined (ACTIVE_PART_TYPE) || defined (SINK_PART_TYPE)
    use Particles_data
#endif

use base_grid_interface
!use base_particle_interface

implicit none

contains

FUNCTION get_time(value)
  ! Needed for worker_code.F90 - SCL 2021oct14
  DOUBLE PRECISION :: value
  INTEGER :: get_time
  call Driver_getSimTime(value)
  get_time=0
END FUNCTION

FUNCTION get_1blk_cell_coords(axis, blockID, limits, coords, nparts)

integer :: get_1blk_cell_coords
integer :: axis, blockID, limits
integer :: nparts
real*8  :: coords(nparts)


call Grid_getCellCoords(axis,blockID,CENTER,.false.,coords, nparts)

get_1blk_cell_coords=0
END FUNCTION get_1blk_cell_coords

FUNCTION evolve_model(value)

  DOUBLE PRECISION, INTENT(IN) :: value
  INTEGER :: evolve_model, num_procs, myID, ierr

  !call RuntimeParameters_set('tmax',value)
  dr_tmax = value
  call Driver_evolveFlash()
!  if (restart) then
!    call Driver_evolveFlash()
!  else
!    call RuntimeParameters_set('tmax',value)
!    call Driver_evolveFlash()
!    restart=.true.
!    !call RuntimeParameters_set('restart',restart)
!    !call RuntimeParameters_get('restart',dr_restart)
!  end if

  dr_nbegin = dr_nstep + 1
  dr_initialSimTime = dr_simTime

  evolve_model=0
END FUNCTION evolve_model

FUNCTION get_cell_volume(block, i, j, k, vol)

  INTEGER :: block, i, j, k, get_cell_volume
  REAL*8  :: vol

  call Grid_getSingleCellVol(block, INTERIOR, [i,j,k], vol)

get_cell_volume=0
END FUNCTION get_cell_volume

! Get the number of processors which have blocks on them.
FUNCTION get_number_of_procs(n)

  integer :: n, get_number_of_procs
  ! Note MESH_COMM is a Flash defined constant.
   call Driver_getNumProcs(MESH_COMM, n)
  ! SCL 09/19/22 Testing global_comm Flash const.
  ! Gets numprocs from Driver_initParallel rather than
  ! Driver_setupParallelEnv
  !call Driver_getNumProcs(GLOBAL_COMM, n)
  get_number_of_procs=0
END FUNCTION get_number_of_procs

FUNCTION get_all_local_num_grids(num_grid_array, nprocs)

  integer :: nprocs, get_all_local_num_grids
  integer :: communicator, myProc, ierr, i
  integer, dimension(nprocs) :: num_grid_array

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  do i=0, nprocs-1

    if (myProc == i) &

      call Grid_getLocalNumBlks(num_grid_array(i+1))

  end do

  if (myProc == 0) then

    call MPI_REDUCE(MPI_IN_PLACE, num_grid_array, nprocs, MPI_INTEGER, &
               MPI_SUM, 0, communicator, ierr)
  else

    call MPI_REDUCE(num_grid_array, num_grid_array, nprocs, MPI_INTEGER, &
                   MPI_SUM, 0, communicator, ierr)
  end if

  get_all_local_num_grids=0
END FUNCTION get_all_local_num_grids

FUNCTION get_number_of_grids(nproc, n)

  INTEGER :: n, local_n, nproc, myProc, ierr, communicator
  INTEGER :: get_number_of_grids
  INTEGER, DIMENSION(MAXBLOCKS) :: list_of_blocks

local_n = 0

call Driver_getComm(GLOBAL_COMM, communicator)
call Driver_getMype(GLOBAL_COMM, myProc)

if (myProc == nproc) then

  call Grid_getLocalNumBlks(local_n)

end if

  call MPI_REDUCE(local_n, n, 1, MPI_INTEGER, MPI_SUM, 0, communicator,ierr)

  get_number_of_grids=0
END FUNCTION get_number_of_grids

FUNCTION grid_update_refinement(gridChanged)

use Driver_data, only : dr_nstep, dr_simTime, dr_simGeneration

integer :: grid_update_refinement
logical, intent(out) :: gridChanged

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")
     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

grid_update_refinement=0
END FUNCTION grid_update_refinement

FUNCTION initialize_grid()

  INTEGER :: initialize_grid
  initialize_grid=0
END FUNCTION initialize_grid

FUNCTION initialize_code()

  INTEGER :: initialize_code, part_init
  print*, "Entered interface initialize_code()"
  print*, "Calling Driver_initParallel()"
  call flush()
  call Driver_initParallel()
  
  print*, "Calling Driver_initFlash()"
  call flush()
  call Driver_initFlash()
  ! Make sure that when we exit an evolve step, Flash actually only
  ! evolves to the end time given.
  
  print*, "Calling RuntimeParameters_set() and RuntimeParameters_get()"
  call flush()
  call RuntimeParameters_set("dr_shortenLastStepBeforeTMax",.true.)
  call RuntimeParameters_get("dr_shortenLastStepBeforeTMax",dr_shortenLastStepBeforeTMax)
  
  print*, "Exiting initialize_code()"
  call flush()
  !call Driver_init()
  !restart = .false.
  initialize_code=0
END FUNCTION initialize_code

FUNCTION cleanup_code()

  INTEGER :: cleanup_code
  call Driver_finalizeFlash()
  cleanup_code=0
END FUNCTION cleanup_code

FUNCTION recommit_parameters()

  INTEGER :: recommit_parameters
  recommit_parameters=0
END FUNCTION recommit_parameters

FUNCTION commit_parameters()

  INTEGER :: commit_parameters
  commit_parameters=0
END FUNCTION commit_parameters

FUNCTION get_global_grid_index_limits(global_indices)

  INTEGER :: global_indices(MDIM)
  INTEGER :: get_global_grid_index_limits
  call Grid_getGlobalIndexLimits(global_indices)
  get_global_grid_index_limits=0
END FUNCTION get_global_grid_index_limits

FUNCTION write_chpt()

integer write_chpt

call IO_writeCheckpoint()

write_chpt=0
END FUNCTION write_chpt

FUNCTION IO_out(output_type, fileNumber)

integer :: IO_out, fileNumber, ierr
character(8) :: output_type
logical :: endrun

! NOTE: io_writeparticles jacks up the particles sorting
! by calling Grid_sortParticles always assuming one type
! of particle. So here we resort after the proper way.

if (trim(output_type)=='chk') then

  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin, &
                 endRun, CHECKPOINT_FILE_ONLY)
  !call Particles_moveAndSort(.true.)

  fileNumber = mod(io_checkpointFileNumber, io_rollingCheckpoint)

else if (trim(output_type)=='pltpart') then

  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin, &
                 endRun, PLOTFILE_AND_PARTICLEFILE)
  !call Particles_moveAndSort(.true.)

  fileNumber = io_plotFileNumber

!else if (trim(output_type)=='all') then

!  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin, &
!                 endRun)

else

  print*, "interface:IO_output :: output filetype ", trim(output_type), " not recognized!"

  fileNumber = -1

end if

call MPI_Barrier(dr_globalComm, ierr)

IO_out=0
END FUNCTION IO_out

FUNCTION IO_num(output_type, fileNumber)

integer :: IO_num, fileNumber
character(8) :: output_type

if (trim(output_type)=='chk') then
  fileNumber = mod(io_checkpointFileNumber, io_rollingCheckpoint)
else if (trim(output_type)=='pltpart') then
  fileNumber = io_plotFileNumber
else
  print*, "interface:IO_num :: filetype ", trim(output_type), " not recognized!"
  fileNumber = -1
end if

IO_num=0
END FUNCTION IO_num

FUNCTION get_output_dir_wrapped(output_dir)
implicit none
integer :: get_output_dir_wrapped
character(len=40), intent(out) :: output_dir

call RuntimeParameters_get("output_directory", output_dir)
output_dir = trim(output_dir)

get_output_dir_wrapped=0
END FUNCTION get_output_dir_wrapped

FUNCTION get_runtime_parameter(rt_name, rt_value)
implicit none
integer :: get_runtime_parameter
character(len=*), intent(in) :: rt_name
real*8, intent(out)          :: rt_value

call RuntimeParameters_get(rt_name, rt_value)

get_runtime_parameter=0
END FUNCTION get_runtime_parameter

FUNCTION timer_summary()
integer :: timer_summary

call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))

timer_summary=0
END FUNCTION timer_summary

subroutine Driver_initParallel ()

  use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs, dr_globalComm, &
       dr_mpiThreadSupport
  !$ use omp_lib



  include "Flash_mpi.h"
  integer :: error, iprovided, errcode

#ifdef _OPENMP
#ifdef FLASH_MPI2_OPENMP
  integer, parameter :: MPI_thread_level = MPI_THREAD_SERIALIZED
#endif
#ifdef __INTEL_COMPILER
  integer(kind=kmp_size_t_kind) :: stksize
#endif
#endif
  logical :: mpiThreadSupport
  mpiThreadSupport = .false.

  !We should use MPI_Init_thread rather than MPI_Init when using multiple
  !threads so that we get a guaranteed level of thread support.

#ifdef FLASH_MPI2_OPENMP
  !We have some OpenMP parallel regions spanning MPI calls - any such
  !MPI calls are currently contained in $omp single sections and so
  !we use MPI_THREAD_SERIALIZED to give us exactly the thread support we need
  !to operate safely.  I print a warning message to the screen when your
  !MPI installation is not providing this level of thread support - it
  !is up to you whether you are happy with this risk.

  !Support Levels                     Description
  !MPI_THREAD_SINGLE     Only one thread will execute.
  !MPI_THREAD_FUNNELED   Process may be multi-threaded, but only main
  !                      thread will make MPI calls (calls are funneled to
  !                      main thread). "Default"
  !MPI_THREAD_SERIALIZED Process may be multi-threaded, any thread can
  !                      make MPI calls, but threads cannot execute MPI
  !                      calls concurrently (MPI calls are serialized).
  !MPI_THREAD_MULTIPLE   Multiple threads may call MPI, no restrictions.

  !The MPI standard says that "a call to MPI_INIT has the same effect as
  !a call to MPI_INIT_THREAD with a required = MPI_THREAD_SINGLE".

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are commented out by Josh since AMUSE already loaded this! - Josh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call MPI_Init_thread(MPI_thread_level, iprovided, error)
!  if (error /= MPI_SUCCESS) then
!     print *, "Error from MPI_Init_thread"
!     stop
!  end if
!#else
!  call MPI_Init (error)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dr_globalComm=FLASH_COMM
  call MPI_Comm_Rank (dr_globalComm, dr_globalMe, error)
  call MPI_Comm_Size (dr_globalComm, dr_globalNumProcs, error)


#ifdef _OPENMP
  if (dr_globalMe == 0) then

# ifdef FLASH_MPI2_OPENMP
     !The default thread support in Open-MPI (in the versions I have used) is
     !MPI_THREAD_SINGLE unless you configure Open-MPI with --enable-mpi-threads.

     !On Cray systems the MPI environment is limited to MPI_THREAD_SINGLE
     !by default.  This can be changed with the environmental variable
     !MPICH_MAX_THREAD_SAFETY - it has possible values of "single", "funneled",
     !"serialized" or "multiple".  To obtain MPI_THREAD_MULTIPLE thread level:
     !1) Set MPICH_MAX_THREAD_SAFETY to multiple in job submission script:
     !   export MPICH_MAX_THREAD_SAFETY="multiple"
     !2) link FLASH against a special MPI library:
     !   -lmpich_threadm.
     write(6,'(a,i3,a,i3)') " [Driver_initParallel]: "//&
          "Called MPI_Init_thread - requested level ", MPI_thread_level, &
          ", given level ", iprovided
     mpiThreadSupport = (iprovided >= MPI_thread_level);
# endif

     if (.not.mpiThreadSupport) then
        write(6,"(/ a /)") " [Driver_initParallel]: WARNING! We do not have "//&
             "a safe level of MPI thread support! (see Driver_initParallel.F90)"
        !write(6,*) "[Driver_initParalllel]: ERROR! MPI thread support too limited"
        !call MPI_Abort (dr_globalComm, errcode, error)
        !stop
     end if
  end if

  !$omp parallel
  if (dr_globalMe == 0) then
     if (omp_get_thread_num() == 0) then
        write(6,'(a,i3)') " [Driver_initParallel]: "//&
             "Number of OpenMP threads in each parallel region", &
             omp_get_num_threads()

        !Add Intel compiler specific code.  It is possible to overflow the
        !stack of the spawned OpenMP threads (e.g. WD_def 3d with block list
        !threading).  The default value for intel software stack on
        !code.uchicago.edu is 4MB (it is useful to print this information).
        !I recommend increasing this to 16MB:
        !export OMP_STACKSIZE="16M".
# ifdef __INTEL_COMPILER
        stksize = kmp_get_stacksize_s() / (1024*1024)
        write(6,'(a,i8,a)') " OpenMP thread stack size:", stksize, " MB"
# endif

        !Add Absoft compiler specific code.  The same loop iteration is
        !executed by multiple threads in parallel do loops that have 1 loop
        !iteration!  This bug happens when compiling the following test
        !problem with Absoft 64-bit Pro Fortran 11.1.4 on code.uchicago.edu.
        !
        !./setup unitTest/Multipole -auto -geometry=cartesian -3d -maxblocks=1 \
        !  +newMpole +noio threadBlockList=True -nxb=64 -nyb=64 -nzb=64
        !
        ! Set lrefine_min = lrefine_max = 1 in the flash.par.
# ifdef __ABSOFT__
        print *, ""
        print *, "WARNING!!!! Absoft compiler OpenMP bug!!!!"
        print *, "A parallel do loop with 1 loop iteration will be executed incorrectly"
        print *, ""
# endif

     end if
  end if
# ifdef DEBUG_THREADING
  write(6,'(a,i3,a,i3)') " [Driver_initParallel]: MPI rank ", dr_globalMe, &
       " has a team that includes thread ", omp_get_thread_num()
# endif
  !$omp end parallel
#endif

  dr_mpiThreadSupport = mpiThreadSupport

end subroutine Driver_initParallel

end module flash_run
