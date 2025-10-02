!!****if* source/Simulation/SimulationMain/Cube/pt_initVoronoiPositions
!!
!! NAME
!!    pt_initVoronoiPositions
!!
!! SYNOPSIS
!!    pt_initVoronoiPositions( logical, INTENT(out) :: success,
!!                             logical, INTENT(out) :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize zero-mass particles with distribution drawn from input
!!    dataset. Arrays for particle positions are dynamically allocated 
!!    and de-allocated here too. 
!!
!! ARGUMENTS
!! success/partPosInitialized   : boolean indicating whether particle positions were
!!                successfully initialized. This is not really relavant
!!                for this version of the routine.
!! updateRefine : is true if the routine wishes to retain the already 
!!                initialized particles instead of reinitializaing them
!!                as the grid refines.
!!
!! NOTES
!! This routine is called by source/Particles/Particles_Initialization/Particles_initPositions.F90
!! Particles_initPositions contains the compiler flag #define VORAMR. This flag must be
!! uncommented in order for this routine to be called.
!!
!! If pt_maxPerProc in the flash.par is not set high enough, this script should abort.
!! The abort may appear as a "double free or corruption (!prev)" C memory error due to
!! IO not being initialized yet at the time of particle placement.
!!
!! Developed by Sean C. Lewis (Drexel University) as part of the VorAMR project.
!!***

!!REORDER(4): solnData
subroutine pt_initVoronoiPositions (partPosInitialized,updateRefine)

  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc, &
       pt_posInitialized, pt_meshMe, pt_globalMe, pt_meshComm, &
       pt_indexList, pt_indexCount
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_getCellCoords,&
       Grid_getBlkIndexLimits, Grid_getBlkData, &
       Grid_getBlkIDFromPos
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get 
  use tree, only: grid_xmin, grid_xmax, grid_ymin, grid_ymax, &
       grid_zmin, grid_zmax
  use Particles_interface, ONLY: Particles_moveAndSort
  !use Simulation_data, ONLY : sim_ptMass, sim_densityThreshold, sim_print

  USE HDF5
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Particles.h"

  logical, intent(INOUT) :: partPosInitialized
  logical, intent(OUT) :: updateRefine

  integer       :: i, j, k, b
  integer       :: xmin, xmax, ymin, ymax, zmin, zmax
  integer       :: blockID, blkCount, ps, pe, part_PE, ierr
  real :: boundBox(2, MDIM)
  integer :: blkList(MAXBLOCKS)
  real,dimension(10) :: x, y, z
  real,dimension(:,:,:,:),pointer :: solnData
  real :: mass, total_mass
  integer :: globalNumParticles

  !Input dataset read-in (EDIT THIS FOR YOUR OWN FILE STRUCTURE)
  CHARACTER(LEN=24), PARAMETER :: filename   = "voramr_input.hdf5"
  CHARACTER(LEN=16), PARAMETER :: coords_dsetname   = "Coordinates"
  CHARACTER(LEN=16), PARAMETER :: group1name = "PartType0"
  ! -----------------------------------------------------------
  ! -----------------------------------------------------------
  INTEGER(HID_T)  :: file_id   ! File identifier
  INTEGER(HID_T)  :: group1_id ! Group1 identifier
  INTEGER(HID_T)  :: coords_dset_id, masses_dset_id, velocs_dset_id
  INTEGER(HID_T)  :: coords_space_id, masses_space_id, velocs_space_id
  INTEGER(HID_T)  :: dtype_id  ! Dataspace type
  INTEGER(KIND=4) :: error      ! error flag
  INTEGER         :: cols, rows ! number of columns and rows of dataset
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords_dset_data ! Dataset memory
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: coords_buffer ! 1D for MPI_Bcast
  INTEGER(HSIZE_T), DIMENSION(2) :: coords_data_dims ! dataset used dimensions
  INTEGER(HSIZE_T), DIMENSION(2) :: coords_max_dims

  real :: xcoord, ycoord, zcoord ! Imported gas particle coordinates
!--------------------------------------------------------------



  ! Get block data 
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)


  call RuntimeParameters_get('use_voramr', use_voramr)
  if (.not.use_voramr) then
     print*, "Skipping VorAMR"
     call Logfile_stamp('Skipping VorAMR', "[pt_initVoronoiPositions]")
     return
  endif
    
  if (pt_globalMe.EQ.MASTER_PE) then
     call Logfile_stamp('VorAMR. Building refined grid', "[pt_initVoronoiPositions]")
  endif
  
  updateRefine=.true. !.false.
  if(partPosInitialized) return

  call RuntimeParameters_get('voramr_input', voramr_input)
  
  pt_numLocal=0
  globalNumParticles=0
  if (pt_globalMe .EQ. MASTER_PE) then
     print *, 'Starting HDF5 Fortran Read'
     CALL h5open_f(error)
     print*, 'Open the file: ', voramr_input !filename
     CALL h5fopen_f (voramr_input, H5F_ACC_RDWR_F, file_id, error)
     print*, 'Open the group: ', group1name
     CALL h5gopen_f (file_id, group1name, group1_id, error)
     print*, 'Open the dataset: ', coords_dsetname
     CALL h5dopen_f(group1_id, coords_dsetname, coords_dset_id, error)
     print*, 'Get dataspace ID and dimensions.'
     CALL h5dget_space_f(coords_dset_id, coords_space_id, error)
     CALL h5sget_simple_extent_dims_f(coords_space_id, coords_data_dims, coords_max_dims, error)
     cols = coords_data_dims(1)
     rows = coords_data_dims(2)
     print*, 'Number of columns and rows:', cols, rows
     print*, ' Dynamically allocate dimensions to coords_dset_data for reading.'
     ALLOCATE(coords_dset_data(cols, rows))
     print*, 'Reading data...'
     CALL h5dread_f(coords_dset_id, H5T_NATIVE_DOUBLE, coords_dset_data, coords_data_dims, error)

     print*, 'Closing file...'
     CALL flush()
     CALL h5gclose_f(group1_id,error)
     CALL h5fclose_f(file_id,error)
     CALL h5close_f(error)
  endif

  ! Broadcast full dataset to all processors

  ! === Step 2: Broadcast metadata ===
  CALL MPI_Bcast(cols, 1, FLASH_INTEGER, MASTER_PE, FLASH_COMM, ierr)
  CALL MPI_Bcast(rows, 1, FLASH_INTEGER, MASTER_PE, FLASH_COMM, ierr)

    print*, "allocating buffer on all ranks"
  CALL flush()
  ! === Step 3: Allocate buffer on all ranks ===
  ALLOCATE(coords_buffer(cols * rows))  ! <== NEW: Allocate 1D buffer on all ranks

    print*, "reshaping buffer"
  CALL flush()
  ! === Step 3.5: On MASTER_PE, copy 2D data into 1D buffer ===
  if (pt_globalMe .EQ. MASTER_PE) then
     coords_buffer = reshape(coords_dset_data, [cols * rows])  ! <== NEW
  endif

    print*, "bcasting 1d buffer"
  CALL flush()
  ! === Step 4: Broadcast 1D buffer ===
  CALL MPI_Bcast(coords_buffer, cols*rows, FLASH_REAL, MASTER_PE, FLASH_COMM, ierr)  ! <== FIXED

      print*, "allocating coords_dset on all ranks"
  CALL flush()
  ! === Step 5: On non-master ranks, reshape back to 2D ===
  if (.NOT. ALLOCATED(coords_dset_data)) then
     ALLOCATE(coords_dset_data(cols, rows))                  ! <== MOVED: Allocation happens here
  endif
  coords_dset_data = reshape(coords_buffer, [cols, rows])    ! <== NEW: Rebuild original 2D array

  print*, 'Placing particles...'
  CALL flush()

  ! loop over blocks on current processor, get limits, 
  ! then loop over all particles and add the ones that are 
  ! inside the current block.

  do b = 1, blkCount
      blockID = blkList(b)
      call Grid_getBlkBoundBox(blockID, boundBox)

      ! get physical domain boundaries on block
      xmin = boundBox(1, IAXIS) 
      xmax = boundBox(2, IAXIS)
      ymin = boundBox(1, JAXIS)
      ymax = boundBox(2, JAXIS)
      zmin = boundBox(1, KAXIS) 
      zmax = boundBox(2, KAXIS) 

      do i=1, rows
         !if (pt_meshMe .eq. pt_globalMe) then
            xcoord = coords_dset_data(1,i)
            ycoord = coords_dset_data(2,i)
            zcoord = coords_dset_data(3,i)

            ! Check if the point is within the bounding box
           if (xcoord >= xmin .and. xcoord <= xmax .and. &
               ycoord >= ymin .and. ycoord <= ymax .and. &
               zcoord >= zmin .and. zcoord <= zmax) then
               
                pt_numLocal = pt_numLocal + 1
                if (pt_numLocal .gt. pt_maxPerProc) then
                        print*,"ABOVE PT_MAXPERPROC!!!!!"
                     call Driver_abortFlash('Particles_initPositions: Particle number exceeds pt_maxPerProc. Increase.')
                endif
                ps = pt_numLocal
                particles(:,ps) = 0
                particles(PROC_PART_PROP,ps) = pt_globalMe
                particles(BLK_PART_PROP,ps)  = blockID
                particles(POSX_PART_PROP,ps) = xcoord
                particles(POSY_PART_PROP,ps) = ycoord
                particles(POSZ_PART_PROP,ps) = zcoord
                particles(TYPE_PART_PROP,ps) = NONEXISTENT
                !NONEXISTENT particle type allows VorAMR particles and their memory to be cleared by Grid_sortParticles() and reused.

            endif ! bounding box check

      enddo
    enddo


  print*, "Out of placing loop, getting global num particles. pt_numLocal=",pt_numLocal
  call flush()
  !call pt_gatherGlobal()
  call Particles_getGlobalNum(globalNumParticles)

  print*, "Out of placing loop, getting global num particles = ",globalNumParticles
  print*, "Deallocating particle coordinate dataset."
  call flush()
  if (ALLOCATED(coords_dset_data)) then
        DEALLOCATE(coords_dset_data)
  endif

  ! === Free temporary buffer ===
  if (ALLOCATED(coords_buffer)) then
        DEALLOCATE(coords_buffer)
  endif

  print*, "Done placing. Exiting Particles_initPositions"
  CALL flush()


  partPosInitialized = .true.
  pt_posInitialized = partPosInitialized
  return

end subroutine pt_initVoronoiPositions
