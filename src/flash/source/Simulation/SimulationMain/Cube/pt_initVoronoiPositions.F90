!!****if* source/Simulation/SimulationMain/Cube/pt_initVoronoiPositions
!!
!! NAME
!!    pt_initVoronoiPositions
!!
!! SYNOPSIS
!!
!!    pt_initVoronoiPositions( logical, INTENT(out) :: success,
!!                             logical, INTENT(out) :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize zero-mass particles with distribution drawn from input
!!    dataset. Arrays for particle positions are dynamically allocated 
!!    and de-allocated here too. 
!! ARGUMENTS
!!
!! success/partPosInitialized   : boolean indicating whether particle positions were
!!                successfully initialized. This is not really relavant
!!                for this version of the routine.
!! updateRefine : is true if the routine wishes to retain the already 
!!                initialized particles instead of reinitializaing them
!!                as the grid refines.
!!
!!***

!!REORDER(4): solnData
subroutine pt_initVoronoiPositions (partPosInitialized,updateRefine)

  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc, &
       pt_posInitialized, pt_meshMe, pt_globalMe, pt_meshComm

  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_getCellCoords,&
       Grid_getBlkIndexLimits, Grid_getBlkData, &
       Grid_getBlkIDFromPos
  
  !use Simulation_data, ONLY : sim_ptMass, sim_densityThreshold, sim_print

  USE HDF5
  implicit none

#include "constants.h"
#include "Flash.h"

  logical, intent(INOUT) :: partPosInitialized
  logical, intent(OUT) :: updateRefine

  integer       :: i, j, k
  integer       :: blockID, blkCount, ps, pe
  real,dimension(10) :: x, y, z
  real,dimension(:,:,:,:),pointer :: solnData
  real :: mass, total_mass
  double precision :: globalNumParticles

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
  INTEGER(HSIZE_T), DIMENSION(2) :: coords_data_dims ! dataset used dimensions
  INTEGER(HSIZE_T), DIMENSION(2) :: coords_max_dims

  real :: xcoord, ycoord, zcoord
!--------------------------------------------------------------

  if (pt_meshMe.EQ.MASTER_PE) print*,'Entering Particles_initPositions.F90'
  updateRefine=.true. !.false.
  if(partPosInitialized) return

  pt_numLocal=0
  if (pt_meshMe .EQ. MASTER_PE) then
     print *, 'Starting HDF5 Fortran Read'
     CALL h5open_f(error)
     print*, 'Open the file: ', filename
     CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
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
  
  print*, 'Placing particles...'
  CALL flush()
  do i=1, rows
     if (pt_meshMe .eq. pt_globalMe) then
        xcoord = coords_dset_data(1,i)
        ycoord = coords_dset_data(2,i)
        zcoord = coords_dset_data(3,i)
        call Grid_getBlkIDFromPos([xcoord, ycoord, zcoord], blockID, pt_meshMe, pt_meshComm)
        pt_numLocal = pt_numLocal + 1
        if (pt_numLocal .gt. pt_maxPerProc) &
             call Driver_abortFlash('Particles_initPositions: Particle number exceeds pt_maxPerProc. Increase.')
        ps = pt_numLocal
        particles(:,ps) = 0
        particles(PROC_PART_PROP,ps) = pt_meshMe
        particles(BLK_PART_PROP,ps)  = blockID
        particles(POSX_PART_PROP,ps) = xcoord
        particles(POSY_PART_PROP,ps) = ycoord
        particles(POSZ_PART_PROP,ps) = zcoord
     endif
  enddo
  call Particles_getGlobalNum(globalNumParticles)
  print*, "Done placing. Exiting Particles_initPositions"
  CALL flush()

  partPosInitialized = .true.
  pt_posInitialized = partPosInitialized
  return

end subroutine pt_initVoronoiPositions
