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
       pt_posInitialized, pt_meshMe, pt_globalMe, pt_meshComm
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_getCellCoords,&
       Grid_getBlkIndexLimits, Grid_getBlkData, &
       Grid_getBlkIDFromPos
  
  !use Simulation_data, ONLY : sim_ptMass, sim_densityThreshold, sim_print

  !USE HDF5
  implicit none

#include "constants.h"
#include "Flash.h"

  logical, intent(INOUT) :: partPosInitialized
  logical, intent(OUT) :: updateRefine

  integer       :: i, j, k, istat
  integer       :: procID, blockID, blkCount, ps, pe
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
  INTEGER(KIND=4) :: error      ! error flag
  INTEGER         :: cols, rows, nrows ! number of columns and rows of dataset
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords_dset_data ! Dataset memory

  real :: xcoord, ycoord, zcoord ! Imported gas particle coordinates
!--------------------------------------------------------------

  call Driver_getComm(GLOBAL_COMM, pt_meshComm)
  call Driver_getMype(GLOBAL_COMM, pt_meshMe)
  if (pt_meshMe .eq. pt_globalMe) then 
     open(unit = 55, file = "test-txt.txt", status = 'old')
     do
        read(55,*,end=999)
        nrows = nrows+1
     end do
999  continue
     close(55)
  endif
  if (pt_meshMe .eq. pt_globalMe) then
     open(unit=55, file="test-txt.txt", status='old')
     ALLOCATE(coords_dset_data(3,nrows))
     print*, "VorAMR input particles: ", nrows
     read(55,*,end=9999) coords_dset_data(:,:)
9999 continue
     close(55)
  endif
  
  if (pt_meshMe.EQ.MASTER_PE) then
     call Logfile_stamp('VorAMR. Building refined grid', "[pt_initVoronoiPositions]")
  endif
  
  updateRefine=.true. !.false.
  if(partPosInitialized) return

  pt_numLocal=0
  
  print*, 'Placing particles...'
  CALL flush()
  do i=1, nrows
     if (pt_meshMe .eq. pt_globalMe) then
        xcoord = coords_dset_data(1,i)
        ycoord = coords_dset_data(2,i)
        zcoord = coords_dset_data(3,i)
        call Grid_getBlkIDFromPos([xcoord, ycoord, zcoord], blockID, procID, pt_meshComm)
        pt_numLocal = pt_numLocal + 1
        if (pt_numLocal .gt. pt_maxPerProc) &
             call Driver_abortFlash('Particles_initPositions: Particle number exceeds pt_maxPerProc. Increase.')
        ps = pt_numLocal
        particles(:,ps) = 0
        particles(PROC_PART_PROP,ps) = procID
        particles(BLK_PART_PROP,ps)  = blockID
        particles(POSX_PART_PROP,ps) = xcoord
        particles(POSY_PART_PROP,ps) = ycoord
        particles(POSZ_PART_PROP,ps) = zcoord
        particles(TYPE_PART_PROP,ps) = NONEXISTENT
        !NONEXISTENT particle type allows VorAMR particles and their memory to be cleared by Grid_sortParticles() and reused.
     endif
  enddo
  !print*, "Broadcasting to other procs."
  !call flush()
  !call MPI_Bcast(particles, MDIM, FLASH_INTEGER, pt_globalMe, pt_meshComm, istat)
  print*, "Out of placing loop, getting global num particles."
  call flush()
  call Particles_getGlobalNum(globalNumParticles)

  print*, "Deallocating particle coordinate dataset."
  call flush()
  DEALLOCATE(coords_dset_data)

  print*, "Done placing. Exiting Particles_initPositions"
  CALL flush()

  partPosInitialized = .true.
  pt_posInitialized = partPosInitialized
  return

end subroutine pt_initVoronoiPositions
