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

use Grid_data, ONLY : gr_eosMode

#ifdef GRAVITY
use Gravity_interface, ONLY: Gravity_accelOneBlock, Gravity_getAccelAtPoint, &
    Gravity_getPotentialAtPoint, Gravity_potentialListOfBlocksPDENonly, &
    Gravity_potentialListOfBlocks, Gravity_accelListOfBlocks

use Gravity_data, ONLY : useGravity
#endif

#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
use Particles_interface, ONLY : Particles_sinkSyncWithParticles, &
    Particles_getGlobalNum, Particles_mapFromMesh, Particles_sinkMoveParticles, &
    Particles_moveAndSort, Particles_longRangeForce, &
    Particles_sinkSortParticles, Particles_addNew
#endif

use Timers_interface, ONLY : Timers_start, Timers_stop

use IO_interface, ONLY : IO_output

use IO_data, ONLY: io_checkpointFileNumber, io_plotFileNumber, &
                   io_rollingCheckpoint

#ifdef ENERGY_INJ
use Particles_interface, ONLY : Particles_energyInjection
#endif

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

! Requiste header files from Flash. Note we use the mpi.h compiled
! with Flash here. This assures that we don't have an MPI mismatch.


!#ifndef MPI_INCLUDED
!#include "Flash_mpi.h"
!#define MPI_INCLUDED
!#endif

logical, save :: restart
character(len=16), save :: data_type, data_var

real*8, save :: force_GoS_x, force_GoS_y, force_GoS_z
real*8, save :: force_SoG_x, force_SoG_y, force_SoG_z


! Pointer to particles array which can be massive or sink particles in Flash

real*8, pointer, save, dimension(:,:) :: particles_pointer
integer, pointer :: num_part_local_ptr

! Pointer for array that captures when new particles are made in Flash to
! pass their tags on to AMUSE.
integer*8, pointer, save, dimension(:) :: new_particles_tags
integer, pointer :: number_new_particles

character(len=4), save :: part_type


contains


FUNCTION set_particle_pointers(part_type_in)

use Driver_data, only : dr_globalComm
#ifdef ACTIVE_PART_TYPE
    use Particles_data, only : particles, pt_numLocal, new_massive_tags, number_new_massive
#endif
#ifdef SINK_PART_TYPE
    use Particles_sinkData, only : particles_local, localnp, new_sink_tags, number_new_sinks
#endif

integer   :: set_particle_pointers, ierr
character(len=4), intent(in) :: part_type_in

#if defined (ACTIVE_PART_TYPE) || defined (SINK_PART_TYPE)
! Here we check if the simulation contains massive particles.
! If yes, it assumes if we are also using sinks they are to gather gas
! and make massive particles and that the massive particles are the ones
! that contribute feedback.

! Set an MPI barrier here to prevent some processors from pointing at
! different particle arrays at the same time.

call MPI_Barrier(dr_globalComm, ierr)

! Point at active particles.
if (part_type_in == 'mass') then

#if defined (ACTIVE_PART_TYPE)

    particles_pointer => particles
    num_part_local_ptr => pt_numLocal

    new_particles_tags => new_massive_tags
    number_new_particles => number_new_massive

    part_type = "mass"

    !if (dr_globalMe == 0) print*, "[set_particle_pointers]: Particle pointers are set to ACTIVE_PART_TYPE!"
#else

    print*, "[set_particle_pointers]: Tried to set mass type but massive &
             particles are not compiled in!"
    call Driver_abortFlash("Tried to set invalid particle pointer.")
#endif

else if (part_type_in == 'sink') then

#if defined (SINK_PART_TYPE)
! Point at sink particles.
    particles_pointer => particles_local
    num_part_local_ptr => localnp

    new_particles_tags => new_sink_tags
    number_new_particles => number_new_sinks

    part_type = "sink"

    !if (dr_globalMe == 0) print*, "[set_particle_pointers]: Particle pointers are set to SINK_PART_TYPE!"

#else

    print*, "[set_particle_pointers]: Tried to set sink type but sink &
             particles are not compiled in!"
    call Driver_abortFlash("Tried to set invalid particle pointer.")

#endif

else

    print*, "[set_particle_pointers]: Invalid part_type_in passed."
    call Driver_abortFlash("Tried to set invalid particle pointer.")

end if
#endif
set_particle_pointers = 0

call flush(6)
END FUNCTION

FUNCTION internal_particle_integration_off()

#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
use Particles_data, only : pt_typeInfo

integer internal_particle_integration_off

! Switch off all integrators for the particles array.
if (dr_globalMe == 0) print*, "[initialize_particle_pointers]: Warning, &
                      switching off all Flash internal integrators for &
                      the particles array."

pt_typeInfo(PART_ADVMETHOD,:) = PT_ADVMETH_NONE
#endif
internal_particle_integration_off=0
END FUNCTION

subroutine get_particle_type_bounds(part_type_in, type_begin, type_end, type_count)
! Get the indicies boundaries of type of particle requested in the particles
! array.

! NOTE: This subroutine assumes the user properly updated the needed array
! using whatever is appropriate, i.e. pt_updateTypeDS for multi-type sims with
! massive or Grid_sinksSortParticles() for sinks, etc.

#if defined (ACTIVE_PART_TYPE) || (SINK_PART_TYPE)
    use Particles_data, only : pt_numLocal, pt_typeInfo
#endif
#ifdef SINK_PART_TYPE
    use Particles_sinkData, only : localnp
#endif


character(len=4), intent(in) :: part_type_in
integer, intent(inout)       :: type_begin, type_end, type_count

!if (dr_globalMe .eq. 0) print*, "part_type_in ", part_type_in

if (part_type_in .eq. 'mass') then

#if defined (ACTIVE_PART_TYPE) && defined (SINK_PART_TYPE)
type_begin = pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE)
type_count = pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE)
type_end   = type_count + type_begin - 1

!print*, "type_begin", type_begin, dr_globalMe
!print*, "type_end", type_end, dr_globalMe
!print*, "type_count", type_count, dr_globalMe

#elif defined (ACTIVE_PART_TYPE) && (NPART_TYPES == 1)

type_begin = 1
type_end   = pt_numLocal
type_count = pt_numLocal

#else

call Driver_abortFlash('[get_particle_type_bounds]: massive particles not in simulaton!')

#endif
else if (part_type_in .eq. 'sink') then

#ifdef SINK_PART_TYPE
  type_begin = 1
  type_end   = localnp
  type_count = localnp

#else

 call Driver_abortFlash('[get_particle_type_bounds]: sink particles not in simulation!')

#endif

else

  call Driver_abortFlash('[get_particle_type_bounds]: part_type not found!')

endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Runtime Parameters Settings
!!!!   These must be set before Flash
!!!!   is initialized (which is when Flash
!!!!   grabs the RT parameters info.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION initialize_restart()

  INTEGER :: initialize_restart
  restart = .false.
  initialize_restart=0
END FUNCTION

FUNCTION get_restart(value)

  INTEGER :: get_restart
  LOGICAL :: value
  call RuntimeParameters_get('restart',value)
  get_restart=0
END FUNCTION

FUNCTION set_restart(value)

  INTEGER :: set_restart
  LOGICAL :: value
  call RuntimeParameters_set('restart',value)
  restart=value
  set_restart=0
END FUNCTION

FUNCTION get_begin_iter_step(value)

  INTEGER :: value
  INTEGER :: get_begin_iter_step
  value = dr_nbegin
  get_begin_iter_step=0
END FUNCTION

FUNCTION set_begin_iter_step(value)

  INTEGER :: value
  INTEGER :: set_begin_iter_step
  dr_nbegin = value
  set_begin_iter_step=0
END FUNCTION

FUNCTION get_max_num_steps(value)

  INTEGER :: value
  INTEGER :: get_max_num_steps
  value = dr_nend
  get_max_num_steps=0
END FUNCTION

FUNCTION set_max_num_steps(value)

  INTEGER :: value
  INTEGER :: set_max_num_steps
  dr_nend = value
  set_max_num_steps=0
END FUNCTION

FUNCTION get_current_step(value)

  INTEGER :: value
  INTEGER :: get_current_step
  value = dr_nstep
  get_current_step=0
END FUNCTION

FUNCTION get_time(value)

  DOUBLE PRECISION :: value
  INTEGER :: get_time
  call Driver_getSimTime(value)
  get_time=0
END FUNCTION

FUNCTION get_end_time(value)

  DOUBLE PRECISION :: value
  INTEGER :: get_end_time
  value = dr_tmax
  get_end_time=0
END FUNCTION

FUNCTION get_timestep(value)

  DOUBLE PRECISION :: value
  INTEGER :: get_timestep
!  call RuntimeParameters_get("dtinit",value)
!!! Here it makes more sense to get it in Driver
!!! also, since a sim is already likely running.

  value = dr_dt  ! changed from dr_dtAdvect - AT, 2019 nov 26
  !value = dr_dtAdvect
  !value = min(dr_dtMin,dr_dtAdvect)
  get_timestep=0
END FUNCTION

FUNCTION set_timestep(value)

  DOUBLE PRECISION :: value
  INTEGER :: set_timestep
  !call RuntimeParameters_set("dtinit",value)
  dr_dt = value  ! This isn't working as intended currently. - JW
!!! Here it makes more sense to set it in Driver
!!! also, since a sim is already likely running.

  dr_dtInit = value
  set_timestep=0
END FUNCTION

FUNCTION set_end_time(value)

  DOUBLE PRECISION :: value
  INTEGER :: set_end_time
  call RuntimeParameters_set('tmax',value)
  dr_tmax = value
!  call Driver_init()
  set_end_time=0
END FUNCTION




! This function returns the coords along a single dimension of the
! block. Limits is the number of cells in that dimension, and should
! be the same as nparts. Note that axis must be an array of length
! nparts for python to give back an array. Axis is i=1, j=2, k=3.

! Currently implemented for 1 proc only.

FUNCTION get_1blk_cell_coords(axis, blockID, procID, limits, coords, nparts)

integer :: get_1blk_cell_coords
integer :: axis, blockID, limits, procID
integer :: nparts
real*8  :: coords(nparts)
  integer :: communicator, myProc, ierr, i

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  if (myProc .eq. procID) then
          call Grid_getCellCoords(axis,blockID,CENTER,.false.,coords, nparts)
  endif

  ! Gather coords from procID to all other processes
  call MPI_Bcast(coords, nparts, FLASH_REAL, procID, communicator, ierr)

get_1blk_cell_coords=0
END FUNCTION

!FUNCTION get_all_1axis_cell_coors(coords,nparts)
!integer :: get_all_1axis_cell_coors


!get_all_1axis_cell_coors=0
!END FUNCTION

!!! Here for evolve model we set tmax to the evolve time,
!!! intialize Flash, hijack the main evolution loop
!!! in Flash then set it for restart such that any evolve call
!!! after is a restart for the loop.

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
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Grid non-variable operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION get_cell_volume(block, i, j, k, vol)

  INTEGER :: block, i, j, k, get_cell_volume
  REAL*8  :: vol

  call Grid_getSingleCellVol(block, INTERIOR, [i,j,k], vol)

get_cell_volume=0
END FUNCTION

! Get the number of processors which have blocks on them.
FUNCTION get_number_of_procs(n)

  integer :: n, get_number_of_procs
  ! Note MESH_COMM is a Flash defined constant.
  call Driver_getNumProcs(MESH_COMM, n)
  get_number_of_procs=0
END FUNCTION

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
END FUNCTION

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
END FUNCTION

FUNCTION set_data_type(type_in)

integer :: set_data_type
character(len=16) :: type_in
data_type=type_in
set_data_type=0
END FUNCTION

FUNCTION set_data_var(var_in)

integer :: set_data_var
character(len=16) :: var_in
data_var=var_in
set_data_var=0
END FUNCTION

FUNCTION grid_update_refinement(gridChanged)

use Driver_data, only : dr_nstep, dr_simTime, dr_simGeneration

integer :: grid_update_refinement
logical, intent(out) :: gridChanged

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")
     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

grid_update_refinement=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get Hydro State interpolates the data
!!! at any point. Currently broken :(
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION get_hydro_state_at_point(x, y, z, vx, vy, vz, &
                         rho, rhovx, rhovy, rhovz, rhoen)

  REAL*8 :: x, y, z, vx, vy, vz
  REAL*8 :: x1(3), x2(3), x3(3)
  REAL*8, DIMENSION(MDIM) :: delta, blockCenter, &
                             local_pos, loc, pos_in_cell, blockSize
  REAL*8 :: weight(3,3,3)
  REAL*8 :: rho, rhovx, rhovy, rhovz, rhoen
  REAL*8 :: rho_cell, rhovx_cell, rhovy_cell, rhovz_cell, en_cell
  INTEGER :: get_hydro_state_at_point
  INTEGER :: i, j, k, blockID, blkLimits(2,MDIM), &
             blkLimitsGC(2,MDIM), Proc_ID, communicator
  INTEGER :: ii, jj, kk
  rho=0.0; rhovx=0.0; rhovy=0.0; rhovz=0.0; rhoen = 0.0
  loc(1)=x; loc(2)=y; loc(3)=z
  x1=0.0; x2=0.0; x3=0.0
  weight=1.0; local_pos=0.0; pos_in_cell=0.0


! Find the i,j,k index of the cell for the queried point x,y,z.

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Grid_getBlkIDFromPos([x,y,z], blockID, Proc_ID, communicator)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getDeltas(blockID, delta)

! Make sure that the guard cells in the 1st layer around the whole block
! are filled with centered values that have been checked against the EOS
! to be self consistent. Note eosMode is set to be the default mode
! already being used from Grid_data.

  call Grid_fillGuardCells(CENTER, ALLDIR, minlayers=1, eosMode=gr_eosMode, doEos=.true.)

do ii=1, MDIM
  local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0 ! Position in the block.
end do

! Index of the cell at the local position in the block, including guard cells.

  i = ceiling(local_pos(1)/delta(1)) + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))

  if (MDIM .gt. 1) then
      j = ceiling(local_pos(2)/delta(2)) + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
  else
      j = 0
  end if

  if (MDIM .gt. 2) then
      k = ceiling(local_pos(3)/delta(3)) + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
  else
      k = 0
  end if

! Where the query point x,y,z is in the local cell at i,j,k where cell size
! is normed to 1 in each dimension and the origin is in the center.

  pos_in_cell =  mod(local_pos,delta)/delta - 0.5

! Weighting towards the inner and outer boundaries of the cell on each
! axis for every cell that borders the cell containing the point. Note that
! if the point is closer to one side, the cells on the opposite side give no
! weight to the interpolation. 1 is cell "to the left" of cell with point,
! 2 is the cell with the point, and 3 is cell "to the right" of cell with point.

  do ii=1,MDIM

    if (pos_in_cell(ii) > 0.0) then
      x1(ii) = 0.0
      x2(ii) = 1.0 - pos_in_cell(ii)
      x3(ii) = pos_in_cell(ii)

    else if (pos_in_cell(ii) < 0.0) then
      x1(ii) = -pos_in_cell(ii)
      x2(ii) = 1.0 + pos_in_cell(ii)
      x3(ii) = 0.0

    else
      x1(ii) = 0.0
      x2(ii) = 1.0
      x3(ii) = 0.0

    end if

  end do

  weight(1,:,:) = weight(1,:,:)*x1(1)
  weight(2,:,:) = weight(2,:,:)*x2(1)
  weight(3,:,:) = weight(3,:,:)*x3(1)
  weight(:,1,:) = weight(:,1,:)*x1(2)
  weight(:,2,:) = weight(:,2,:)*x2(2)
  weight(:,3,:) = weight(:,3,:)*x3(2)
  weight(:,:,1) = weight(:,:,1)*x1(3)
  weight(:,:,2) = weight(:,:,2)*x2(3)
  weight(:,:,3) = weight(:,:,3)*x3(3)


!  do ii=1,3
!    do jj=1,3
!      do kk=1,3

!        weight(ii,jj,kk) = x1(ii)*x2(jj)*x3(kk)

!      end do
!    end do
!  end do


  do ii=1, 3
    do jj=1, 3
      do kk=1, 3

        call Grid_getPointData(blockID,CENTER,DENS_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rho_cell)
        rho  = rho + weight(ii,jj,kk)*rho_cell
        call Grid_getPointData(blockID,CENTER,VELX_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovx_cell)
        rhovx = rhovx + weight(ii,jj,kk)*rhovx_cell
        call Grid_getPointData(blockID,CENTER,VELY_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovy_cell)
        rhovy = rhovy + weight(ii,jj,kk)*rhovy_cell
        call Grid_getPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],rhovz_cell)
        rhovz = rhovz + weight(ii,jj,kk)*rhovz_cell
        call Grid_getPointData(blockID,CENTER,ENER_VAR,EXTERIOR,[i-2+ii,j-2+jj,k-2+kk],en_cell)
        rhoen = rhoen + weight(ii,jj,kk)*en_cell

      end do
    end do
  end do
  get_hydro_state_at_point = 0
END FUNCTION

FUNCTION get_potential(i, j, k, index_of_grid, potential)

  integer :: i, j, k, index_of_grid, get_potential
  real*8  :: potential
#ifdef GRAVITY
  call Grid_getPointData(index_of_grid, CENTER, GPOT_VAR, INTERIOR, [i, j, k], potential)
#endif
  get_potential=0
END FUNCTION get_potential

!!! Currently implementing!!! -JW 12-15-14
!!! Seems to be working properly now !!!  -JW 12-19-14

!FUNCTION get_potential_at_point(eps, x, y, z, potential)
!
!integer   :: get_potential_at_point
!integer   :: ProcID, communicator, blockID
!integer, parameter :: n_attrib=1
!real*8, dimension(LOW:HIGH,MDIM) :: bndbox
!integer, dimension(2,n_attrib) :: attrib
!real*8    :: eps, x, y, z
!real*8, dimension(n_attrib) :: potential
!real*8, dimension(MDIM)  :: loc, deltas
!real*8, pointer, dimension(:,:,:,:) :: solndata

!  loc = (/x, y, z/)
!  attrib(1,1) = GPOT_PART_PROP
!  attrib(2,1) = GPOT_VAR

!!print*, "Getting potential at ", loc

!  call Driver_getComm(GLOBAL_COMM, communicator)
!  call Grid_getBlkIDFromPos(loc, blockID, ProcID, communicator)
!  call Grid_getBlkPtr(blockID,solndata,CENTER)
!  call Grid_getBlkBoundBox(blockID,bndbox)
!  call Grid_getDeltas(blockID,deltas)
!  call Particles_mapFromMesh(1, 1, attrib, loc, bndbox, &
!                               deltas, solndata, potential)
!  call Grid_releaseBlkPtr(blockID,solndata)


!get_potential_at_point=0
!END FUNCTION

!!! This version of get_gravity matches the "spirit" of AMUSE's
!!! idea of this call. However since we want to do as few calls
!!! as possible in our production code, we use a verison that
!!! also updates the particle locations in Flash, then uses
!!! the tree to calculate this by direct sum (more acc than
!!! the finite differencing used below).

!FUNCTION get_gravity_at_point(eps, x, y, z, gax, gay, gaz, nparts)
!
!integer   :: get_gravity_at_point, i, nparts, blkLimits(2,MDIM), &
!             blkLimitsGC(2,MDIM), ii, jj, kk, iii
!integer   :: ProcID, communicator, blockID, myProc, ierr, sts
!integer, parameter :: n_attrib=3
!real*8, dimension(LOW:HIGH,MDIM) :: bndbox
!integer, dimension(2,n_attrib) :: attrib
!real*8, dimension(nparts)   :: x, y, z, gax, gay, gaz
!real*8    :: eps, error
!real*8, dimension(nparts,MDIM) :: gravity
!real*8, dimension(MDIM)  :: loc, deltas, blockSize, blockCenter, local_pos
!real*8, dimension(MDIM)  :: cell_loc
!real*8, dimension(MDIM,2):: grav_cell, cell_locs
!real*8, dimension(2)     :: x_cell, y_cell, z_cell
!real*8, dimension(MDIM, GRID_IHI_GC, GRID_JHI_GC, GRID_KHI_GC) :: gvec
!real*8, pointer, dimension(:,:,:,:) :: solndata
!integer, save :: numcalls=0
!
!!!!!! NOTE: This function REQUIRES the user to define the variable for
!!!!!! graviational acceleration in each direction on the grid as
!!!!!! GACX, GACY, and GACZ in the Config file. Note this is done if using
!!!!!! BHTree gravity and you do bhtreeAcc=1 on setup line. See Flash
!!!!!! users guide for how to define variables for other gravity solvers.
!
!!!!!! ALSO NOTE: The Multigrid gravity solver will give you nonsense back
!!!!!! if the number of guard cells is less than four. I found out the
!!!!!! hard way. Don't be like me.
!
!!!!!! FINAL NOTE: You MUST include the unit Particles/ParticlesMapping/
!!!!!! Quadratic in your setup for this function to work. If you don't
!!!!!! this function will return zero gravity!
!
!!print*, "Now we're in get_gravity_at_point."
!
!call Grid_fillGuardCells(CENTER, ALLDIR, minLayers=NGUARD, &
!                         eosMode=gr_eosMode, doEos=.false.)
!
!attrib(1,1) = ACCX_PART_PROP
!attrib(2,1) = GACX_VAR
!attrib(1,2) = ACCY_PART_PROP
!attrib(2,2) = GACY_VAR
!attrib(1,3) = ACCZ_PART_PROP
!attrib(2,3) = GACZ_VAR
!
!gravity = 0.0
!
! do  i=1, nparts
!
!
!
!  loc = (/x(i), y(i), z(i)/)
!
!!  numcalls = numcalls + 1
!
!!  print*, "Call # ", numcalls
!!  print*, "Getting gravity at ", loc
!
!  call Driver_getMype(GLOBAL_COMM, myProc)
!  call Driver_getComm(GLOBAL_COMM, communicator)
!  call Grid_getBlkIDFromPos(loc, blockID, ProcID, communicator)
!
!  if (myProc .eq. ProcID) then
!
!!!!! This method uses finite differencing on the grid potential.
!
!!    print*, "blockID = ", blockID
!
!!!!! This line likely only needs to be called for the Multigrid solver.
!!    call Gravity_accelOneBlock(blockID,NGUARD,gvec)
!
!
!!!    call Grid_getBlkPtr(blockID,solndata,CENTER)
!
!!!    solndata(GRAX_VAR,:,:,:) = gvec(1,:,:,:)
!!!    solndata(GRAY_VAR,:,:,:) = gvec(2,:,:,:)
!!!    solndata(GRAZ_VAR,:,:,:) = gvec(3,:,:,:)
!
!!!!! This method assumes we already calculated g accel when we
!!!!! called Grid_solvePoisson sometime earlier. Note we might have to
!!!!! now call Grid_solvePoisson during driver_init. - JW
!
!    call Grid_getBlkPtr(blockID,solndata,CENTER)
!    call Grid_getBlkBoundBox(blockID,bndbox)
!    call Grid_getDeltas(blockID,deltas)
!
!
!!!    print*, "Cell gravity for block is: ", solndata(GACX_VAR,:,:,:) &
!!!            , solndata(GACY_VAR,:,:,:), solndata(GACZ_VAR,:,:,:)
!
!    call Particles_mapFromMesh(1, 3, attrib, loc, bndbox, &
!                                 deltas, solndata, gravity(i,:))
!
!
!!!    Grid_mapMeshToParticles (particles, part_props,&
!!!                                    numParticles,posAttrib,&
!!!                                    numAttrib, attrib,&
!!!                                    mapType,gridDataStruct)
!
!
!    call Grid_releaseBlkPtr(blockID,solndata)
!
!!    call Grid_getBlkCenterCoords(blockID, blockCenter)
!!    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
!!    call Grid_getBlkPhysicalSize(blockID, blockSize)
!!    call Grid_getDeltas(blockID, deltas)
!
!
!
!
!
!
!!    ! Position in the block.
!!    do iii=1, MDIM
!!      local_pos(iii) = loc(iii) - blockCenter(iii) + blockSize(iii)/2.0
!!    end do
!
!!  ! Find the cell the loc is in. We use the local pos in the block including
!!  ! guard cells.
!
!!    ii = ceiling(local_pos(1)/deltas(1))  &
!!       + (blkLimits(LOW,IAXIS)-blkLimitsGC(LOW,IAXIS))
!
!!    if (MDIM .gt. 1) then
!!        jj = ceiling(local_pos(2)/deltas(2))  &
!!           + (blkLimits(LOW,JAXIS)-blkLimitsGC(LOW,JAXIS))
!!    else
!!        jj = 0
!!    end if
!
!!    if (MDIM .gt. 2) then
!!        kk = ceiling(local_pos(3)/deltas(3))  &
!!           + (blkLimits(LOW,KAXIS)-blkLimitsGC(LOW,KAXIS))
!!    else
!!        kk = 0
!!    end if
!
!!    ! Now we determine if the point is to the left or right of the
!!    ! cell center of the cell it resides in. Then we get the cell centers
!!    ! to the left and right of the point for linear interpolation and
!!    ! the gravity to the left and right of the point in each
!!    ! dimension.
!
!!    if (nint(mod(local_pos(1),deltas(1))) .eq. 1) then
!
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(1,1) = cell_loc(1)
!!        call Grid_getSingleCellCoords([ii+1, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(1,2) = cell_loc(1)
!
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(1,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii+1, jj, kk],grav_cell(1,2))
!!    else
!
!!        call Grid_getSingleCellCoords([ii-1, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(1,1) = cell_loc(1)
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(1,2) = cell_loc(1)
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii-1, jj, kk],grav_cell(1,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(1,2))
!!    end if
!
!!    if (nint(mod(local_pos(2),deltas(2))) .eq. 1) then
!
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(2,1) = cell_loc(2)
!!        call Grid_getSingleCellCoords([ii, jj+1, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(2,2) = cell_loc(2)
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(2,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj+1, kk],grav_cell(2,2))
!!    else
!
!!        call Grid_getSingleCellCoords([ii, jj-1, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(2,1) = cell_loc(2)
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(2,2) = cell_loc(2)
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj-1, kk],grav_cell(2,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(2,2))
!!    end if
!
!!    if (nint(mod(local_pos(3),deltas(3))) .eq. 1) then
!
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(3,1) = cell_loc(3)
!!        call Grid_getSingleCellCoords([ii, jj, kk+1],blockID,CENTER, &
!!                                  EXTERIOR, cell_loc)
!!        cell_locs(3,2) = cell_loc(3)
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(3,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk+1],grav_cell(3,2))
!!    else
!
!!        call Grid_getSingleCellCoords([ii, jj, kk-1],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(3,1) = cell_loc(3)
!!        call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!                                      EXTERIOR, cell_loc)
!!        cell_locs(3,2) = cell_loc(3)
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk-1],grav_cell(3,1))
!!        call Grid_getPointData(blockID,CENTER,GACX_VAR,EXTERIOR, &
!!                               [ii, jj, kk],grav_cell(3,2))
!!    end if
!
!!!    call Grid_getBlkPtr(blockID,solndata,CENTER)
!
!!!    grav_cell(1,iii) = solndata(GACX_VAR,ii,jj,kk)
!!!    grav_cell(2,iii) = solndata(GACY_VAR,ii,jj,kk)
!!!    grav_cell(3,iii) = solndata(GACZ_VAR,ii,jj,kk)
!
!!! Get the locations of the cells for the interpolation.
!!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(1,1))
!!!    call Grid_getSingleCellCoords([ii+1, jj, kk],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(1,2))
!!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(2,1))
!!!    call Grid_getSingleCellCoords([ii, jj+1, kk],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(2,2))
!!!    call Grid_getSingleCellCoords([ii, jj, kk],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(3,1))
!!!    call Grid_getSingleCellCoords([ii, jj, kk+1],blockID,CENTER, &
!!!                                  EXTERIOR, cell_locs(3,2))
!
!!!    call Grid_releaseBlkPtr(blockID, solndata)
!
!!    ! Finally we linearly interpolate the gravity between the values.
!
!!    do iii=1,MDIM
!
!!        call ut_polint(cell_locs(iii,:), grav_cell(iii,:), 2, &
!!                       loc(iii), gravity(iii), error)
!
!!    end do
!
!!    gax(i) = gravity(1); gay(i) = gravity(2); gax(i) = gravity(3)
!
!!    print*, "Local cell gravity:"
!!    print*, grav_cell
!!    print*, "Interpolated gravity:"
!!    print*, gravity
!
!!    if (MyProc .ne. 0) then
!
!!      call MPI_SEND(gravity, 3, MPI_DOUBLE_PRECISION, 0, 1, communicator, ierr)
!
!!      print*, "Sent ", gravity
!
!!    end if
!
!  end if
!
!!  if (myProc .ne. ProcID .and. myProc .eq. 0) then
!
!!    call MPI_RECV(gravity, 3, MPI_DOUBLE_PRECISION, ProcID, MPI_ANY_TAG, communicator, sts, ierr)
!
!!    print*, "Recieved ", gravity
!
!!  end if
!
!!  if (myProc .eq. 0) then
!
!!    print*, "I'm zero and grav is ", gravity
!
!!  end if
!
!!  gax(i) = gravity(1); gay(i) = gravity(2); gaz(i) = gravity(3)
!
!end do
!
!!print*, "Out of loop."
!
!call MPI_Reduce(gravity(:,1), gax, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
!call MPI_Reduce(gravity(:,2), gay, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
!call MPI_Reduce(gravity(:,3), gaz, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
!
!!print*, "Info sent."
!
!!print*, "Gravity array = ", gravity, myProc
!!print*, "gax =", gax, myProc
!
!get_gravity_at_point=0
!END FUNCTION


FUNCTION get_potential_at_point(eps, x, y, z, gpot, nparts)
INTEGER :: nparts, i, get_potential_at_point, MyPe, blockID, ProcID
INTEGER :: communicator, ierror
REAL*8, DIMENSION(nparts)  :: eps, x, y, z, gpot, locpot
!LOGICAL :: has_data

#ifdef GRAVITY

!has_data = .false.
locpot = 0.0
gpot = 0.0

#ifdef TREE

call Driver_getMype(GLOBAL_COMM, MyPe)
call Driver_getComm(GLOBAL_COMM, communicator)

do i=1, nparts

  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], blockID, ProcID, communicator)

  if (MyPe .eq. ProcID) then

    call Gravity_getPotentialAtPoint(x(i), y(i), z(i), locpot(i))
!   has_data = .true.
!    print*, "Local potential = ", locpot(i)

  end if

end do

!if (has_data) then

  call MPI_Reduce(locpot, gpot, nparts, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierror)

!end if

#else

do i=1, nparts

  call Gravity_getPotentialAtPoint(x(i), y(i), z(i), gpot(i))

end do

#endif

#else

print*, "ERROR! Gravity unit not compiled in!"
stop
#endif

get_potential_at_point=0
END FUNCTION


FUNCTION get_accel_gas_on_particles(eps, x, y, z, gax, gay, gaz, nparts)
INTEGER :: nparts, i, get_accel_gas_on_particles, MyPe, blockID, ProcID
INTEGER :: communicator, ierror
REAL*8, DIMENSION(nparts)  :: eps, x, y, z, gax, gay, gaz
REAL*8, DIMENSION(nparts)  :: local_gx, local_gy, local_gz

LOGICAL, PARAMETER :: Debug=.false.

#ifdef GRAVITY

if (.not. useGravity) return

!has_data = .false.
local_gx=0.0; local_gy=0.0; local_gz=0.0
gax = 0.0; gay = 0.0; gaz = 0.0
force_GoS_x=0.0; force_GoS_y=0.0; force_GoS_z=0.0


! Note TREE gravity has a local tree of ALL cells on each proc
! so you only run this on the proc with the point in it.

! But multigrid only has local blocks on each proc so its run on all
! and then internally returns the summed gravity.

#ifdef TREE
call Driver_getMype(GLOBAL_COMM, MyPe)
call Driver_getComm(GLOBAL_COMM, communicator)

do i=1, nparts

  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], blockID, ProcID, communicator)

  if (MyPe .eq. ProcID) then

    call Gravity_getAccelAtPoint(x(i), y(i), z(i), local_gx(i), local_gy(i), local_gz(i))
!    has_data = .true.

  end if

end do



!print*, "On proc ", MyPe, "grav is", local_gx, local_gy, local_gz

!if (has_data) then

!call MPI_Reduce(local_gx, gax, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!call MPI_Reduce(local_gy, gay, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!call MPI_Reduce(local_gz, gaz, nparts, MPI_DOUBLE_PRECISION, &
!                  MPI_SUM, 0, communicator, ierror)
!end if

#else
do i=1, nparts

    call Gravity_getAccelAtPoint(x(i), y(i), z(i), gax(i), gay(i), gaz(i))

end do
#endif

! Calculate the forces acting on the particles from the gas.
! For debugging purposes.
if (Debug .and. (dr_globalMe .eq. MASTER_PE)) then
do i=1, nparts

    force_GoS_x = force_GoS_x + gax(i)*particles_global(MASS_PART_PROP,i)
    force_GoS_y = force_GoS_y + gay(i)*particles_global(MASS_PART_PROP,i)
    force_GoS_z = force_GoS_z + gaz(i)*particles_global(MASS_PART_PROP,i)

end do

     if (dr_globalMe .eq. MASTER_PE) &
       & write(*,'(A,4(1X,E17.10))') 'Particles_sinkAccelGasOnSinks: Total force GAS->SINKS (time, x,y,z) = ', &
       & dr_simTime, force_GoS_x, force_GoS_y, force_GoS_z
end if
! No gravity
#else
     print*, "ERROR! Gravity unit not compiled in!"
     stop
#endif


get_accel_gas_on_particles=0
END FUNCTION

!!! These two functions represent the latest bridge work and should
!!! be faster than the direct summation bridge. JW 12-20-15.

!!! This magical routine Particles_longRangeForce calculates
!!! the gravitational acceleration and potential on the particles
!!! by CIC mapping of these quantities to the particles and putting
!!! the values in GPOT_PART_PROP,ACCX_PART_PROP,etc. Note the
!!! acceleration comes from a finite differencing of the potential using
!!! Gravity_accelListOfBlocks. After calling this function, we can just
!!! get the accelerations by calling the tags.

!!! Note WEIGHTED map type gives the traditional CIC linear weighted mapping.

FUNCTION get_gravity_gas_on_particles(dt, kick_number)
!#define debug_force
integer :: get_gravity_gas_on_particles, p, blockCount
integer :: blocks(MAXBLOCKS), kick_number, numAttrib, attrib(2,3)
integer :: type_begin, type_end, type_count, ierr
real*8 :: dt, sum_force_norm
real*8, allocatable :: accx(:), accy(:), accz(:)
!logical, parameter :: debug= .true. !.false.
character(len=20) :: file_name
logical :: sink_separate = .false.

#ifdef GRAVITY

if (.not. useGravity) return
!print*, "Kicking particles with bridge."
#ifdef debug_force
force_GoS_x=0.0; force_GoS_y=0.0; force_GoS_z=0.0
#endif
!if (localnp .ne. 0) then

!print*, "In gas->particles, Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(BLK_PART_PROP,1:localnp)

!end if
! May not need to call this here if we call it in particles_on_gas.
!call Particles_advance(1.0, 1.0)
!call pt_sinkGatherGlobal()
!call Particles_sinkSyncWithParticles(sink_to_part=.true.)

if (kick_number == 2) then

  ! Make sure the accelerations for gravity are filled in the guard cells
  ! since the longRangeForce subroutine interpolates gravity.
    call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.false., &
            selectBlockType=LEAF,unitReadsMeshDataOnly=.true.)

end if
  !print*, "In gas->sinks, Particle block prop = ", particles_local(BLK_PART_PROP,1:localnp)
! Just sinks.
!#ifdef SINK_PART_TYPE
!#ifndef ACTIVE_PART_TYPE

!#ifdef debug_force
!    if (dr_globalMe == 0) print*, "Inside just gas on sinks kick."
!#endif

!    call Particles_longRangeForce(particles_local, localnp, WEIGHTED)

!    do p=1, localnp

!    !  print*, particles_local(ACCX_PART_PROP,p), particles_local(ACCY_PART_PROP,p), &
!    !          particles_local(ACCZ_PART_PROP,p)

!      particles_local(VELX_PART_PROP,p) = particles_local(VELX_PART_PROP,p) &
!                                + particles_local(ACCX_PART_PROP,p)*dt

!      particles_local(VELY_PART_PROP,p) = particles_local(VELY_PART_PROP,p) &
!                                + particles_local(ACCY_PART_PROP,p)*dt

!      particles_local(VELZ_PART_PROP,p) = particles_local(VELZ_PART_PROP,p) &
!                                + particles_local(ACCZ_PART_PROP,p)*dt


!    !  print*, particles_local(ACCX_PART_PROP,p), particles_local(ACCY_PART_PROP,p), &
!    !          particles_local(ACCZ_PART_PROP,p)
!    end do
!    call pt_sinkGatherGlobal()
!#endif
!#endif

!! Just massive.
!#ifdef ACTIVE_PART_TYPE
!#ifndef SINK_PART_TYPE

!#ifdef debug_force
!    if (dr_globalMe == 0) print*, "Inside just gas on massive kick."
!#endif

!    call Particles_longRangeForce(particles, pt_numLocal, WEIGHTED)

!    do p=1, pt_numLocal


!      particles(VELX_PART_PROP,p) = particles(VELX_PART_PROP,p) &
!                                + particles(ACCX_PART_PROP,p)*dt

!      particles(VELY_PART_PROP,p) = particles(VELY_PART_PROP,p) &
!                                + particles(ACCY_PART_PROP,p)*dt

!      particles(VELZ_PART_PROP,p) = particles(VELZ_PART_PROP,p) &
!                                + particles(ACCZ_PART_PROP,p)*dt

!    end do
!#endif
!#endif

!  !print*, "In gas->sinks, Particle block prop = ", particles_local(BLK_PART_PROP,1:localnp)

!  !print*, "The particle masses are ", particles_local(MASS_PART_PROP,1:localnp)

!  !print*, "The accelerations of the particles are"

!! Both sinks and massive particles, meaning that the bridge is between the massive
!! and gas+sinks.

!#if defined(ACTIVE_PART_TYPE) && defined(SINK_PART_TYPE)

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
    print*, "Inside gas on both types kick."
#endif

  ! This is the acceleration of the gas on the massive particles.
  ! Also now includes the sinks acceleration, since I added that
  ! back in Grid_getAccelOneRow in
  ! Gravity/Poisson/BHTree43/Couple_AMUSE_Sinks_and_Stars.

  ! Now the particles array has both sinks and massive particles, but
  ! we only want to kick with the massive particles. So have to pass
  ! only that info. Note that we must have updated the pt_typeInfo structure
  ! here before we call to get the numbers for the arrays.
  ! This should be done at the end of Particles_advance. - JW

    call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

    !print*, "Particles count in gravity_gas_on_particles = ", type_count

    call Particles_longRangeForce(particles_pointer(:,type_begin:type_end), type_count, WEIGHTED)

if (sink_separate) then

    allocate(accx(type_count))
    allocate(accy(type_count))
    allocate(accz(type_count))


    accx = particles_pointer(ACCX_PART_PROP,type_begin:type_end)
    accy = particles_pointer(ACCY_PART_PROP,type_begin:type_end)
    accz = particles_pointer(ACCZ_PART_PROP,type_begin:type_end)

    !print*, "accx(1) = ", accx(1), dr_globalMe

! Map the gravitational acceleration of the sinks accelerating
! the massive particles.

    numAttrib=3
    particles(ACCX_PART_PROP,:) = 0.
    particles(ACCY_PART_PROP,:) = 0.
    particles(ACCZ_PART_PROP,:) = 0.

    attrib(GRID_DS_IND,1)=SGAX_VAR
    attrib(PART_DS_IND,1)=ACCX_PART_PROP

    !call Grid_mapMeshToParticles(particles_pointer(:,type_begin:type_end), &
    !                                     NPART_PROPS,BLK_PART_PROP, &
    !   type_count,[POSX_PART_PROP, POSY_PART_PROP,POSZ_PART_PROP], &
    !                              numAttrib,attrib,WEIGHTED,CENTER)

    !print*, "particles accx(1) =", particles(ACCX_PART_PROP,1), dr_globalMe

    attrib(GRID_DS_IND,2)=SGAY_VAR
    attrib(PART_DS_IND,2)=ACCY_PART_PROP

    !call Grid_mapMeshToParticles(particles_pointer(:,type_begin:type_end), &
    !                                     NPART_PROPS,BLK_PART_PROP, &
    !   type_count,[POSX_PART_PROP, POSY_PART_PROP,POSZ_PART_PROP], &
    !                              numAttrib,attrib,WEIGHTED,CENTER)

    attrib(GRID_DS_IND,3)=SGAZ_VAR
    attrib(PART_DS_IND,3)=ACCZ_PART_PROP

    call Grid_mapMeshToParticles(particles_pointer(:,type_begin:type_end), &
                                         NPART_PROPS,BLK_PART_PROP, &
       type_count,[POSX_PART_PROP, POSY_PART_PROP,POSZ_PART_PROP], &
                                  numAttrib,attrib,WEIGHTED,CENTER)

    !print*, "particles accx(1) =", particles(ACCX_PART_PROP,1), dr_globalMe

    do p=type_begin, type_end
    ! Kick massive particles with the gravity of the gas+sinks.
      particles_pointer(VELX_PART_PROP,p) = particles_pointer(VELX_PART_PROP,p) &
                                + (accx(p)+particles_pointer(ACCX_PART_PROP,p))*dt

      particles_pointer(VELY_PART_PROP,p) = particles_pointer(VELY_PART_PROP,p) &
                                + (accy(p)+particles_pointer(ACCY_PART_PROP,p))*dt

      particles_pointer(VELZ_PART_PROP,p) = particles_pointer(VELZ_PART_PROP,p) &
                                + (accz(p)+particles_pointer(ACCZ_PART_PROP,p))*dt
    end do

    !print*, "particles accx(1) =", particles(ACCX_PART_PROP,1), dr_globalMe

    deallocate(accx)
    deallocate(accy)
    deallocate(accz)

else

    do p=type_begin, type_end

      !print*, "For", p, "particles accx =", particles(ACCX_PART_PROP,p), dr_globalMe
    ! Kick massive particles with the gravity of the gas+sinks.
      particles_pointer(VELX_PART_PROP,p) = particles_pointer(VELX_PART_PROP,p) &
                                + particles_pointer(ACCX_PART_PROP,p)*dt

      particles_pointer(VELY_PART_PROP,p) = particles_pointer(VELY_PART_PROP,p) &
                                + particles_pointer(ACCY_PART_PROP,p)*dt

      particles_pointer(VELZ_PART_PROP,p) = particles_pointer(VELZ_PART_PROP,p) &
                                + particles_pointer(ACCZ_PART_PROP,p)*dt
    end do
end if

!#endif

!if (localnp .ne. 0) then

!print*, "In gas->particles, after calc, Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)

!end if

#ifdef debug_force
! Calculate the forces acting on the particles from the gas.
! For debugging purposes.

!if (part_type == 'sink') &
!    call pt_sinkGatherGlobal()

!if (dr_globalMe .eq. MASTER_PE) then

do p=type_begin, type_end

    force_GoS_x = force_GoS_x + particles_pointer(ACCX_PART_PROP,p)*particles_pointer(MASS_PART_PROP,p)
    force_GoS_y = force_GoS_y + particles_pointer(ACCY_PART_PROP,p)*particles_pointer(MASS_PART_PROP,p)
    force_GoS_z = force_GoS_z + particles_pointer(ACCZ_PART_PROP,p)*particles_pointer(MASS_PART_PROP,p)

end do

call MPI_ALLREDUCE(MPI_IN_PLACE, force_GoS_x, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, force_GoS_y, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, force_GoS_z, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)


if (dr_globalMe .eq. MASTER_PE) &
    & write(*,'(A,4(1X,E17.10))') 'Particles_sinkAccelGasOnSinks: Total force GAS->SINKS (time, x,y,z) = ', &
    & dr_simTime, force_GoS_x, force_GoS_y, force_GoS_z
!end if

if (dr_globalMe .eq. MASTER_PE) then

if (dr_SimTime == 0.0) then

file_name = "forceGoS.dat"
open(unit=12, file=trim(file_name))
file_name = "forceSoG.dat"
open(unit=13, file=trim(file_name))
file_name = "force_error.dat"
open(unit=14, file=trim(file_name))

else

file_name = "forceGoS.dat"
open(unit=12, file=trim(file_name), position='append')
file_name = "forceSoG.dat"
open(unit=13, file=trim(file_name), position='append')
file_name = "force_error.dat"
open(unit=14, file=trim(file_name), position='append')

end if

sum_force_norm = sqrt((abs(force_GoS_x) + abs(force_SoG_x)/ 2.0)**2.0 + &
                      (abs(force_GoS_y) + abs(force_SoG_y)/ 2.0)**2.0 + &
                      (abs(force_GoS_z) + abs(force_SoG_z)/ 2.0)**2.0)

write(12,'(4(1X,E22.15))') dr_simTime, force_GoS_x,force_GoS_y,force_GoS_z
write(13,'(4(1X,E22.15))') dr_simTime, force_SoG_x,force_SoG_y,force_SoG_z
write(14,'(4(1X,E22.15))') dr_simTime, &
                             (force_GoS_x+force_SoG_x) / sum_force_norm, &
                             (force_GoS_y+force_SoG_y) / sum_force_norm, &
                             (force_GoS_z+force_SoG_z) / sum_force_norm

close(12)
close(13)
close(14)

end if
#endif

! No gravity
#else
     print*, "ERROR! Gravity unit not compiled in!"
     stop
#endif


get_gravity_gas_on_particles=0
END FUNCTION get_gravity_gas_on_particles

!!! Here we call Gravity_potentialListOfBlocks as in the main Driver loop.
!!! But we only include the mapping of the SINK particle densities to the grid
!!! variable PDEN (which is usually only used for MASSIVE not SINK particles)
!!! so that the gravity solver uses this density on the grid for solving
!!! potential, which is stored in SGPT_VAR.

!!! It then solves for the proper potential on the grid. We then call
!!! Gravity_accelListOfBlocks (as in gas_on_particles) to get the accelerations
!!! on the gas. Then finally we kick the gas velocity with these accelerations
!!! over the bridge timestep. Note that this should be the half timestep
!!! which we are assuming the user passed properly (dt here should be 1/2 bridge total dt).

FUNCTION get_gravity_particles_on_gas(dt, kick_number)

use gr_ptInterface, ONLY : gr_ptVerifyBlock


integer :: get_gravity_particles_on_gas, blockID, blockCount, ierr
integer :: blocks(MAXBLOCKS), blkLimits(2,NDIM), blkLimitsGC(2,NDIM)
integer :: ii,jj,kk, particlesPerBlk(MAXBLOCKS, NPART_PROPS), kick_number, p
real*8    :: dt
real*8, pointer, dimension(:,:,:,:) :: solndata
logical, save :: regrid = .true.
real*8 :: del(3), dVol
!logical, parameter :: debug= .true. !.false.

#ifdef GRAVITY
if (.not. useGravity) return
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
   print*, "Kicking gas with bridge."
force_SoG_x=0.0; force_SoG_y=0.0; force_SoG_z=0.0
#endif

! First map the sink densities to the PDEN variable.
! Nevermind, this is done in the routine below!
!call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR)

!if (localnp .ne. 0) then

!print*, "Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)

!end if

! Note this MUST be called here, as this updates the particles array
! with the particles_local masses and positions. Then the particles array
! is used to calculate the potential from the particles gravity for the
! bridge (in the SGPT variable). Note doing it this way means we can use
! either particles, sinks or both in the bridge. - JW.

if (part_type == 'sink') then
    call Particles_sinkSyncWithParticles(sink_to_part=.true.)
!if (localnp .ne. 0) then
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
    print*, "Called SyncWithPart."
#endif
end if

!print*, "Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)

!end if

!call pt_sinkGatherGlobal()

! For right now, always assume regrid. To make it faster in the future
! we should make this depend on if this is the 1st or 2nd kick.

! Also note, it appears that in Particles_advance regrid was set as a
! parameter which was always false. So it would never look to see if
! particles have moved do to grid refinement.
call Particles_moveAndSort(regrid)

!if (localnp .ne. 0) then
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
print*, "Called Particles_moveAndSort."
#endif
!print*, "Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(BLK_PART_PROP,1:localnp)

!end if

!call Particles_moveAndSort(regrid)

!print*, "Particle block prop = ", particles_local(BLK_PART_PROP,1:localnp)

!call gr_ptVerifyBlock(particles_local, NPART_PROPS, localnp, sink_maxSinks)

call Grid_getListOfBlocks(LEAF, blocks, blockCount)


!#ifdef debug
!    call Grid_getBlkPtr(blocks(1),solndata,CENTER)

!    print*, "Grav potential from particles is:"
!    print*, solndata(SGPT_VAR,5,5,5)
!    print*, "Max grav accelx from particles on gas is:"
!    print*, maxval(solndata(SGAX_VAR,:,:,:))

!    call Grid_releaseBlkPtr(blocks(1),solndata)
!#endif
! Now solve the potential everywhere on the Grid.
! Also, for safety we should probably pass a different potential index
! to represent the potential of the particles->gas on the grid.
! Remember, at all other times FLASH should not see the particles to
! calculate gravity on the gas.
! Lets call it "Sinks on Gas PoTential" -> SGPT_VAR.

! This is handing back positive potentials in some places... something
! definitely seems wrong here.... maybe they define some reference pot?

! Calling Particles_advance before a potential solve might fix the bug
! in the gravity bridge... (this is the order its done in Driver_evolve)- JW
!call Particles_advance(1.0, 1.0)

!print*, "Particle block prop = ", particles_local(BLK_PART_PROP,1:localnp)

! Note that only during the codes evolve do the particles or gas actually
! move. So the second kick needs updated info on the gas and particle locations.
! The first kick only needs it on the very first call to bridge.
if (kick_number == 2) then
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "In kick 2"
#endif
  ! Fill guard cells because the hydro evolve step has moved the gas density
  ! field since the last time guard cells were updated.
  ! I had this commented out, perhaps this is done somewhere else?! 4/30/17.
  call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.false., selectBlockType=ACTIVE_BLKS) !LEAF

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Guard cells filled."
#endif

  ! First get the potential of just the gas.
  ! NOTE: This is important to do if any stars were made from sink particles,
  ! because now the sink masses are much lower (close to zero likely) and so
  ! the SGAX, SGAY, SGAZ accels on the grid are much, much higher than they
  ! should be. Since these are updated in this routine, calling it here makes
  ! sense if we formed star particles.
  ! It may make sense to not call the grid solver itself here, and just call
  ! the sink routines only. We can try this later.

  call Gravity_potentialListOfBlocks(blockCount, blocks, GPOT_VAR)

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Called Gravity_potGas."
#endif

! Then calculate the potential of the gas + particles, and subtract the gas
! potential.

! What are we kicking with here? If just sinks around then...

#ifdef SINK_PART_TYPE
#ifndef ACTIVE_PART_TYPE

  call Gravity_potentialListOfBlocksPDENonly(blockCount, blocks, BGPT_VAR, SINK_PART_TYPE)

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Called sink only version of Gravity_potentialListOfBlocksPDENonly."
#endif

#else
! There's both types and we are only kicking with the massive particles.
  call Gravity_potentialListOfBlocksPDENonly(blockCount, blocks, BGPT_VAR, ACTIVE_PART_TYPE)

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Called sink and massive version of Gravity_potentialListOfBlocksPDENonly."
#endif

#endif

#elif defined(ACTIVE_PART_TYPE)
! There's only massive particles and we are still kicking with them.
  call Gravity_potentialListOfBlocksPDENonly(blockCount, blocks, BGPT_VAR, ACTIVE_PART_TYPE)

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Called massive only version of Gravity_potentialListOfBlocksPDENonly."
#endif

#else
! Something here broke.
  print*, "[get_gravity_particles_on_gas]: No proper particle type found to kick."
  call Driver_abortFlash("No proper particle type found to kick.")
#endif


!if (localnp .ne. 0) then
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Called Gravity_PDENonly."
#endif
!print*, "Particles x, y, z positions are:"
!write(*,'(E10.4)') particles_local(POSX_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSY_PART_PROP,1:localnp)
!write(*,'(E10.4)') particles_local(POSZ_PART_PROP,1:localnp)

!end if
!print*, "Particle block prop = ", particles_local(BLK_PART_PROP,1:localnp)

! Don't forget to notify for this (including sink accels in case we took mass
! from them to make stars just before the bridge kick)!

call Grid_notifySolnDataUpdate( (/ GPOT_VAR, BGPT_VAR, SGAX_VAR, SGAY_VAR, SGAZ_VAR /) )

! Filling of guard cells could be important here... we did just change the
! potential variables on the grid.

  call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.false., selectBlockType=ACTIVE_BLKS) !LEAF
#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
  print*, "Guard cells filled."
#endif
! Now calculate the graviational acceleration from the sinks in each direction.
! Store this in the sink->gas acceleration variable.

  call Gravity_accelListOfBlocks(blockCount, blocks, IAXIS, BGAX_VAR, BGPT_VAR)
  call Gravity_accelListOfBlocks(blockCount, blocks, JAXIS, BGAY_VAR, BGPT_VAR)
  call Gravity_accelListOfBlocks(blockCount, blocks, KAXIS, BGAZ_VAR, BGPT_VAR)

end if

! Finally kick the gas velocity based on this acceleration.
! We may have to map this acceleration using CIC the same way we map it
! to the particles...

do blockID=1, blockCount

    call Grid_getBlkIndexLimits(blocks(blockID), blkLimits, blkLimitsGC, CENTER)
    call Grid_getBlkPtr(blocks(blockID),solndata,CENTER)

!#ifdef debug
!    if (blockID==1) then
!    print*, "Grav potential from particles is:"
!    print*, solndata(SGPT_VAR,5,5,5)
!    print*, "Max accelx from particles on gas is:"
!    print*, maxval(abs(solndata(SGAX_VAR,:,:,:)))
!    end if
!#endif

    ! Now lets kick the interior cells with the acceleration from the sinks.

    solndata(VELX_VAR,:,:,:) = solndata(VELX_VAR,:,:,:) + solndata(BGAX_VAR,:,:,:)*dt
    solndata(VELY_VAR,:,:,:) = solndata(VELY_VAR,:,:,:) + solndata(BGAY_VAR,:,:,:)*dt
    solndata(VELZ_VAR,:,:,:) = solndata(VELZ_VAR,:,:,:) + solndata(BGAZ_VAR,:,:,:)*dt

    ! For debugging purposes only.
#ifdef debug_force
    call Grid_getDeltas(blocks(blockID), del)
    dVol = del(1)*del(2)*del(3)

    !if (dr_globalMe .eq. MASTER_PE) then

        do ii = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
            do jj = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                do kk = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)

              force_SoG_x = force_SoG_x + solndata(BGAX_VAR,ii,jj,kk)*solndata(DENS_VAR,ii,jj,kk)*dVol
              force_SoG_y = force_SoG_y + solndata(BGAY_VAR,ii,jj,kk)*solndata(DENS_VAR,ii,jj,kk)*dVol
              force_SoG_z = force_SoG_z + solndata(BGAZ_VAR,ii,jj,kk)*solndata(DENS_VAR,ii,jj,kk)*dVol

                end do
            end do
        end do
    !end if

#endif
    ! Zero out the accelerations so they don't get added during the
    ! normal hydro solution.

    if (kick_number == 1) then
        solndata(BGAX_VAR,:,:,:) = 0.0
        solndata(BGAY_VAR,:,:,:) = 0.0
        solndata(BGAZ_VAR,:,:,:) = 0.0
    end if


    call Grid_releaseBlkPtr(blocks(blockID),solndata)

end do

#if defined(ACTIVE_PART_TYPE) && defined(SINK_PART_TYPE)
! Now kick the sink particles with the gravity of the massive particles.

! First map the acceleration from the massive particles onto the sinks.

!!! NOTE !!! This assumes the BGPT_VAR potential is up to date for this
!            kick, so make sure you call the kick for particles on gas FIRST.
!            Also note that here the BGPT_VAR is the potential from the
!            massive particles that kicks the sinks and the gas.

#ifdef debug_force
if (dr_globalMe .eq. MASTER_PE) &
print*, "Kicking sinks with massive gravity."
#endif

call Particles_longRangeForce(particles_local, localnp, WEIGHTED, BGPT_VAR)

! Now kick the sinks.

do p=1, localnp

!  print*, particles_local(ACCX_PART_PROP,p), particles_local(ACCY_PART_PROP,p), &
!          particles_local(ACCZ_PART_PROP,p)

  particles_local(VELX_PART_PROP,p) = particles_local(VELX_PART_PROP,p) &
                            + particles_local(ACCX_PART_PROP,p)*dt

  particles_local(VELY_PART_PROP,p) = particles_local(VELY_PART_PROP,p) &
                            + particles_local(ACCY_PART_PROP,p)*dt

  particles_local(VELZ_PART_PROP,p) = particles_local(VELZ_PART_PROP,p) &
                            + particles_local(ACCZ_PART_PROP,p)*dt


!  print*, particles_local(ACCX_PART_PROP,p), particles_local(ACCY_PART_PROP,p), &
!          particles_local(ACCZ_PART_PROP,p)
end do
call pt_sinkGatherGlobal()
#endif

#ifdef debug_force

#if defined(ACTIVE_PART_TYPE) && defined(SINK_PART_TYPE)

        do p=1, localnp

              force_SoG_x = force_SoG_x + particles_local(ACCX_PART_PROP,p)*particles_local(MASS_PART_PROP,p)
              force_SoG_y = force_SoG_y + particles_local(ACCY_PART_PROP,p)*particles_local(MASS_PART_PROP,p)
              force_SoG_z = force_SoG_z + particles_local(ACCZ_PART_PROP,p)*particles_local(MASS_PART_PROP,p)

        end do
#endif

call MPI_ALLREDUCE(MPI_IN_PLACE, force_SoG_x, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, force_SoG_y, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, force_SoG_z, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dr_globalComm, ierr)

if (dr_globalMe .eq. MASTER_PE) &
print*, "Updated the cell velocities."

     if (dr_globalMe .eq. MASTER_PE) &
       & write(*,'(A,4(1X,E17.10))') 'Particles_sinkAccelSinksOnGas: Total force SINKS->GAS (time, x,y,z) = ', &
       & dr_simTime, force_SoG_x, force_SoG_y, force_SoG_z

#endif
call Grid_notifySolnDataUpdate( (/ VELX_VAR, VELY_VAR, VELZ_VAR /) )

! No gravity
#else
     print*, "ERROR! Gravity unit not compiled in!"
     stop
#endif


get_gravity_particles_on_gas=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Other Grid Maintenence. Many are not in use now.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION initialize_grid()

  INTEGER :: initialize_grid
  initialize_grid=0
END FUNCTION

FUNCTION initialize_code()

  INTEGER :: initialize_code, part_init
  call Driver_initParallel()
  call Driver_initFlash()
  ! Make sure that when we exit an evolve step, Flash actually only
  ! evolves to the end time given.

  call RuntimeParameters_set("dr_shortenLastStepBeforeTMax",.true.)
  call RuntimeParameters_get("dr_shortenLastStepBeforeTMax",dr_shortenLastStepBeforeTMax)
  !call Driver_init()
  restart = .false.
#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
  ! Initialize the particle list pointer for the Flash particles array.
  part_init = set_particle_pointers('mass')
  ! Switch off internal integration of the particles array in Flash.
  part_init = internal_particle_integration_off()
#endif
  initialize_code=0
END FUNCTION

FUNCTION cleanup_code()

  INTEGER :: cleanup_code
  call Driver_finalizeFlash()
  cleanup_code=0
END FUNCTION

FUNCTION recommit_parameters()

  INTEGER :: recommit_parameters
  recommit_parameters=0
END FUNCTION

FUNCTION commit_parameters()

  INTEGER :: commit_parameters
  commit_parameters=0
END FUNCTION

FUNCTION get_global_grid_index_limits(global_indices)

  INTEGER :: global_indices(MDIM)
  INTEGER :: get_global_grid_index_limits
  call Grid_getGlobalIndexLimits(global_indices)
  get_global_grid_index_limits=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PARTICLE FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION get_number_of_particles(n)

  INTEGER :: n, get_number_of_particles, ierr
  INTEGER :: type_begin, type_end
  n = 0
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, n)


!print*, "part_type =", part_type, dr_globalMe
!print*, "n =", n, dr_globalMe

! I'm suspicious that the other MPI_REDUCE was returning too quickly on the
! root node. - JW

call MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, dr_globalComm, ierr)

! MPI_IN_PLACE occurs only on the root process. On other processors it errors.

!if (dr_globalMe == 0) then
!    call MPI_REDUCE(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, 0, dr_globalComm, ierr)
!else
!    call MPI_REDUCE(n, n, 1, MPI_INTEGER, MPI_SUM, 0, dr_globalComm, ierr)
!end if
!  n = localnp
!  n = localnpf

!!! Note, according to Klaus this function will sync sinks to the particle
!!! file. If needed in the future (would we run with a mix of sink and
!!! non-sink?).

!  call Particles_sinkSyncWithParticles(sink_to_part=.true.)

#endif
!#ifdef particle_exist
!  call Particles_getGlobalNum(n)
!#endif
get_number_of_particles=0
END FUNCTION

FUNCTION get_particle_position_array(tags, x, y, z, nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: x, y, z, tags
  integer :: get_particle_position_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

x = 0.0
y = 0.0
z = 0.0

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            x(i) = particles_pointer(POSX_PART_PROP, QSindex(j))
            y(i) = particles_pointer(POSY_PART_PROP, QSindex(j))
            z(i) = particles_pointer(POSZ_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            x(i) = 0.0
            y(i) = 0.0
            z(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    x = 0.0
    y = 0.0
    z = 0.0
endif

call MPI_AllReduce(MPI_IN_PLACE, x, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, y, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, z, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

  ! Sort by particle tag. Note that input positions should also be
  ! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)

  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then ! Yes, then lets do them all at once.

    do i=1, localnpf

      x(i) = particles_global(POSX_PART_PROP, QSindex(i))
      y(i) = particles_global(POSY_PART_PROP, QSindex(i))
      z(i) = particles_global(POSZ_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

      do j=1, nparts

        do i=1, localnpf

          if (particles_global(iptag,i) .eq. tags(j)) then
!          if (id_sorted(QSindex(i)) .eq. tags(j)) then ! Check for matching tag.

!            x(j) = particles_global(POSX_PART_PROP, QSindex(i))
!            y(j) = particles_global(POSY_PART_PROP, QSindex(i))
!            z(j) = particles_global(POSZ_PART_PROP, QSindex(i))
            x(j) = particles_global(POSX_PART_PROP, i)
            y(j) = particles_global(POSY_PART_PROP, i)
            z(j) = particles_global(POSZ_PART_PROP, i)

          end if

        end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if
#endif
#endif
get_particle_position_array=0
END FUNCTION

FUNCTION get_particle_velocity_array(tags,vx,vy,vz,nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: get_particle_velocity_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

vx = 0.0
vy = 0.0
vz = 0.0

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            vx(i) = particles_pointer(VELX_PART_PROP, QSindex(j))
            vy(i) = particles_pointer(VELY_PART_PROP, QSindex(j))
            vz(i) = particles_pointer(VELZ_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            vx(i) = 0.0
            vy(i) = 0.0
            vz(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    vx = 0.0
    vy = 0.0
    vz = 0.0
endif

call MPI_AllReduce(MPI_IN_PLACE, vx, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vy, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vz, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      vx(i) = particles_global(VELX_PART_PROP, QSindex(i))
      vy(i) = particles_global(VELY_PART_PROP, QSindex(i))
      vz(i) = particles_global(VELZ_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then

!        if (id_sorted(QSindex(i)) .eq. tags(j)) then

!          vx(j) = particles_global(VELX_PART_PROP, QSindex(i))
!          vy(j) = particles_global(VELY_PART_PROP, QSindex(i))
!          vz(j) = particles_global(VELZ_PART_PROP, QSindex(i))
          vx(j) = particles_global(VELX_PART_PROP, i)
          vy(j) = particles_global(VELY_PART_PROP, i)
          vz(j) = particles_global(VELZ_PART_PROP, i)

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if
#endif
#endif
get_particle_velocity_array=0
END FUNCTION

FUNCTION get_particle_acceleration_array(tags, ax, ay, az, nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: ax, ay, az, tags
  integer :: get_particle_acceleration_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

ax = 0.0
ay = 0.0
az = 0.0

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            ax(i) = particles_pointer(ACCX_PART_PROP, QSindex(j))
            ay(i) = particles_pointer(ACCY_PART_PROP, QSindex(j))
            az(i) = particles_pointer(ACCZ_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            ax(i) = 0.0
            ay(i) = 0.0
            az(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    ax = 0.0
    ay = 0.0
    az = 0.0
endif

call MPI_AllReduce(MPI_IN_PLACE, ax, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, ay, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, az, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

  ! Sort by particle tag. Note that input positions should also be
  ! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)

  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then ! Yes, then lets do them all at once.

    do i=1, localnpf

      ax(i) = particles_global(ACCX_PART_PROP, QSindex(i))
      ay(i) = particles_global(ACCY_PART_PROP, QSindex(i))
      az(i) = particles_global(ACCZ_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

      do j=1, nparts

        do i=1, localnpf

          if (particles_global(iptag,i) .eq. tags(j)) then
!          if (id_sorted(QSindex(i)) .eq. tags(j)) then ! Check for matching tag.

!            x(j) = particles_global(POSX_PART_PROP, QSindex(i))
!            y(j) = particles_global(POSY_PART_PROP, QSindex(i))
!            z(j) = particles_global(POSZ_PART_PROP, QSindex(i))
            ax(j) = particles_global(ACCX_PART_PROP, i)
            ay(j) = particles_global(ACCY_PART_PROP, i)
            az(j) = particles_global(ACCZ_PART_PROP, i)

          end if

        end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if
#endif
#endif
get_particle_acceleration_array=0
END FUNCTION

FUNCTION get_particle_mass(tags,mass,nparts)

  integer :: nparts, MyPe
  integer :: get_particle_mass, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: mass, tags
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

mass = 0.0

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            mass(i) = particles_pointer(MASS_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            mass(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    mass = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, mass, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#endif
get_particle_mass=0
END FUNCTION

FUNCTION get_particle_oldmass(tags,mass,nparts)

  integer :: nparts, MyPe
  integer :: get_particle_oldmass, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: mass, tags
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

mass = 0.0

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            mass(i) = particles_pointer(OLD_PMASS_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            mass(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    mass = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, mass, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#endif
get_particle_oldmass=0
END FUNCTION

FUNCTION get_particle_creation_time(tags,creation_time,nparts)

  integer :: nparts, MyPe
  integer :: get_particle_creation_time, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: creation_time, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            creation_time(i) = particles_pointer(CREATION_TIME_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            creation_time(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    creation_time = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, creation_time, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#endif
get_particle_creation_time=0
END FUNCTION

FUNCTION set_particle_creation_time(tags,creation_time,nparts)

  integer :: nparts
  double precision :: creation_time(nparts), tags(nparts)
  integer :: set_particle_creation_time, i, p, j, myProc, local_index, oldj
  integer*8 :: local_tag
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(CREATION_TIME_PART_PROP, QSindex(j)) = creation_time(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif

#endif

set_particle_creation_time=0
END FUNCTION

! This function should be called on any restart with shared tags
! between the stars and sinks.
FUNCTION set_starting_local_tag_numbers()
use Particles_interface, ONLY : pt_gatherGlobal
use Particles_data, ONLY : allproc_particles, pt_numGlobal
use Particles_sinkData, ONLY : local_tag_number

integer set_starting_local_tag_numbers, lp

    call pt_gatherGlobal()

!    In our version where all particles share this tag method, we can't
!    call this only on sinks, since for restarts or starts with particles it messes
!    up the tag tracking for all particles. Making this called on all particles array. - JW
!    Here this is called on particles_global... why? Is the fix to just make a global array
!    for regular particles array?
    local_tag_number = 0
    do lp = 1, pt_numGlobal
        if (get_ppe(int(allproc_particles(TAG_PART_PROP,lp),8)) .EQ. dr_globalMe) then
           local_tag_number = max(local_tag_number, get_pno(int(allproc_particles(TAG_PART_PROP,lp),8)))
        endif
    end do

set_starting_local_tag_numbers=0
END FUNCTION

FUNCTION get_sink_mean_cs(tags,cs,nparts)

  integer :: nparts, MyPe
  integer :: get_sink_mean_cs, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: cs, tags
  integer :: type_begin, type_end, type_count, offset
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

cs = 0.0

#if defined (SINK_PART_TYPE)

call get_particle_type_bounds('sink', type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    ! Offset for particles location in particles array possibly not starting at
    ! first index.
    offset = 0
    if (type_begin /= 1) offset = type_begin - 1

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p-offset) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            cs(i) = particles_pointer(CSGM_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            cs(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    cs = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, cs, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
#endif
get_sink_mean_cs = 0
return
END FUNCTION get_sink_mean_cs

FUNCTION get_sink_mean_vel_array(tags,vx,vy,vz,nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: get_sink_mean_vel_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINKS_AND_STARS)

#ifdef bisect

call get_particle_type_bounds('sink', type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_local(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            vx(i) = particles_local(VXGM_PART_PROP, QSindex(j))
            vy(i) = particles_local(VYGM_PART_PROP, QSindex(j))
            vz(i) = particles_local(VZGM_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            vx(i) = 0.0
            vy(i) = 0.0
            vz(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    vx = 0.0
    vy = 0.0
    vz = 0.0
endif

call MPI_AllReduce(MPI_IN_PLACE, vx, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vy, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vz, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else
call Driver_abortFlash('[get_sink_mean_vel_array]: Not using bisect!')
#endif
#else
call Driver_abortFlash('[get_sink_mean_vel_array]: Not using sinks!')
#endif

get_sink_mean_vel_array=0
END FUNCTION

FUNCTION get_sink_var_vel_array(tags,vx,vy,vz,nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: get_sink_var_vel_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINKS_AND_STARS)

#ifdef bisect

call get_particle_type_bounds('sink', type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_local(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            vx(i) = particles_local(VXGV_PART_PROP, QSindex(j))
            vy(i) = particles_local(VYGV_PART_PROP, QSindex(j))
            vz(i) = particles_local(VZGV_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            vx(i) = 0.0
            vy(i) = 0.0
            vz(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    vx = 0.0
    vy = 0.0
    vz = 0.0
endif

!print*, "gas var in first sink for x=", vx(1), "on", dr_globalMe

call MPI_AllReduce(MPI_IN_PLACE, vx, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vy, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, vz, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else
call Driver_abortFlash('[get_sink_var_vel_array]: Not using bisect!')
#endif
#else
call Driver_abortFlash('[get_sink_var_vel_array]: Not using sinks!')
#endif

get_sink_var_vel_array=0
END FUNCTION

FUNCTION get_sink_ang_mom_array(tags,lx,ly,lz,nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: lx, ly, lz, tags
  integer :: get_sink_ang_mom_array, i, j, p, oldj, ierr
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds('sink', type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_local(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.
        ! Note that since the inputs are sorted by tag, the last tag found
        ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            lx(i) = particles_local(X_ANG_PART_PROP, QSindex(j))
            ly(i) = particles_local(Y_ANG_PART_PROP, QSindex(j))
            lz(i) = particles_local(Z_ANG_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            lx(i) = 0.0
            ly(i) = 0.0
            lz(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else
    lx = 0.0
    ly = 0.0
    lz = 0.0
endif

call MPI_AllReduce(MPI_IN_PLACE, lx, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, ly, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)
call MPI_AllReduce(MPI_IN_PLACE, lz, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else
call Driver_abortFlash('[get_sink_ang_mom_array]: Not using bisect!')
#endif
#else
call Driver_abortFlash('[get_sink_ang_mom_array]: Not using sinks!')
#endif

!print*, "lx =", lx, dr_globalMe
!print*, "ly =", ly, dr_globalMe
!print*, "lz =", lz, dr_globalMe

get_sink_ang_mom_array=0
END FUNCTION


FUNCTION get_particle_nion(tags,nion,nparts)
! Get the number of ionizing photons from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_nion, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: nion, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            nion(i) = particles_pointer(MASS_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            nion(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    nion = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, nion, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      nion(i) = particles_global(NION_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          nion(j) = particles_global(NION_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if
#endif
get_particle_nion=0
END FUNCTION

FUNCTION get_particle_eion(tags,eion,nparts)
! Get the energy of ionizing photons OVER 13.6 eV
! (how much actually heats the gas) from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_eion, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: eion, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            eion(i) = particles_pointer(MASS_PART_PROP, QSindex(j))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            eion(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    eion = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, eion, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else

call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      eion(i) = particles_global(EION_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          eion(j) = particles_global(EION_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if
#endif
get_particle_eion=0
END FUNCTION

FUNCTION get_particle_sigh(tags,sigh,nparts)
! Get the Lyman cross section of ionizing photons from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_sigh, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: sigh, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            sigh(i) = particles_pointer(SIGH_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            sigh(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    sigh = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, sigh, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      sigh(i) = particles_global(SIGH_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          sigh(j) = particles_global(SIGH_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

get_particle_sigh=0
END FUNCTION

FUNCTION get_particle_rel_mass(tags,rel_mass,nparts)
! Get the SeBa relative mass from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_rel_mass, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: rel_mass, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            rel_mass(i) = particles_pointer(REL_MASS_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            rel_mass(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    rel_mass = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, rel_mass, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      rel_mass(i) = particles_global(REL_MASS_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          rel_mass(j) = particles_global(REL_MASS_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_rel_mass=0
END FUNCTION







FUNCTION get_particle_rel_age(tags,rel_age,nparts)
! Get the SeBa relative age from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_rel_age, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: rel_age, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            rel_age(i) = particles_pointer(REL_AGE_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            rel_age(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    rel_age = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, rel_age, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      rel_age(i) = particles_global(REL_AGE_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          rel_age(j) = particles_global(REL_AGE_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_rel_age=0
END FUNCTION



FUNCTION get_particle_corem(tags,corem,nparts)
! Get the SeBa core mass from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_corem, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: corem, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            corem(i) = particles_pointer(COREM_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            corem(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    corem = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, corem, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      corem(i) = particles_global(COREM_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          corem(j) = particles_global(COREM_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_corem=0
END FUNCTION



FUNCTION get_particle_co_corem(tags,co_corem,nparts)
! Get the SeBa CO core mass from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_co_corem, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: co_corem, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            co_corem(i) = particles_pointer(CO_COREM_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            co_corem(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    co_corem = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, co_corem, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      co_corem(i) = particles_global(CO_COREM_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          co_corem(j) = particles_global(CO_COREM_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_co_corem=0
END FUNCTION


FUNCTION get_particle_stype(tags,stype,nparts)
! Get the SeBa stellar type from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_stype, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: stype, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            stype(i) = particles_pointer(STYPE_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            stype(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    stype = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, stype, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      stype(i) = particles_global(STYPE_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          stype(j) = particles_global(STYPE_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_stype=0
END FUNCTION


FUNCTION get_particle_radius(tags,radius,nparts)
! Get the SeBa radius from the particle by tag.
  integer :: nparts, MyPe
  integer :: get_particle_radius, i, j, p, oldj, ierr
  double precision, dimension(nparts) :: radius, tags
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            radius(i) = particles_pointer(RADIUS_PART_PROP, QSindex(j))
            !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            radius(i) = 0.0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    radius = 0.0

endif

call MPI_AllReduce(MPI_IN_PLACE, radius, nparts, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, dr_globalcomm, ierr)

#else


call Driver_getMype(GLOBAL_COMM, MyPe)
call pt_sinkGatherGlobal()

! I only need one proc to do this.

if (MyPe .eq. 0) then

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

  allocate(QSindex(localnpf))
  allocate(id_sorted(localnpf))

  do p = 1, localnpf
     id_sorted(p) = int(particles_global(iptag,p),8)
  end do

  call NewQsort_IN(id_sorted, QSindex)


  ! Are we updating every particle in the simulation?

  if (nparts .eq. localnpf) then

  ! Yes, then lets do them all at once.

    do i=1, localnpf

      radius(i) = particles_global(RADIUS_PART_PROP, QSindex(i))

    end do

  else

  ! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then
!        if (id_sorted(QSindex(i)) .eq. int(tags(j))) then

          radius(j) = particles_global(RADIUS_PART_PROP, i)
!          mass(j) = particles_global(ipm, QSindex(i))

        end if

      end do

    end do

  end if

  deallocate(QSindex)
  deallocate(id_sorted)

end if

#endif

!#endif
get_particle_radius=0
END FUNCTION



FUNCTION set_particle_position(tags,x,y,z,nparts)

  integer :: nparts, MyPe
  double precision, dimension(nparts) :: x, y, z, tags
  integer :: set_particle_position, i, j, p, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(POSX_PART_PROP, QSindex(j)) = x(i)
            particles_pointer(POSY_PART_PROP, QSindex(j)) = y(i)
            particles_pointer(POSZ_PART_PROP, QSindex(j)) = z(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif

#else


call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then ! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(POSX_PART_PROP, QSindex(i)) = x(i)
    particles_global(POSY_PART_PROP, QSindex(i)) = y(i)
    particles_global(POSZ_PART_PROP, QSindex(i)) = z(i)

  end do

else

! If not doing them all, have to do it by tag number.

    do j=1, nparts

      do i=1, localnpf

        if (particles_global(iptag,i) .eq. tags(j)) then

!        if (id_sorted(QSindex(i)) .eq. tags(j)) then ! Check for matching tag.

!          particles_global(POSX_PART_PROP, QSindex(i)) = x(j)
!          particles_global(POSY_PART_PROP, QSindex(i)) = y(j)
!          particles_global(POSZ_PART_PROP, QSindex(i)) = z(j)
          particles_global(POSX_PART_PROP, i) = x(j)
          particles_global(POSY_PART_PROP, i) = y(j)
          particles_global(POSZ_PART_PROP, i) = z(j)

        end if

      end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(POSX_PART_PROP,i) = particles_global(POSX_PART_PROP, i)
    particles_local(POSY_PART_PROP,i) = particles_global(POSY_PART_PROP, i)
    particles_local(POSZ_PART_PROP,i) = particles_global(POSZ_PART_PROP, i)

deallocate(QSindex)
deallocate(id_sorted)

end do
#endif

#endif
set_particle_position=0
END FUNCTION

FUNCTION set_particle_velocity(tags,vx,vy,vz,nparts)

  integer :: nparts
  double precision, dimension(nparts) :: vx, vy, vz, tags
  integer :: set_particle_velocity, i, p, j, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(VELX_PART_PROP, QSindex(j)) = vx(i)
            particles_pointer(VELY_PART_PROP, QSindex(j)) = vy(i)
            particles_pointer(VELZ_PART_PROP, QSindex(j)) = vz(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif

#else

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)


! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(VELX_PART_PROP, QSindex(i)) = vx(i)
    particles_global(VELY_PART_PROP, QSindex(i)) = vy(i)
    particles_global(VELZ_PART_PROP, QSindex(i)) = vz(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

       if (particles_global(iptag,i) .eq. tags(j)) then

!      if (id_sorted(QSindex(i)) .eq. tags(j)) then

!        particles_global(VELX_PART_PROP, QSindex(i)) = vx(j)
!        particles_global(VELY_PART_PROP, QSindex(i)) = vy(j)
!        particles_global(VELZ_PART_PROP, QSindex(i)) = vz(j)
        particles_global(VELX_PART_PROP, i) = vx(j)
        particles_global(VELY_PART_PROP, i) = vy(j)
        particles_global(VELZ_PART_PROP, i) = vz(j)

      end if

    end do

  end do

end if

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(VELX_PART_PROP,i) = particles_global(VELX_PART_PROP, i)
    particles_local(VELY_PART_PROP,i) = particles_global(VELY_PART_PROP, i)
    particles_local(VELZ_PART_PROP,i) = particles_global(VELZ_PART_PROP, i)

end do

deallocate(QSindex)
deallocate(id_sorted)
#endif

#endif
set_particle_velocity=0
END FUNCTION

FUNCTION set_particle_mass(tags,mass, nparts)

  integer :: nparts
  double precision :: mass(nparts), tags(nparts)
  integer :: set_particle_mass, i, p, j, myProc, local_index, oldj
  integer*8 :: local_tag
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(MASS_PART_PROP, QSindex(j)) = mass(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif
#endif
set_particle_mass=0
END FUNCTION

FUNCTION set_particle_oldmass(tags,mass, nparts)

  integer :: nparts
  double precision :: mass(nparts), tags(nparts)
  integer :: set_particle_oldmass, i, p, j, myProc, local_index, oldj
  integer*8 :: local_tag
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(OLD_PMASS_PART_PROP, QSindex(j)) = mass(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif
#endif
set_particle_oldmass=0
END FUNCTION

FUNCTION set_particle_ang_mom(tags,lx,ly,lz,nparts)

  integer :: nparts
  double precision, dimension(nparts) :: lx, ly, lz, tags
  integer :: set_particle_ang_mom, i, p, j, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.
    ! Note that since the inputs are sorted by tag, the last tag found
    ! becomes the new lower bound for the next search.

        if (j .ne. -1) then
            particles_pointer(X_ANG_PART_PROP, QSindex(j)) = lx(i)
            particles_pointer(Y_ANG_PART_PROP, QSindex(j)) = ly(i)
            particles_pointer(Z_ANG_PART_PROP, QSindex(j)) = lz(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif

#else

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)


! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(X_ANG_PART_PROP, QSindex(i)) = lx(i)
    particles_global(Y_ANG_PART_PROP, QSindex(i)) = ly(i)
    particles_global(Z_ANG_PART_PROP, QSindex(i)) = lz(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

       if (particles_global(iptag,i) .eq. tags(j)) then

!      if (id_sorted(QSindex(i)) .eq. tags(j)) then

!        particles_global(VELX_PART_PROP, QSindex(i)) = vx(j)
!        particles_global(VELY_PART_PROP, QSindex(i)) = vy(j)
!        particles_global(VELZ_PART_PROP, QSindex(i)) = vz(j)
        particles_global(X_ANG_PART_PROP, i) = lx(j)
        particles_global(Y_ANG_PART_PROP, i) = ly(j)
        particles_global(Z_ANG_PART_PROP, i) = lz(j)

      end if

    end do

  end do

end if

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(X_ANG_PART_PROP,i) = particles_global(X_ANG_PART_PROP, i)
    particles_local(Y_ANG_PART_PROP,i) = particles_global(Y_ANG_PART_PROP, i)
    particles_local(Z_ANG_PART_PROP,i) = particles_global(Z_ANG_PART_PROP, i)

end do

deallocate(QSindex)
deallocate(id_sorted)
#endif

#endif
set_particle_ang_mom=0
END FUNCTION

FUNCTION set_particle_nion(tags, nion, nparts)
! Set the particle's number of ionizing photons by tag number.
  integer :: nparts
  double precision :: nion(nparts), tags(nparts)
  integer :: set_particle_nion, i, p, j, myProc, local_index, oldj
  integer*8 :: local_tag
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(NION_PART_PROP, QSindex(j)) = nion(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(NION_PART_PROP, QSindex(i)) = nion(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(NION_PART_PROP, i) = nion(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(NION_PART_PROP,i) = particles_global(NION_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

set_particle_nion=0
!stop

!do i=1, localnp

!    if (particles_local(NION_PART_PROP, i) .ne. 0.0) print*, particles_local(NION_PART_PROP, i)

!end do

END FUNCTION

FUNCTION set_particle_eion(tags, eion, nparts)
! Set the particle's energy of ionizing photons OVER 13.6 eV
! (the energy heating the gas) by tag number.
  integer :: nparts
  double precision :: eion(nparts), tags(nparts)
  integer :: set_particle_eion, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(EION_PART_PROP, QSindex(j)) = eion(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else
! Are we updating every particle in the simulation?

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)


if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(EION_PART_PROP, QSindex(i)) = eion(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(EION_PART_PROP, i) = eion(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(EION_PART_PROP,i) = particles_global(EION_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)

#endif

set_particle_eion=0
!stop

END FUNCTION

FUNCTION set_particle_sigh(tags, sigh, nparts)
! Set the particle's Lyman cross section for ionizing photons by tag number.
  integer :: nparts
  double precision :: sigh(nparts), tags(nparts)
  integer :: set_particle_sigh, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(SIGH_PART_PROP, QSindex(j)) = sigh(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(SIGH_PART_PROP, QSindex(i)) = sigh(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(SIGH_PART_PROP, i) = sigh(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(SIGH_PART_PROP,i) = particles_global(SIGH_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif


set_particle_sigh=0
END FUNCTION


FUNCTION set_particle_npep(tags, nion, nparts)
! Set the particle's number of photoelectric photons by tag number.
  integer :: nparts
  double precision :: nion(nparts), tags(nparts)
  integer :: set_particle_npep, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
! Are we using radiation?
! Are we using photoelectric heating?
#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

    if (j .ne. -1) then
        particles_pointer(NPEP_PART_PROP, QSindex(j)) = nion(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
    else ! If not found (j=-1), the particle is not on this proc. Skip.
        j = oldj
    end if

end do
deallocate(QSindex)
deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(NPEP_PART_PROP, QSindex(i)) = nion(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(NPEP_PART_PROP, i) = nion(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(NPEP_PART_PROP,i) = particles_global(NPEP_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)

#endif

set_particle_npep=0

END FUNCTION

FUNCTION set_particle_epep(tags, eion, nparts)
! Set the particle's average energy of photoelectric photons from 5.6-13.6 eV
! (the energy heating the gas) by tag number.
  integer :: nparts
  double precision :: eion(nparts), tags(nparts)
  integer :: set_particle_epep, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
! Are we using radiation?
! Are we using photoelectric heating?
#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(EPEP_PART_PROP, QSindex(j)) = eion(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(EPEP_PART_PROP, QSindex(i)) = eion(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(EPEP_PART_PROP, i) = eion(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(EPEP_PART_PROP,i) = particles_global(EPEP_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif


set_particle_epep=0

END FUNCTION

FUNCTION set_particle_sigd(tags, sigh, nparts)
! Set the particle's Lyman cross section for dust for photoelectric photons by tag number.
  integer :: nparts
  double precision :: sigh(nparts), tags(nparts)
  integer :: set_particle_sigd, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
! Are we using radiation?
! Are we using photoelectric heating?
#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(SPEP_PART_PROP, QSindex(j)) = sigh(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(SPEP_PART_PROP, QSindex(i)) = sigh(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(SPEP_PART_PROP, i) = sigh(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(SPEP_PART_PROP,i) = particles_global(SPEP_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif


set_particle_sigd=0

END FUNCTION

FUNCTION set_particle_rel_mass(tags, rel_mass, nparts)
! Set the SeBa relative mass of a particle by tag number.
  integer :: nparts
  double precision :: rel_mass(nparts), tags(nparts)
  integer :: set_particle_rel_mass, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(REL_MASS_PART_PROP, QSindex(j)) = rel_mass(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(REL_MASS_PART_PROP, QSindex(i)) = rel_mass(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(REL_MASS_PART_PROP, i) = rel_mass(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(REL_MASS_PART_PROP,i) = particles_global(REL_MASS_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_rel_mass=0
END FUNCTION



FUNCTION set_particle_rel_age(tags, rel_age, nparts)
! Set the SeBa relative age of a particle by tag number.
  integer :: nparts
  double precision :: rel_age(nparts), tags(nparts)
  integer :: set_particle_rel_age, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(REL_AGE_PART_PROP, QSindex(j)) = rel_age(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(REL_AGE_PART_PROP, QSindex(i)) = rel_age(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(REL_AGE_PART_PROP, i) = rel_age(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(REL_AGE_PART_PROP,i) = particles_global(REL_AGE_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_rel_age=0
END FUNCTION



FUNCTION set_particle_corem(tags, corem, nparts)
! Set the SeBa core mass of a particle by tag number.
  integer :: nparts
  double precision :: corem(nparts), tags(nparts)
  integer :: set_particle_corem, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(COREM_PART_PROP, QSindex(j)) = corem(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(COREM_PART_PROP, QSindex(i)) = corem(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(COREM_PART_PROP, i) = corem(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(COREM_PART_PROP,i) = particles_global(COREM_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_corem=0
END FUNCTION


FUNCTION set_particle_co_corem(tags, co_corem, nparts)
! Set the SeBa CO core mass of a particle by tag number.
  integer :: nparts
  double precision :: co_corem(nparts), tags(nparts)
  integer :: set_particle_co_corem, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(CO_COREM_PART_PROP, QSindex(j)) = co_corem(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(CO_COREM_PART_PROP, QSindex(i)) = co_corem(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(CO_COREM_PART_PROP, i) = co_corem(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(CO_COREM_PART_PROP,i) = particles_global(CO_COREM_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_co_corem=0
END FUNCTION


FUNCTION set_particle_stype(tags, stype, nparts)
! Set the SeBa stellar type of a particle by tag number.
  integer :: nparts
  double precision :: stype(nparts), tags(nparts)
  integer :: set_particle_stype, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(STYPE_PART_PROP, QSindex(j)) = stype(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(STYPE_PART_PROP, QSindex(i)) = stype(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(STYPE_PART_PROP, i) = stype(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(STYPE_PART_PROP,i) = particles_global(STYPE_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_stype=0
END FUNCTION

FUNCTION set_particle_radius(tags, radius, nparts)
! Set the SeBa radius of a particle by tag number.
  integer :: nparts
  double precision :: radius(nparts), tags(nparts)
  integer :: set_particle_radius, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
!#ifdef FERVENT

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            particles_pointer(RADIUS_PART_PROP, QSindex(j)) = radius(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
        end if

    end do
    deallocate(QSindex)
    deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(RADIUS_PART_PROP, QSindex(i)) = radius(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(RADIUS_PART_PROP, i) = radius(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(RADIUS_PART_PROP,i) = particles_global(RADIUS_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

!#endif

set_particle_radius=0
END FUNCTION


FUNCTION set_particle_wind_mass(tags, dmdt, nparts)
! Set the particle's wind dM/dt by tag number.
  integer :: nparts
  double precision :: dmdt(nparts), tags(nparts)
  integer :: set_particle_wind_mass, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

! Are we using winds?
#ifdef WIND_INJ

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

    if (j .ne. -1) then
        particles_pointer(DMDT_PART_PROP, QSindex(j)) = dmdt(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
    else ! If not found (j=-1), the particle is not on this proc. Skip.
        j = oldj
    end if

end do
deallocate(QSindex)
deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(DMDT_PART_PROP, QSindex(i)) = dmdt(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(DMDT_PART_PROP, i) = dmdt(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(DMDT_PART_PROP,i) = particles_global(DMDT_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

#else

print*, "[interface.F90]: WARNING! Called particles_set_wind_dmdt but &
         winds are not compiled into Flash!"
call Driver_abortFlash("No winds compiled in, but attempting to set wind dmdt.")

#endif

set_particle_wind_mass=0

END FUNCTION

FUNCTION set_particle_wind_vel(tags, velw, nparts)
! Set the particle's wind dM/dt by tag number.
  integer :: nparts
  double precision :: velw(nparts), tags(nparts)
  integer :: set_particle_wind_vel, i, p, j, myProc, local_index, local_tag, oldj
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

! Are we using winds?
#ifdef WIND_INJ

#ifdef bisect

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

    ! If found, set particle attribute accordingly.

    if (j .ne. -1) then
        particles_pointer(VELW_PART_PROP, QSindex(j)) = velw(i)
        !print*, "Set a local particle attrib on proc ", dr_globalMe
    else ! If not found (j=-1), the particle is not on this proc. Skip.
        j = oldj
    end if

end do
deallocate(QSindex)
deallocate(id_sorted)
endif
#else

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(VELW_PART_PROP, QSindex(i)) = velw(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
       if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global mass", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(ipm, QSindex(i)) , myProc
!        print*, "tags and submitted mass", tags(j), mass(j), myProc

        !particles_global(ipm, QSindex(i)) = mass(j)
        particles_global(VELW_PART_PROP, i) = velw(j)
!        print*, particles_global(ipm, i)

       end if

    end do

  end do

end if


! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(VELW_PART_PROP,i) = particles_global(VELW_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do
deallocate(QSindex)
deallocate(id_sorted)
#endif

#else

print*, "[interface.F90]: WARNING! Called particles_set_wind_dmdt but &
         winds are not compiled into Flash!"
call Driver_abortFlash("No winds compiled in, but attempting to set wind dmdt.")

#endif

set_particle_wind_vel=0

END FUNCTION

FUNCTION set_particle_gpot(tags,gpot,nparts)

  integer :: nparts
  double precision :: gpot(nparts), tags(nparts)
  integer :: set_particle_gpot, i, p, j, myProc, local_index, local_tag
  integer*8, dimension(:), allocatable :: QSindex, id_sorted

#ifdef SINK_PART_TYPE

call Driver_getMype(GLOBAL_COMM, myProc)

call pt_sinkGatherGlobal()

! Sort by particle tag. Note that input positions should also be
! ordered by tag number then.

allocate(QSindex(localnpf))
allocate(id_sorted(localnpf))

do p = 1, localnpf
   id_sorted(p) = int(particles_global(iptag,p),8)
end do

call NewQsort_IN(id_sorted, QSindex)

! Are we updating every particle in the simulation?

if (nparts .eq. localnpf) then

! Yes, then lets do them all at once.

  do i=1, localnpf

    particles_global(GPOT_PART_PROP, QSindex(i)) = gpot(i)

  end do

else

! If not doing them all, have to do it by tag number.

  do j=1, nparts

    do i=1, localnpf

!      local_index = id_sorted(QSindex(i))
!      local_tag   = int(tags(j))

      !if (local_index == local_tag) then
!      if (id_sorted(QSindex(i)) == int(tags(j))) then
      if (particles_global(iptag,i) .eq. tags(j)) then

!        print*, "id_sorted and global gpot", id_sorted(QSindex(i)), particles_global(iptag, QSindex(i)), &
!                             & particles_global(GPOT_PART_PROP, QSindex(i)) , myProc
!        print*, "tags and submitted gpot", tags(j), gpot(j), myProc

        !particles_global(GPOT_PART_PROP, QSindex(i)) = gpot(j)
        particles_global(GPOT_PART_PROP, i) = gpot(j)
        !print*, particles_global(GPOT_PART_PROP, i)

      end if

    end do

  end do

end if

deallocate(QSindex)
deallocate(id_sorted)

! The first localnp particles on a processor in particles_global are
! the local particle arrays in order.

do i=1, localnp

    particles_local(GPOT_PART_PROP,i) = particles_global(GPOT_PART_PROP, i)
!    print*, "Final local and global mass."
!    print*, particles_local(ipm,i), myProc
!    print*, particles_global(ipm,i), myProc

end do


!do i=1, nparts

!    particles_local(GPOT_PART_PROP,n(i)) = gpot(i)

!end do
call pt_sinkGatherGlobal()
#endif
set_particle_gpot=0
END FUNCTION

FUNCTION get_particle_gpot(n,gpot,nparts)

  integer :: nparts
  double precision, dimension(nparts) :: gpot
  integer :: n(nparts), get_particle_gpot, i

#ifdef SINK_PART_TYPE

do i=1, nparts

    gpot(i) = particles_local(GPOT_PART_PROP,n(i))

end do

#endif
get_particle_gpot=0
END FUNCTION

!FUNCTION set_particle_prop(tags, prop, nparts)
!! Set any particular particle property by tag number.
!! TAGS MUST BE PASSED IN ASCENDING ORDER.
!  integer :: nparts
!  double precision :: prop(nparts), tags(nparts)
!  integer :: set_particle_prop, i, p, j, oldj
!  integer*8, dimension(:), allocatable :: QSindex, id_sorted

!! Are we using radiation?
!!#ifdef FERVENT

!!#define bisect_test

!!#ifdef bisect

!! Sort by particle tag. Note that input array should also be
!! ordered by tag number then.

!allocate(QSindex(localnp))
!allocate(id_sorted(localnp))

!do p = 1, localnp
!   id_sorted(p) = int(particles_local(iptag,p),8)
!end do

!call NewQsort_IN(id_sorted, QSindex)

!! Initial lower bound is 1.
!j = 1

!do i=1, nparts

!    oldj = j
!    j = bisect_search(tags(i), j, localnp, localnp, real(id_sorted,8))

!    ! If found, set particle attribute accordingly.

!    if (j .ne. -1) then
!        particles_local(EION_PART_PROP, QSindex(j)) = eion(i)
!        !print*, "Set a local particle attrib on proc ", dr_globalMe
!    else ! If not found (j=-1), the particle is not on this proc. Skip.
!        j = oldj
!    end if

!end do

!END FUNCTION

!!! This gives the acceleration on the particles due to the gas.
!!! We will use this to add to the acceleration on the particles in
!!! PH4 during the gravity bridge.

! Note: Do we need to pass an array of the particle tags here?
! This is going to be the new get_gravity_at_point.

!!! NOTE RIGHT NOW I'M ASSUMING YOU ARE UPDATING THE PARTICLES 1 thru n
!!! IN PROPER ORDER. So x,y,z need to be in proper order.

!FUNCTION get_gravity_at_point(eps, x, y, z, gax, gay, gaz, nparts)
!!FUNCTION get_accel_gas_on_particles(eps, x, y, z, gax, gay, gaz, nparts)
!
!  INTEGER :: nparts
!  DOUBLE PRECISION :: eps
!  DOUBLE PRECISION, DIMENSION(nparts) :: x, y, z, gax, gay, gaz
!!  LOGICAL :: usePart, useSinkPart
!  LOGICAL :: correct_location
!  INTEGER :: i, n(nparts) !, get_accel_gas_on_particles
!  INTEGER :: get_gravity_at_point


!#ifdef SINK_PART_TYPE
!!call Particles_sinkAccelGasOnSinks() ! Calculate the accel from the gas on sinks.

!call pt_sinkGatherGlobal()

!gax=0.0; gay=0.0; gaz=0.0

!correct_location = .true.

!check_pos: do i=1, nparts

!  if (x(i) .ne. particles_local(POSX_PART_PROP, i) .or. &
!      y(i) .ne. particles_local(POSY_PART_PROP, i) .or. &
!      z(i) .ne. particles_local(POSZ_PART_PROP, i)) then

!      correct_location = .false.
!      exit check_pos

!  end if

!end do check_pos

!if (correct_location) then

!  call Particles_sinkAccelGasOnSinks()

!  do i=1, nparts

!!    gax(i)=0.0; gay(i)=0.0; gaz(i)=0.0

!    gax(i) = particles_local(ACCX_PART_PROP,i)
!    gay(i) = particles_local(ACCY_PART_PROP,i)
!    gaz(i) = particles_local(ACCZ_PART_PROP,i)

!  end do

!!  return

!else

!  write(*,*) "Updating particle positions."

!   n = (/ (i, i = 1, nparts) /)

!  i = set_particle_position(n, x, y, z, nparts)

!  call Particles_sinkAccelGasOnSinks()

!  do i=1, nparts

!!    gax(i)=0.0; gay(i)=0.0; gaz(i)=0.0

!    gax(i) = particles_local(ACCX_PART_PROP,i)
!    gay(i) = particles_local(ACCY_PART_PROP,i)
!    gaz(i) = particles_local(ACCZ_PART_PROP,i)

!  end do

!end if

!!  else if (usePart) then
!!#elif defined particle_exist
!!    call Particles_longRangeForce()
!!    ax = particles_global(ACCX_PART_PROP,n)
!!    ay = particles_global(ACCY_PART_PROP,n)
!!    az = particles_global(ACCZ_PART_PROP,n)
!!  else
!#else
!    write(*,*) "No particles in this simulation found."
!#endif
!!  end if
!  get_gravity_at_point=0
!!  get_accel_gas_on_particles=0
!END FUNCTION

FUNCTION get_num_part_prop(n)

  integer :: n, get_num_part_prop

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

  n = NPART_PROPS
#endif
  get_num_part_prop=0
END FUNCTION


FUNCTION make_sink(x, y, z, tags, nparts)

integer :: nparts, n_created
real*8, dimension(nparts)   :: x, y, z, tags, local_tags
real*8, dimension(NPART_PROPS, nparts) :: part_copy
integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk
real*8   :: time
integer :: block, Proc_ID, myProc, make_sink, part_num, &
           communicator, i, ierror

#ifdef SINK_PART_TYPE

! This makes the particles in parallel, but how do you know the order
! is the same as the order for AMUSE? You don't, so you'll have to order
! them correctly at the end of the call. - JW

n_created = 0
tags = 0.0
local_tags = 0.0

!write(*,*) nparts

do i=1, nparts

  Proc_ID = 0
  block = 0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)
  call Grid_getBlkIDFromPos([x(i),y(i),z(i)], block, Proc_ID, communicator)
  call Driver_getSimTime(time)

  if (myProc .eq. Proc_ID) then

    part_num = pt_sinkCreateParticle(x(i), y(i), z(i), time, block, Proc_ID)

    local_tags(i)  = particles_local(iptag, part_num) ! should be i not n_created?

  end if

end do


call MPI_Reduce(local_tags, tags, nparts, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierror)


call pt_sinkGatherGlobal()
call Particles_sinkSyncWithParticles(.true.)
print*, "total number of particles =", localnp, pt_numLocal
#else

call Driver_abortFlash("Called make_sink but sinks not compiled in!")

#endif
make_sink=0
END FUNCTION

! Add massive/active type particles, not sinks. - JW
FUNCTION add_particles(x, y, z, tags, nparts)

integer :: nparts, n_created
real*8, dimension(nparts)   :: x, y, z, tags, local_tags
real*8, dimension(NPART_PROPS, nparts) :: part_copy
integer, dimension(MAXBLOCKS, NPART_TYPES) :: particlesPerBlk
real*8   :: time(nparts)
integer :: block(nparts), Proc_ID(nparts), myProc, add_particles, &
           communicator, i, ierror, init_num_parts, num_parts, &
           type_begin, type_end

#ifdef ACTIVE_PART_TYPE

!print*, "In add_particles!"

! This makes the particles in parallel, but how do you know the order
! is the same as the order for AMUSE? You don't, so you'll have to order
! them correctly at the end of the call. - JW

n_created = 0
tags = 0.0
local_tags = 0.0

call get_particle_type_bounds(part_type, type_begin, type_end, init_num_parts)


    call Particles_addNew(nparts, x, y, z, n_created, local_tags, ptype_in=real(ACTIVE_PART_TYPE,8))

call Particles_moveAndSort(.true.)

call get_particle_type_bounds(part_type, type_begin, type_end, num_parts)


call MPI_Reduce(local_tags, tags, nparts, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, dr_globalComm, ierror)

#else

call Driver_abortFlash("Called add_particles but active particles not compiled in!")

#endif

!print*, "total number of particles =", pt_numLocal

add_particles=0
END FUNCTION

FUNCTION new_source_flag(flag)
#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
use Particles_data, ONLY : new_source
#endif
logical, intent(in) :: flag
integer new_source_flag, ierr

#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
call MPI_ALLREDUCE(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, dr_globalComm, ierr)

new_source = flag
#endif
new_source_flag = 0
END FUNCTION

FUNCTION get_particle_tags(n, tags, nparts)

integer :: get_particle_tags, nparts, i, j, p, MyPe, n(nparts)
real*8  :: tags(nparts)
integer*8, dimension(:), allocatable :: QSindex, id_sorted

integer   :: communicator, ierr
integer   :: num_array(dr_globalNumProcs)
integer*8 :: new_tags_array(nparts)
integer   :: disp(dr_globalNumProcs), rank_minus_one
integer   :: type_begin, type_end, type_count
real*8    :: real_tags(nparts)

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

!call pt_sinkGatherGlobal()

!nparts = localnpf

! Get all the tags. It doesn't make sense to get "some" tags,
! since tags are how we tell different particles apart. How would you
! distigush "which tags" you wanted?

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

!print*, "type count", type_count, dr_globalMe
!print*, "type_begin", type_begin, dr_globalMe
!print*, "type_end", type_end, dr_globalMe
!print*, "nparts", nparts, dr_globalMe

!call Driver_getMype(GLOBAL_COMM, MyPe)

  disp = 0
  rank_minus_one = dr_globalNumProcs - 1
  call Driver_getComm(GLOBAL_COMM, communicator)
  !call Driver_getMype(GLOBAL_COMM, MyPe)

  ! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array. This can be
  ! gotten from get_number_of_new_tags.

  ! Make an array of the # of incoming particles from each processor.
  call MPI_Gather(type_count, 1, MPI_INTEGER, &
                  num_array, 1, MPI_INTEGER, &
                  0, communicator, ierr)

  ! Set the displacement for the incoming data based on how many
  ! particles are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.

  do i=1, dr_globalNumProcs-1

    disp(i+1) = disp(i) + num_array(i)

  end do

!print*, int(particles_pointer(TAG_PART_PROP,type_begin:type_begin+1)), dr_globalMe
!print*, type_count, dr_globalMe
!print*, new_tags_array(1), dr_globalMe
!print*, num_array(1), dr_globalMe
!print*, disp, dr_globalMe
!print*, communicator, dr_globalMe
!print*, ierr, dr_globalMe
!call flush(6)

  ! Now actually gather the tags using the variable length array
  ! gather command in MPI.
  
  call MPI_Gatherv(int(particles_pointer(TAG_PART_PROP,type_begin:type_end),8), &
                   type_count, MPI_LONG, new_tags_array, num_array, &
                   disp, MPI_LONG, 0, communicator, ierr)

! I only need one proc to do this.

if (dr_globalMe .eq. 0) then

! Sort by particle tag # so that the returned tags are always in ascending order.

allocate(QSindex(nparts))
allocate(id_sorted(nparts))

!  do p = 1, localnpf
!     id_sorted(p) = int(particles_global(iptag,p),8)
!  end do

  call NewQsort_IN(new_tags_array, QSindex)

  do i=1, nparts
! implicit type conversion from int to real here!
    tags(i) = real(new_tags_array(QSindex(i)),8)

  end do

  deallocate(QSindex)
  deallocate(id_sorted)

  !print*, "Tags on MasterPE =", tags

end if

#endif

get_particle_tags=0

END FUNCTION

FUNCTION get_particle_proc(tags, procs, nparts)

  integer :: nparts, MyPe
  integer :: get_particle_proc, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: tags
  integer, dimension(nparts) :: procs
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            procs(i) = int(particles_pointer(PROC_PART_PROP, QSindex(j)))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            procs(i) = 0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    procs = 0

endif

call MPI_AllReduce(MPI_IN_PLACE, procs, nparts, MPI_INTEGER, &
                   MPI_SUM, dr_globalcomm, ierr)

#endif

get_particle_proc=0

END FUNCTION

FUNCTION get_particle_block(tags, blocks, nparts)

  integer :: nparts, MyPe
  integer :: get_particle_block, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: tags
  integer, dimension(nparts) :: blocks
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, set particle attribute accordingly.

        if (j .ne. -1) then
            blocks(i) = int(particles_pointer(MASS_PART_PROP, QSindex(j)))

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            blocks(i) = 0
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

else

    blocks = 0

endif

call MPI_AllReduce(MPI_IN_PLACE, blocks, nparts, MPI_INTEGER, &
                   MPI_SUM, dr_globalcomm, ierr)

#endif
get_particle_block=0

END FUNCTION

! Note! You have to pass an array of the number of particles for
! this function due to the interface requiring the same length
! input array as the output array. - JW

FUNCTION get_new_tags(new_tags_length, tags, nparts)

integer   :: nparts
real*8    :: tags(nparts)
integer   :: new_tags_length(nparts), get_new_tags
integer   :: communicator, ierr, i
integer*8 :: new_tags_array(nparts)
integer   :: num_array(dr_globalNumProcs)
integer   :: disp(dr_globalNumProcs), MyPe, rank_minus_one
integer*8, dimension(:), allocatable :: QSindex, id_sorted

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)

  !print*, "In get_new_tags", dr_globalMe
  !print*, "number_new_particles=", number_new_particles, dr_globalMe
  !print*, "new_particles_tags=", new_particles_tags(1:number_new_particles), dr_globalMe

  disp = 0
  rank_minus_one = dr_globalNumProcs - 1
  call Driver_getComm(GLOBAL_COMM, communicator)
  !call Driver_getMype(GLOBAL_COMM, MyPe)

  ! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array. This can be
  ! gotten from get_number_of_new_tags.

  ! Make an array of the # of incoming particles from each processor.
  call MPI_Gather(number_new_particles, 1, MPI_INTEGER, &
                  num_array, 1, MPI_INTEGER, &
                  0, communicator, ierr)

  ! Set the displacement for the incoming data based on how many
  ! particles are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.

  do i=1, dr_globalNumProcs-1

    disp(i+1) = disp(i) + num_array(i)

  end do

  ! Now actually gather the tags using the variable length array
  ! gather command in MPI.
  call MPI_Gatherv(new_particles_tags, number_new_particles, MPI_LONG, &
                  new_tags_array, num_array, disp, MPI_LONG, &
                  0, communicator, ierr)

! I only need one proc to do this.

if (dr_globalMe .eq. 0) then

! Sort by particle tag # so that the returned tags are always in ascending order.

allocate(QSindex(nparts))
allocate(id_sorted(nparts))

!  do p = 1, localnpf
!     id_sorted(p) = int(particles_global(iptag,p),8)
!  end do

  call NewQsort_IN(new_tags_array, QSindex)

  do i=1, nparts
! implicit type conversion from int to real here!
    tags(i) = real(new_tags_array(QSindex(i)),8)

  end do

  deallocate(QSindex)
  deallocate(id_sorted)

  !print*, "[get_new_tags]: Tags on MasterPE =", tags

end if


  !print*, "Reals are ", real_tags

  !new_tags_array = floor(real_tags)

  !print*, "Gathered tags are ", new_tags_array(1:sum(num_array))

  !print*, "Ierr =", ierr

!new_tags_array  = new_tags(1:new_tags_length)
#endif
get_new_tags = 0
END FUNCTION

FUNCTION get_number_of_new_tags(new_tag_num)

integer :: new_tag_num, get_number_of_new_tags
integer :: communicator, ierr

#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)
  call Driver_getComm(GLOBAL_COMM, communicator)

  ! Number_new_sinks is the local processor number that have
  ! been created. We combine these to give back to the interface.

  call MPI_Reduce(number_new_particles, new_tag_num, 1, MPI_INTEGER, &
                MPI_SUM, 0, communicator, ierr)

!new_tag_num = number_new_sinks
#endif
get_number_of_new_tags = 0
END FUNCTION

FUNCTION clear_new_tags()

integer :: clear_new_tags
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)
number_new_particles = 0
new_particles_tags = 0.0
#endif
clear_new_tags = 0
END FUNCTION

FUNCTION particles_gather()
integer :: particles_gather
#ifdef SINK_PART_TYPE
call pt_sinkGatherGlobal()
#else
print*, "particles_gather: Warning, sink particles not in this simulation!"
#endif
particles_gather=0
END FUNCTION

FUNCTION particles_sort()

!use Particles_interface, ONLY : Particles_advanceNoCreate
!
integer :: particles_sort
!real*8 :: dummy_dtOld, dummy_dtNew
#if defined (SINK_PART_TYPE) || defined (ACTIVE_PART_TYPE)
#ifdef SINK_PART_TYPE
call pt_sinkGatherGlobal()
#endif
!particles(:,1:pt_numLocal) = particles_local(:,1:localnp)
call Particles_moveAndSort(regrid=.true.)
!
!call Particles_sinkSyncWithParticles(sink_to_part=.true.)
!

!call Particles_sinkSyncWithParticles(sink_to_part=.true.)

!print*, "particles_gather: Warning, sink particles not in this simulation!"
#endif
particles_sort=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! subroutine Particles_sinkRemove
!!!
!!! Remove sink particles from the simulation by tag number.
!!! Note it takes the particle out of: particles_local,
!!! particles_global and particles arrays.
!!! Also note it assumes tags are passed in ascending order!!!!!!**!!*!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION remove_particles(tags, nparts)

  integer :: nparts, MyPe
  integer :: remove_particles, i, j, p, oldj, ierr, counter
  double precision, dimension(nparts) :: tags
  integer, dimension(nparts) :: remove_index
  integer :: type_begin, type_end, type_count
  integer*8, dimension(:), allocatable :: QSindex, id_sorted
#if defined(SINK_PART_TYPE) || defined(ACTIVE_PART_TYPE)
! Initialize the remove_index array.
remove_index = -1

call get_particle_type_bounds(part_type, type_begin, type_end, type_count)

if (type_count .ge. 1) then

! Sort by particle tag. Note that input array should also be
! ordered by tag number then.

    allocate(QSindex(type_count))
    allocate(id_sorted(type_count))

    do p = type_begin, type_end
       id_sorted(p) = int(particles_pointer(iptag,p),8)
    end do

    if (type_count .eq. 1) then
        QSindex = 1
    else
        call NewQsort_IN(id_sorted, QSindex)
    endif

    ! Initial lower bound is 1.
    j = 1

    do i=1, nparts

        oldj = j
        j = bisect_search(tags(i), j, type_count, type_count, real(id_sorted,8))

        ! If found, add the index to the remove array.

        if (j .ne. -1) then
            remove_index(i) = QSindex(j)

        else ! If not found (j=-1), the particle is not on this proc. Skip.
            j = oldj
            !remove_index(i) = -1
        end if

    end do

    deallocate(QSindex)
    deallocate(id_sorted)

endif

! Sort the removal indicies from largest to smallest, since we remove
! particles by shifting the higher indicies down one, overwriting the
! removed particle.
call ut_qsort(remove_index, nparts, .false.)

do i=1, nparts

    if (remove_index(i) .gt. 0) then

! If both types are defined, we have to deal with massive and sink
! differently.
#if defined (SINK_PART_TYPE) && defined (ACTIVE_PART_TYPE)
        if (part_type .eq. 'mass') then

            ! Remove from the particles array.
            particles_pointer(:,remove_index(i):pt_numLocal) = particles_pointer(:,remove_index(i)+1:pt_numLocal+1)
            ! Correct the counts of all particles in this array and the count for massive type.
            pt_numLocal = pt_numLocal - 1
            pt_typeInfo(PART_LOCAL, ACTIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL, ACTIVE_PART_TYPE) - 1

        else if (part_type .eq. 'sink') then

            ! Remove from the particles_local array and the particles array.
            particles_pointer(:,remove_index(i):localnp) = particles_pointer(:,remove_index(i)+1:localnp+1)
            type_begin = pt_typeInfo(PART_TYPE_BEGIN,SINK_PART_TYPE)
            particles(:,type_begin+remove_index(i):pt_numLocal) = particles(:,type_begin+remove_index(i)+1:pt_numLocal+1)
            ! Correct the counts of all particles in both arrays and the count for sink type.
            pt_numLocal = pt_numLocal - 1
            localnp     = localnp - 1
            pt_typeInfo(PART_LOCAL, SINK_PART_TYPE) = pt_typeInfo(PART_LOCAL, SINK_PART_TYPE) - 1

        else

            print*, "[remove_particles]: Unrecognized part type given."

        end if

#elif defined (SINK_PART_TYPE) && !defined (ACTIVE_PART_TYPE)

        if (part_type .eq. 'sink') then
            ! Remove from the particles_local array and the particles array.
            particles_pointer(:,remove_index(i):localnp) = particles_pointer(:,remove_index(i)+1:localnp+1)
            particles(:,type_begin+remove_index(i):pt_numLocal) = particles(:,type_begin+remove_index(i)+1:pt_numLocal+1)
            ! Correct the counts of all particles in both arrays.
            pt_numLocal = pt_numLocal - 1
            localnp     = localnp - 1

        else

            print*, "[remove_particles]: Part type not sink and should be."

        end if

#elif !defined (SINK_PART_TYPE) && defined (ACTIVE_PART_TYPE)

        if (part_type .eq. 'mass') then

            ! Remove from the particles array.
            particles_pointer(:,remove_index(i):pt_numLocal) = particles_pointer(:,remove_index(i)+1:pt_numLocal+1)
            ! Correct the counts of all particles in this array.
            pt_numLocal = pt_numLocal - 1

        else

            print*, "[remove_particles]: Part type not mass and should be."

        end if
#endif
    end if

end do

!call Particles_moveAndSort(.false.)
#endif
remove_particles=0
return
END FUNCTION

FUNCTION remove_all_particles(clear_local_tag)
use Particles_sinkData, only : local_tag_number
implicit none
integer :: remove_all_particles
logical, intent(in) :: clear_local_tag

if (dr_globalMe == 0) then
    print*, "[interface.F90:remove_particles]: WARNING! Flash resetting all particle ", &
        "arrays of type ", part_type, " to NONEXISTENT and number of particles to zero on all procs!"
end if

particles_pointer  = NONEXISTENT
num_part_local_ptr = 0

! NOTE: Only on the first clearing of all particles would you want to reset the
! global tags. I.E. if you wanted to reset all the sink and massive particles,
! you'd clear the global tags for both on the first call where you removed all
! the sinks (and cleared tags for both) but not when you removed the massive,
! since both particle types share the same tags and tag generation method. - JW

if (clear_local_tag) then
    ! Also reset all tags.
    local_tag_number = 0
end if

remove_all_particles = 0
return
END FUNCTION
!!! Give a velocity kick to one block in Flash over a time step dt.

! Here we pass accel arrays for more than one block at a time, along with
! an array with the number of blocks on each proc in ascending order (the
! same order as the accel arrays). We then loop over the accel matching
! the correct blocks with the right processor.

FUNCTION kick_block(accel_x, accel_y, accel_z, blockID, block_arr, limits, dt, nparts)

integer :: kick_block, nparts, numBlocks, start_ind, end_ind
integer, dimension(nparts) :: blockID, block_arr
integer, dimension(dr_globalNumProcs) :: allProcs
integer :: limits, myProc, Proc_ID, communicator, ierr, ii, jj
integer, dimension(2,MDIM)  :: blkLimits, blkLimitsGC
integer :: l1,l2,l3,h1,h2,h3
real*8, dimension(nparts) :: accel_x, accel_y, accel_z
!real*8, dimension(numBlocks, limits, limits, limits) &
!                          :: acc3x, acc3y, acc3z
real*8, dimension(:,:,:,:), allocatable :: acc3x, acc3y, acc3z
real*8, dimension(limits, limits, limits) :: ovx,ovy,ovz
integer, dimension(4)     :: array_shape, order
real*8                    :: dt
real*8, dimension(MDIM)   :: blockCenter
real*8, pointer, dimension(:,:,:,:) :: solndata

!#ifdef GRAVITY

!call Driver_getMype(GLOBAL_COMM, myProc)
!call Driver_getComm(GLOBAL_COMM, communicator)

!call MPI_AllGather(myProc, 1, MPI_INTEGER, allProcs, 1, &
!                MPI_INTEGER, communicator, ierr)


!! Calculate the total number of blocks. This is how long the accel arrays
!! have data for.
!numBlocks = sum(block_arr(:dr_globalNumProcs))

!!print*, "num blocks =", numBlocks

!allocate(acc3x(numBlocks, limits, limits, limits))
!allocate(acc3y(numBlocks, limits, limits, limits))
!allocate(acc3z(numBlocks, limits, limits, limits))

!! We passed the accel arrays flattened. Reform them in the proper shape
!! with struct: blockID, i, j, k (cell coords).
!array_shape = (/ numBlocks, limits, limits, limits /)
!order = (/ 4, 3, 2, 1 /) ! Reconstruct using C ordering.
!!order = (/ 1, 2, 3, 4 /)  ! Reconstruct using Fortran ordering.
!! Reconstruct the arrays to match what we passed before it was flattened
!! in Python.
!acc3x = reshape(accel_x, array_shape, ORDER=order)
!acc3y = reshape(accel_y, array_shape, ORDER=order)
!acc3z = reshape(accel_z, array_shape, ORDER=order)

!!print*, shape(acc3x)

!!print*, acc3x(1,:,:,:)
!!print*, acc3x(2,5,2,3)
!!return

!start_ind = 0
!end_ind   = 0

!! Calculate the lower and upper indices for a block on this proc.

!!print*, "All procs =", allprocs

!do jj=1, dr_globalNumProcs

!    start_ind = 1 + end_ind
!    end_ind   = sum(block_arr(:jj))

!    if (myProc .eq. allProcs(jj)) then

!    !print*, "start_ind = ", start_ind, myProc
!    !print*, "end_in = ", end_ind, myProc

!    ! Loop over all the blocks in the arrays.
!        do ii=start_ind, end_ind !numBlocks

!        ! Figure out which proc we are on and which one the block is on.
!        ! This is clunky, but not sure how to do it better (yet).

!        !call Grid_getBlkCenterCoords(blockID(ii), blockCenter)
!        !call Grid_getBlkIDFromPos(blockCenter, blockID(ii), Proc_ID, communicator)

!        !if (myProc .eq. Proc_ID) then

!            ! Verify the limits are correct.
!            call Grid_getBlkIndexLimits(blockID(ii), blkLimits, blkLimitsGC, CENTER)

!            if (.NOT. (limits /= blkLimits(2,1))) then
!              print*, "kick_block: Limits given don't match actual block limits. Aborting!"
!              stop
!            end if
!            !print*, "Updating soln in block ", blockID(ii), myProc

!            l1 = blkLimitsGC(LOW,IAXIS)+NGUARD
!            h1 = blkLimitsGC(HIGH,IAXIS)-NGUARD
!            l2 = blkLimitsGC(LOW,JAXIS)+NGUARD
!            h2 = blkLimitsGC(HIGH,JAXIS)- NGUARD
!            l3 = blkLimitsGC(LOW,KAXIS)+NGUARD
!            h3 = blkLimitsGC(HIGH,KAXIS)- NGUARD

!        !    start_ind = (ii-1)*(nparts/numBlocks)+1
!        !    end_ind   = ii*(nparts/numBlocks)

!        !    print*, "Start_ind = ", start_ind, "End_ind =", end_ind

!        !    print*, acc3x(1,2,3)

!            !stop
!            ! Now actually update the velocity of the gas from the kick.
!            ! Note that we only kick the interior cells, not the guard cells.
!            call Grid_getBlkPtr(blockID(ii),solndata,CENTER)

!            !print*, "Random soln velx before = ", solndata(VELX_VAR,6,6,6)

!            !print*, "acc3x =", acc3x(ii,2,2,2)*dt
!            !print*, shape(acc3x(ii,:,:,:))
!            !print*, shape(solndata(VELX_VAR,l1:h1,l2:h2,l3:h3))
!            !print*, shape(ovx)
!            !print*, "dt = ", dt

!            !ovx = solndata(VELX_VAR,l1:h1,l2:h2,l3:h3)
!            !ovy = solndata(VELY_VAR,l1:h1,l2:h2,l3:h3)
!            !ovz = solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3)


!            solndata(VELX_VAR,l1:h1,l2:h2,l3:h3) = &
!            solndata(VELX_VAR,l1:h1,l2:h2,l3:h3) + acc3x(ii,:,:,:)*dt
!            solndata(VELY_VAR,l1:h1,l2:h2,l3:h3) = &
!            solndata(VELY_VAR,l1:h1,l2:h2,l3:h3) + acc3y(ii,:,:,:)*dt
!            solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3) = &
!            solndata(VELZ_VAR,l1:h1,l2:h2,l3:h3) + acc3z(ii,:,:,:)*dt

!            !print*, "Random soln velx after = ", solndata(VELX_VAR,6,6,6)

!            call Grid_releaseBlkPtr(blockID(ii),solndata)

!        end do

!    end if

!end do

!! This makes sure guard cells are refilled if Flash needs to do
!! interpolation or averaging after we do this.

!call Grid_notifySolnDataUpdate( (/VELX_VAR,VELY_VAR,VELZ_VAR/) )

!deallocate(acc3x)
!deallocate(acc3y)
!deallocate(acc3z)
!#endif
kick_block=0
END FUNCTION

FUNCTION kick_grid(dt)

!! Give a velocity kick from gravity of sinks for time step dt.
!! ??? Then make sure the acceleration from the sinks is zeroed. ???

integer :: kick_grid, num_blks, blk_list(MAXBLOCKS), blockID
real*8  :: dt
character(len=20) :: file_name
real*8 :: sum_force_norm
logical, parameter :: debug=.true.

!real*8, pointer, dimension(:,:,:,:) :: solndata
#ifdef GRAVITY
force_SoG_x=0.0; force_SoG_y=0.0; force_SoG_z=0.0

!call pt_sinkGatherGlobal()

!call Grid_notifySolnDataUpdate( (/VELX_VAR,VELY_VAR,VELZ_VAR/) )

!call Grid_getListOfBlocks(LEAF, blk_list, num_blks)

!call Particles_sinkKickGas(num_blks, blk_list, dt, force_SoG_x,force_SoG_y,force_SoG_z)

print*, "[kick_grid]: I do nothing now!"

!call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.false., selectBlockType=LEAF)

!!! This is now done in Particles_sinkKickGas!!!

!do blockID=1, num_blks

!    call Grid_getBlkPtr(blockID,solndata,CENTER)

!    solndata(VELX_VAR,:,:,:) = solndata(VELX_VAR,:,:,:) + solndata(SGAX_VAR,:,:,:)*dt
!    solndata(VELY_VAR,:,:,:) = solndata(VELY_VAR,:,:,:) + solndata(SGAY_VAR,:,:,:)*dt
!    solndata(VELZ_VAR,:,:,:) = solndata(VELZ_VAR,:,:,:) + solndata(SGAZ_VAR,:,:,:)*dt

!    ! Zero out the accelerations so they don't get added during the
!    ! normal hydro solution.

!    solndata(SGAX_VAR,:,:,:) = 0.0
!    solndata(SGAY_VAR,:,:,:) = 0.0
!    solndata(SGAZ_VAR,:,:,:) = 0.0

!    solndata(SGXO_VAR,:,:,:) = 0.0
!    solndata(SGYO_VAR,:,:,:) = 0.0
!    solndata(SGZO_VAR,:,:,:) = 0.0

!    call Grid_releaseBlkPtr(blockID,solndata)



!end do

!call Grid_notifySolnDataUpdate( (/SGAX_VAR,SGAY_VAR,SGAZ_VAR,SGXO_VAR,SGYO_VAR,SGZO_VAR/) )

!if (debug .and. (dr_globalMe .eq. MASTER_PE)) then

!if (dr_SimTime == 0.0) then

!file_name = "forceGoS.dat"
!open(unit=12, file=trim(file_name))
!file_name = "forceSoG.dat"
!open(unit=13, file=trim(file_name))
!file_name = "force_error.dat"
!open(unit=14, file=trim(file_name))

!else

!file_name = "forceGoS.dat"
!open(unit=12, file=trim(file_name), position='append')
!file_name = "forceSoG.dat"
!open(unit=13, file=trim(file_name), position='append')
!file_name = "force_error.dat"
!open(unit=14, file=trim(file_name), position='append')

!end if

!sum_force_norm = sqrt((abs(force_GoS_x) + abs(force_SoG_x)/ 2.0)**2.0 + &
!                      (abs(force_GoS_y) + abs(force_SoG_y)/ 2.0)**2.0 + &
!                      (abs(force_GoS_z) + abs(force_SoG_z)/ 2.0)**2.0)

!write(12,'(4(1X,E22.15))') dr_simTime, force_GoS_x,force_GoS_y,force_GoS_z
!write(13,'(4(1X,E22.15))') dr_simTime, force_SoG_x,force_SoG_y,force_SoG_z
!write(14,'(4(1X,E22.15))') dr_simTime, &
!                             abs(force_GoS_x+force_SoG_x) / sum_force_norm, &
!                             abs(force_GoS_y+force_SoG_y) / sum_force_norm, &
!                             abs(force_GoS_z+force_SoG_z) / sum_force_norm

!close(12)
!close(13)
!close(14)

!end if

!#ifdef TREE
!if (associated(loc_t)) then
!  call free_tree(loc_t)
!  end if
!call tree_build_grav(DENS_VAR)
!#endif
#endif
kick_grid=0
END FUNCTION
#ifdef TREE
FUNCTION make_particle_tree()

integer :: make_particle_tree

call create_particle_tree(particles_local(POSX_PART_PROP,1:localnp), &
                          particles_local(POSY_PART_PROP,1:localnp), &
                          particles_local(POSZ_PART_PROP,1:localnp), &
                          particles_local(MASS_PART_PROP,1:localnp), &
                          localnp)

make_particle_tree=0
END FUNCTION
#endif

FUNCTION make_stars(dt)

logical :: made_stars
real*8, intent(in) :: dt
integer make_stars

#ifdef AMUSE_STARS

!call Particles_sinkMakeStars(dt)
!call dens_removal(made_stars)
call Particles_starForm()

#endif

make_stars=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MORE GRID STUFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Here's a function that injects energy into the grid.
! It injects energy as a cube around the injection point, using
! the cloud in cell methods for particle mapping in Flash.
#ifdef ENERGY_INJ
FUNCTION energy_injection(energy, fracKin, mass, xloc, yloc, zloc, dt)

integer :: energy_injection
real*8, intent(in)  :: energy, mass, xloc, yloc, zloc
real*8, intent(inout) :: fracKin, dt

dt = 0.0

!print*, "[energy_injection]: Proc", dr_globalMe, "reports", energy, fracKin, mass, xloc, yloc, zloc, dt

call Particles_energyInjection(energy, fracKin, mass, xloc, yloc, zloc, dt)

energy_injection=0
END FUNCTION

#else
FUNCTION energy_injection(energy, fracKin, mass, xloc, yloc, zloc, dt)

integer :: energy_injection
real*8, intent(in)  :: energy, mass, xloc, yloc, zloc
real*8, intent(inout) :: fracKin, dt

print*, "WARNING: Energy injection called but not defined in interface.F90!"

energy_injection=0
return
END FUNCTION
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IO stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION write_chpt()

integer write_chpt

call IO_writeCheckpoint()

write_chpt=0
END FUNCTION

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
  call Particles_moveAndSort(.true.)

  fileNumber = mod(io_checkpointFileNumber, io_rollingCheckpoint)

else if (trim(output_type)=='pltpart') then

  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin, &
                 endRun, PLOTFILE_AND_PARTICLEFILE)
  call Particles_moveAndSort(.true.)

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
END FUNCTION

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
END FUNCTION

FUNCTION get_output_dir_wrapped(output_dir)
implicit none
integer :: get_output_dir_wrapped
character(len=40), intent(out) :: output_dir

call RuntimeParameters_get("output_directory", output_dir)
output_dir = trim(output_dir)

get_output_dir_wrapped=0
END FUNCTION

FUNCTION get_runtime_parameter(rt_name, rt_value)
implicit none
integer :: get_runtime_parameter
character(len=*), intent(in) :: rt_name
real*8, intent(out)          :: rt_value

call RuntimeParameters_get(rt_name, rt_value)

get_runtime_parameter=0
END FUNCTION

FUNCTION timer_summary()
integer :: timer_summary

call Timers_getSummary( max(0,dr_nstep-dr_nbegin+1))

timer_summary=0
END FUNCTION

FUNCTION bisect_search(x, indl, indh, n, array)
! Bisection search from Numerical Reciepes, Bill Press et al.
real*8, intent(in) :: x
integer, intent(in) :: indl, indh, n ! low, high indices and length of the array
real*8, dimension(:), intent(in) :: array
integer :: jl, jh, jm ! low, high and mid index
integer :: bisect_search
logical :: ascend

bisect_search = -1

! First just check the end, since bisection looks at the n-1 elements.

if (x == array(n)) then
    bisect_search = n
    return
end if

! If the array has no elements, return immediately.
! Also return if the array has only 1 element, since we checked
! that element above.
if (n <= 1) return

jl = indl
jh = indh

! Is the array sorted by ascending or descending order?

ascend = array(n-1) >= array(1)

! Now search array.
do while ((jh - jl) > 1)

    jm = SHIFTR( (jh+jl), 1) ! Midpoint

    if ( (x >= array(jm)) .eqv. ascend) then
        jl = jm
    else
        jh = jm
    end if

end do

! If the actual value is not found, return -1.
if (x == array(jl)) then
    bisect_search = jl
else
    bisect_search = -1
end if

return

END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parallel initialization copied from FLASH -Josh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!****if* source/Driver/DriverMain/Driver_initParallel
!!
!! NAME
!!
!!  Driver_initParallel
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing interface,
!!  the number of processors in a run and each processing
!!  element
!!
!!
!!  ARGUMENTS
!!    myPE : current processor
!!    numProcs : number of processors
!!
!!
!!
!!***

!In mpif.h MPI_VERSION is an integer (and thus can't be used to conditionally
!compile code) and so we allow the user to define FLASH_MPI1 to indicate
!that they have an MPI-1 implementation.

#ifdef _OPENMP
#ifndef FLASH_MPI1
#define FLASH_MPI2_OPENMP
#endif
#endif

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
