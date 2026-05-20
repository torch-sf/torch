!===============================================================================
!
!  Subroutine: calc_local_block_contributions_sink
!
!===============================================================================
!
! Calculates local column densities for all centres of all zones of
! all leaf blocks. This routine is specific to ray-tracing from point sources.
! Uses the 'Fast Voxel Traversal Algorithm for Ray Tracing' by
! J. Amanatides and Andrew Woo, Proc. Eurographics, Aug. 1987
! (http://www.cs.yorku.ca/~amana/research/grid.pdf)
! See also Abel, Norman & Madau  ApJ 523, 66.
! Author: Manuel Jung with mods by Shyam Menon (2022)
!===============================================================================

MODULE calc_local_sink_mod
  use tree, only: lrefine, grid_xmin, grid_xmax, grid_ymin, grid_ymax, &
       grid_zmin, grid_zmax
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords, &
    Grid_getDeltas
  use RadTrans_hybridCharModule, ONLY: equal, qdr_bezier, check_incoords, NB, KD, obeybox, &
    valid, obeyintervall, precision
  use raytrace_data
  use fvt, ONLY : FastVoxelTraversal, intersection, &
    fvt_init, fvt_check, fvt_valid, fvt_get, fvt_dstep, fvt_checkLT

IMPLICIT NONE
PRIVATE
  INTEGER, PARAMETER :: NDIM = 3
  real                              :: opac
  real, dimension(NDIM)             :: x0, x1
  integer, dimension(NDIM)          :: b0, b1
  integer                           :: state
  real, dimension(NDIM)             :: rescale, rescale_inv
  real, dimension(NDIM)             :: direction
  real                              :: dtold, told, tend, tbegin
  real                              :: taud
  real, dimension(NDIM) :: rtmin, rtmax, dx, bmin, bmax, regFact
  real :: levelFact
  real, DIMENSION(NDIM) :: grid_min, grid_max
  integer :: face

  PUBLIC :: calc_local_block_contributions_sink

CONTAINS

subroutine calc_local_block_contributions_sink( &
     regFact_, blk,  &
     src, point, direction_,               &
     maxLevel,                             &
     taud_,face_, dtaucell_,dlcontrib_,     &
     debug_rt)

#undef VERBOSE_RT

#include "Flash.h"
#include "constants.h"

!
  implicit none
!
  integer, intent(in)                       :: maxLevel
  real,    intent(in), dimension(NDIM)      :: regFact_, point, src, direction_
  real,    intent(out)                      :: taud_, dtaucell_, dlcontrib_
  integer, intent(in)                       :: face_ ! <> 0 -> face values (1: x-face, 2: y-face, 3: z-face). Else: centers
  integer, intent(in) :: blk
  logical, intent(in) :: debug_rt

!
  integer, parameter    :: q = MAXCELLS
  real, dimension(q,NDIM) :: x
  real :: norm
  integer               :: n
  REAL :: t
  INTEGER, DIMENSION(NDIM) :: voxel, dir
  INTEGER :: cut, stopcnt, i
  REAL :: ttemp
  REAL, DIMENSION(NDIM) :: b0c, b1c
  real :: length, dtaucell, lcontrib, dlcontrib, sinkpos(NDIM)
  logical :: printval
  REAL,DIMENSION(NDIM) :: target_point
  !-------------------------------------------------------------------------------

  grid_min = (/grid_xmin, grid_ymin, grid_zmin/)
  grid_max = (/grid_xmax, grid_ymax, grid_zmax/)
  regFact = regFact_
  face = face_

  norm = SQRT(SUM(direction_**2))
  IF(norm.le.0.) THEN
    taud_ = 0.
    dtaucell_ = 0.
    dlcontrib_ = 0.
    return
  END IF

  direction = direction_ / norm

  !Reset values to zero
  dtold = 0.
  told = 0.
  state = 0
  opac = 0.

!
!  Calculate the correction factor when the current block
!  is not at the finest level of refinement.
!
  levelFact = 2.**(lrefine(blk)-maxLevel)

!  Get the physical coordinates and calculate the min and max
!  regular coordinates of this block.
!
  call Grid_getDeltas(blk, dx)
  do n=1,NDIM
    call Grid_getCellCoords(n, blk, CENTER, .TRUE., x(:,n), MAXCELLS)

!    Get the physical coordinates and calculate the min and max
!    coordinates of this block. We use the corner of the first ghost cell as the starting point.
!    This are then the boundaries of the staggered RT mesh.
!
    rtmin(n) = x(NGUARD-1,n) !This is cell centre of second guard cell in the -ve dir
    rtmax(n) = x(NGUARD*KD(n)+NB(n)+2,n) !This is cell centre of second guard cell in the +ve dir
  end do

  !rt min and max at the corner of the first guard cell in both dirs
  rtmin = rtmin+0.5*dx 
  rtmax = rtmax-0.5*dx

  ! Bounding box of 'active' domain
  bmin = rtmin + 1*dx
  bmax = rtmax - 1*dx

  !TODO: Write this correctly
  if(debug_rt) then 
    print *, 'Debug mode'
  endif

    ! Reset all traced quantities

  taud     = 0.0
  dtaucell = 0.0
  dlcontrib = 0.0
  length = 0.0 
  lcontrib = 0.0

  !Rescale physical length to FVT voxel length
  rescale = (regFact * levelFact)
  !Inverse scaling
  rescale_inv = 1./rescale
  ! Add a small t value, so intersection does not find the point itself.
  ! -> intersection searches for t>=0
  x1 = (point - rtmin) * rescale !Final destination of FVT
  !Defining the voxel indices of the FVT grid
  b0 = 0 
  b1 = NINT((rtmax-rtmin) * rescale)

  call obeybox(b0,b1,x1) !Ensure that the dest is bounded in FVT grid

  !Get the point, x0, where the ray enters the FVT voxel grid.
  ttemp  = intersection(-direction, x1, REAL(b0), REAL(b1))
  !This x0 would set the limits of the FVT grid
  x0 = ttemp * (-direction) + x1 !Ray enters x0

  !Safety checks
  WHERE(equal(x0,REAL(b0)+1.0)) x0 = REAL(b0)+1.0
  WHERE(equal(x0,REAL(b1)-1.0)) x0 = REAL(b1)-1.0

  sinkpos = src
  !In case sink is at a common face/corner of blocks, precision errors can lead to
  ! it being missed when determining whether the sink is in the block or not.
  ! The below prevents this
  if((abs(src(1)-bmin(1)) * rescale(1)) .lt. precision) sinkpos(1) = bmin(1)
  if((abs(src(2)-bmin(2)) * rescale(2)) .lt. precision) sinkpos(2) = bmin(2)
  if((abs(src(3)-bmin(3)) * rescale(3)) .lt. precision) sinkpos(3) = bmin(3)
  if((abs(src(1)-bmax(1)) * rescale(1)) .lt. precision) sinkpos(1) = bmax(1)
  if((abs(src(2)-bmax(2)) * rescale(2)) .lt. precision) sinkpos(2) = bmax(2)
  if((abs(src(3)-bmax(3)) * rescale(3)) .lt. precision) sinkpos(3) = bmax(3)

  !Check if sink lies within block
  !TODO: Maybe there is a better way to check if the sink lies within the block?
  if(check_incoords(bmin, bmax, sinkpos)) then
    ! The sink is located inside this block.
    ! We still start at the boundary of this block as usual
    ! since this is the easiest, but only start integrating
    ! at the source itself.
    ! This way we do more cuts than necessary, but
    ! the other blocks do a similar amount of cuts, so
    ! it does not matter.
    tbegin = SQRT(SUM(((sinkpos-rtmin)*rescale-x0)**2))
  else
    !Sink not in block. Start FVT from the point where the ray enters the 'active' block
    tbegin  = intersection(direction, x0, REAL(b0)+1.0, REAL(b1)-1.0)
  endif

  ! This is necessary, since otherwise there may be small deviations to
  ! the internal tend of FastVoxelTraversalLT ('distance')
  tend = SQRT(SUM((x1-x0)**2)) !Final destination t value

  cut = 0
  state = 0

  !Start FVT only if tend>tbegin
  if(tend.gt.tbegin.and..not.equal(tend,tbegin)) then
    call fvt_init(x0, x1, b0, b1) !Initialise FVT
    stopcnt = 0
    !CheckLT checks if current spearhead of ray has crossed the destination point
    if(.not.fvt_checkLT()) &
      stopcnt = stopcnt + 1

    !The loop through the FVT grid starts here
    do while(stopcnt.le.1)
      !Note that the voxel traversal and optical depth accumulation occurs in a delayed fashion
      !First a step is taken from t0 to t1. Then the contribution is accumulated for the step from t-1 to t0.
      !By the end of the loop, the accumulation upto the previous voxel has been done.
      !This means that the loop would end 2 steps after the ray spearhead crosses the dest.
      !This is why the loop stops if stopcnt goes greater than 1; this allows 2 extra steps

      call fvt_dstep() !Take an FVT step, i.e. t0 to t1
      call fvt_get(voxel,dir,t) !Get the voxel corresponding to the current spearhead
      !Get the optical depth traversed in the last step, i.e. from t-1 to t0
      call onNextVoxelSink(voxel,dir,t,taud,opac,told,dtold,state,cut,length,dtaucell,lcontrib,debug_rt)
      !Check if crossed
      if(.not.fvt_checkLT()) &
        stopcnt = stopcnt + 1

      !Do the following - 
      ! i) take optical depth to last cell centre as optical depth 'upto' cell
      ! ii) Compute the optical depth step taken from the last to this cell as the absorption optical depth
      ! iii) Use the segment length from the last step of the FVT as the cell length intersection
      if(stopcnt .gt. 1 .and. face .eq. 0) then
        !Move back to the last cell centre with the FVT
        taud = taud - dtaucell
        taud = MAX(taud,0.0)
        dlcontrib = lcontrib * rescale_inv(1)
      endif

      if(debug_rt) & 
        print *,'voxel,t,tend,taud,length,stopcnt,dtaucell,dlcontrib,lcontrib', voxel, t,tend, taud, dtold,stopcnt,dtaucell, &
        & dlcontrib,lcontrib

    end do
    !End voxel traversal loop
  end if

  ! In case of face destinations we normalize the optical depth to the ray
  ! length. Therefore we don't need to calculate the length again, when we
  ! recombine this ray parts with interpolation in create_cut_block_list,
  ! where it is rescaled to the actual ray length after/while interpolation.
  if(face.ne.0.and. exp(-taud).ne. 1.) then
    ! normalize by ray length in units of number of  cells at highest refinement
    ! level
    taud = taud*levelFact/length
  end if

  !Set final returned values
  taud_ = taud
  if(face .eq. 0) then 
    dtaucell_ = dtaucell
    dlcontrib_ = dlcontrib
  else
    dtaucell_ = 0.0
    dlcontrib_ = 0.0
  endif
end subroutine calc_local_block_contributions_sink

SUBROUTINE onNextVoxelSink(voxel,dir,t,taud,opac,told,dtold,state,cut,length,dtau,dlength,debug_rt)

  IMPLICIT NONE
  INTEGER, DIMENSION(NDIM), INTENT(IN) :: voxel, dir
  REAL, INTENT(IN) :: t
  REAL, INTENT(INOUT) :: taud
  REAL, INTENT(INOUT) :: opac
  REAL, INTENT(INOUT) :: told, dtold, length, dtau, dlength
  INTEGER, INTENT(INOUT) :: state
  INTEGER, INTENT(INOUT) :: cut
  LOGICAL, INTENT(IN) :: debug_rt
  REAL, DIMENSION(NDIM) :: ray
  REAL :: dt
  REAL :: a,b

  !Update voxel cut counter
  cut = cut + 1

  ! This is our current spearhead of the ray.
  ray = t * direction + x0
  !Step size of latest step
  dt = t - told

  if(debug_rt) print *, 'OPAC_VAR of Voxel,',voxel,' = ',solnData(OPAC_VAR,voxel(1)+NGUARD,voxel(2)+NGUARD,voxel(3)+NGUARD)
  if(debug_rt) print *, 'OPAC_VAR of last voxel, resP = ',opac

  !Only add contribution if the starting point has been crossed
  IF(tbegin.le.told) THEN
    !Catch cases where the first step would be partial, i.e. not = a full voxel step
    !This would almost always be the case, as long as the sink is not directly at the corner
    if(state.eq.0) then
      !a is the fraction of voxel length beyond which the sink lies
      a = 1.0 - (told - tbegin) / dtold
      state = 1
    else
      a = 0.0
    end if
    !Only do this
    if(tend.le.told.and.face.ne.0) then
      b = (told - tend) / dtold
      state = 2
    else
      b = 1.0
    end if


    !Calculate optical depth traversed in last step
    dtau = opac * (b-a)* dtold * rescale_inv(1)
    !The effective length in the last step
    dlength = (b-a)*dtold
    !Add this to cumulative counters
    length = length + dlength
    taud = taud + dtau

    if(debug_rt) print *,'opac,a,b,dtau', solnData(OPAC_VAR,voxel(1)+NGUARD,voxel(2)+NGUARD,voxel(3)+NGUARD), a, b, dtau

  END IF

  !TODO: Confirm this
  opac = solnData(OPAC_VAR,voxel(1)+NGUARD,voxel(2)+NGUARD,voxel(3)+NGUARD)
  dtold = dt
  told = t

END SUBROUTINE onNextVoxelSink

END MODULE calc_local_sink_mod
!===============================================================================







