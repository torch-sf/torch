!===============================================================================
!
!  Subroutine: calc_local_block_contributions_3DRT
!
!===============================================================================
!
! Calculates local column densities for all centres of all zones of
! all leaf blocks.
! Uses the 'Fast Voxel Traversal Algorithm for Ray Tracing' by
! J. Amanatides and Andrew Woo, Proc. Eurographics, Aug. 1987
! (http://www.cs.yorku.ca/~amana/research/grid.pdf)
! See also Abel, Norman & Madau  ApJ 523, 66.
!
!===============================================================================

MODULE calc_local_mod
  use tree, only: lrefine, grid_xmin, grid_xmax, grid_ymin, grid_ymax, &
       grid_zmin, grid_zmax
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords, &
    Grid_getDeltas
  use RadTrans_hybridCharModule, ONLY: equal, qdr_bezier, check_incoords, NB, KD, obeybox, &
    valid, obeyintervall
  use raytrace_data
  use fvt, ONLY : FastVoxelTraversal, intersection, &
    fvt_init, fvt_check, fvt_valid, fvt_get, fvt_dstep, fvt_checkLT

IMPLICIT NONE
PRIVATE
  INTEGER, PARAMETER :: NDIM = 3
  real, dimension(2)                :: resP, res, resN
  real, dimension(NDIM)             :: x0, x1
  integer, dimension(NDIM)          :: b0, b1
  integer                           :: state
  real, dimension(NDIM)             :: rescale, rescale_inv
  real                              :: dtau1,dtau2
  real, dimension(NDIM)             :: direction
  real                              :: dtold, told, tend, tbegin
  real                              :: u,v,w
  real                              :: taud, intensity
  real, dimension(NDIM) :: rtmin, rtmax, dx, bmin, bmax, regFact
  real :: levelFact
  real, DIMENSION(NDIM) :: grid_min, grid_max
  logical :: sink
  integer :: face
  real :: lambda_akku

  PUBLIC :: calc_local_block_contributions_3DRT

CONTAINS

subroutine calc_local_block_contributions_3DRT( &
     regFact_, blk,  &
     src, point, direction_,               &
     maxLevel,                             &
     taud_, intensity_,                    &
     lambda_akku_,                         &
     face_, &
     sink_,dtaucontrib,dlcontrib)

#undef VERBOSE_RT
!PPD for debugging
#undef DEBUG_RT

#include "Flash.h"
#include "constants.h"

!
  implicit none
!
  integer, intent(in)                       :: maxLevel
  real,    intent(in), dimension(NDIM)      :: regFact_, point, src, direction_
  real,    intent(out)                      :: taud_, intensity_, dtaucontrib, dlcontrib
  real,    intent(inout)                    :: lambda_akku_
  integer, intent(in)                       :: face_ ! <> 0 -> face values (1: x-face, 2: y-face, 3: z-face). Else: centers
  integer, intent(in) :: blk
  logical, intent(in) :: sink_

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
  real :: length, dtaucell, lcontrib
#ifdef DEBUG_RT
  logical :: printval
  REAL,DIMENSION(NDIM) :: target_point
#endif
  !-------------------------------------------------------------------------------

! print all inputs:
!  print *, "regFact:   ",regFact
!  print *, "levelFact: ", levelFact
!  print *, "dx:        ", dx
!  print *, "x:         ", x
!  print *, "blk:       ", blk
!  print *, "point:     ", point
!  print *, "direction: ", direction_
!  print *, "maxLevel:  ", maxLevel
!  print *, "dIdt_fact: ", dIdt_fact
!  print *, "face:      ", face

  grid_min = (/grid_xmin, grid_ymin, grid_zmin/)
  grid_max = (/grid_xmax, grid_ymax, grid_zmax/)
  regFact = regFact_
  sink = sink_
  face = face_

!  IF(.not.ALL(point(1:2).le.0.008)) return

  norm = SQRT(SUM(direction_**2))
  IF(norm.le.0.) THEN
    taud_ = 0.
    intensity_ = 0.
    lambda_akku_ = 0.
    return
  END IF

  direction = direction_ / norm

  dtold = 0.
  told = 0.
  state = 0
  res = 0.
  resP = 0.
  resN = 0.

!
!  Calculate the correction factor when the current block
!  is not at the finest level of refinement.
!
  levelFact = 2.**(lrefine(blk)-maxLevel)
!
!  Get the physical coordinates and calculate the min and max
!  regular coordinates of this block.
!
  call Grid_getDeltas(blk, dx)
  do n=1,NDIM
    call Grid_getCellCoords(n, blk, CENTER, .TRUE., x(:,n), MAXCELLS)

!    Get the physical coordinates and calculate the min and max
!    coordinates of this block. We use two ghost cells.
!    This are then the boundaries of the staggered RT mesh.
!
    rtmin(n) = x(NGUARD-1,n)
    rtmax(n) = x(NGUARD*KD(n)+NB(n)+2,n)
  end do
  bmin = rtmin + 1.5*dx
  bmax = rtmax - 1.5*dx

!Debugging option to print values for a given point
#ifdef DEBUG_RT
  target_point = (/4.3999999999999974E+017, 0.0, 0.0/)
  if(point(1) .eq. target_point(1) .and. point(2) .eq. target_point(2) .and. point(3) .eq. target_point(3)) then 
    printval = .true.
    print *, 'Debugging and printing info for ',target_point
  else
    printval = .false.
  endif
#endif

#ifdef VERBOSE_RT
     write(6,*) '----------------------------------------------'
     write(6,*) 'starting raytracing on block', blk
     write(6,*) 'lrefine',lrefine(blk)
     write(6,*) '----------------------------------------------'
#endif
  ! Reset all traced quantities

  taud        = 0.0
  dtau1       = 0.0
  dtau2       = 0.0

  u           = 0.0
  v           = 0.0
  w           = 0.0

  intensity   = 0.0
  lambda_akku = 0.0
  dtaucontrib = 0.0
  dlcontrib = 0.0


  length = 0.

  rescale = (regFact * levelFact)
  rescale_inv = 1./rescale
  ! Add a small t value, so intersection does not find the point itself.
  ! -> intersection searches for t>=0
  x1 = (point - rtmin) * rescale
  b0 = 0
  b1 = NINT((rtmax-rtmin) * rescale)

  call obeybox(b0,b1,x1)
  ttemp  = intersection(-direction, x1, REAL(b0)+0.5, REAL(b1)-0.5)
  !x0 = (tend+0.0001) * (-direction) + x1
  x0 = ttemp * (-direction) + x1
  !WHERE(equal(x0,NINT(x0))) x0 = NINT(x0)
  WHERE(equal(x0,REAL(b0)+1.5)) x0 = REAL(b0)+1.5
  WHERE(equal(x0,REAL(b1)-1.5)) x0 = REAL(b1)-1.5

 ! call obeybox(REAL(b0)+0.5,REAL(b1)-0.5,x0)
  if(sink.and.check_incoords(bmin, bmax, src)) then
    ! The sink is located inside this block.
    ! We still start at the boundary of this block as usual
    ! since this is the easiest, but only start integrating
    ! at the source itself.
    ! This way we do more cuts than necessary, but
    ! the other blocks do a similar amount of cuts, so
    ! it does not matter.
    tbegin = SQRT(SUM(((src-rtmin)*rescale-x0)**2))
  else
    tbegin  = intersection(direction, x0, REAL(b0)+1.5, REAL(b1)-1.5)
  end if

  ! This is necessary, since otherwise there may be small deviations to
  ! the internal tend of FastVoxelTraversalLT ('distance')
  tend = SQRT(SUM((x1-x0)**2))

  ! If we calculate a face destination point and
  ! tbegin==tend, this is a corner destination point. It is not
  ! reached by any part of the ray actually cutting the block.
  ! We still need an approximation of the values at this point.
  ! Therefore we include one ghost cell to calculate a ray, which
  ! has the same length as the first ray, which actually cuts the
  ! block.
  ! For this we enlarge the intersection box in all directions, but the
  ! face direction.
!  if(face.ne.0.and.(tbegin.ge.tend*0.999999)) then
!    b0c = REAL(b0)
!    b1c = REAL(b1)
!    ttemp  = intersection(-direction, x1, b0c-0.5, b1c+0.5)
!    x0 = ttemp * (-direction) + x1
!    DO i=1,NDIM
!      IF(i.eq.face) THEN
!        b0c(i) = b0c(i) + 1.5
!        b1c(i) = b1c(i) - 1.5
!      ELSE
!        b0c(i) = b0c(i) + 0.5
!        b1c(i) = b1c(i) - 0.5
!      END IF
!    END DO
!    WHERE(equal(x0,b0c)) x0 = b0c
!    WHERE(equal(x0,b1c)) x0 = b1c
!    tbegin  = intersection(direction, x0, b0c, b1c)
!    tend = SQRT(SUM((x1-x0)**2))
!  end if

!  print *, "New ray --------------------------------------------"
!  print *, "b0: ", b0
!  print *, "b1: ", b1
!  print *, "regFact, levelFact: ", regFact, levelFact
!  print *, "rescale, rescale_inv: ", rescale, rescale_inv
!  print *, "rtmin: ", rtmin
!  print *, "x0: ", x0
!  print *, "x0 rescaled: ",x0*rescale_inv + rtmin
!  print *, "x1: ", x1
!  print *, "x1 rescaled: ",x1*rescale_inv + rtmin
!  print *, "dirFact: ", dirFact
!  print *, "tbegin: ", tbegin
!  print *, "tend:   ", tend

!  if(tend.gt.tbegin) &
!    CALL FastVoxelTraversal(x0, x1, b0, b1, .TRUE., onNextVoxel)
  cut = 0
  state = 0
  if(tend.gt.tbegin.and..not.equal(tend,tbegin)) then
    call fvt_init(x0, x1, b0, b1)
    stopcnt = 0
    if(.not.fvt_checkLT()) &
      stopcnt = stopcnt + 1
    do while(stopcnt.le.1)
      call fvt_dstep()
      call fvt_get(voxel,dir,t)
!      print *,"dir, t, voxel", dir(3), t, voxel(3)
#ifdef CHECKS_RT
      IF(ALL(dir.eq.0)) THEN
        print *,"dir: ", dir
        print *,"voxel: ", voxel
        print *,"t: ", t
        print *, "tbegin: ", tbegin
        print *, "tend:   ", tend
        call Driver_abortFlash('calc_local_block_contributions_3DRT: dir is zero, this should not happen!')
      END IF
#endif
      call onNextVoxel(voxel,dir,t,taud,intensity,resP,res,told,dtau1,lambda_akku,dtold,state,cut,length,dtaucell,lcontrib)
      if(.not.fvt_checkLT()) &
        stopcnt = stopcnt + 1

      !Do the following - 
      ! i) take optical depth to last cell centre as optical depth 'upto' cell
      ! ii) Compute the optical depth step taken from the last to this cell as the absorption optical depth
      ! iii) Use the segment length from the last step of the FVT as the cell length intersection
      if(stopcnt .gt. 1 .and. face .eq. 0 .and. sink) then
        !Move back to the last cell centre with the FVT
        taud = taud - dtaucell
        taud = MAX(taud,0.0)
        dlcontrib = lcontrib * rescale_inv(1)
      endif

      !Debugging option
#ifdef DEBUG_RT
      if(printval) &
        print *,'2:', voxel, t, taud, dtold, dtau1, length,stopcnt,dtaucell, dlcontrib
#endif
    end do
  end if

!#ifdef CHECKS_RT
!  if(taud.le.0.) then
!    print *, "taud: ", taud
!    print *, "length: ", length
!    print *, "face: ", face
!    print *, "b0: ", b0
!    print *, "b1: ", b1
!    print *, "regFact, levelFact: ", regFact, levelFact
!    print *, "rescale, rescale_inv: ", rescale, rescale_inv
!    print *, "rtmin: ", rtmin
!    print *, "x0: ", x0
!    print *, "x0 rescaled: ",x0*rescale_inv + rtmin
!    print *, "x1: ", x1
!    print *, "x1 rescaled: ",x1*rescale_inv + rtmin
!    print *, "dirFact: ", direction
!    print *, "tbegin: ", tbegin
!    print *, "tend:   ", tend
!    call Driver_abortFlash('calc_local_block_contributions_3DRT: taud <=0, this should not happen!')
!  end if
!#endif

  ! In case of face destinations we normalize the optical depth to the ray
  ! length. Therefore we don't need to calculate the length again, when we
  ! recombine this ray parts with interpolation in create_cut_block_list,
  ! where it is rescaled to the actual ray length after/while interpolation.
  if(face.ne.0.and. exp(-taud).ne. 1.) then
    intensity = intensity / (1.-exp(-taud))
    ! normalize by ray length in units of number of  cells at highest refinement
    ! level
    taud = taud*levelFact/length
  end if

!  if(.not.face.and.dir(3).gt.0.) then
!    print *, "---------------------------"
!    print *, "length: ", length
!    print *, "tend-tbegin: ", tend-tbegin
!    print *, "tbegin: ", tbegin
!    print *, "tend:   ", tend
!    print *,"voxel: ", voxel
!    print *, "x0: ", x0
!    print *, "x1: ", x1
!    print *, "dirFact: ", dirFact
!    print *, "t:", t
!    print *, "cut: ", cut
!    print *, "taud: ", taud
!    print *, "face: ", face
!  end if

#ifdef CHECKS_RT
  IF(cut.gt.0.and..not.((face.ne.0.and.state.eq.2).or.(face.eq.0.and.state.lt.2))) THEN
    print *, "state: ", state
    call Driver_abortFlash('calc_local_block_contributions_3DRT: wrong state!')
  END IF
  if(tend.gt.tbegin.and..not.equal(tend,tbegin).and..not.equal(length,tend-tbegin)) THEN
    print *, "length: ", length
    print *, "tend-tbegin: ", tend-tbegin
    print *, "tbegin: ", tbegin
    print *, "tend:   ", tend
    print *,"voxel: ", voxel
    print *, "x0: ", x0
    print *, "x1: ", x1
    print *, "dirFact: ", direction
    print *, "t:", t
    print *, "cut: ", cut
    print *, "taud: ", taud
    print *, "face: ", face
    call Driver_abortFlash('calc_local_block_contributions_3DRT: Integrated ray length is wrong!')
  END IF
!  IF(.not.sink.and.intensity.LE.0.AND.tend.gt.tbegin.and..not.equal(tend,tbegin)) THEN
!    print *, "cutno: ", cut
!    print *, "b0: ", b0
!    print *, "b1: ", b1
!    print *, "x0: ", x0
!    print *, x0*rescale_inv + rtmin
!    print *, "x1: ", x1
!    print *, x1*rescale_inv + rtmin
!    print *, "tbegin: ", tbegin
!    print *, "tend:   ", tend
!    print *, "t:      ", t
!    print *, "sink:   ", sink
!    call Driver_abortFlash('calc_local_block_contributions_3DRT: intensity zero!')
!  END IF
#endif
#ifdef CHECKS_RT
    if(v.lt.0..or.v.ge.1.) then
      print *, "voxel: ", voxel
      print *, "dir: ", dir
      print *, "t: ", t
!      print *, "dt: ", dt
      print *,"lambda: ", v
      print *,"res(2): ", res(2)
      print *,"resN(2): ", resN(2)
      print *,"dtau_prev: ", dtau2
      call Driver_abortFlash('calc_local_block_contributions_3DRT: lambda < 0.')
    end if

#endif


!MJ: Not sure about this part. This is necessary, because otherwise there
!    is no contribution in the 0th/1th cell.
!   if(face) then
!     taud = taud + dtau1
!
!     if(dtau1>0.0+precision) then
!        !dI = qdr_olsonkunasz(sourceP,source,0.0,dtau1,0.0,2,u,v,w)
!        !dI = qdr_olsonkunasz(sourceP,source,source,dtau1,dtau1,2,u,v,w)
!        dI = qdr_olsonkunasz(resP(2),res(2),res(2),dtau1,dtau1,2,u,v,w)
!     else
!        dI = 0.0
!     endif
!
!     intensity = intensity * exp(-dtau1) + dI
!
!     IF(taud.EQ.0.0) THEN
!       !print *,"taud zero, zone: ", zone
!       call Driver_abortFlash('calc_local_block_contributions_3DRT: '&
!         // 'taud zero!')
!     END IF
!   end if

  taud_ = taud
  intensity_ = intensity
  lambda_akku_ = lambda_akku
  if(face .eq. 0) then 
    dtaucontrib = dtaucell
  else
    dtaucontrib = 0.0
  endif
end subroutine calc_local_block_contributions_3DRT

#ifndef CHECKS_RT
pure &
#endif
SUBROUTINE onNextVoxel(voxel,dir,t,taud,intensity,resP,res,told,dtau1,lambda_akku,dtold,state,cut,length,dtau,dlength)
#ifdef ERAD_VAR
  use RadTrans_data, ONLY: rt_hydro_type
#endif
  IMPLICIT NONE
  INTEGER, DIMENSION(NDIM), INTENT(IN) :: voxel, dir
  REAL, INTENT(IN) :: t
  REAL, INTENT(INOUT) :: taud, intensity
  REAL, DIMENSION(2), INTENT(INOUT) :: resP, res
  REAL, DIMENSION(2) :: resN
  REAL, INTENT(INOUT) :: told, dtau1, lambda_akku, dtold, length, dtau, dlength
  INTEGER, INTENT(INOUT) :: state
  INTEGER, INTENT(INOUT) :: cut
  REAL :: dtau2, u, v, w, dtautmp
  REAL, DIMENSION(NDIM) :: ray
  REAL :: dI, dImax
  INTEGER, DIMENSION(2) :: vars
  REAL :: dt
  REAL :: a,b

  cut = cut + 1

  if(sink) then 
   vars= (/ OPAC_VAR, SOUR_VAR /)
  else 
#ifdef ERAD_VAR
   if(rt_hydro_type .eq. 1) then 
    vars = (/ TAUR_VAR, SOUR_VAR /)
   else
    vars = (/ TAUP_VAR, SOUR_VAR /)
   endif
#else
   vars = (/ TAUP_VAR, SOUR_VAR /)
#endif
  endif

  !Use stellar opacities if sink else use planck mean opacity
  

#ifdef CHECKS_RT
  IF(.not.ALL(-1.LE.voxel.AND.voxel.LE.10)) THEN
    print *,"Voxel: ", voxel
    print *, "x0: ", x0
    print *, "x0 rescaled: ",x0*rescale_inv + rtmin
    print *, "x1: ", x1
    print *, "x1 rescaled: ",x1*rescale_inv + rtmin
    print *, "dir: ", dir
    print *, "tbegin: ", tbegin
    print *, "t: ", t
    print *, "dt: ", dt
    print *, "tend:   ", tend
    call Driver_abortFlash('calc_local_block_contributions_3DRT: voxel out of range [0,9]')
  END IF
#endif

  ! This is our current spearhead of the ray.
  ray = t * direction + x0
  dt = t - told

!if(dir(3).gt.0.and..not.face) then
!  print *,"Voxel: ", voxel
!  print *,"dir:   ", dir
!  print *,"t:     ", t
!  print *,"dt:    ", dt
!  print *,"ray    ", ray
!  print *, "dt,tbegin,told,t,tend: ", dt,tbegin,told,t,tend
!end if

  call interpolate(voxel, dir, ray, vars, resN)

  dtau2 = 0.5 * (resN(1) + res(1)) * dt * rescale_inv(1)

  IF(tbegin.le.told) THEN
#ifdef CHECKS_RT
    if(dtold.le.0..or.cut.lt.3) then
      print *,"cut: ", cut
      print *,"told, t: ", told, t
      print *,"dtold, dt: ", dtold, dt
      print *,"tbegin, tend: ", tbegin, tend
      print *,"x0: ", x0
      print *,"x1: ", x1
      print *,"direction: ", direction
      print *,"bmin: ", bmin
      print *,"bmax: ", bmax
      call Driver_abortFlash('calc_local_block_contributions_3DRT: dtold <= 0. or cut < 3')
    end if
#endif
    if(state.eq.0) then
      a = 1.0 - (told - tbegin) / dtold
    else
      a = 0.0
    end if
    if(tend.le.told.and.face.ne.0) then
#ifdef CHECKS_RT
      IF(state.eq.2) THEN
        call Driver_abortFlash('calc_local_block_contributions_3DRT: state is already 2. Should not happen.')
      END IF
#endif
      b = (told - tend) / dtold
      state = 2
    else
      b = 1.0
    end if
#ifdef CHECKS_RT
    if(.not.((0..le.a).and.(a.le.b).and.(b.le.1.))) then
      print *,"state: ", state
      print *,"a, b: ", a, b
      print *,"tbegin, tend, told, dtold: ", tbegin, tend, told, dtold
      call Driver_abortFlash('calc_local_block_contributions_3DRT: 0 <= a <= b <= 1 is not fullfilled.')
    end if
#endif
    !write(*,"(6F6.2)") tbegin, tend, told, dtold, a, b
    !Second order scheme
    dtau = (resP(1) * (b-a) + 0.5*(res(1)-resP(1))*(b**2-a**2)) * dtold*rescale_inv(1)
    !First order scheme
    ! dtau = resP(1) * (b-a)* dtold * rescale_inv(1)
    dlength = (b-a)*dtold
    length = length + dlength
    taud = taud + dtau
!    dtau = dtau * rescale_inv(1)
!    if(dir(3).gt.0) &
!      print *,"taud: ",taud

    ! We have valid resP, res and resN values => 3rd order integration
#ifdef CHECKS_RT
    IF(dtau1.LE.0.OR.dtau2.LE.0) THEN
      print *, "cutno: ", cut
      print *,"voxel: ", voxel
      print *,"dtau1: ", dtau1
      print *,"dtau2: ", dtau2
      print *, "resN(1), res(1): ", resN(1), res(1)
      print *, "rescale_inv(1): ", rescale_inv(1)
      print *, "x0: ", x0
      print *, "x1: ", x1
      print *, "t:", t
      print *, "tbegin: ", tbegin
      print *, "tend: ", tend
      print *, "told: ", told
      print *, "dt:   ", dt
      print *, "dtold:", dtold
      print *, "sour: ", resP(2), res(2), resN(2)
      call Driver_abortFlash('calc_local_block_contributions_3DRT: dtau1 or dtau2 <= 0.')
    END IF
#endif

    ! For the quadrature we have to use integration limits
    ! based on optical depths and not based on distance.
    ! Therefore we modify them, if we have to.
    SELECT CASE(state)
    CASE(0)
      state = 1
      a = 1.-dtau/dtau1
    CASE(2)
      if(a.gt.0.) then
        ! cut dtau1: dtau1 = dtautemp + dtau + dtaurest
        dtautmp = (resP(1) * (a-0.) + 0.5*(res(1)-resP(1))*(a**2-0.**2)) * dtold*rescale_inv(1)
        a = dtautmp/dtau1
        b = (dtautmp + dtau)/dtau1
      else
        b = dtau/dtau1
      end if
    END SELECT

    if(.not.sink) then
      ! Sinks only need the optical depth

      !call qdr_olsonkunasz(resP(2),res(2),resN(2),dtau1,dtau2,4,a,b,u,v,w,dI)
      call qdr_bezier(resP(2),res(2),resN(2),dtau1,dtau2,a,b,u,v,w,dI)

      ! quadrature limiter in case the gradients of the source function and optical depth
      ! have opposite signs (use for order=2 or 3 of qdr_olsonkunasz
      !if(PRODUCT(res-resP).lt.0.) then
      dImax = 0.5d0*(PRODUCT(resP)+PRODUCT(res)) * dtold * rescale_inv(1)
      dI = min(dI,dImax)
      !end if

      ! finally, compute local intensity contribution for the current characteristic
      intensity = intensity * exp(-dtau) + dI
      lambda_akku = v
    end if
  END IF

  dtold = dt
  told = t

  resP = res
  res  = resN
  dtau1 = dtau2

END SUBROUTINE onNextVoxel


#ifndef CHECKS_RT
pure &
#endif
SUBROUTINE interpolate(voxel, dir, ray_, vars, res)
IMPLICIT NONE
  integer, dimension(NDIM), intent(in) :: voxel, dir
  real, dimension(NDIM), intent(in) :: ray_
  integer, intent(IN), dimension(:) :: vars
  real, intent(OUT),dimension(SIZE(vars)) :: res
  real, dimension(NDIM) :: ray
  integer :: i, v
#define SPEEDOPTIMIZED
#ifndef SPEEDOPTIMIZED
  integer j,k,m,n
  integer, dimension(2,2,NDIM) :: pos
  real, dimension(2,2) :: q
#endif
  integer, dimension(2,2), parameter :: jmask = RESHAPE((/ 0, 1, 0, 1 /), (/ 2, 2 /))
  integer, dimension(2,2), parameter :: kmask = RESHAPE((/ 0, 0, 1, 1 /), (/ 2, 2 /))
  real, dimension(2,2) :: xpos
  integer, dimension(NDIM) :: offset, vox
  real :: q11,q21,q12,q22
  INTERFACE
#ifndef CHECKS_RT
    pure &
#endif
    function interpolate_bilinear(q11,q21,q12,q22,           &
                                  x1,x2,y1,y2,               &
                                  x,y)
    real :: interpolate_bilinear
    real, intent(in) :: q11,q21,q12,q22
    real, intent(in) :: x1,x2,y1,y2
    real, intent(in) :: x,y
    end function
  END INTERFACE


  ray = ray_

!  if(face) then
!  print *,"dir: ", dir
!  print *,"ray: ", ray
!  print *,"voxel: ", voxel
!  end if

  offset = MAX(0,dir)
!  offset = ABS(MIN(0,dir))

!  SELECT CASE(SUM(ABS(dir)))
!  CASE(1,2,3)
    ! regular step accross a cell face
    !leaveI = voxel !+ offset!+ KD*NGUARD

    do i=1,NDIM
      if(dir(i).NE.0) then
#ifdef VERBOSE_RT
        write(6,*) 'interpolation using x-face'
#endif

#ifdef SPEEDOPTIMIZED
! This is a speed optimized version, where all constants are written out.
        select case(i)
        case(1)
          xpos(:,1) = REAL(voxel(2)+jmask(:,1))
          xpos(:,2) = REAL(voxel(3)+kmask(1,:))
          vox(1) = voxel(1) + NGUARD - 1 + offset(1)
          vox(2) = voxel(2) + NGUARD - 1
          vox(3) = voxel(3) + NGUARD - 1
          do v = 1,SIZE(vars)
            q11 = solnData(vars(v), vox(1), vox(2)    , vox(3))
            q21 = solnData(vars(v), vox(1), vox(2) + 1, vox(3))
            q12 = solnData(vars(v), vox(1), vox(2)    , vox(3) + 1)
            q22 = solnData(vars(v), vox(1), vox(2) + 1, vox(3) + 1)
            res(v) = interpolate_bilinear(&
              q11,q21,q12,q22,&
              xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
              ray(2),ray(3))
          end do
        case(2)
          xpos(:,1) = REAL(voxel(3)+jmask(:,1))
          xpos(:,2) = REAL(voxel(1)+kmask(1,:))
          vox(1) = voxel(1) + NGUARD - 1
          vox(2) = voxel(2) + NGUARD - 1 + offset(2)
          vox(3) = voxel(3) + NGUARD - 1
          do v = 1,SIZE(vars)
            q11 = solnData(vars(v), vox(1)    , vox(2), vox(3))
            q21 = solnData(vars(v), vox(1)    , vox(2), vox(3) + 1)
            q12 = solnData(vars(v), vox(1) + 1, vox(2), vox(3))
            q22 = solnData(vars(v), vox(1) + 1, vox(2), vox(3) + 1)
            res(v) = interpolate_bilinear(&
              q11,q21,q12,q22,&
              xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
              ray(3),ray(1))
          end do
        case(3)
          xpos(:,1) = REAL(voxel(1)+jmask(:,1))
          xpos(:,2) = REAL(voxel(2)+kmask(1,:))
          vox(1) = voxel(1) + NGUARD - 1
          vox(2) = voxel(2) + NGUARD - 1
          vox(3) = voxel(3) + NGUARD - 1 + offset(3)
          do v = 1,SIZE(vars)
            q11 = solnData(vars(v), vox(1)    , vox(2)    , vox(3))
            q21 = solnData(vars(v), vox(1) + 1, vox(2)    , vox(3))
            q12 = solnData(vars(v), vox(1)    , vox(2) + 1, vox(3))
            q22 = solnData(vars(v), vox(1) + 1, vox(2) + 1, vox(3))
            res(v) = interpolate_bilinear(&
              q11,q21,q12,q22,&
              xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
              ray(1),ray(2))
          end do
        end select

#else
! This is the logical variant, where all indezies are hidden
! in the pos array. This leads to less efficient code than
! the written out variant above.
        j = MOD(i,3) + 1
        k = MOD(i+1,3) + 1

        pos(:,:,i) = voxel(i) + NGUARD - 1 + offset(i)
        pos(:,:,j) = voxel(j) + NGUARD - 1 + jmask
        pos(:,:,k) = voxel(k) + NGUARD - 1 + kmask

        xpos(:,1) = REAL(voxel(j) + jmask(:,1))
        xpos(:,2) = REAL(voxel(k) + kmask(1,:))

        ! correct for precision errors
        call obeyintervall(xpos(:,1),ray(j))
        call obeyintervall(xpos(:,2),ray(k))

#ifdef CHECKS_RT
        if(.not.(valid(xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),ray(j),ray(k)))) then
           write(*,*) 'ERROR, initial coordinates for interpolation'
           write(*,*) 'are invalid in step, aborting...  '
           write(*,*) 'i,j,k',i,j,k
           write(*,*) 'pos:  ', pos
           write(*,*) 'xpos: ', xpos
           write(*,*) 'dir',direction
           write(*,*) 'ray(j),ray(k)',ray(j),ray(k)
           call Driver_abortFlash('calc_local_block_contributions_3DRT: ' &
                                  //'invalid coordinates for bilin. interpolation')
        endif
#endif

        do v = 1,SIZE(vars)
          do n=1,2
            do m=1,2
              q(m,n) = solnData(vars(v), &
                                pos(m,n,1), &
                                pos(m,n,2), &
                                pos(m,n,3))
            end do
          end do
          res(v) = interpolate_bilinear(&
            q(1,1),q(2,1),q(1,2),q(2,2),&
            xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
            ray(j),ray(k))
! This is a little bit faster than the one directly above:
!          q11 = solnData(vars(v),pos(1,1,1),pos(1,1,2),pos(1,1,3))
!          q21 = solnData(vars(v),pos(2,1,1),pos(2,1,2),pos(2,1,3))
!          q12 = solnData(vars(v),pos(1,2,1),pos(1,2,2),pos(1,2,3))
!          q22 = solnData(vars(v),pos(2,2,1),pos(2,2,2),pos(2,2,3))
!          res(v) = interpolate_bilinear(&
!            q11,q21,q12,q22,&
!            xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
!            ray(j),ray(k))
        end do
#endif

        return
      endif

    end do

#ifdef CHECKS_RT
    if(ANY(isnan(res))) then
       write(6,*) '============================='
       write(6,*) 'NaN error after interpolation'
       write(6,*) 'vars: ', vars
       write(6,*) 'res: ', res
       write(6,*) '============================='
       call Driver_abortFlash('NaN error in subroutine calc_local_block_contributions_3DRT') 
    endif
#endif

!  CASE(4)
!    ! diagonal step, interpolation not necessary
!    leaveI = voxel+KD*NGUARD-1
!    do v=1,SIZE(vars)
!      res(v) = solnData(vars(v), leaveI(1), leaveI(2), leaveI(3))
!    end do
!#ifdef CHECKS_RT
!  CASE DEFAULT
!    print *, "dir: ", dir
!    call Driver_abortFlash('calc_local_block_contributions_3DRT: Invalid direction sum.')
!#endif
!  END SELECT
END SUBROUTINE interpolate

END MODULE calc_local_mod
!===============================================================================
