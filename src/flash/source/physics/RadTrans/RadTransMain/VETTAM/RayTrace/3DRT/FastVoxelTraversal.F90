! Manuel Jung

MODULE fvt
  use RadTrans_hybridCharModule, ONLY: equal
PRIVATE
  INTEGER, PARAMETER :: NDIM = 3
  REAL, PARAMETER :: EPS=1.e-6
  REAL, DIMENSION(NDIM) :: tDelta,tMax,n
  INTEGER, DIMENSION(NDIM) :: voxel,dx,dir,oldvoxel,dims
  REAL, DIMENSION(NDIM) :: x1
  REAL :: t, distance, told
  REAL, PARAMETER :: gridmax = HUGE(1.e+10)

INTERFACE inGrid
  MODULE PROCEDURE inGridReal, inGridInt
END INTERFACE

INTERFACE intersection
  MODULE PROCEDURE intersection_real, intersection_int
END INTERFACE

PUBLIC :: FastVoxelTraversal, FastVoxelTraversalList, intersection, FastVoxelTraversalLT, &
  fvt_init, fvt_dstep, fvt_step, fvt_check, fvt_valid, fvt_get, fvt_checkLT, SIGN3, distance

CONTAINS

ELEMENTAL FUNCTION FRAC(x) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x
  REAL :: res
  res = (x-FLOOR(x))
END FUNCTION FRAC

#ifndef CHECKS_RT
PURE &
#endif
FUNCTION intersection_int(n,x0,b0,b1) RESULT(t)
  IMPLICIT NONE
  REAL,DIMENSION(:),INTENT(IN) :: n,x0
  INTEGER,DIMENSION(:),INTENT(IN) :: b0,b1
  REAL :: t
  t = intersection_real(n,x0,REAL(b0),REAL(b1))
END FUNCTION intersection_int

#ifndef CHECKS_RT
PURE &
#endif
FUNCTION intersection_real(n,x0,b0,b1) RESULT(t)
  IMPLICIT NONE
  REAL,DIMENSION(:),INTENT(IN) :: n,x0
  REAL,DIMENSION(:),INTENT(IN) :: b0,b1
  REAL :: t
  REAL :: tmin,tmax,t1,t2
  INTEGER :: i
  tmin = -HUGE(tmin)
  tmax =  HUGE(tmax)
  DO i=1,SIZE(n)
    IF(n(i).NE.0.) THEN
      t1 = (b0(i)-x0(i))/n(i)
      t2 = (b1(i)-x0(i))/n(i)

      tmin = MAX(tmin, MIN(t1, t2))
      tmax = MIN(tmax, MAX(t1, t2))
#ifdef CHECKS_RT
    ELSE IF(x0(i) .LT. b0(i) .OR. x0(i) .GT. b1(i)) THEN
      print *,"n:  ", n
      print *,"x0: ", x0
      print *,"b0: ", b0
      print *,"b1: ", b1
      print *,"ray does not intersect with grid/box."
#endif
    END IF
  END DO

  !t=tmin, unless tmin<0, in which case you're inside the box and t=tmax.
  IF(tmin.LT.0) THEN
    t = tmax
  ELSE
    t = tmin
  END IF
END FUNCTION intersection_real

ELEMENTAL FUNCTION SIGN3(x) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x
  REAL :: res
  IF(x.NE.0.) THEN
    res = SIGN(1.,x)
  ELSE
    res = 0.
  END IF
END FUNCTION

FUNCTION inGridReal(x,b0,b1) RESULT(res)
  IMPLICIT NONE
  REAL,DIMENSION(:),INTENT(IN) :: x
  INTEGER,DIMENSION(:),INTENT(IN) :: b0,b1
  LOGICAL :: res
  res = ALL(REAL(b0) .LE. x) .AND. ALL(x .LE. REAL(b1))
END FUNCTION

FUNCTION inGridInt(x,dims) RESULT(res)
  IMPLICIT NONE
  INTEGER,DIMENSION(:),INTENT(IN) :: x,dims
  LOGICAL :: res
  res = ALL(x .LT. dims)
END FUNCTION

SUBROUTINE fvt_init(x1p,x2,b0,b1)
  IMPLICIT NONE
  REAL, DIMENSION(:),INTENT(IN) :: x1p,x2
  INTEGER, DIMENSION(:),INTENT(IN) :: b0,b1
  REAL, DIMENSION(SIZE(x1p)) :: diff

  INTEGER :: i

  dims = b1-b0

  diff = x2 - x1p
  distance = SQRT(SUM(diff*diff))
  n = diff/distance
  dx = INT(SIGN3(diff))

  ! Assume that x1p and x2 always are inside the box marked
  ! by the b0 and b1 corners.
!  IF(.NOT.inGrid(x1p,b0,b1)) THEN
!    t = intersection(n, x1p, b0, b1)
!    x1 = Prec(n * t + x1p)
!  ELSE
    t = 0.
    x1 = x1p
   told = 0.
!  END IF

  DO i=1,SIZE(x1)
    IF(dx(i).NE.0) THEN
      tDelta(i) = MIN(dx(i)/n(i), gridmax)
    ELSE
      tDelta(i) = gridmax
    END IF
    IF(dx(i).GT.0) THEN
      tMax(i) = tDelta(i) * (1.-FRAC(x1(i)))
    ELSE IF(dx(i).LT.0) THEN
      tMax(i) = tDelta(i) * FRAC(x1(i))
    ELSE
      tMax(i) = gridmax
    END IF
    voxel(i) = FLOOR(x1(i))
  END DO
  ! Constant offset so that t==distance if x2 is reached.
  tMax = tMax + t
  dir = 0

END SUBROUTINE fvt_init

! check if we should loop again (using .lt.)
! also break, if we are nearly at distance.
! Otherwise floating point precision errors lead
! sometimes to very short steps of an illegal voxel.
PURE FUNCTION fvt_checkLT() RESULT(res)
  IMPLICIT NONE
  LOGICAL :: res
  res = t .LE. distance*(0.999999) !.and. .not. equal(t,distance)
END FUNCTION fvt_checkLT

! check if we should loop again (using .le.)
PURE FUNCTION fvt_check() RESULT(res)
  IMPLICIT NONE
  LOGICAL :: res
  res = t .LE. distance
END FUNCTION fvt_check

! take a single normal step
SUBROUTINE fvt_step()
  IMPLICIT NONE
  INTEGER :: i
  oldvoxel = voxel

  CALL MINVALLOC(tMax,t,i)
  dir = 0
  dir(i) = dx(i)
  voxel(i) = voxel(i) + dx(i)
  tMax(i) = tMax(i) + tDelta(i)
END SUBROUTINE fvt_step

! take a single step and allow
! diagonal steps
SUBROUTINE fvt_dstep()
  IMPLICIT NONE
  INTEGER :: i

  oldvoxel = voxel

  t = MINVAL(tMax)
  dir = 0

  ! also diagonal
  DO i=1,NDIM
    IF(ABS(t-tMax(i)) .LT. EPS) THEN
      dir(i) = dx(i)
      voxel(i) = voxel(i) + dx(i)
      tMax(i) = tMax(i) + tDelta(i)
    END IF
  END DO
END SUBROUTINE fvt_dstep

SUBROUTINE fvt_get(v,d,localt)
  IMPLICIT NONE
  INTEGER, DIMENSION(NDIM),INTENT(OUT) :: v, d
  REAL, INTENT(OUT) :: localt
  v = oldvoxel
  d = dir
  localt = t
END SUBROUTINE

! check if was a valid step
! should not be needed and is bad
! for performance
FUNCTION fvt_valid() RESULT(res)
  IMPLICIT NONE
  LOGICAL :: res
  res = inGrid(oldvoxel, dims)
END FUNCTION fvt_valid

SUBROUTINE FastVoxelTraversalList(x1,x2,b0,b1,allowDiagonal,list,n)
  IMPLICIT NONE
  REAL, DIMENSION(:),INTENT(IN) :: x1,x2
  INTEGER, DIMENSION(:),INTENT(IN) :: b0,b1
  LOGICAL, INTENT(IN) :: allowDiagonal
  INTEGER, DIMENSION(:,:), INTENT(INOUT) :: list
  INTEGER, INTENT(INOUT) :: n
  INTEGER :: maxlen
  maxlen = SIZE(list,2)
  n = 0
  CALL FastVoxelTraversal(x1,x2,b0,b1,allowDiagonal,fn)

  CONTAINS
  SUBROUTINE fn(x,dir,t)
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: x, dir
    REAL, INTENT(IN) :: t
    IF(n.LT.maxlen) THEN
      n = n + 1
      list(:,n) = x
    END IF
  END SUBROUTINE fn
END SUBROUTINE FastVoxelTraversalList


! Cutoff e.g. -1.e-7 to 0.
ELEMENTAL REAL FUNCTION Prec(x)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x
  IF(ABS(x-NINT(x)).LT.EPS) THEN
    Prec = NINT(x)
  ELSE
    Prec = x
  END IF
END FUNCTION Prec

SUBROUTINE FastVoxelTraversal(x1p,x2,b0,b1,allowDiagonal,fn)
  IMPLICIT NONE
  REAL, DIMENSION(:),INTENT(IN) :: x1p,x2
  INTEGER, DIMENSION(:),INTENT(IN) :: b0,b1
  LOGICAL, INTENT(IN) :: allowDiagonal
  INTERFACE
    SUBROUTINE fn(x,dir,t)
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(IN) :: x,dir
      REAL, INTENT(IN) :: t
    END SUBROUTINE fn
  END INTERFACE
  REAL, DIMENSION(SIZE(x1p)) :: tDelta,tMax,x1
  REAL, DIMENSION(SIZE(x1p)) :: n,diff
  INTEGER, DIMENSION(SIZE(x1p)) :: voxel,dx,dims,dir,oldvoxel
  REAL, PARAMETER :: gridmax = HUGE(1.e+10)
  REAL :: t, distance
  INTEGER :: i

  dims = b1-b0

  diff = x2 - x1p
  distance = SQRT(SUM(diff*diff))
  n = diff/distance
  dx = INT(SIGN3(diff))
  IF(.NOT.inGrid(x1p,b0,b1)) THEN
    t = intersection(n, x1p, b0, b1)
    x1 = Prec(n * t + x1p)
  ELSE
    t = 0.
    x1 = x1p
  END IF

  DO i=1,SIZE(x1)
    IF(dx(i).NE.0) THEN
      tDelta(i) = MIN(dx(i)/n(i), gridmax)
    ELSE
      tDelta(i) = gridmax
    END IF
    IF(dx(i).GT.0) THEN
      tMax(i) = tDelta(i) * (1.-FRAC(x1(i)))
    ELSE IF(dx(i).LT.0) THEN
      tMax(i) = tDelta(i) * FRAC(x1(i))
    ELSE
      tMax(i) = gridmax
    END IF
    voxel(i) = FLOOR(x1(i))
  END DO
  ! Constant offset so that t==distance if x2 is reached.
  tMax = tMax + t
  dir = 0
  !print *,voxel,x1

  DO WHILE (inGrid(REAL(voxel),b0,b1).AND.t.LE.distance)
    oldvoxel = voxel

    t = MINVAL(tMax)
    dir = 0
    IF(.NOT.allowDiagonal) THEN
      ! no diagonal
      i = MINLOC(tMax,1)
      dir(i) = dx(i)
      voxel(i) = voxel(i) + dx(i)
      tMax(i) = tMax(i) + tDelta(i)
    ELSE
      ! also diagonal
      DO i=1,SIZE(x1p)
        !print *, i, ABS(t-tMax(i)), EPS
        IF(ABS(t-tMax(i)) .LT. EPS) THEN
          dir(i) = dx(i)
          voxel(i) = voxel(i) + dx(i)
          tMax(i) = tMax(i) + tDelta(i)
        END IF
      END DO
    END IF
    ! Check first if integer voxel coordinates are valid, since
    ! for negative directions they start one cell too early
    IF(inGrid(oldvoxel,dims)) &
      CALL fn(oldvoxel,dir,t)
  END DO
END SUBROUTINE FastVoxelTraversal

SUBROUTINE FastVoxelTraversalLT(x1p,x2,b0,b1,allowDiagonal,fn)
  IMPLICIT NONE
  REAL, DIMENSION(:),INTENT(IN) :: x1p,x2
  INTEGER, DIMENSION(:),INTENT(IN) :: b0,b1
  LOGICAL, INTENT(IN) :: allowDiagonal
  INTERFACE
    SUBROUTINE fn(x,dir,t)
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(IN) :: x,dir
      REAL, INTENT(IN) :: t
    END SUBROUTINE fn
  END INTERFACE
  REAL, DIMENSION(SIZE(x1p)) :: tDelta,tMax,x1
  REAL, DIMENSION(SIZE(x1p)) :: n,diff
  INTEGER, DIMENSION(SIZE(x1p)) :: voxel,dx,dims,dir,oldvoxel
  REAL, PARAMETER :: gridmax = HUGE(1.e+10)
  REAL :: t, distance
  INTEGER :: i

  dims = b1-b0

  diff = x2 - x1p
  distance = SQRT(SUM(diff*diff))
  n = diff/distance
  dx = INT(SIGN3(diff))
  IF(.NOT.inGrid(x1p,b0,b1)) THEN
    t = intersection(n, x1p, b0, b1)
    x1 = Prec(n * t + x1p)
  ELSE
    t = 0.
    x1 = x1p
  END IF

  DO i=1,SIZE(x1)
    IF(dx(i).NE.0) THEN
      tDelta(i) = MIN(dx(i)/n(i), gridmax)
    ELSE
      tDelta(i) = gridmax
    END IF
    IF(dx(i).GT.0) THEN
      tMax(i) = tDelta(i) * (1.-FRAC(x1(i)))
    ELSE IF(dx(i).LT.0) THEN
      tMax(i) = tDelta(i) * FRAC(x1(i))
    ELSE
      tMax(i) = gridmax
    END IF
    voxel(i) = FLOOR(x1(i))
  END DO
  ! Constant offset so that t==distance if x2 is reached.
  tMax = tMax + t
  dir = 0
  !print *,voxel,x1

  DO WHILE (inGrid(REAL(voxel),b0,b1).AND.t.LT.distance)
    oldvoxel = voxel

    t = MINVAL(tMax)
    dir = 0
    IF(.NOT.allowDiagonal) THEN
      ! no diagonal
      i = MINLOC(tMax,1)
      dir(i) = dx(i)
      voxel(i) = voxel(i) + dx(i)
      tMax(i) = tMax(i) + tDelta(i)
    ELSE
      ! also diagonal
      DO i=1,SIZE(x1p)
        !print *, i, ABS(t-tMax(i)), EPS
        IF(ABS(t-tMax(i)) .LT. EPS) THEN
          dir(i) = dx(i)
          voxel(i) = voxel(i) + dx(i)
          tMax(i) = tMax(i) + tDelta(i)
        END IF
      END DO
    END IF
    ! Check first if integer voxel coordinates are valid, since
    ! for negative directions they start one cell too early
    IF(inGrid(oldvoxel,dims)) &
      CALL fn(oldvoxel,dir,t)
  END DO
END SUBROUTINE FastVoxelTraversalLT

! Mix of MINVAL and MINLOC, so only one
! search has to be done
PURE SUBROUTINE MINVALLOC(arr, val, i)
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: arr
  REAL, INTENT(OUT) :: val
  INTEGER, INTENT(OUT) :: i
  i = MINLOC(arr,1)
  val = arr(i)
END SUBROUTINE MINVALLOC

END MODULE fvt
