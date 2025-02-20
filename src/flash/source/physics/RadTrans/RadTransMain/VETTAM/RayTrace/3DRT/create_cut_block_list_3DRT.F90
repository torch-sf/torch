!===============================================================================
!
!  Subroutine: create_cut_block_list
!
!===============================================================================
!
! Create a list of blocks cut by the ray connecting source and destination.
! Uses the 'Fast Voxel Traversal Algorithm for Ray Tracing' by
! J. Amanatides and Andrew Woo, Proc. Eurographics, Aug. 1987
! (http://www.cs.yorku.ca/~amana/research/grid.pdf)
! See also Abel, Norman & Madau  ApJ 523, 66.
!
! NOTE: in all comments in this subroutine the term 'block' is used for
! the blocks at the highest refinement level, whether they exist or not,
! i.e. they are the ones stored in the blockMapping array.
!
!===============================================================================

#undef VERBOSE_RT
!#define VERBOSE_RT

MODULE create_cut_block_mod
#include "Flash.h"
#include "constants.h"

  use RadTrans_TreeCommData, ONLY : g_lrefine, g_bndbox

  use RadTrans_HybridCharModule

  use Driver_interface, ONLY : Driver_abortFlash
  use fvt, ONLY : FastVoxelTraversal, intersection, &
    fvt_init, fvt_check, fvt_valid, fvt_get, fvt_dstep, SIGN3
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY: gr_oneBlock, gr_globalDomain
  use raytrace_data

IMPLICIT NONE
  include 'mpif.h'
PRIVATE
  real                :: distance
  integer             :: lastProcID,lastBlockID, lastPatchID
  real                :: told
  real, dimension(NDIM) :: x0, x1, rescale, rescale_inv, src
  integer, dimension(NDIM) :: b0, b1
  real, dimension(NDIM) :: globalDomainMin, globalDomainMax, dirFact
  logical :: first
  real :: contributiond, contributioni
  integer :: iAindex, maxLevel, nrOfPe
  integer :: maxNrOfLeafBlocks
  logical :: sink
  integer, parameter :: cblen = 1000
  real, dimension(3,cblen) :: cb
  integer :: cutno
  real :: t, tend

PUBLIC :: create_cut_block_list_3DRT

CONTAINS

subroutine create_cut_block_list_3DRT(     &
     myPe, myB, &
     contributiond_, contributioni_,    &
     src_, dest,             &
     regFact,    &
     dirFact_,                &
     iAindex_, &
     maxNrOfLeafBlocks_, &
     nrOfPe_, maxLevel_, sink_)
  !
  !
  implicit none
!
!
  integer, intent(in)                                :: myPe, myB
  real,    intent(out)                               :: contributiond_
  real,    intent(inout)                             :: contributioni_
  real, dimension(NDIM), intent(in)                  :: dest, regFact, dirFact_, src_
  integer, intent(in)                                :: maxNrOfLeafBlocks_
  integer, intent(in)                                :: nrOfPe_
  integer, intent(in)                               :: iAindex_
  integer, intent(in)                                :: maxLevel_
  logical, intent(in)                                :: sink_

  integer, dimension(NDIM) :: voxel, dir_
  !real :: t, tend
  real :: tdelta, dtau
  integer, dimension(NDIM) :: oldvoxel, dx, dxbnd, blkbnd
  real, dimension(NDIM) :: tMax, pos, dirFact_inv
  integer :: i
  !logical :: res
  real, dimension(2) :: res


  dirFact = dirFact_
  globalDomainMin = gr_globalDomain(LOW,IAXIS:KAXIS)
  globalDomainMax = gr_globalDomain(HIGH,IAXIS:KAXIS)
  iAindex = iAindex_
  maxLevel = maxLevel_
  nrOfPe = nrOfPe_
  maxNrOfLeafBlocks = maxNrOfLeafBlocks_
  sink = sink_
  src = src_

  contributiond = 0.0
  contributioni = contributioni_
  cb = 0.
  cutno = 0
  told = 0.0
  first = .true.


! NOTE: rays are traced BACKWARDS from destination to the edges of the domain.

  ! traverse always trough complete blocks, e.g. the step sizes can vary from block
  ! to block. This is way more effecient, otherwise there are to many steps.
  dirFact = -dirFact

  distance  = SQRT(SUM(dirFact**2))

  if (distance.le.0.) then
    contributiond_ = 0.
    contributioni_ = 0.
    return
  end if

  dirFact = dirFact / distance

  WHERE(dirFact.ne.0)
    dirFact_inv = distance/dirFact
  ELSE WHERE
    dirfact_inv = 1.e+100
  END WHERE

  ! Initialization for the new (generic) FastVoxelTraversal

  rescale = 1./regFact
  rescale_inv = regFact
  ! translate the coordinates to begin at 0 and rescale them,
  ! so a cell as maxlevel has the size of unity (1,1,1)
  b0 = 0
  b1 = NINT((globalDomainMax - globalDomainMin) * regFact)
  x0 = (dest - globalDomainMin) * regFact

  ! create_cut_block_list is called with a ray (startpoint 'dest', direction 'dirFact'),
  ! but the Voxel Traversal algorithm needs a line segment. Therefore we calculate the
  ! point, where the ray (which destination is our starting point, so we follow it
 ! backwards) intersects with the domain boundaries.
  call obeybox(b0,b1,x0)
  if(.not.sink) then
    t = intersection(dirFact, x0, b0, b1)
    x1 = dirFact * t + x0
  else
    if(.not.check_incoords(real(b0), real(b1), (src-globalDomainMin)*regFact)) then
!      print *,"source not in this block"
!      print *,"src: ", (src-globalDomainMin)*regFact
      t = intersection(dirFact, x0, b0, b1)
      x1 = dirFact * t + x0
    else
!      print *,"source in this block"
     x1 = (src - globalDomainMin)*regFact
    end if
  end if
  WHERE(equal(x1,b0)) x1=b0
  WHERE(equal(x1,b1)) x1=b1
  WHERE(ABS(x1-NINT(x1)).lt.precision) x1 = NINT(x1)
#ifdef CHECKS_RT
  IF(.not.check_incoords(real(b0), real(b1), x0).OR.&
     .not.check_incoords(real(b0), real(b1), x1)) THEN
    print *,"x0: ", x0
    print *,"x1: ", x1
    print *,"dirFact: ", dirFact
    print *,"b0: ", b0
    print *,"b1: ", b1
    call Driver_abortFlash("ERROR: create_cut_block_list: x0 or x1 not in b0 -> b1 box.")
  END IF
#endif

!  print *, "-- Raytrace --"
!  print *, "sink: ", sink, src
!  print *, "from: ",x0
!  print *, "to:   ",x1
!  print *, "dirFact: ", dirFact
!  print *, "gmin: ",b0
!  print *, "gmax: ",b1


  !print *,"start ========================="
  voxel = INT(x0)
  tend = SQRT(SUM((x1-x0)**2))
  t=0.
  dx = INT(SIGN3(dirFact))
  pos = (x0+precision*dx)*rescale+globalDomainMin

  lastProcID = myPe
  lastBlockID = myB
  lastPatchID = myB + lastProcID * maxNrOfLeafBlocks
!    else

!  print *,"blkmin: ", (g_bndbox(1,:,lastBlockID, lastProcID)-globalDomainMin)*rescale_inv
!  print *,"blkmax: ", (g_bndbox(2,:,lastBlockID, lastProcID)-globalDomainMin)*rescale_inv

  ! dx=-1 -> 1
  ! dx= 1 -> 2
  ! dx= 0 -> 1?
  dxbnd = MAX(0,dx)+1
  do while(t.lt.tend)
!    print *,"loop ======>"
!    oldvoxel = voxel

    !print *, "levelFact: ", levelFact*8
!    print *, "levelCheck: ",  NINT((g_bndbox(2,:,lastBlockID,lastProcID) -g_bndbox(1,:,lastBlockID,lastProcID))* rescale_inv)
    do i=1,NDIM
      IF(dx(i).ne.0) then
        blkbnd(i) = NINT((g_bndbox(dxbnd(i),i,lastBlockID,lastProcID) - globalDomainMin(i)) * rescale_inv(i))
      ELSE
        blkbnd(i) = HUGE(blkbnd(i))
      END IF
    end do
!    print *,"next boundary at: ", blkbnd

    !WHERE(dirFact.ne.0.)
    !  tMax = (REAL(blkbnd)-(REAL(x0)+t*dirFact))/dirFact
    !ELSE WHERE
    !  tMax = HUGE(0.)
    !END WHERE
    tMax = (REAL(blkbnd)-(REAL(x0)+t*dirFact))*dirFact_inv
    !tdelta = intersection(dirFact, x0+t*dirFact, blkmin, blkmax)

!    print *,"tMax: ", tMax
    tdelta = MINVAL(tMax)
    tDelta = MIN(tDelta,tend-t)
    t = t + tdelta
    if(tend-t.lt.5.e-2) t = tend
    dir_ = 0
    !print *,"t, tend: ", t, tend
    cutno=cutno+1
    IF(t.lt.tend) THEN
    ! also diagonal
    DO i=1,NDIM
      IF(ABS(tdelta-tMax(i)) .LT. precision) THEN
        dir_(i) = dx(i)
        exit
!        voxel(i) = voxel(i) + dx(i)
!        voxel(i) = blkbnd(i) + MIN(dx(i),0)
!        blkbnd(i) = blkbnd(i) + dx(i)*
      END IF
    END DO
    pos = x0+t*dirFact
    WHERE(ABS(pos-NINT(pos)).lt.precision) pos = NINT(pos)

    oldvoxel = FLOOR(pos)-MAX(dir_,0)
    voxel = FLOOR(pos+precision*dirFact)

!    print *,"voxel: ", oldvoxel
!    print *,"dir: ", dir_
!    print *,"Nextvoxel: ", voxel
!    print *,"t: ", t
!    print *,"x: ", pos
    call onNextVoxel2(oldvoxel,voxel,dir_,t,&
      lastPatchID, lastBlockID, lastProcID, res)
    cb(1,cutno) = res(1)
    cb(2,cutno) = res(2)
    !else
    else
      cb(1:2,cutno) = 0.
    END IF
    if(cutno.gt.1) &
      cb(3,cutno-1) = tDelta
  end do
!  print *,"stop  ========================="

  IF(cutno.ge.cblen) THEN
    ! cutno is equal to cblen. This means cblen is too small and has to be increased.
    print *,"cutno: ", cutno
    print *,"cblen: ", cblen
    call Driver_abortFlash('create_cut_block_list_3DRT: Post-mortem: cutno>=cblen. Too many cuts, increase cblen.')
  END IF

  DO i=cutno,1,-1
    dtau = cb(1,i)*cb(3,i)*distance
    contributiond = contributiond + dtau
    if(.not.sink) &
      contributioni = contributioni * exp(-dtau) + cb(2,i) * (1. - exp(-dtau))
  END DO

  contributiond_ = contributiond
  contributioni_ = contributioni

end subroutine create_cut_block_list_3DRT

subroutine onNextVoxel2(voxel,NextVoxel,dir_,t,&
  lastPatchID, lastBlockID, lastProcID, res)
  implicit none
  integer, dimension(:), intent(in) :: voxel, dir_, NextVoxel
  real, intent(in) :: t
  integer, intent(inout) :: lastPatchID, lastBlockID, lastProcID
  real, dimension(2), intent(out) :: res
  integer, dimension(SIZE(voxel)) :: block_dir
  real, dimension(SIZE(voxel)) :: x, ray
#ifdef CHECKS_RT
  real, dimension(SIZE(voxel)) :: tmp
#endif
  integer :: BlockID, ProcID, PatchID, i
  integer, dimension(NDIM) :: idx
  INTEGER, DIMENSION(2), PARAMETER :: vars = (/ OPAC_VAR, SOUR_VAR /)
  !print *,"new voxel: ", voxel

  ! get actual coordinates from voxel index inside the voxel
  ! the exact position inside the cell does not matter (in this case)
  x = (NextVoxel+0.5) * rescale + globalDomainMin
  !print *,"NextX: ", x/rescale-globalDomainMin
  DO i=1,3
    IF(dir_(i).eq.0) THEN
      IF(x(i) .lt. g_bndbox(1,i,lastBlockID,lastProcID)) THEN
        x(i) = g_bndbox(1,i,lastBlockID,lastProcID)
      ELSE IF(x(i).gt.g_bndbox(2,i,lastBlockID,lastProcID)) THEN
        x(i) = g_bndbox(2,i,lastBlockID,lastProcID)
      END IF
    END IF
  END DO

#ifdef CHECKS_RT
  tmp = x
  WHERE(dir_.ne.0) tmp=g_bndbox(1,:,lastBlockID,lastProcID)
  IF(.NOT.check_inblock(lastBlockID,lastProcID,tmp)) THEN
    print *, "block_dir:     ", dir_
    print *, "dirfact:       ", dirFact
    print *, "x:             ", x
    print *, "bmin:          ", g_bndbox(1,:,lastBlockID,lastProcID)
    print *, "bmax:          ", g_bndbox(2,:,lastBlockID,lastProcID)
    print *, "voxel: ", voxel
    print *, "NextVoxel: ", NextVoxel
    print *, "x as int: ", x/rescale-globalDomainMin
    print *, "bmin as int: ", g_bndbox(1,:,lastBlockID,lastProcID)/rescale-globalDomainMin
    print *, "bmax as int: ", g_bndbox(2,:,lastBlockID,lastProcID)/rescale-globalDomainMin
    call Driver_abortFlash("onNextVoxel2 ERROR: create_cut_block_list: Blk does not include the point of interest.")
  END IF
#endif

  block_dir = dir_

  CALL getNeighbour2(lastBlockID,lastProcID,block_dir,x,BlockID,ProcID)
!#ifdef CHECKS_RT
!  CALL checkNeighbour(lastBlockID,lastProcID,BlockID,ProcID,block_dir)
!#endif
  PatchID = find_leafindex(nrOfPe, maxNrOfLeafBlocks, nrOfLeafBlocksList, leafListAll, BlockID, ProcID)
!  PatchID = PatchID + ProcID * maxNrOfLeafBlocks

  ray = dirFact*t + x0
  WHERE(ABS(ray-NINT(ray)).lt.precision) ray=NINT(ray)
#ifdef CHECKS_RT
  IF(.not.check_inblock(lastBlockID,lastProcID,ray*rescale+globalDomainMin)) THEN
    print *, "bndmin:      ", g_bndbox(1,1:3,lastBlockID,lastProcID)
    print *, "bndmax:      ", g_bndbox(2,1:3,lastBlockID,lastProcID)
    print *, "ray rescale: ", ray*rescale+globalDomainMin
    print *, "ray:         ", ray
    call Driver_abortFlash("onNextVoxel2 ERROR: create_cut_block_list: ray not in block.")
  END IF
#endif

  idx = get_blockcoord_int(NextVoxel,BlockID,ProcID)
  if(.not.sink) then
    call interpolate(idx, -block_dir, get_blockcoord_real(ray,BlockID,ProcID), &
      PatchID, ProcID, vars, res)
  else
    call interpolate(idx, -block_dir, get_blockcoord_real(ray,BlockID,ProcID), &
      PatchID, ProcID, vars(1:1), res(1:1))
  end if

  lastPatchID = PatchID
  lastBlockID = BlockID
  lastProcID = ProcID
end subroutine onNextVoxel2

#ifdef CHECKS_RT
  pure function check_inblock(blk, p, idx)
!-------------------------------------------------------------------------------
    implicit none
    logical :: check_inblock
    integer, intent(in) :: blk, p
    real,    intent(in) :: idx(NDIM)
!-------------------------------------------------------------------------------
    check_inblock = .NOT.(     ANY(g_bndbox(1,1:3,blk,p).GT.idx) &
                          .OR. ANY(g_bndbox(2,1:3,blk,p).LT.idx))
  end function check_inblock
#endif

  pure function get_blockcoord_int(voxel, blk, p)
    integer, dimension(NDIM) :: get_blockcoord_int
    integer, dimension(NDIM), intent(in) :: voxel
    integer, intent(in) :: blk, p
    real, dimension(NDIM) :: bndmin
    integer, dimension(NDIM) :: bndminI
    real :: levelFact

    bndmin = g_bndbox(1,:,blk,p)
    bndminI = NINT((bndmin - globalDomainMin) * rescale_inv)

    levelFact = 2.0**(REAL(g_lrefine(blk,p)-maxLevel))
    get_blockcoord_int = INT((voxel - bndminI)*levelfact)
  end function get_blockcoord_int

  pure function get_blockcoord_real(x, blk, p)
    real, dimension(NDIM) :: get_blockcoord_real
    real, dimension(NDIM), intent(in) :: x
    integer, intent(in) :: blk, p
    real, dimension(NDIM) :: bndmin
    real :: levelFact

    bndmin = g_bndbox(1,:,blk,p)
    bndmin = (bndmin - globalDomainMin) * rescale_inv

    levelFact = 2.0**(REAL(g_lrefine(blk,p)-maxLevel))
    get_blockcoord_real = (x - bndmin)*levelfact

!    print *,"x: ", x
!    print *,"get_blockcoord: ", get_blockcoord_real
  end function get_blockcoord_real

  pure function get_blockcoord_src(x, blk, p)
    real, dimension(NDIM) :: get_blockcoord_src
    real, dimension(NDIM), intent(in) :: x
    integer, intent(in) :: blk, p
    real, dimension(NDIM) :: bndmin
    real :: levelFact

    bndmin = g_bndbox(1,:,blk,p)
    !bndmin = (bndmin - globalDomainMin) * rescale_inv

    levelFact = 2.0**(REAL(g_lrefine(blk,p)-maxLevel))
    get_blockcoord_src = (x - bndmin)*levelfact

  end function get_blockcoord_src

  function left_block(blk, p, idx)
!-------------------------------------------------------------------------------
    implicit none
    integer, dimension(NDIM) :: left_block
    integer, intent(in) :: blk, p
    real,    intent(in) :: idx(NDIM)
!-------------------------------------------------------------------------------
    real, dimension(NDIM) :: bndmin, bndmax
    bndmin(:) = g_bndbox(1,:,blk,p)
    bndmax(:) = g_bndbox(2,:,blk,p)
    WHERE(idx.LT.bndmin)
      left_block = -1
    ELSE WHERE(idx.GT.bndmax)
      left_block = 1
    ELSE WHERE
      left_block = 0
    END WHERE
  end function left_block

#ifndef CHECKS_RT
  pure &
#endif
  function find_leafindex(nrOfPe, maxNrOfLeafBlocks, nrOfLeafBlocksList, leafListAll, blk, p)
!-------------------------------------------------------------------------------
    implicit none
    integer :: find_leafindex
    integer, intent(in)                                :: nrOfPe
    integer, intent(in), dimension(0:nrOfPe-1)         :: nrOfLeafBlocksList
    integer, intent(in)                                :: maxNrOfLeafBlocks
    integer, intent(in), dimension(MaxNrOfLeafBlocks,0:nrOfPe-1)   :: leafListAll
    integer, intent(in)                                :: blk, p
    integer :: i, l, r
!-------------------------------------------------------------------------------
    find_leafindex = -1

! linear search
!    do i = 1, nrOfLeafBlocksList(p)
!      if (blk .eq. leafListAll(i, p)) THEN
!        find_leafindex = i
!        return
!      end if
!    enddo

    ! binary search
    l = 1
    r = nrOfLeafBlocksList(p)
    do while (l.le.r)
      i = (l+r)/2
      if(blk.eq.leafListAll(i, p)) THEN
        find_leafindex = i
        return
      else if(blk .lt. leafListAll(i,p)) then
        r = i - 1
      else
        l = i + 1
      end if
    end do
#ifdef CHECKS_RT
    if(find_leafindex.eq.-1) then
      print *,"searching for blk: ", blk
      print *,"p: ", p
      print *,"leafListAll(:,p): ", leafListAll(:,p)
      call Driver_abortFlash("ERROR: create_cut_block_list:find_leafindex: Could not find leaf index.")
    end if
#endif

    return
  end function find_leafindex

!#ifdef CHECKS_RT
!  subroutine checkNeighbour(b,p,b_new,p_new,block_dir)
!!-------------------------------------------------------------------------------
!    implicit none
!    integer,intent(in) :: b,p,b_new,p_new
!    integer,dimension(NDIM),intent(in) :: block_dir
!    real,dimension(3) :: distance
!!-------------------------------------------------------------------------------
!
!    distance = ABS(g_coord(:,b_new,p_new)-g_coord(:,b,p)) &
!               -0.5*(g_size(:,b_new,p_new) +g_size(:,b,p))
!    IF(.NOT.ALL(distance.LT.precision*g_size(1,b,p))) THEN
!      write(*,*) "dir_block:        ", block_dir
!      write(*,*) "oldBlock:         ", b,p
!      write(*,*) "oldBlock center:  ", g_coord(:,b,p)
!      write(*,*) "oldBlock size:    ", g_size(:,b,p)
!      write(*,*) "oldBlock lrefine: ", g_lrefine(b,p)
!      write(*,*) "oldBlock bmin     ", g_bndbox(1,:,b,p)
!      write(*,*) "oldBlock bmax     ", g_bndbox(2,:,b,p)
!      write(*,*) "newBlock:         ", b_new,p_new
!      write(*,*) "newBlock center:  ", g_coord(:,b_new,p_new)
!      write(*,*) "newBlock size:    ", g_size(:,b_new,p_new)
!      write(*,*) "newBlock lrefine: ", g_lrefine(b_new,p_new)
!      write(*,*) "newBlock bmin     ", g_bndbox(1,:,b_new,p_new)
!      write(*,*) "newBlock bmax     ", g_bndbox(2,:,b_new,p_new)
!      write(*,*) "abs diff centers: ", ABS(g_coord(:,b_new,p_new)-g_coord(:,b,p))
!      write(*,*) "<= sum sizes/2 ?: ", 0.5*(g_size(:,b_new,p_new)+g_size(:,b,p))
!      write(*,*) "distances:        ", distance
!      call Driver_abortFlash("ERROR: create_cut_block_list: Blocks are not neighbours.")
!    END IF
!  end subroutine checkNeighbour
!#endif

  SUBROUTINE getNeighbour2(b,p,block_dir,x,b_new,p_new)
!-------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: b, p
    integer, intent(in), dimension(NDIM) :: block_dir
    real, dimension(3), intent(in) :: x
    integer, intent(out) :: b_new, p_new
    integer :: dir, n
    integer :: numNegh
    integer, dimension(BLKNO:PROCNO) :: neghBlkProc
    INTEGER, DIMENSION(MDIM) :: edge
    integer :: leafid, i, j, k
!-------------------------------------------------------------------------------

    edge(:) = (/ 1+K1D, 1+K2D, 1+K3D /) ! For 3D this initializes all directions to CENTER==2
    DO n=1,NDIM
      IF(block_dir(n).GT.0) THEN
        edge(n) = RIGHT_EDGE
        dir = 2*n
        EXIT
      ELSE IF(block_dir(n).LT.0) THEN
        edge(n) = LEFT_EDGE
        dir = 2*n-1
        EXIT
      END IF
    END DO
    leafid = rblocklistAll(b,p)+1
    numNegh = SurrBlkSumNumNegh(dir,leafid)
    IF(numNegh.eq.1) THEN ! Might be 0, if it is a physical boundary
      neghBlkProc(BLKNO:PROCNO) = SurrBlkSumArray(:,numNegh,dir,leafid)
    ELSE IF(numNegh.eq.0) THEN
      ! do nothing - physical boundary
      print *, "numNegh: ", numNegh
      call Driver_abortFlash("ERROR: Physical boundary should never be hit. One possible solution is to increase the &
      & precision in RadTrans_hybridCharModule.F90 so that we do not overshoot the ray &
      & to the domain boundary.")
    ELSE IF(numNegh.eq.4) THEN
      SELECT CASE(n)
      CASE(1)
        j=2
        k=3
      CASE(2)
        j=1
        k=3
      CASE(3)
        j=1
        k=2
      END SELECT
      i = 1
      IF(x(j).GT.(g_bndbox(1,j,b,p) + 0.5 * (g_bndbox(2,j,b,p)-g_bndbox(1,j,b,p)))) THEN
        i = 2
      END IF
      IF(x(k).GT.(g_bndbox(1,k,b,p) + 0.5 * (g_bndbox(2,k,b,p)-g_bndbox(1,k,b,p)))) THEN
        i = i + 2
      END IF
      neghBlkProc(BLKNO:PROCNO) = SurrBlkSumArray(:,i,dir,leafid)
      ! do nothing for now
    ELSE
      print *, "numNegh: ", numNegh
      call Driver_abortFlash("ERROR: numNegh number unexpected/unknown.")
    END IF
    b_new = neghBlkProc(BLKNO)
    p_new = neghBlkProc(PROCNO)
  END SUBROUTINE getNeighbour2

#ifndef CHECKS_RT
PURE &
#endif
SUBROUTINE interpolate(voxel, dir, ray_, BlockID, ProcID, vars, res)
IMPLICIT NONE
  integer, dimension(NDIM), intent(in) :: voxel, dir
  real, dimension(NDIM), intent(in) :: ray_
  integer, INTENT(IN) :: BlockID,ProcID
  integer, intent(IN), dimension(:) :: vars
  real, intent(OUT),dimension(SIZE(vars)) :: res
  real, dimension(2,2) :: q
!  integer, dimension(2,2,NDIM) :: pos
  real, dimension(NDIM) :: ray
  integer :: i,j,k,v,m,n,a,b
  integer, dimension(2,2), parameter :: jmask = RESHAPE((/ 0, 1, 0, 1 /), (/ 2, 2 /))
  integer, dimension(2,2), parameter :: kmask = RESHAPE((/ 0, 0, 1, 1 /), (/ 2, 2 /))
  real, dimension(2,2) :: xpos
  integer, dimension(NDIM) :: offset
  integer :: facecode
#undef SPEEDOPTIMIZED
!#define SPEEDOPTIMIZED
#ifdef SPEEDOPTIMIZED
  real :: q11,q12,q21,q22
#endif
  integer, dimension(2) :: vox
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

  offset = MAX(0,dir)
!  offset = ABS(MIN(0,dir))

  facecode = getfacecode(dir)

  SELECT CASE(SUM(ABS(dir)))
  CASE(1,2,3)
    ! regular step accross a cell face

    do i=1,NDIM
      if(dir(i).NE.0) then

#ifdef SPEEDOPTIMIZED
! This is a speed optimized version, where all constants are written out.
        SELECT CASE(i)
        CASE(1)
          xpos(:,1) = REAL(voxel(2) + jmask(:,1))
          xpos(:,2) = REAL(voxel(3) + kmask(1,:))
          vox(1) = voxel(2) + NGUARD
          vox(2) = voxel(3) + NGUARD
          IF(sink) THEN
            v = offset(1)+1
            q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
            q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              res(1) = interpolate_bilinear(&
                q11,q21,q12,q22,&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(2),ray(3))
          ELSE
            do v = 1,SIZE(vars)
              q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              res(v) = interpolate_bilinear(&
                q11,q21,q12,q22,&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(2),ray(3))
            end do
          END IF
        CASE(2)
          xpos(:,1) = REAL(voxel(1) + jmask(:,1))
          xpos(:,2) = REAL(voxel(3) + kmask(1,:))
          vox(1) = voxel(1) + NGUARD
          vox(2) = voxel(3) + NGUARD
          IF(sink) THEN
            v = offset(2)+1
            q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
            q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              res(1) = interpolate_bilinear(&
                q11,q21,q12,q22,&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(1),ray(3))
          ELSE
            do v = 1,SIZE(vars)
              q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              res(v) = interpolate_bilinear(&
                q11,q21,q12,q22,&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(1),ray(3))
            end do
          END IF
        CASE(3)
          xpos(:,1) = REAL(voxel(1) + jmask(:,1))
          xpos(:,2) = REAL(voxel(2) + kmask(1,:))
          vox(1) = voxel(1) + NGUARD
          vox(2) = voxel(2) + NGUARD
          IF(sink) THEN
            v = offset(3)+1
            q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
            q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
            q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
            res(1) = interpolate_bilinear(&
              q11,q21,q12,q22,&
              xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
              ray(1),ray(2))
          ELSE
            do v = 1,SIZE(vars)
              q11 = faceValueAll(vox(1)    , vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q21 = faceValueAll(vox(1) + 1, vox(2)    , facecode,v,BlockID,iAindex,ProcID)
              q12 = faceValueAll(vox(1)    , vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              q22 = faceValueAll(vox(1) + 1, vox(2) + 1, facecode,v,BlockID,iAindex,ProcID)
              res(v) = interpolate_bilinear(&
                q11,q21,q12,q22,&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(1),ray(2))
            end do
          END IF
        END SELECT

#else
! This is the logical variant, where all indezies are hidden
! in the pos array. This leads to less efficient code than
! the written out variant above.
        j = MOD(i,3) + 1
        k = MOD(i+1,3) + 1

        !pos(:,:,i) = voxel(i) + offset(i)
        !pos(:,:,j) = voxel(j) + NGUARD + jmask
        !pos(:,:,k) = voxel(k) + NGUARD + kmask
        vox(1) = voxel(j) + NGUARD
        vox(2) = voxel(k) + NGUARD

        ! This can happen, if we hit a edge of the block.
        ! The interpolation is than happening effeciently only between
        ! two corners, so it makes no difference if we take the last
        ! cell face.
        !IF(ANY(pos(:,:,j).gt.8)) pos(:,:,j) = pos(:,:,j) - 1
        !IF(ANY(pos(:,:,k).gt.8)) pos(:,:,k) = pos(:,:,k) - 1
        !IF(ANY(pos(:,:,j).lt.0)) pos(:,:,j) = pos(:,:,j) + 1
        !IF(ANY(pos(:,:,k).lt.0)) pos(:,:,k) = pos(:,:,k) + 1
        WHERE(vox.lt.NGUARD) vox=vox+1

        !WHERE(pos.gt.8) pos = pos - 1
        !WHERE(pos.lt.0) pos = pos + 1

#ifdef CHECKS_RT
        IF(vox(1).gt.8+NGUARD.or.ANY(vox(1)+jmask(:,1).lt.0+NGUARD).OR. &
           vox(2).gt.8+NGUARD.or.ANY(vox(2)+kmask(:,2).lt.0+NGUARD)) THEN
           print *, "vox", vox(:) 
           call Driver_abortFlash('create_cut_block_list_3DRT: ' &
                                  //'invalid corner index for bilin. interpolation')
        END IF
#endif

        xpos(:,1) = REAL(voxel(j) + jmask(:,1))
        xpos(:,2) = REAL(voxel(k) + kmask(1,:))

        ! correct for precision errors
        !call obeyintervall(xpos(:,1),ray(j))
        !call obeyintervall(xpos(:,2),ray(k))

#ifdef CHECKS_RT
        if(.not.(valid(xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),ray(j),ray(k)))) then
           write(*,*) 'ERROR, initial coordinates for interpolation'
           write(*,*) 'are invalid in step, aborting...  '
           write(*,*) 'i,j,k',i,j,k
           write(*,*) 'vox:  ', vox
           write(*,*) 'xpos: ', xpos
           write(*,*) 'ray(j),ray(k)',ray(j),ray(k)
           call Driver_abortFlash('create_cut_block_list_3DRT: ' &
                                  //'invalid coordinates for bilin. interpolation')
        endif
#endif
          IF(sink) THEN
            v = offset(i)+1
            do n=1,2
              do m=1,2
                a = MAX(MIN(vox(1)+jmask(m,n),NXB+NGUARD), NGUARD)
                b = MAX(MIN(vox(2)+kmask(m,n),NXB+NGUARD), NGUARD)
                q(m,n) = faceValueAll(a,b,facecode,v,BlockID,iAindex,ProcID)
              end do
            end do
            res(1) = interpolate_bilinear(&
              q(1,1),q(2,1),q(1,2),q(2,2),&
              xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
              ray(j),ray(k))
          ELSE
            do v = 1,SIZE(vars)
              do n=1,2
                do m=1,2
                  a = MAX(MIN(vox(1)+jmask(m,n),NXB+NGUARD), NGUARD)
                  b = MAX(MIN(vox(2)+kmask(m,n),NXB+NGUARD), NGUARD)
                  q(m,n) = faceValueAll(a,b,facecode,v,BlockID,iAindex,ProcID)
                end do
              end do
              res(v) = interpolate_bilinear(&
                q(1,1),q(2,1),q(1,2),q(2,2),&
                xpos(1,1),xpos(2,1),xpos(1,2),xpos(2,2),&
                ray(j),ray(k))
            end do
          END IF
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
!  CASE(3)
!    ! diagonal step, interpolation not necessary
!    if(sink) then
!      vox(1:2) = voxel(1:2) + offset(1:2) + NGUARD
!      v = offset(1)+1
!      res(1) = faceValueAll(vox(1),vox(2),facecode,v,BlockID,iAindex,ProcID)
!    else
!    vox(1:2) = voxel(1:2) + offset(1:2) + NGUARD
!    do v=1,SIZE(vars)
!      res(v) = faceValueAll(vox(1),vox(2),facecode,v,BlockID,iAindex,ProcID)
!    end do
!    end if
#ifdef CHECKS_RT
  CASE DEFAULT
    call Driver_abortFlash('create_cut_block_list_3DRT: Invalid direction sum.')
#endif
  END SELECT
END SUBROUTINE interpolate

#ifndef CHECKS_RT
pure &
#endif
integer function getfacecode(block_dir)
  implicit none
  integer, dimension(NDIM), intent(in) :: block_dir
  integer :: n
  do n=NDIM,1,-1
    IF(block_dir(n).NE.0) THEN
      getfacecode = n
    END IF
  end do
#ifdef CHECKS_RT
  if(getfacecode.LT.1.or.getfacecode.gt.3) then
    print *,"block_dir: ", block_dir
    print *,"getfacecode: ", getfacecode
    call Driver_abortFlash("ERROR: create_cut_block_list: Wrong block_dir in getfacecode().")
  end if
#endif
end function getfacecode

END MODULE create_cut_block_mod
