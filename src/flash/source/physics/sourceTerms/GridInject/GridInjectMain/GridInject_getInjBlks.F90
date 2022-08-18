!!****f* source/physics/sourceTerms/GridInject/GridInject_getInjBlks
!!
!! NAME
!!
!!  GridInject_getInjBlks
!!
!! SYNOPSIS
!!
!!  call GridInject_getInjBlks(
!!      real(IN)      :: xloc
!!      real(IN)      :: yloc
!!      real(IN)      :: zloc
!!      real(IN)      :: radius
!!      integer(OUT)  :: injBlks(MAXBLOCKS)
!!      integer(OUT)  :: injBlkNum
!!  )
!!
!! DESCRIPTION
!!
!!  Find all leaf blocks ON THIS PROC that overlap injection sphere
!!  specified by (xloc, yloc, zloc, radius).
!!
!!  Only works for 3-D.
!!
!!  WARNING - does not work with periodic boundary conditions
!!
!! ARGUMENTS
!!
!!  xloc      : x-coordinate of injection sphere center
!!  yloc      : y-coordinate of injection sphere center
!!  zloc      : z-coordinate of injection sphere center
!!  radius    : injection sphere radius
!!  injBlks   : array of overlapping blockIDs (stored in 1:injBlkNum)
!!  injBlkNum : number of overlapping blocks
!!
!!***
subroutine GridInject_getInjBlks (xloc, yloc, zloc, radius, injBlks, InjBlkNum)

  use Grid_interface, ONLY : Grid_getBlkCenterCoords, &
                             Grid_getBlkPhysicalSize, &
                             Grid_getListOfBlocks

  implicit none

#include "Flash.h"
#include "constants.h"

  real, intent(IN)      :: xloc, yloc, zloc, radius
  integer, intent(OUT)  :: injBlks(MAXBLOCKS)
  integer, intent(OUT)  :: injBlkNum

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: n, blockID
  real :: blkCtr(MDIM), blkSize(MDIM)
  real :: xcoll, ycoll, zcoll

  injBlkNum = 0
  injBlks = 0

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do n = 1, blockCount
    blockID = blockList(n)

    call Grid_getBlkCenterCoords(blockID,blkCtr)
    call Grid_getBlkPhysicalSize(blockID,blkSize)

    ! WARNING - does not work with periodic boundary conditions

    ! If all three directions of the block have been pierced by the injection
    ! sphere then this is an injection block.
    ! Logic is not exact; if inject location lies diagonally outside a
    ! block corner, the block can be marked overlapping even if sphere
    ! radius does NOT intersect block.
    !if (abs(blkCtr(1) - xloc) .le. (0.5*blkSize(1)+radius) .and. &
    !    abs(blkCtr(2) - yloc) .le. (0.5*blkSize(2)+radius) .and. &
    !    abs(blkCtr(3) - zloc) .le. (0.5*blkSize(3)+radius) &
    !) then

    ! exact collision detection for sphere and rectangular prism
    ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
    xcoll = max(blkCtr(1)-0.5*blkSize(1),min(xloc,blkCtr(1)+0.5*blkSize(1)))
    ycoll = max(blkCtr(2)-0.5*blkSize(2),min(yloc,blkCtr(2)+0.5*blkSize(2)))
    zcoll = max(blkCtr(3)-0.5*blkSize(3),min(zloc,blkCtr(3)+0.5*blkSize(3)))
    if ((xcoll-xloc)**2+(ycoll-yloc)**2+(zcoll-zloc)**2 < radius**2) then
      injBlkNum = injBlkNum + 1
      injBlks(injBlkNum) = blockID
    end if
  end do

  return
end subroutine GridInject_getInjBlks
