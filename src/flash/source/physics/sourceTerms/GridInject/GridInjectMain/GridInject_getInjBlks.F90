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
!!  Find all leaf blocks ON THIS PROC that overlap injection region.
!!  Injection region is a cube of width 2*radius centered on (xloc,yloc,zloc)
!!  for fast collision detection.
!!
!!  Only works for 3-D.
!!
!! ARGUMENTS
!!
!!  xloc      : x-coordinate of injection cube center
!!  yloc      : y-coordinate of injection cube center
!!  zloc      : z-coordinate of injection cube center
!!  radius    : half-width of injection cube
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
  real :: blockCenter(MDIM), blockSize(MDIM)
  real :: loc(3)

  loc = [xloc, yloc, zloc]

  injBlkNum = 0
  injBlks = 0

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do n = 1, blockCount
    blockID = blockList(n)

    call Grid_getBlkCenterCoords(blockID,blockCenter)
    call Grid_getBlkPhysicalSize(blockID,blockSize)

    ! If all three directions of the block have been pierced by the injection
    ! sphere then this is an injection block.
    !
    ! Logic is not exact; if inject location lies diagonally outside a
    ! block corner, the block can be marked overlapping even if sphere
    ! radius does NOT intersect block.
    if (abs(blockCenter(1) - xloc) .le. (0.5*blockSize(1)+radius) .and. &
        abs(blockCenter(2) - yloc) .le. (0.5*blockSize(2)+radius) .and. &
        abs(blockCenter(3) - zloc) .le. (0.5*blockSize(3)+radius) &
    ) then
      injBlkNum = injBlkNum + 1
      injBlks(injBlkNum) = blockID
    end if
  end do

  return
end subroutine GridInject_getInjBlks
