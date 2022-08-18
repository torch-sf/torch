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

  implicit none

#include "Flash.h"

  real, intent(IN)      :: xloc, yloc, zloc, radius
  integer, intent(OUT)  :: injBlks(MAXBLOCKS)
  integer, intent(OUT)  :: injBlkNum

  return
end subroutine GridInject_getInjBlks
