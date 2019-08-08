!!****if* source/Particles/localAPI/pt_mapFromMeshQuadratic
!!
!! NAME
!!
!!  pt_mapFromMeshQuadratic
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!  
!!
!!***

subroutine pt_mapFromMeshInZone (numAttrib, attrib, pos, bndBox,&
     deltaCell,solnVec, partAttribVec)
  

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: numAttrib
  integer, dimension(2, numAttrib),intent(IN) :: attrib
  real,dimension(MDIM), INTENT(in) :: pos,deltaCell
  real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox
  real, pointer :: solnVec(:,:,:,:)
  real,dimension(numAttrib), intent(OUT) :: partAttribVec

  partAttribVec = 0.0
  return
  
end subroutine pt_mapFromMeshInZone

!===============================================================================

