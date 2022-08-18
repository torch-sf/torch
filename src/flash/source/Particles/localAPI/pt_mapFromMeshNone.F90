!!****if* source/Particles/localAPI/pt_mapFromMeshQuadratic
!!
!! NAME
!!
!!  pt_mapFromMeshQuadratic
!!
!! SYNOPSIS
!!
!!  pt_mapFromMeshQuadratic(integer, INTENT(in)    :: numAttrib,
!!                        integer, INTENT(in)    :: attrib(2,numAttrib),
!!                           real, INTENT(in)    :: pos(MDIM),
!!                           real, INTENT(in)    :: bndBox(2,MDIM),
!!                           real, INTENT(in)    :: deltaCell(MDIM),
!!                           real, pointer       :: solnVec(:,:,:,:),
!!                           real, INTENT(OUT)   :: partAttribVec(numAttrib))
!!
!! DESCRIPTION
!!  
!! returns given values, does not map grid gas varialbles to velocities
!!
!!***

subroutine pt_mapFromMeshNone (numAttrib, attrib, pos, bndBox,&
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
  
end subroutine pt_mapFromMeshNone

!===============================================================================

