!!****f* source/physics/RadTrans/RadTransMain/VETTAM/RadTrans_finalize
!!
!! NAME
!!
!!  RadTrans_finalize
!!
!!
!! SYNOPSIS
!!
!!  RadTrans_finalize()
!!  
!!
!! DESCRIPTION
!! 
!!  Finalize unit scope variables which might mean deallocation, etc.
!!
!!***


#include "petsc/finclude/petscksp.h"

subroutine RadTrans_finalize()
  use RadTrans_data
  implicit none
  type(PetscErrorCode) :: ierr

  if (.not. rt_useRadTrans) return

  if (alloced) call Petsc_dealloc()

  call PetscFinalize(ierr)

  return
end subroutine RadTrans_finalize
