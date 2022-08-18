!!****if* source/Particles/localAPI/pt_initPositionsLattice
!!
!! NAME
!! SYNOPSIS
!!
!! DESCRIPTION
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!! PARAMETERS
!!***


subroutine pt_initPositionsNParPerZone (blockID,success)


  implicit none

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  success = .FALSE.
  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsNParPerZone


