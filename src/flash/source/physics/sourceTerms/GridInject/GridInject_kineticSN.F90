!!****f* source/physics/sourceTerms/GridInject/GridInject_kineticSN
!!
!! NAME
!!
!!  GridInject_kineticSN
!!
!! SYNOPSIS
!!
!!  call GridInject_kineticSN(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: energy
!!      real(IN)          :: mass
!!      optional,logical(IN) :: snap_to_grid
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit SN onto grid with both kinetic and thermal energy
!!
!! ARGUMENTS
!!
!!  xloc    : where to inject
!!  yloc    : where to inject
!!  zloc    : where to inject
!!  energy  : amount of energy to inject
!!  mass    : amount of mass to inject
!!  snap_to_grid : for testing/debugging
!!
!!***
subroutine GridInject_kineticSN (xloc, yloc, zloc, energy, mass, snap_to_grid)

  implicit none

  real, intent(IN)              :: xloc, yloc, zloc, energy, mass
  logical, optional, intent(IN) :: snap_to_grid

  return
end subroutine GridInject_kineticSN
