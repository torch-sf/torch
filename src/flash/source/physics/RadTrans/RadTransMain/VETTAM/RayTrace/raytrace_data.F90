!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/raytrace_data
!!
!!  NAME 
!!    raytrace_data
!!
!!  SYNOPSIS
!!    use raytrace_data
!!
!!  DESCRIPTION 
!!    Stores data for HybridChar3DRT Ray Tracer
!!
!!***

#include "Flash.h"

module raytrace_data
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
  implicit none

  logical, save :: rt_useRayTrace
  integer, save :: rt_nPhi
  integer, save :: rt_nTheta

  integer, save :: rt_nrOfAngleGroups
  integer, save :: rt_maxNrOfBoundIter
  integer, save :: rt_ALI
  real, save    :: rt_epsilon

  integer, save :: rt_healpix_nSide
  integer, save :: rt_healpix_randomize
  
  real, save :: rt_dirX
  real, save :: rt_dirY
  real, save :: rt_dirZ

  real, save :: irradiation
  !
  ! Internal reference at which dr_simGeneration the memory allocation
  ! has been done. If dr_simGeneration is raised, a refinement or
  ! derefinement has happend and we have to reallocate memory as well
  ! as initialize it again.
  INTEGER, save :: allocGeneration

  real, dimension(:,:,:,:,:,:,:), ALLOCATABLE :: faceValueAll
  integer, dimension(:), ALLOCATABLE        :: nrOfLeafBlocksList
  integer, dimension(:,:), ALLOCATABLE      :: leafListAll

  real, dimension(:,:,:,:), pointer :: solnData

  integer, dimension(:,:,:,:), allocatable :: SurrBlkSumArray
  integer, dimension(:,:), allocatable :: SurrBlkSumNumNegh
  type (AllBlockRegions_t), allocatable :: SurrBlkSum(:)
  integer, DIMENSION(:,:), allocatable :: rblockListAll

end module raytrace_data
