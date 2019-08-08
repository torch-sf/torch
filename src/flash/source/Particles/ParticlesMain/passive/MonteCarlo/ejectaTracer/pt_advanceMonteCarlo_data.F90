!===============================================================================

module pt_advanceMonteCarlo_data
!===============================================================================

  implicit none

!-------------------------------------------------------------------------------

!! random number seed
  integer*4, save :: jran
  real, save      :: ej_TempUnfreeze, ej_TimeUnfreeze
  logical ,save   :: ej_shockUnfreeze

!! order of traversal matrix
! 1,2,3|1,3,2|2,1,3|2,3,1|3,1,2|3,2,1|
  integer*2, save, dimension(3,6)  :: orderM = RESHAPE((/1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/),shape(orderM))
! (/(/1,2,3/),(/1,3,2/),(/2,1,3/),(/2,3,1/),(/3,1,2/),(/3,2,1/)/)

end module pt_advanceMonteCarlo_data
