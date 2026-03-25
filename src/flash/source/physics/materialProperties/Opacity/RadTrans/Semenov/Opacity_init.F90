!!****f* source/physics/materialProperties/Opacity/RadTrans/Semenov
!!
!! NAME
!!
!!  Opacity_init
!!
!!
!! SYNOPSIS
!!
!!  call Opacity_init()
!!
!! DESCRIPTION
!!
!! Initialiazed data for the Semonov opacity model using run time
!! parameters.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!***
subroutine Opacity_init()
  use SemenovOpacities, ONLY: init_opacities
  implicit none

  call init_opacities()

end subroutine Opacity_init
