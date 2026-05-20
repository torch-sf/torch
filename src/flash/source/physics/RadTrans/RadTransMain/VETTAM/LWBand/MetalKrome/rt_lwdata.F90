!!****if* source/physics/sourceTerms/Cool/CoolMain/RadTransVET/EUVIonise/rt_lwdata.F90
!!
!! NAME
!!
!!  rt_lwdata
!!
!! SYNOPSIS
!!  rt_lwdata()
!!
!! DESCRIPTION
!!  Stores the local data for the FUV dissociation routines
!!
!!  AUTHOR
!!  Shyam Harimohan Menon(2023) 
!!
!!
!!***
module rt_lwdata
  !==============================================================================
  
  implicit none

#include "Flash.h"
#include "constants.h"

  ! Runtime parameters
  real, save    :: bfive, fpump, energyPerDissociation, avgEnergyLW
  logical, save :: useH2Dissociate, useHIshield
  character(len=MAX_STRING_LENGTH) :: lwdiss_type

#ifdef C_SPEC
  real, save :: cA
#endif

#ifdef CO_SPEC
  real, save :: coA
#endif

!==============================================================================

end module rt_lwdata