!!****if* source/physics/sourceTerms/GridInject/GridInjectMain/GridInject_data
!!
!! NAME
!!
!!    GridInject_data
!!
!! SYNOPSIS
!!
!!    use GridInject_data
!!
!! DESCRIPTION
!!
!!    Holds variables within GridInject unit scope
!!
!!***

module GridInject_data

#include "Flash.h"
#include "constants.h"

  implicit none

  integer, save :: gi_maxref

  character(len=MAX_STRING_LENGTH), save :: wind_cons_quant
  logical, save :: wind_add_therm_e
  logical, save :: wind_mass_load
  real, save    :: wind_target_temp
  real, save    :: wind_min_radius
  logical, save :: wind_var_radius
  logical, save :: wind_perturb
  real, save    :: wind_perturb_stdev

end module GridInject_data
