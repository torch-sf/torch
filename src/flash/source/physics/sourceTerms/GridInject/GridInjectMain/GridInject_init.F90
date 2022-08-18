!!****f* source/physics/sourceTerms/GridInject/GridInject_init
!!
!! NAME
!!
!!  GridInject_init
!!
!! SYNOPSIS
!!
!!  call GridInject_init()
!!
!! DESCRIPTION
!!
!!  Perform various initializations for the GridInject unit.
!!
!! ARGUMENTS
!!
!!  None
!!
!!***
subroutine GridInject_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use GridInject_data

  implicit none

  call RuntimeParameters_get("lrefine_max", gi_maxref) ! from paramesh unit

  call RuntimeParameters_get("wind_cons_quant",    wind_cons_quant)
  call RuntimeParameters_get("wind_add_therm_e",   wind_add_therm_e)
  call RuntimeParameters_get("wind_mass_load",     wind_mass_load)
  call RuntimeParameters_get("wind_target_temp",   wind_target_temp)
  call RuntimeParameters_get("wind_min_radius",    wind_min_radius)
  call RuntimeParameters_get("wind_var_radius",    wind_var_radius)
  call RuntimeParameters_get("wind_perturb",       wind_perturb)
  call RuntimeParameters_get("wind_perturb_stdev", wind_perturb_stdev)

  return
end subroutine GridInject_init
