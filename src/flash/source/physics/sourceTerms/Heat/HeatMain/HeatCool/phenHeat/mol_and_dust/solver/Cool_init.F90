!!****f* source/physics/sourceTerms/Cool/molecules+dust/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***



subroutine Cool_init()

#include "constants.h"
#include "Flash.h"

  use Cool_data
  use Driver_data, ONLY : dr_globalMe
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  integer :: cpts

  call RuntimeParameters_get('T_cool_min', T_cool_min)
  call RuntimeParameters_get('nd_cool_min', nd_cool_min)
  call RuntimeParameters_get('nd_cool_max', nd_cool_max)
  call RuntimeParameters_get('tstep_cool_factor', tstep_cool_factor)
  call RuntimeParameters_get('useDustCool', useDustCool)
  !call RuntimeParameters_get('initialize_dust_restart', initialize_dust_restart)
  !call RuntimeParameters_get('tdust_restart', tdust_restart)
  call RuntimeParameters_get('he_int_method', he_int_method)

  he_int_method = trim(he_int_method)


  call RuntimeParameters_get('T_max', T_max)
  call PhysicalConstants_get('ideal gas constant', gasConstant)
  !call RuntimeParameters_get('eos_singleSpeciesA', mu_mol)

  call PhysicalConstants_get('pi', pi)
  call PhysicalConstants_get('Newton', newton)
  call PhysicalConstants_get('Stefan-Boltzmann', sigma)
  call PhysicalConstants_get('Boltzmann',kB)
  call RuntimeParameters_get('gamma', gamma)
  call PhysicalConstants_get('proton mass',mp)

  gammam1 = gamma - 1.

  call RuntimeParameters_get('T_max_core', T_max_core)
  call RuntimeParameters_get('T_max_core_radius', T_max_core_radius)

  cpts = get_cooling_data(cool_dat,DENS_PTS,TEMP_PTS)

! minimal temperature and number density available in cooling table
  T_min    = cool_dat(1,1,1)
  nd_min   = cool_dat(1,1,2)

if (dr_globalMe == 0) then

  print*, "[Dust/Mol Cooling]: Dust/Mol cooling on =", useDustCool, dr_globalMe
  print*, "[Dust/Mol Cooling]: T_cool_min, T_max, table T_min =", T_cool_min, T_max, T_min, dr_globalMe
  print*, "[Dust/Mol Cooling]: nd_cool_min, nd_cool_max, table nd_min =", nd_cool_min, nd_cool_max, nd_min, dr_globalMe
  print*, "[Dust/Mol Cooling]: Integration method for heating and cooling is: ", he_int_method
  
end if

!   call RuntimeParamters_get("density_floor", density_floor)
!   call RuntimeParamters_get("smlrho1", smlrho1)
!   call RuntimeParamters_get("smlrho1_dist", smlrho1_dist)
!   call RuntimeParamters_get("smlrho2", smlrho2)
!   call RuntimeParamters_get("smlrho2_dist", smlrho2_dist)
!   call RuntimeParamters_get("smlrho3", smlrho3)
!   call RuntimeParamters_get("smlrho3_dist", smlrho3_dist)
!   call RuntimeParamters_get("smlrho4", smlrho4)
!   call RuntimeParamters_get("smlrho4_dist", smlrho4_dist)
!   call RuntimeParamters_get("small_tstep", small_tstep)
!   call RuntimeParamters_get("lref_dens_floor", lref_dens_floor)
!   call RuntimeParamters_get("cfl", cfl)

  return
end subroutine Cool_init
