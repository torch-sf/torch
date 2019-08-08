!!****if* source/Simulation/SimulationMain/SNSBBox/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  Set neutral and ionized medium properties in the Multispecies unit.
!!
!!***

subroutine Simulation_initSpecies()

#include "Multispecies.h"
#include "Flash.h"

#ifndef IHP_SPEC
  return
#else

  use Simulation_data
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  ! Simulation_initSpecies() is called by Multispecies_init(),
  ! which precedes Simulation_init()
  ! Must grab parameters here, cannot get in Simulation_init()
  ! or else Multispecies parameters will NOT be set!
  call RuntimeParameters_get('sim_A_n', sim_A_n)
  call RuntimeParameters_get('sim_gamma_n', sim_gamma_n)
  call RuntimeParameters_get('sim_A_i', sim_A_i)
  call RuntimeParameters_get('sim_gamma_i', sim_gamma_i)
  call RuntimeParameters_get('he_abundM', sim_abundM)  !  from Heat/HeatMain/HeatCool/phenHeat
  call RuntimeParameters_get('he_metal', sim_metal)  !  from Heat/HeatMain/HeatCool/phenHeat

  ! effective weight of atomic hydrogen
  sim_abar = 1.0 + sim_abundM*sim_metal

  ! total number of constituents for neutral medium
  if (sim_A_n .lt. 0) then
    sim_A_n = sim_abar / (1. + sim_abundM )
  endif

  ! total number of constituents for ionized medium
  if (sim_A_i .lt. 0) then
    sim_A_i = sim_abar / (2. + sim_abundM )
  endif

  call Multispecies_setProperty(IHP_SPEC, A, sim_A_i)
  call Multispecies_setProperty(IHP_SPEC, Z, sim_A_i)
  call Multispecies_setProperty(IHP_SPEC, GAMMA, sim_gamma_i)

  call Multispecies_setProperty(IHA_SPEC, A, sim_A_n)
  call Multispecies_setProperty(IHA_SPEC, Z, sim_A_n)
  call Multispecies_setProperty(IHA_SPEC, GAMMA, sim_gamma_n)

#endif
! ifndef IHP_SPEC ... else

end subroutine Simulation_initSpecies
