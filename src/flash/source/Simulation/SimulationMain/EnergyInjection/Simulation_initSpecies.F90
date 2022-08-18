!!****if* source/Simulation/SimulationMain/EnergyInjection/Simulation_initSpecies
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
!!  This routine will initialize the species and species values needed for
!!  the EnergyInjection setup, which advects two fluids with different
!!  Gamma values
!!
!!***

subroutine Simulation_initSpecies()

  use Simulation_data
  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
!#include "RuntimeParameters.h"
#include "Multispecies.h"
#include "Flash.h"
!#include "Multispecies_interface.h"

  call RuntimeParameters_get('sim_A_n', sim_A_n)
  call RuntimeParameters_get('sim_gamma_n', sim_gamma_n)

  call RuntimeParameters_get('sim_A_i', sim_A_i)
  call RuntimeParameters_get('sim_gamma_i', sim_gamma_i)

  call RuntimeParameters_get('he_abundM', sim_abundM)
  call RuntimeParameters_get('he_metal', sim_metal)

! effective weight of atomic hydrogen
  if ((sim_A_i .lt. 0) .or. (sim_A_n .lt. 0)) then
    sim_abar = 1.0 + sim_abundM*sim_metal
  endif

! total number of constituents for neutral medium
  if (sim_A_n .lt. 0) then
    sim_A_n = sim_abar / (1. + sim_abundM )
  endif

! total number of constituents for ionized medium
  if (sim_A_i .lt. 0) then
    sim_A_i = sim_abar / (2. + sim_abundM )
  endif

!  print*,'For IHA_SPEC we use sim_A_n=',sim_A_n
!  print*,'For IHP_SPEC we use sim_A_i=',sim_A_i

  call Multispecies_setProperty(IHP_SPEC, A, sim_A_i)
  call Multispecies_setProperty(IHP_SPEC, Z, sim_A_i)
  call Multispecies_setProperty(IHP_SPEC, GAMMA, sim_gamma_i)

  call Multispecies_setProperty(IHA_SPEC, A, sim_A_n)
  call Multispecies_setProperty(IHA_SPEC, Z, sim_A_n)
  call Multispecies_setProperty(IHA_SPEC, GAMMA, sim_gamma_n)

end subroutine Simulation_initSpecies
