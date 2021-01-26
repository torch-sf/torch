!!****if* source/Simulation/SimulationMain/EnergyInjection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for EnergyInjection test problem
!!
!! ARGUMENTS
!!
!!   
!!
!!
!!***
subroutine Simulation_init() 
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  ! get the runtime parameters relevant for this problem

  call RuntimeParameters_get('gamma', sim_gamma)

  call PhysicalConstants_get('ideal gas constant', sim_gasconstant)
  call PhysicalConstants_get('proton mass', sim_protonMass)

  call RuntimeParameters_get('amTemp', sim_amTemp)
  call RuntimeParameters_get('amNumDens', sim_amNumDens)

  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)
  
  call RuntimeParameters_get( 'bx0', sim_magx)
  call RuntimeParameters_get( 'by0', sim_magy)
  call RuntimeParameters_get( 'bz0', sim_magz)

  call RuntimeParameters_get( 'sim_tdust', sim_tdust)

  call RuntimeParameters_get('killdivb', sim_killdivb)

  sim_abar = 1.0 + sim_abundM*sim_metal

end
