!!****if* source/Simulation/SimulationMain/SNSBBox/Simulation_init
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
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!!  Reference:  Gardiner & Stone JCP 205(2005),509-539
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Eos_interface, ONLY : Eos

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('Bx0',     sim_Bx0)
  call RuntimeParameters_get('By0',     sim_By0)
  call RuntimeParameters_get('Bz0',     sim_Bz0)
  call RuntimeParameters_get('rho',     sim_rho)
  call RuntimeParameters_get('pres',    sim_pres)
  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('smallp',  sim_smallP)  ! from another unit

  call RuntimeParameters_get('sim_tdust', sim_tdust)
  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)
  call RuntimeParameters_get('killdivb', sim_killdivb)

end subroutine Simulation_init
