!!****if* source/Simulation/SimulationMain/EnergyInjection/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: EnergyInjection
!!
!!
!!***
module Simulation_data
  implicit none
#include "constants.h"
#include "Eos.h"
! single fluid stuff
  real, save :: sim_gamma, sim_abar
  real, save :: sim_gasconstant, sim_protonMass

!  ambient stuff
  real, save :: sim_amTemp, sim_amNumDens, sim_tdust
  real, save :: sim_magx, sim_magy, sim_magz

  !! chemistry parameters
  real, save :: sim_init_Hp
! for chemistry: more gammas more fun
  real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
  real, save :: sim_abundM, sim_metal
  logical, save :: sim_killdivb

end module Simulation_data
