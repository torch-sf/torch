!!****if* source/Simulation/SimulationMain/Cube_Velocity/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!   use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Stratbox_AT
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  real, save :: sim_smallp, sim_smlrho, sim_gamma

  !! stratbox
  logical, save :: sim_useStrat
  real, save :: sim_p, sim_rho, sim_pIGM, sim_rhoIGM
  real, save :: sim_aParm1, sim_aParm2, sim_aParm3, sim_aParm4

  !! chemistry
  real, save :: sim_tdust
  real, save :: sim_init_Hp ! sim_init_H2, sim_init_CO
  real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
  real, save :: sim_abar, sim_abundM, sim_metal

  !! B field
  real, save :: sim_magx, sim_magy, sim_magz

  !! 3D arrays with initial conditions - first the runtime parameter
  character(len=255), save :: sim_velcubeFile
  real, save :: sim_machTurb  ! to help set the correct velocity magnitude
  !! 3D arrays with initial conditions
  integer, save :: sim_nCD(MDIM)
  real, save :: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin
  real, allocatable, save :: sim_velxArr(:,:,:), sim_velyArr(:,:,:)
  real, allocatable, save :: sim_velzArr(:,:,:)

end module Simulation_data
