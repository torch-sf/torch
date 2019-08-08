!!****if* source/Simulation/SimulationMain/cube/Simulation_data
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
!!  Stores the local data for Simulation setup: cube
!!
!!
!!***


module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  character(len=255),save :: sim_cubeFile

  integer, save :: sim_comm, sim_myPE
  real, save :: smallp, smlrho, smallX, sim_boltz, sim_mH, sim_pi, sim_gamma

  !! 3D arrays with initial conditions
  integer, save :: sim_nCD(MDIM)
  real, save :: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin
  real, allocatable, save :: sim_densArr(:,:,:), sim_presArr(:,:,:)
  real, allocatable, save :: sim_gpotArr(:,:,:), sim_velxArr(:,:,:)
  real, allocatable, save :: sim_velyArr(:,:,:), sim_velzArr(:,:,:)


end module Simulation_data


