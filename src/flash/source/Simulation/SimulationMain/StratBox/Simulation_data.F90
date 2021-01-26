!!****if* source/Simulation/SimulationMain/StratBox/Simulation_data
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
!!  Stores the local data for Simulation setup: StratBox
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
  !! stratbox, nested refinement
  real, save :: sim_zrefcut1, sim_zrefcut2, sim_zrefcut3, sim_zrefcut4
  real, save :: sim_zrefcut5, sim_zrefcut6, sim_zrefcut7
  integer, save :: sim_targetRef
  logical, save :: sim_useNestedRef

  !! chemistry
  real, save :: sim_tdust
  real, save :: sim_init_Hp ! sim_init_H2, sim_init_CO
  real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
  real, save :: sim_abar, sim_abundM, sim_metal

  !! B field
  real, save :: sim_magx, sim_magy, sim_magz

  logical, save :: sim_withStaticGrav

<<<<<<< HEAD
  logical, save :: sim_stirLayer
  real, save :: sim_stirH
=======
  logical, save :: sim_killdivb
>>>>>>> f7dbe88... Bfield fix for StratBox

end module Simulation_data
