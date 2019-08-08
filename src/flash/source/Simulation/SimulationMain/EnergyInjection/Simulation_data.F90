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
   real, save :: sim_gamma, sim_smallX
   real, save :: sim_molarMass, sim_gasconstant
   real, save :: sim_ptMass, sim_protonMass

   real, save :: sim_p1x, sim_p1y, sim_p1z, sim_Eph, sim_Nph, sim_tdust
   integer, save :: sim_meshMe, sim_nPtot

!  ambient stuff
   real, save :: sim_amTemp, sim_amNumDens, sim_amDens, sim_abar
	
  !! *** EOS Parameters *** !!
!  real, save, dimension(EOS_NUM) :: sim_eosArr
!  integer, save :: sim_vecLen, sim_mode

  !! chemistry parameters
  real, save :: sim_init_H2, sim_init_Hp, sim_init_CO
! for chemistry: more gammas more fun
  real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
  real, save :: sim_abundM, sim_metal
  real, save :: sim_magx, sim_magy, sim_magz

end module Simulation_data
