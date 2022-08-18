!!****if* source/Simulation/SimulationMain/SNSBBox/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: SNSBBox
!!
!!
!!***

module Simulation_data

  implicit none

  real, save :: sim_gamma
  real, save :: sim_smallP
  real, save :: sim_Bx0, sim_By0, sim_Bz0
  real, save :: sim_rho, sim_pres

  real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
  real, save :: sim_init_Hp, sim_tdust
  real, save :: sim_abundM, sim_metal
  real, save :: sim_abar  ! required by rt_init.F90
  logical, save :: sim_killdivb

end module Simulation_data
