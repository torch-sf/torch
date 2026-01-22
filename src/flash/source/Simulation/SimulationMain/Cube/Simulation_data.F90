!!****if* source/Simulation/SimulationMain/Cube/Simulation_data
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
!!  Stores the local data for Simulation setup: Cube
!!
!!
!!***


module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

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
  
#ifdef TRACER_FIELDS
  !! Arrays and variables for elements.
  integer, save :: sim_nelements
  real, allocatable, save :: sim_elemArr(:,:,:,:)
#endif

! single fluid stuff

   real, save :: sim_molarMass, sim_gasconstant
   real, save :: sim_protonMass

   real, save :: sim_tdust
   integer, save :: sim_meshMe

   !  ambient stuff
   real, save :: sim_abar
                    
   !! chemistry parameters
   real, save :: sim_init_H2, sim_init_Hp, sim_init_CO
   ! for chemistry: more gammas more fun
   real, save :: sim_A_n, sim_gamma_n, sim_A_i, sim_gamma_i
   real, save :: sim_abundM, sim_metal
   
   real, save :: sim_magx, sim_magy, sim_magz
   logical, save :: sim_killdivb
! New Parameterized heating and cooling parameters
!  real, save :: sim_Z, sim_G0, sim_pe_h, sim_cr_h
!  logical, save :: sim_stratify_heating, sim_constant_heating

   !! VorAMR stuff
   logical, save :: use_voramr, use_localRef, refPartCount, center_localRef
   character(len=255),save :: voramr_source, voramr_input
   real, save :: localRef_x, localRef_y, localRef_z, localRef_r
   !! Derefinement outside rectangular region of interest
   logical, save :: use_deref
   real, save :: deref_xl, deref_xr, deref_yl, deref_yr, deref_zl, deref_zr
   integer, save :: deref_lref
   
   !! static grav field parameters
   logical, save :: sim_withStaticGrav
   real, save :: sim_aParm1, sim_aParm2, sim_aParm3, sim_aParm4

 end module Simulation_data


