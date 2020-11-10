!!****if* source/Simulation/SimulationMain/StratBox/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer myPE)
!!
!! ARGUMENTS
!!
!!    myPE      Current Processor Number
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: i,j,k, istat, ii,jj,kk
  character(len=255),save :: dumstr

  real :: eint, ek ! for sim_cubeFile velocities
  integer :: sim_comm, sim_myPE


  !AT 20190221 - use the multispecies gamma instead.  We are always using 5/3 so no effect
  !call RuntimeParameters_get('gamma', sim_gamma) ! overlap w/ Eos_init.F90
  call RuntimeParameters_get('smlrho', sim_smlrho) ! declared by Hydro/HydroMain/unsplit
  call RuntimeParameters_get('smallp', sim_smallp) ! declared by Hydro/HydroMain/unsplit
  !call RuntimeParameters_get('smallX', sim_smallX) ! declared by Hydro/HydroMain/unsplit

  call RuntimeParameters_get('sim_bx0', sim_magx)
  call RuntimeParameters_get('sim_by0', sim_magy)
  call RuntimeParameters_get('sim_bz0', sim_magz)

  ! stratbox initialization
  call RuntimeParameters_get('sim_useStrat',  sim_useStrat)
  call RuntimeParameters_get('sim_p',  sim_p)
  call RuntimeParameters_get('sim_rho',  sim_rho)
  call RuntimeParameters_get('sim_pIGM',  sim_pIGM)
  call RuntimeParameters_get('sim_rhoIGM',  sim_rhoIGM)
  ! stratbox, gravity profile parameters, stellar disk
  call RuntimeParameters_get('sim_aParm1',  sim_aParm1)
  call RuntimeParameters_get('sim_aParm2',  sim_aParm2)
  call RuntimeParameters_get('sim_aParm3',  sim_aParm3)
  call RuntimeParameters_get('sim_aParm4',  sim_aParm4)
  ! stratbox, z-nested refinement
  call RuntimeParameters_get('sim_zrefcut1', sim_zrefcut1)
  call RuntimeParameters_get('sim_zrefcut2', sim_zrefcut2)
  call RuntimeParameters_get('sim_zrefcut3', sim_zrefcut3)
  call RuntimeParameters_get('sim_zrefcut4', sim_zrefcut4)
  call RuntimeParameters_get('sim_zrefcut5', sim_zrefcut5)
  call RuntimeParameters_get('sim_zrefcut6', sim_zrefcut6)
  call RuntimeParameters_get('sim_zrefcut7', sim_zrefcut7)
  call RuntimeParameters_get('sim_useNestedRef', sim_useNestedRef)
  call RuntimeParameters_get('sim_targetRef', sim_targetRef)

  ! chemistry
  call RuntimeParameters_get('sim_tdust', sim_tdust)
  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)

  ! stratbox, toggle staticgrav
  call RuntimeParameters_get('sim_withStaticGrav', sim_withStaticGrav)

  call RuntimeParameters_get('killdivb', sim_killdivb)

end subroutine Simulation_init
