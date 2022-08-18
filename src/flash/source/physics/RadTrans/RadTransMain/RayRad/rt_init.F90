!!****if* source/physics/RadTrans/RadTransMain/MGD/rt_init
!!
!!  NAME 
!!
!!  rt_init
!!
!!  SYNOPSIS
!!
!!  call rt_init()
!!
!!  DESCRIPTION 
!!    Initialize local data for each radiative transfer model
!!
!!***

! c  is rt_speedlt ! Speed of light
! kb is rt_boltz ! Boltzmann constant
! sb is rt_radconst  ! stefan-boltzmann constant
subroutine rt_init
  use rt_data

  use RadTrans_data, ONLY: rt_useRadTrans, rt_meshCopyCount, &
       rt_meshMe, rt_acrossMe, rt_radconst, rt_boltz, rt_speedlt

  use Driver_interface, ONLY : Driver_abortFlash
  use Simulation_data, ONLY : sim_abar
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Eos_interface, ONLY : Eos_getAbarZbar

  implicit none

#include "Flash.h"
#include "constants.h"

  real, allocatable, dimension(:) :: temparr, lumarr, massarr
  real    :: totflux
  integer :: i, s
  real    :: temper
  integer :: ntab

!  call RuntimeParameters_get('rt_nrOfSpectralClasses', rt_nrOfSpectralClasses)
  call RuntimeParameters_get('rt_maxHchange', rt_maxHchange)
  call RuntimeParameters_get("rt_rayTrace", rt_rayTrace)
  call RuntimeParameters_get("abar", rt_abar)

  call RuntimeParameters_get('he_abundM', rt_abundM)
  call RuntimeParameters_get('he_metal', rt_metal)

! some of these initialisations are taken from doric_rad_ini
  call PhysicalConstants_get("Planck", rt_planck)
  call PhysicalConstants_get("ideal gas constant", rt_idealgas)
  call PhysicalConstants_get("proton mass", rt_protonMass)
  call PhysicalConstants_get("Newton", rt_Newton)
  call PhysicalConstants_get("Stefan-Boltzmann", rt_stboltz)

  rt_abar = 1.0 + rt_abundM*rt_metal
  
  tpic2 = 2.0 * PI / (rt_speedlt*rt_speedlt)
  tc2 = 2.0 / (rt_speedlt*rt_speedlt)

! DORIC abundances
!
! Set heavy element abundances (Osterbrock 1989, Table 5.13)
! C abundance is used to ensure non-zero electron densities)
! gives a minimum of electrons to do stuff with
  abu_c = 7.1e-7

! some parameters for the photoionization heating
  call RuntimeParameters_get('gamma', rt_gamma1)
  rt_gamma1 = rt_gamma1 - 1.0

end subroutine rt_init
