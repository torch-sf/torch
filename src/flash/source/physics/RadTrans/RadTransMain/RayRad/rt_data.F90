!! used values
!! alpha = 	2.59e-13 (bh00) case B recombination coefficient of H at 1e4 K 
!! fitted function to alpha(dens) are not used 
!! sigma =  6.30e-18 (ah0) Lyman limit HI photoionization cross section


!!****if* source/physics/RadTrans/RadTransMain/MGD/rt_data
!!
!!  NAME 
!!    rt_data
!!
!!  SYNOPSIS
!!    use rt_data
!!
!!  DESCRIPTION 
!!    Stores data for MGD
!!
!!***
module rt_data
  implicit none
  
! runtime parameters
  real, save    :: rt_planck  
  real, save    :: rt_idealgas
  real, save    :: rt_protonMass
  real, save    :: rt_Newton
  real, save    :: rt_stboltz
  real, save    :: rt_maxHchange

  real, save    :: rt_abar
  real, save    :: rt_gamma1

  real, save    :: rt_abundM, rt_metal

  real, save    :: rt_dt
  real, save, dimension(5) :: rt_dt_pos

  logical,save :: rt_rayTrace

! Physical constants (initialized in doric_get_constants).

  real, save :: tpic2  ! 2*pi/c^2
  real, save :: tc2  ! 2*pi/c^2

! Heavy element abundances.

  real,save  :: abu_c

! Astrophysical quantities.

  real, parameter :: r_solar=6.96e10
  real, parameter :: l_solar=3.862e33

! 13.6 ev in ergs
  real, parameter :: rt_ionHIener = 2.17863955e-11

! Set the values of physical constants:
!   ev2k:   conversion factor between evs and kelvins.
!   ev2erg: conversion factor between evs and ergs.
!   ev2fr:  conversion factor between evs and sclfre hertz.
     
  real, parameter :: ev2k=1.0/8.617e-05
  real, parameter :: ev2erg=1.602e-12
  real, parameter :: ev2fr=0.241838e15

! Hydrogen recombination coefficient (case B)
! (See Osterbrock (1989), table 2.1)
     
  real, parameter :: albpow=-0.7
  real, parameter :: bh00=2.59e-13

! Collisional ionization data for hydrogen (Cox 1970).
     
  real, parameter :: eth0=13.598
  real, parameter :: xih0=1.0
  real, parameter :: fh0=0.83
  real, parameter :: colh0=1.3e-8*fh0*xih0/(eth0*eth0)
  real, parameter :: temph0=eth0*ev2k
  real, parameter :: hionen=eth0*ev2erg

! Photoionization constants for hydrogen (Osterbrock 1989):
!  ah0:   cross section at threshold.
!  frth0: threshold frequency.
!  sh0:   power of frequency dependence of cross section.
     
  real, parameter :: frth0=ev2fr*eth0
  real, parameter :: ah0=6.30e-18
  real, parameter :: betah0=1.34
  real, parameter :: sh0=2.8


end module rt_data
