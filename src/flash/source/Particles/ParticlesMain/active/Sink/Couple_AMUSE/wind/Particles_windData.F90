!!****if* source/Particles/ParticlesMain/active/Sink/Couple_AMUSE/wind/Particles_windData
!!
!! NAME
!!
!!    Particles_windData
!!
!! SYNOPSIS
!!
!!    Particles_windData()
!!
!! DESCRIPTION
!!
!!    Module to hold local variables and data types for wind unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

module Particles_windData

  implicit none

  real*8, save :: min_wind_dt      = 1d99
! Wind injection radius max. Negative means use 3.5*sqrt(3.0)*min_dx
! Default ref_radius updated to 10*delta in inject_direct.F90 
! to account for jet injection region -SA 20230918
  real*8, save :: ref_radius       = -1.0 
  real*8, save :: min_radius       = 0.0d0 ! Wind injection radius min.
  real*8, save :: min_wind_mass    = 0.0d0 ! smallest star that makes a wind (in grams).
  real*8, save :: wind_target_temp = 1d6 ! Target temperature from wind shock.

! Jet parameters - Added 20230728 SA
  real*8, save :: min_jet_mass     = 1.988d35 ! smallest star that makes a jet (value matches Config)
  real*8, save :: max_jet_mass     = 1.988d31 ! largest star that makes a jet (value matches Config)
  real*8, save :: jet_time         = 0.0d0 ! length of time over which to inject jet (in seconds)

  logical, save    :: use_wind_compute_dt =.true.
  logical, save    :: mass_load           =.false. ! Mass load winds?
  logical, save    :: add_therm_e         =.false.
  logical, save    :: var_radius          =.false.
  logical, save    :: perturb_velocity    =.false.
  real*8,  save    :: perturb_std_dev

end module Particles_windData
