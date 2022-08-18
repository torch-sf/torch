!!****if* source/physics/sourceTerms/Heat/SN+Heat+Cool/Heat_data
!!
!! NAME
!!
!!  Hydro_data
!!
!!
!! SYNOPSIS
!!
!!  use Hydro_data
!!
!!
!! DESCRIPTION
!!
!!  Placeholder module containing private Bouchut5 MHD solver-specific data
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

!!TODO READ IN RHO AMBIENT AND PRESSURE AMBIENT FOR HEATING ROUTING
!!he_ should be sn_ but renaming stuff is no fun 
module SN_data

#include "Flash.h"

  implicit none

  real,save		:: he_smallpres, he_smalldens
  !SN parameter 
  real, save		:: he_tsn1, he_tsn2, he_r_init, he_r_exp_max
  !SN number 
  integer, save		::  he_nsndt, he_nSN, he_nms
  !vertical and radial extend of SN 
  real, save		:: he_hstar1, he_hstar2, he_erstar1, he_erstar2
  !set the maximum possible number of superbubbles and use it to declare vars
  integer, save		:: he_TracerPerSN
  ! Everybody should know these!
  integer, save	:: he_meshNumProcs, he_meshMe, he_oldStartTag
  !extends of the simulation box 
  real, save		:: he_imin, he_imax, he_jmin, he_jmax, he_kmin, he_kmax
  ! SN energy
  real, save		:: he_exp_energy, sn_max_temp
  ! distance to edge of SN region, for MC particle injection
  real, save		:: sn_MCTracerShellDis
  ! distance to edge of SN region, for MC particle injection
  real, save		:: sn_shellThicky
  real, save		:: sn_MassRatio
  ! mass to be collected to determine SN blast
  real, save		:: he_Mejc
!  switches for different physics, super-novae, -bubbles and stellar winds, stratification and subcycling 
  logical,save	:: he_stratifySN, he_radialSN
  logical,save	:: he_useSN, he_useSNrandom
  logical,save	:: he_SNmapToGrid, he_useSNTracer
  logical				:: he_exp_flag, he_veltracer
  ! Storage for timestep calculation
  real, save		:: he_SNminstep, he_maxradius
! factor by which cooling time is multiplied for more subcycling
  real					:: he_newDt

  real,allocatable,save, dimension(:,:,:) :: he_tracerPosX
  real,allocatable,save, dimension(:,:,:) :: he_tracerPosY
  real,allocatable,save, dimension(:,:,:) :: he_tracerPosZ

! for log file 
  character(len=80),save  :: sn_outfile = "SN_duds"
  integer, parameter :: sn_funit_evol = 15

end module SN_data
