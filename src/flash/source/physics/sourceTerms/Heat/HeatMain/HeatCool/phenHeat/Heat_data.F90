!!  ++ / 
!!  + / -   Heating/Cooling
!!   / --
!!written by C. Baczynski, 2012-2013

module Heat_data

#include "Flash.h"

  implicit none

  !heating parameter
  real, save		:: he_crIonRate, he_crIonNH, he_crIonExp, he_crIonEnergy
  real, save    :: he_Gzero, he_pe_norm
  character(len=4), save :: he_pe_recipe
  !some tolerances and initial values for root finder in heating calculation
!  real, save		:: he_x1, he_x2, he_tol
  ! radiation variables
  real, save		:: he_tradmin, he_tradmax, he_dradmin, he_dradmax, he_h_uv
  ! heating thresholds
  real, save		:: he_theatmin, he_theatmax, he_absTmin, he_absTmax
  real, save    :: he_dust_sputter_temp ! Temperature at which dust things switch off.
  !physical constants in cgs
  real, save		:: he_boltz, he_protonmass, he_smallpres
  logical,save	:: he_stratifyHeat
  logical,save	:: he_useHeat
  logical, save :: he_use_cr_heating
  ! gas properties
  real, save    :: he_abar, he_abundM, he_metal
  logical,save	:: he_coolOff
! factor by which cooling time is multiplied for more subcycling
  real, save		:: he_subfactor
! 2/3 <m> /kboltz for the implicit solver
!  real, save		:: he_kconst

! maximum number of subcycles
  integer, save	:: he_meshMe

! switch to implicit scheme if cooling rate falls below implicitTol
  real, save		:: hy_eswitch, he_dtThres

! for output
  integer, parameter :: he_funit_log = 22
  character(len=80),save  :: he_outfile = "phenHeat"

end module Heat_data
