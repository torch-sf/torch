!!****if* source/physics/RadTrans/RadTransMain/VETTAM/PhotoIonization/rt_ionisedata.F90
!!
!! NAME
!!
!!  rt_ionisedata
!!
!! SYNOPSIS
!!  rt_ionisedata()
!!
!! DESCRIPTION
!!  Stores the local data for the EUV ionising routines
!!
!!  AUTHOR
!!  Shyam Harimohan Menon(2022) 
!!
!!
!!***
module rt_ionisedata
  !==============================================================================
  
  implicit none

#include "Flash.h"
#include "constants.h"

  ! Runtime parameters
  real, save    :: alpha_A, alpha_B, alpha_rec_constant, alpha_ground_constant, C_tss, C_therm, ion_conv_rtol, ion_conv_atol
  integer, save :: ion_type, ion_conv_maxits
  logical, save :: useEUVIonize, useH2Ionize, ion_implicit, ion_ots,  multiple_ionbands
#ifdef IHA_SPEC
  real, save    :: hA, hpA, elecA, h2A
#endif
  real, save    :: hnu, energyPerIonH, energyPerIonH2, energyPerIonH_13p6_15p2
  real, save    :: ion_sigmaH, ion_sigmaH2, ion_sigmaH2_15p2_infty, ion_sigmaH_13p6_15p2, &
                    ion_sigmaH_15p2_infty
  character(len=MAX_STRING_LENGTH) :: alpha_type

!==============================================================================

end module rt_ionisedata