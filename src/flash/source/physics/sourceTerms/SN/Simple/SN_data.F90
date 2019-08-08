!!****if* source/physics/sourceTerms/SN/Simple/SN_data
!!
!! NAME
!!
!!    SN_data
!!
!! SYNOPSIS
!!
!!    use SN_data
!!
!! DESCRIPTION
!!
!!    Holds variable within SN unit scope
!!
!!***

module SN_data

#include "Flash.h"
#include "constants.h"

  implicit none

  !SN parameter 
  logical,save  :: sn_useSN, sn_stratifySN, sn_SNmapToGrid
  real, save    :: sn_tsn1, sn_tsn2, sn_tstop, sn_r_init, sn_r_exp_max
  real, save    :: sn_exp_energy, sn_Mejc, sn_max_temp
  real, save    :: sn_hstar1, sn_hstar2
  !SN number 
  integer, save ::  sn_nsndt, sn_nSN, sn_nms, sn_callRNG
  ! extent of the simulation box 
  real, save    :: sn_imin, sn_imax, sn_jmin, sn_jmax, sn_kmin, sn_kmax
  ! Storage for timestep calculation
  real, save    :: sn_SNminstep
  real          :: sn_newDt
  ! I/O
  integer, save :: sn_meshMe
  character (len=MAX_STRING_LENGTH), save :: sn_outputDir
  character(len=MAX_STRING_LENGTH) :: sn_outfile = "SNfeedback.dat"
  integer, parameter :: sn_funit = 15

  ! modify how to deposit SN energy/momentum on grid
  logical, save :: sn_kinetic

  ! how to place SNe on grid
  character(len=MAX_STRING_LENGTH), save :: sn_fieldMode
  real, save :: sn_single_x, sn_single_y, sn_single_z

end module SN_data
