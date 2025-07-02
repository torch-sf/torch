!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RadTrans_data
!!
!! NAME
!!
!!  RadTrans_data
!!
!! SYNOPSIS
!!  RadTrans_data()
!!
!! DESCRIPTION
!!  Stores the local data for the VETTAM module
!!
!!  AUTHOR
!!  Shyam Harimohan Menon(2020-2024)
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in the VETTAM unit.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have over written these values with the flash.par values
!!   for your specific run.  
!!
!!
!!***

module RadTrans_data

#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscvec.h"

  use petscksp
  use petscvec
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t

!==============================================================================
  use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
  implicit none

#include "Flash.h"
#include "constants.h"
!include "mpif.h"

  ! Runtime parameters
  logical, save :: rt_useRadTrans, rt_compute_Dt, rt_update_hydro,rt_debug, &
                  use_hlle_correction, rt_sink_implicit, rt_picard_correction, &
                  rt_sinkheat, rt_sinkmom, rt_ovcterms, noSink_VETTAM
  integer, save :: rt_Dt_type, rt_hydro_type, rt_temp_type

  real, save    :: rt_dtFactor, rt_smalle, rt_smallf
  real, save    :: rt_T_min

  ! Physical constants:
  real, save :: rt_radconst, rt_speedlt, rt_boltz ! Radiation constant, Speed of light, Boltzmann constant

  !Parameters to control frequency of VETTAM (only use rt every dt)
  real, save :: rt_useeverydt

  !Boundary condition stuff
  INTEGER, save :: xl_rad_BC, xr_rad_BC, yl_rad_BC, &
                   yr_rad_BC, zl_rad_BC, zr_rad_BC
  INTEGER, save, DIMENSION(2,NDIM) :: rad_BC
  real, save    :: rt_boundary_T, &
                    fradx_boundary_value, frady_boundary_value, &
                    fradz_boundary_value, rt_stellar_opacity

  !Tolerance/iterations and solver specifics of the implicit Matrix solver
  real, save    :: rt_rtol, rt_abstol, rt_dtol, rt_picard_rtol, &
                   rt_picard_abstol,rt_picard_theta,rt_nr_rtol, rt_nr_dtol
  integer, save :: rt_maxits, rt_minits,rt_picard_maxits,rt_nr_maxits 
  character(len=MAX_STRING_LENGTH), save :: rt_preconditioner_str, &
    rt_linearsolver_str, rt_stellarop_type, rt_etens

  !Multiple frequency bands (by default 3)
  integer, save :: rt_freqbands
  !RT bands; by default this is : IR, FUV/Optical, EUV, PE and Lyman Werner respectively in order; 
  ! this has to be in order of decreasing energy from 5 to 1
  character(len=MAX_STRING_LENGTH), save :: rt_band1, rt_band2, rt_band3, rt_band4, rt_band5
  !dust opacity types in each band
  character(len=MAX_STRING_LENGTH), save :: dust_opacity_ir, dust_opacity_fuv, dust_opacity_euv
  !Opacity in different bands if opacity type is constant
  real, save :: dusttoGasRatio, dust_iropac_value,dust_fuvopac_value, dust_euvopac_value, dust_lwopac_value, dust_peopac_value
  !Parameter for current band being updated. This is just for internal bookkeeping for updates; not set by user.
  character(len=MAX_STRING_LENGTH), save :: current_band

  !Log stuff (for all bands)
  integer, save:: cumulative_iter_band1, cumulative_iter_band2, cumulative_iter_band3, max_iter_band1 &
                    ,max_iter_band2, max_iter_band3, min_iter_band1, min_iter_band2, min_iter_band3, &
                    cumulative_iter_band4, cumulative_iter_band5, max_iter_band4, max_iter_band5, &
                    min_iter_band4, min_iter_band5

  ! petsc
  INTEGER, save :: blockCount = 0
  INTEGER, save, DIMENSION(MAXBLOCKS) :: blockList 
  type(tKSP) :: ksp
  !A separate Xvec for each band to store the previous time guess
  type(tVec) :: Xvec_Band1, Xvec_Band2, Xvec_Band3,Xvec_Band4,Xvec_Band5,Yvec
  type(tMat) :: Amat
  PetscFortranAddr ctx
  type(PetscInt) :: ncells


  ! AMR Related Stuff.
  type (AllBlockRegions_t), allocatable :: SurrBlkSum(:)
  integer, DIMENSION(:,:), allocatable :: rblockListAll
  integer, allocatable ::  NeghLevels(:,:,:,:)

  INTEGER, PARAMETER :: comm = MPI_COMM_WORLD
  LOGICAL, save :: alloced
  logical, parameter :: Petsc_Log = .true.

  !Legacy stuff in FLASH Radtrans default MGD solver that (might) be used in other parts of the code
  integer, save :: rt_meshMe ! Process rank
  ! The number of replicated meshes active in this simulation
  integer, save :: rt_meshCopyCount
  ! The mesh number of this process
  integer, save :: rt_acrossMe
  ! The across communicator
  integer, save :: rt_acrossComm
  ! Global communicator for all processes
  integer, save :: rt_globalComm

  ! Structure that holds context information on the current operation,
  ! for debugging
  type(RadTrans_dbgContext_t),save,target :: rt_dbgContext

!==============================================================================

end module RadTrans_data
