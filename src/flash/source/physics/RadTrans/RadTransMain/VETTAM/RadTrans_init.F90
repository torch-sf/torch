!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RadTrans_init
!!
!! NAME
!!
!!  RadTrans_init
!!
!! SYNOPSIS
!!
!!  RadTrans_init()
!!
!! DESCRIPTION
!! 
!! Initialise the VETTAM module.
!!
!! ARGUMENTS
!!
!!   
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine RadTrans_init()

  use RadTrans_data
  use rt_interface, ONLY: rt_init
  use Driver_interface, ONLY            : Driver_abortFlash, Driver_getMype, Driver_getComm
  use Driver_data, ONLY                 : dr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"
#include "Flash.h"
#include "petsc/finclude/petscksp.h"

  implicit none
  PetscErrorCode :: ierr
  INTEGER        :: bc
!==============================================================================

  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)
  call RuntimeParameters_get("rt_debug",rt_debug)

  
  call RuntimeParameters_get("rt_T_min",rt_T_min)

  !Torch stuff 
  call rt_init

  ! Uniform UV background
  call RuntimeParameters_get("use_uv_bkgd", use_uv_bkgd)
  call RuntimeParameters_get("uv_bkgd_ion_rate", uv_bkgd_ion_rate)
  call RuntimeParameters_get("uv_bkgd_heat_rate", uv_bkgd_heat_rate)

  !Flag for HLLE Correction (see Jiang et al 2013)
  call RuntimeParameters_get("rt_hlle_correction",use_hlle_correction)

  ! these are not used in the iterative implementation (no timestep based on cooling)
  call RuntimeParameters_get("rt_compute_Dt",rt_compute_Dt)
  call RuntimeParameters_get("rt_Dt_type",rt_Dt_type)
  call RuntimeParameters_get("rt_dtFactor",rt_dtFactor)

  ! Limits to erad/flux
  call RuntimeParameters_get("rt_smalle",rt_smalle)
  call RuntimeParameters_get("rt_smallf",rt_smallf)

  ! PETSc Stuff
  call RuntimeParameters_get("rt_rtol", rt_rtol)
  call RuntimeParameters_get("rt_abstol", rt_abstol)
  call RuntimeParameters_get("rt_dtol", rt_dtol)
  call RuntimeParameters_get("rt_maxits", rt_maxits)
  call RuntimeParameters_get("rt_minits", rt_minits)
  call RuntimeParameters_get("rt_preconditioner", rt_preconditioner_str)
  call RuntimeParameters_get("rt_linearsolver", rt_linearsolver_str)

  !Multiband stuff
  call RuntimeParameters_get("rt_freqbands",rt_freqbands)
  call RuntimeParameters_get("rt_band1",rt_band1)
  call RuntimeParameters_get("rt_band2",rt_band2)
  call RuntimeParameters_get("rt_band3",rt_band3)
  call RuntimeParameters_get("rt_band4",rt_band4)
  call RuntimeParameters_get("rt_band5",rt_band5)

   ! Opacity stuff
  !Opacity types for each band
  call RuntimeParameters_get("dusttoGasRatio", dusttoGasRatio)
  call RuntimeParameters_get("dust_opacity_ir",dust_opacity_ir)
  call RuntimeParameters_get("dust_opacity_fuv",dust_opacity_fuv)
  call RuntimeParameters_get("dust_opacity_euv",dust_opacity_euv)
  !Constant opacity values if fixed opacity used
  call RuntimeParameters_get("dust_iropac_value",dust_iropac_value)
  call RuntimeParameters_get("dust_fuvopac_value",dust_fuvopac_value)
  call RuntimeParameters_get("dust_euvopac_value",dust_euvopac_value)
  call RuntimeParameters_get("dust_lwopac_value",dust_lwopac_value)
  call RuntimeParameters_get("dust_peopac_value",dust_peopac_value)
  call RuntimeParameters_get("rt_stellarop_type",rt_stellarop_type)
  call RuntimeParameters_get("rt_stellar_opacity",rt_stellar_opacity)

  !Eddington Tensor computation flag
  call RuntimeParameters_get("rt_etens",rt_etens)

  !O(v/c) terms correction flag
  call RuntimeParameters_get("rt_ovcterms",rt_ovcterms)

  !Hydro
  call RuntimeParameters_get("rt_update_hydro",rt_update_hydro)  

  !Flags for sink heating/momentum contributions
  call RuntimeParameters_get("rt_sinkheat",rt_sinkheat)
  call RuntimeParameters_get("rt_sinkmom",rt_sinkmom)

  !Implicit sink heating or not
  call RuntimeParameters_get("rt_sink_implicit",rt_sink_implicit)
  
  !Picard temperature correction
  call RuntimeParameters_get("rt_picard_correction",rt_picard_correction)

  !Flag for whether to use VETTAM if there are no sinks
  call RuntimeParameters_get("noSink_VETTAM",noSink_VETTAM)

  !Flag for whether to use VETTAM if there are no stars
  call RuntimeParameters_get("noStar_VETTAM",noStar_VETTAM)

  !Parameters to control frequency of VETTAM (only use rt every dt; default 0)
  call RuntimeParameters_get("rt_useeverydt",rt_useeverydt)

  if(rt_update_hydro) then
    !Picard iteration stuff
    call RuntimeParameters_get("rt_hydro_type",rt_hydro_type)  
    call RuntimeParameters_get("rt_temp_type",rt_temp_type)  
    if(rt_hydro_type .gt. 1) then
      call RuntimeParameters_get("rt_picard_rtol",rt_picard_rtol)
      call RuntimeParameters_get("rt_picard_abstol",rt_picard_abstol)
      call RuntimeParameters_get("rt_picard_theta",rt_picard_theta)  
      call RuntimeParameters_get("rt_picard_maxits",rt_picard_maxits) 
      if(rt_temp_type .eq. 1) then
        call RuntimeParameters_get("rt_nr_rtol",rt_nr_rtol)
        call RuntimeParameters_get("rt_nr_dtol",rt_nr_dtol)
        call RuntimeParameters_get("rt_nr_maxits",rt_nr_maxits)
      endif
    endif
  endif

  

  ! Boundary conditions
  call RuntimeParameters_get("xl_rad_BC", xl_rad_BC)
  rad_BC(LOW,IAXIS) = xl_rad_BC
  call RuntimeParameters_get("xr_rad_BC", xr_rad_BC)
  rad_BC(HIGH,IAXIS) = xr_rad_BC
!#if NDIM> CONSTANT_ONE
#if K2D
  call RuntimeParameters_get("yl_rad_BC", yl_rad_BC)
  rad_BC(LOW,JAXIS) = yl_rad_BC
  call RuntimeParameters_get("yr_rad_BC", yr_rad_BC)
  rad_BC(HIGH,JAXIS) = yr_rad_BC
#if NDIM>2  
  call RuntimeParameters_get("zl_rad_BC", zl_rad_BC)
  rad_BC(LOW,KAXIS) = zl_rad_BC
  call RuntimeParameters_get("zr_rad_BC", zr_rad_BC)
  rad_BC(HIGH,KAXIS) = zr_rad_BC
#endif
#endif 
  call RuntimeParameters_get("rt_boundary_T", rt_boundary_T)

  call RuntimeParameters_get("fradx_boundary_value", fradx_boundary_value)
#if NDIM>1  
  call RuntimeParameters_get("frady_boundary_value", frady_boundary_value)
#if NDIM>2
  call RuntimeParameters_get("fradz_boundary_value", fradz_boundary_value)
#endif
#endif  

  ! Obtain physical constants
  call PhysicalConstants_get("speed of light",rt_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",rt_boltz)
  ! Radiation constant (often: aR)
  rt_radconst = 4. * rt_boltz/rt_speedlt 

  if (.not. rt_useRadTrans .and. dr_globalMe .eq. MASTER_PE) then
     write(6,*) 'WARNING:  You have included the RadTrans unit but have set '
     write(6,*) '   the runtime parameter useRadTrans to FALSE'
     write(6,*) '   No cooling will occur but RadTrans_init will continue.'
  end if
  
  PETSC_COMM_WORLD=MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  IF(ierr.NE.0) call Driver_abortFlash("Error while initializing Petsc.")

  blockCount = 0
  alloced = .FALSE.

  !Legacy variables
  call Driver_getMype(MESH_COMM,rt_meshMe)
  call RuntimeParameters_get("meshCopyCount", rt_meshCopyCount)  
  call Driver_getMype(MESH_ACROSS_COMM, rt_acrossMe)
  call Driver_getComm(MESH_ACROSS_COMM, rt_acrossComm)
  call Driver_getComm(GLOBAL_COMM,rt_globalComm)

#ifdef RAYTRACE_3DRT
  !Initialise variables for ray-tracer
  call raytrace_init()
#endif

  return

end subroutine RadTrans_init
