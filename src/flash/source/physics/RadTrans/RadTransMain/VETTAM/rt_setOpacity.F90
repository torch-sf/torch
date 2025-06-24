!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_setOpacity
!!
!! NAME
!!  
!!  rt_setOpacity
!!
!!
!! SYNOPSIS
!! 
!!  rt_setOpacity()
!!  
!! DESCRIPTION
!!
!!	Set dust and gas opacities for the VETTAM radiation module
!!
!!***
!!***

#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"
SUBROUTINE rt_setOpacity()

  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe
  use Grid_interface
  use SemenovOpacities, only: getOpacity_planck, getOpacity_rosseland
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Heat_data, ONLY: he_dust_sputter_temp
#if defined(IHA_SPEC) && defined(UEUV_VAR)
  use rt_ionisedata, ONLY : hA, hpA, elecA, h2A, ion_sigmaH, ion_sigmaH2, ion_sigmaH2_15p2_infty, &
                            & ion_sigmaH_13p6_15p2, ion_sigmaH_15p2_infty, useH2Ionize
  use Multispecies_interface, ONLY : Multispecies_getProperty
#elif defined(IONY_MSCALAR)
  use rt_ionisedata, ONLY: ion_sigmaH
  use Eos_data, only: eos_singleSpeciesA
#endif
  implicit none

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, myPE, i, j, k
  real, dimension(:,:,:,:), pointer :: solnData
  real :: rho, temp, kappa, opac_planck, opac_rosseland, iony, opac_gas, nH, nH2, ion_cross_section
  real, save :: mH
  logical, save :: firstCall = .true.

  if(firstCall) then
    if ((dr_globalMe .eq. MASTER_PE)) &
      & print *, 'Dust to Gas ratio adopted = ',dusttoGasRatio
    call PhysicalConstants_get("proton mass",mH)
    call RuntimeParameters_get('sim_A_n', hA)
    call RuntimeParameters_get('sim_A_i', hpA)
#if defined(IHA_SPEC) && defined(UEUV_VAR)
    !Set the mass per atoms here; although rt_ionise will do this, this is called after
    !call Multispecies_getProperty(IHA_SPEC,A,hA)
    !call Multispecies_getProperty(IHP_SPEC,A,hpA)
    !call Multispecies_getProperty(ELEC_SPEC,A,elecA)
    !call Multispecies_getProperty(H2_SPEC,A,h2A)
    !Convert this to mass
    !hA = 1.0
    !hpA = 1.0
    elecA = 0.0
    hA = hA * mH
    hpA = hpA * mH
    elecA = elecA * mH
    h2A = h2A * mH
#endif
    firstCall = .false.
  endif !firstCall
  
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    !Get processor ID
    call Driver_getMype(GLOBAL_COMM,myPE)

    do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

          rho    = solnData(DENS_VAR,i,j,k)
          opac_gas = 0.0 

          !Opacity EUV due to dust
          if(current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
          & current_band .eq. 'EUV_15P2_INFTY') then
            if(dust_opacity_euv .eq. "Fixed") then
              opac_planck = dust_euvopac_value * rho
              opac_rosseland = opac_planck
            else if(dust_opacity_euv .eq. "Scattering") then
              opac_planck = 0.0
              opac_rosseland = dust_euvopac_value*rho
            else
              call Driver_abortFlash("[VETTAM]: Unrecognised opacity type for EUV band.")
            endif

          !Opacity due to atomic hydrogen
#ifdef IONY_MSCALAR
            iony = solnData(IONY_MSCALAR,i,j,k)
            nH = rho/(eos_singleSpeciesA*mH)
            ion_cross_section = ion_sigmaH
            opac_gas = ion_cross_section * nH * iony
          !KROME Case
#elif defined(IHA_SPEC) && defined(UEUV_VAR)
            iony = (solnData(IHA_SPEC,i,j,k)/hA)/(solnData(IHA_SPEC,i,j,k)/hA+ &
              & solnData(IHP_SPEC,i,j,k)/hpA)
            nH = rho*solnData(IHA_SPEC,i,j,k)/hA + rho*solnData(IHP_SPEC,i,j,k)/hpA
            if(current_band .eq. 'EUV_15P2_INFTY') then
              ion_cross_section = ion_sigmaH_15p2_infty
            else if(current_band .eq. 'EUV_13P6_15P2') then
              ion_cross_section = ion_sigmaH_13p6_15p2
            else !i.e. the EUV >13.6 eV band
              ion_cross_section = ion_sigmaH
            endif
            opac_gas = ion_cross_section * nH * iony
            !Contribution to gas opacity from H2
#ifdef H2_SPEC
            nH2 = rho * solnData(H2_SPEC,i,j,k)/h2A
            !EUV band with hnu>15.2
            if(current_band .eq. 'EUV_15P2_INFTY' .and. useH2Ionize) then
              ion_cross_section = ion_sigmaH2_15p2_infty
            !EUV band with hnu>13.6
            else if(current_band .eq. 'EUV' .and. useH2Ionize) then
              ion_cross_section = ion_sigmaH2
            !Any other EUV bands (e.g. 13.6<hnu<15.2) cross section is zero for H2
            else
              ion_cross_section = 0.0
            endif
            opac_gas = opac_gas + ion_cross_section * nH2
#endif
#endif
#ifdef TAUH_VAR
            solnData(TAUH_VAR,i,j,k)  = opac_gas
#endif  
          !Opacity FUV or Lyman Werner Band or Photoelectric heating band
          !Note, although the LW band can also be absorbed by H2, this is taken care of by rt_lw.F90; here only dust
          else if(current_band .eq. 'FUV' .or. current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW' &
               .or. current_band .eq. 'PE') then
            
            if(dust_opacity_fuv .eq. "Fixed") then
              if(current_band .eq. 'FUV') opac_planck = dust_fuvopac_value * rho
              if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') opac_planck = dust_lwopac_value * rho
              if(current_band .eq. 'PE') opac_planck = dust_peopac_value * rho
              opac_rosseland = opac_planck
            else if(dust_opacity_fuv .eq. "Scattering") then 
              opac_planck = 0.0
              if(current_band .eq. 'FUV') opac_rosseland = dust_fuvopac_value * rho
              if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') opac_rosseland = dust_lwopac_value * rho
              if(current_band .eq. 'PE') opac_rosseland = dust_peopac_value * rho  
            else
              call Driver_abortFlash("[VETTAM]: Unrecognised opacity type for FUV band.")
            endif

          else if(current_band .eq. 'IR') then 
            if(dust_opacity_ir .eq. "Fixed") then
              opac_planck = dust_iropac_value * rho
              opac_rosseland = opac_planck
            
            else if(dust_opacity_ir .eq. "Scattering") then
              opac_planck = 0.0
              opac_rosseland = dust_iropac_value*rho
            
            !Semenov opacities
            else if(dust_opacity_ir .eq. "Semenov" .or. dust_opacity_ir .eq. "ScatSemenov" &
                .or. dust_opacity_ir .eq. "SemenovTrad" .or. dust_opacity_ir .eq. "ScatSemenovTrad") then
#ifdef TEMP_VAR
            if(dust_opacity_ir .eq. "Semenov" .or. dust_opacity_ir .eq. "ScatSemenov") then 
              temp = solnData(TEMP_VAR,i,j,k)
            else
              !If using SemenovTrad or ScatSemenov Trad use radiation temperature
              !ERAO is used instead of ERAD as this can be overwritten by other bands
              temp = (solnData(ERAO_VAR,i,j,k)/rt_radconst)**0.25
            endif
#else
            !If Isothermal-like EOS, use trad since temp not available
            temp = (solnData(ERAO_VAR,i,j,k)/rt_radconst)**0.25
#endif
            call getOpacity_planck(temp,rho,kappa)
            opac_planck = 0.0 ! Default to zero
            if(dust_opacity_ir .eq. "Semenov" .or. dust_opacity_ir .eq. "SemenovTrad") opac_planck = kappa * rho
            call getOpacity_rosseland(temp,rho,kappa)
            opac_rosseland = kappa*rho

            else
              call Driver_abortFlash("[VETTAM]: Unrecognised opacity type for IR band.")
            endif
          endif

          !Scale opacities by dust-to-gas ratio, relative to solar, for dust opacities
          opac_planck = opac_planck * dusttoGasRatio
          opac_rosseland = opac_rosseland * dusttoGasRatio
#ifdef TEMP_VAR
          !Destroy dust opacities due to thermal sputtering if gas temperatures > 10^6 K
          !TODO: Implement non-thermal sputtering
          if(solnData(TEMP_VAR,i,j,k) .gt. he_dust_sputter_temp) then
            opac_planck = 1.0e-27 ! set floor to prevent NANs in boundary conditions
            opac_rosseland = 1.0e-27
          endif
#endif

          !Total EUV opacity = dust + hydrogen gas opacity (non-zero only for EUV band)
          solnData(TAUP_VAR,i,j,k) = opac_planck + opac_gas
          solnData(TAUR_VAR,i,j,k) = opac_rosseland + opac_gas

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do


end SUBROUTINE rt_setOpacity
