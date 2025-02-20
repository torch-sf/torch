!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/Photoionization/ThermoChemistry/rt_ionise
!
!!  rt_ionise
!!
!!
!! SYNOPSIS
!! 
!!  rt_ionise(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt,
!!          real(IN) :: time)
!!  
!! DESCRIPTION
!! Photoionization coupled to VETTAM for modelling EUV radiation. This subroutine calculates the ionization rates of H and H2.
!!   and stores it in the solution variables PHIO_VAR and IH2_VAR respectively. 
!! AUTHOR
!!  Shyam Harimohan Menon (2022-2023)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!  time : the current time
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
SUBROUTINE rt_ionise(blockCount_,blockList_,dt,time)
  use rt_ionisedata
  use Grid_interface
  use Grid_data, ONLY : gr_smallx
  use Eos_interface
  use Hydro_data, ONLY: hy_useHydro
  use RadTrans_interface, ONLY: RadTrans
  use Driver_data, ONLY: dr_globalMe, dr_nstep, dr_simTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RadTrans_data, ONLY: current_band, rt_freqbands, rt_band1, rt_band2, rt_band3, rt_band4, rt_band5, dusttoGasRatio

#ifdef SINK_PART_TYPE
  use Particles_sinkData, ONLY : useSinkParticles, localnpf
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
#endif

  use Logfile_interface, ONLY : Logfile_stamp
#ifdef IONY_MSCALAR
  use Eos_data, ONLY : eos_singleSpeciesA
#endif
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
#ifdef IHA_SPEC
  use Multispecies_interface, ONLY : Multispecies_getProperty
#endif
  
  
  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time
  integer                                    :: b,blockID, nAllParticles, myPE, i, j, k, band_count
  integer,dimension(2,MDIM)                  :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer            :: solnData
  logical                                    :: exist, gcmask(NUNK_VARS), converged, gcMaskLogged = .false., ion_band_present
  real                                       :: IonizationRate
  logical, save                              :: first_call = .true., log_ion = .true.
  real, save                                 :: mH, pi, kb
  character(len=MAX_STRING_LENGTH)           :: substep_reason, my_band


  if(first_call) then
#if defined(IONY_MSCALAR) || defined(IHA_SPEC)
    !Read runtime params
    call RuntimeParameters_get("useEUVIonize",useEUVIonize)
    call RuntimeParameters_get("useH2Ionize",useH2Ionize)
    call RuntimeParameters_get("ion_OTS",ion_OTS)
    call RuntimeParameters_get("ion_implicit",ion_implicit)
    call RuntimeParameters_get("ion_sigmaH",ion_sigmaH)
    call RuntimeParameters_get("ion_sigmaH_13p6_15p2",ion_sigmaH_13p6_15p2)
    call RuntimeParameters_get("ion_sigmaH_15p2_infty",ion_sigmaH_15p2_infty)
    call RuntimeParameters_get("ion_sigmaH2",ion_sigmaH2)
    call RuntimeParameters_get("ion_sigmaH2_15p2_infty",ion_sigmaH2_15p2_infty)
    call RuntimeParameters_get("ion_type",ion_type)
    call RuntimeParameters_get("alpha_type",alpha_type)
    call RuntimeParameters_get("alpha_A",alpha_A)
    call RuntimeParameters_get("alpha_B",alpha_B)
    call RuntimeParameters_get("hnu",hnu)
    !Thermal energy per ionisation; user sets this only if THERMPION is not a sink property
    !Otherwise the THERMPION of sink(s) is used for energyPerIon; this is set in RadTrans.F90
#if defined(EION_PART_PROP) || defined(THERMPION_13P6_15P2_PART_PROP) || defined(THERMPION_15P2_INFTY_PART_PROP)
    if ((dr_globalMe == MASTER_PE)) print *, 'energyPerIonH would be set by the sinks'
#else
    !TODO: Only one energyPerIonH for the default case; may need to add one more if this is used outside PopIII
    call RuntimeParameters_get("energyPerIonH",energyPerIonH)
#endif
    !Same for H2
#if defined(THERMPION_15P2_INFTY_H2_PART_PROP) || defined(THERMPION_H2_PART_PROP)
    if ((dr_globalMe == MASTER_PE)) print *, 'energyPerIonH2 would be set by the sinks'
#else
    call RuntimeParameters_get("energyPerIonH2",energyPerIonH2)
#endif

    call PhysicalConstants_get("proton mass",mH)
    call PhysicalConstants_get("pi",pi)
    call PhysicalConstants_get("Boltzmann",kb)
#ifdef IHA_SPEC
    !Save mass per particle in units of mH for species
    call Multispecies_getProperty(IHA_SPEC,A,hA)
    call Multispecies_getProperty(IHP_SPEC,A,hpA)
#ifdef H2_SPEC
    call Multispecies_getProperty(H2_SPEC,A,h2A)
#endif
    !Convert this to mass
    hA = hA * mH
    hpA = hpA * mH
    h2A = h2A * mH
#endif

    !Set recombination coefficient. Note that if alpha_type is "default" the recombination coefficient would be reset based on temperature.
    !i.e. this value is only ever used if alpha_type is "constant"
    !These are not used here, but are still used in the radiation solve for the inclusion of the diffuse recombination radiation (if on)
    if(ion_OTS) then
      alpha_rec_constant = alpha_B !Case-B
      alpha_ground_constant = 0.0 !No net recombinations to ground state
    else
      alpha_rec_constant = alpha_A !Case-A 
      alpha_ground_constant = alpha_A - alpha_B !alpha_A-alpha_B
    endif

    ion_band_present=.false.
    multiple_ionbands = .false.
    !Log file
    do band_count = 1, rt_freqbands
      !Choose current frequency band
      SELECT CASE(band_count)
        CASE(1)
          my_band = rt_band1
        CASE(2)
          my_band = rt_band2
        CASE(3)
          my_band = rt_band3
        CASE(4)
          my_band = rt_band4
        CASE(5)
          my_band = rt_band5
        CASE DEFAULT
          if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,I0)") 'Band Count, No of bands: ', band_count, rt_freqbands
          call Driver_abortFlash("[VETTAM]: Something is wrong; band_count is greater than 5, should not be happening!")
      END SELECT
      if((my_band .eq. 'EUV' .or. my_band .eq. 'EUV_13P6_15P2' &
      & .or. my_band .eq. 'EUV_15P2_INFTY') .and. log_ion) then
        !Checks to see if multiple EUV bands are present
        if(ion_band_present) then
          multiple_ionbands = .true.
        else
          ion_band_present = .true.
        endif
      endif !EUV bands endif
    end do


    first_call = .false.
    !Return to RadTrans.F90 (where the first call will be made from; so the block of code above is basically all for initialising)
    return
#else
    call Driver_abortFlash("[VETIonize]: Either IONY or IHA_SPEC variables have to be present to use VETIonize.")
#endif
  endif

  if (.not. useEUVIonize)then 
    return
  endif

  call Timers_start("Ionize")
  !Get processor ID
  call Driver_getMype(GLOBAL_COMM,myPE)

   !Set old time ionisation state and reset heating rate to zero
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
#ifdef IHA_SPEC
          solnData(IONO_VAR,i,j,k) = (solnData(IHA_SPEC,i,j,k)/hA)/(solnData(IHA_SPEC,i,j,k)/hA+ &
            & solnData(IHP_SPEC,i,j,k)/hpA)
#elif defined(IONY_MSCALAR)
          solnData(IONO_VAR,i,j,k) = solnData(IONY_MSCALAR,i,j,k)
#endif
          solnData(IONP_VAR,i,j,k) = solnData(IONO_VAR,i,j,k)
          !Reset ionisation heating rate to zero
          solnData(PHHE_VAR,i,j,k) = 0.0
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

  !Obtain solutions and impart momentum/dust heating
  if(multiple_ionbands .and. current_band .eq. 'EUV_15P2_INFTY') then
    if(hy_useHydro) then 
      call rt_ionMomentum(blockCount_,blockList_,dt,time)
#ifdef REIR_VAR
      if(dusttoGasRatio .gt. 0) call rt_ionDustRep(blockCount_,blockList_,dt,time)
#endif
    endif
    !Now do the solution for the 13P6 band
    current_band = 'EUV_13P6_15P2'
    !Set opacity
    call rt_setOpacity()
    call rt_sinkInject(blockCount_, blockList_, dt) !This has to be called because the current JSTR_VAR is for the 15P2 band
    call Radiation_Step(dt)
  endif

  !Now momentum and energy: note if there are multiple bands, here it is the 13P6 band for which this done
  !If only a single band it is done here as well. 
  if(hy_useHydro) then 
    call rt_ionMomentum(blockCount_,blockList_,dt,time)
#ifdef REIR_VAR
    if(dusttoGasRatio .gt. 0) call rt_ionDustRep(blockCount_,blockList_,dt,time)
#endif
  endif

  !Set the ionization rates in units of per second to pass to KROME
  call SetIonRates()

  call Timers_stop("Ionize")

CONTAINS 

SUBROUTINE SetIonRates
  use Grid_interface
  integer :: b, blockID
  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: rho, iony, nH, nH0, nHplus, ne, IonizationRate, IonizationRate_H, IonizationRate_H2, nH2
  real :: frac_15p2_infty_H, frac_15p2_infty_H2
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
          
          rho = solnData(DENS_VAR,i,j,k)
          ! Neutral fraction
          !Neutral fraction = nH0/(nH0 + nHplus)
          ! Total hydrogen nuclei number density = nHplus + nH0
          nH0 = rho*max(solnData(IHA_SPEC,i,j,k),gr_smallx)/hA
          nHplus = rho*max(solnData(IHP_SPEC,i,j,k),gr_smallx)/hpA
          !IonizationRate from Radiation solve
          IonizationRate    = solnData(UEUV_VAR,i,j,k)/(hnu)
#if defined(H2_SPEC)
          nH2 = rho * max(solnData(H2_SPEC,i,j,k),gr_smallx)/h2A
#endif
          !If H2 ionization is switched on, compute the fraction of photons absorbed in the hnu>15.2 eV band by H and H2 respectively
          
          !See if the band for 15P2_INFTY exists; this considers both multiple band cases, or if only the higher energy EUV band present
          if(useH2Ionize .and. (multiple_ionbands .or. current_band .eq. 'EUV_15P2_INFTY')) then
            frac_15p2_infty_H2 = ion_sigmaH2_15p2_infty*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH_15p2_infty*nH0/(solnData(TAUH_VAR,i,j,k))
          !Consider case where only one single EUV band; same as above, but use appropriately band-averaged cross sections
          else if(useH2Ionize .and. current_band .eq. "EUV") then
            frac_15p2_infty_H2 = ion_sigmaH2*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH*nH0/(solnData(TAUH_VAR,i,j,k))
          !This considers the case where no H2 ionization, or if only the 13P6_15P2 band is present; both shouldn't ionize H2
          else
            frac_15p2_infty_H2 = 0.0
            frac_15p2_infty_H = 1.0
          endif
          !Ionization rate from the higher energy ionization band
          IonizationRate_H2 = solnData(UEUV_VAR,i,j,k)/(hnu) * frac_15p2_infty_H2
          IonizationRate_H = solnData(UEUV_VAR,i,j,k)/(hnu) * frac_15p2_infty_H
          !Ionization rate from the lower energy ionization band (if it exists)
          !HEUV will be zero if this band doesnt exist
          IonizationRate_H = IonizationRate_H + solnData(HEUV_VAR,i,j,k)/(hnu)
#ifdef PHIO_VAR
          !Ionization rate of H in units of per second -> this goes to KROME
          solnData(PHIO_VAR,i,j,k) = IonizationRate_H/nH0
#endif
#ifdef IH2_VAR
          !Ionization rate of H2 in units of per second -> this goes to KROME
          solnData(IH2_VAR,i,j,k) = IonizationRate_H2/nH2
#endif

#ifdef PHHE_VAR
          !If a cell-by-cell heating rate per ionization exists use that; else use the global value
#ifdef EIH_VAR
          energyPerIonH = solnData(EIH_VAR,i,j,k)
#endif

#ifdef EIH2_VAR
          energyPerIonH2 = solnData(EIH2_VAR,i,j,k)
#endif    
          !Heating rate due to ionizations
          solnData(PHHE_VAR,i,j,k) = IonizationRate_H * energyPerIonH + IonizationRate_H2 * energyPerIonH2
#endif
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

end SUBROUTINE SetIonRates

end SUBROUTINE rt_ionise
