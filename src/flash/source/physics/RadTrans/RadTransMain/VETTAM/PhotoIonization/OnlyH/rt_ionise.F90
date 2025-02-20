!! Photoionization coupled to VETTAM for modelling EUV radiation. This uses the solution for the EUV photon field
!! and combines it with a rate equation solver for the ionization-recombination balance. For this reaction, we use
!! a method analogous to the implementation in Kim et al. 2018 ART version, plus alternative approaches - see RateEquations.F90
!! rt_ionisehydro.F90 supplies the energy and momentum terms on gas due to photoionization.
!! 
!! Two data containers possible thus far for storing the ionization state: i) a bare-basics version that uses only Hydrogen and stores
!! the neutral fraction in a mass-scalar called IONY_MSCALAR. ii) Coupled to the KROME network, using Multispecies container IHA_SPEC. PP directive
!! checks are implemented here to invoke the appropriate implementations, depending on the version compiled in.
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
  use Eos_interface
  use Hydro_data, ONLY: hy_useHydro
  use RadTrans_interface, ONLY: RadTrans
  use Driver_data, ONLY: dr_globalMe, dr_nstep, dr_simTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RadTrans_data, ONLY: current_band, rt_freqbands, rt_band1, rt_band2, rt_band3, dusttoGasRatio

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
  integer                                    :: b,blockID, nAllParticles, myPE, i, j, k
  integer,dimension(2,MDIM)                  :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer            :: solnData
  logical                                    :: exist, gcmask(NUNK_VARS), converged, gcMaskLogged = .false., ion_band_present
  real                                       :: dt_substep, cumulative_dt, IonizationRate
  logical, save                              :: first_call = .true., log_ion = .true.
  real, save                                 :: mH, pi, kb
  character(len=MAX_STRING_LENGTH)           :: substep_reason, my_band
  character(len=80)                          :: outfile
  integer, parameter                         :: funit = 26
  integer                                    :: nr_its, band_count, nr_substeps, nr_nh, nr_th
  real                                       :: dt_substep_min, dt_substep_max


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
    call RuntimeParameters_get("C_tss",C_tss)
    call RuntimeParameters_get("C_therm",C_therm)
    call RuntimeParameters_get("ion_conv_rtol",ion_conv_rtol)
    call RuntimeParameters_get("ion_conv_atol",ion_conv_atol)
    call RuntimeParameters_get("ion_conv_maxits",ion_conv_maxits)
    call RuntimeParameters_get("hnu",hnu)
    !Thermal energy per ionisation; user sets this only if THERMPION is not a sink property
    !Otherwise the THERMPION of sink(s) is used for energyPerIon; this is set in Cool.F90
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
    call Multispecies_getProperty(ELEC_SPEC,A,elecA)
    call Multispecies_getProperty(H2_SPEC,A,h2A)
    !Convert this to mass
    hA = hA * mH
    hpA = hpA * mH
    elecA = elecA * mH
    h2A = h2A * mH
#endif

    !Set recombination coefficient. Note that if alpha_type is "default" the recombination coefficient would be reset based on temperature.
    !i.e. this value is only ever used if alpha_type is "constant"
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
        CASE DEFAULT
          if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,I0)") 'Band Count, No of bands: ', band_count, rt_freqbands
          call Driver_abortFlash("[VETTAM]: Something is wrong; band_count is greater than 3, should not be happening!")
      END SELECT
      if((my_band .eq. 'EUV' .or. my_band .eq. 'EUV_13P6_15P2' &
      & .or. my_band .eq. 'EUV_15P2_INFTY') .and. log_ion) then
        !Checks to see if multiple EUV bands are present
        if(ion_band_present) then
          multiple_ionbands = .true.
        else
          ion_band_present = .true.
        endif

        !Create substep file
        outfile = trim('Substep_') // trim(my_band) // trim('.log')

        if(dr_globalMe .eq. MASTER_PE) then
          write(*,"(A,A)") 'Writing Ionizations substep info to ',trim(outfile)
          inquire(file=outfile,exist=exist)
          !Create the file only if it does not exist
          if(.not. exist) then
            open(funit, file=trim(outfile), status='new')
            write(funit,'(A10,3X,A5,3X,A6,3X,A13,3X,A10,3X,A10,3X,A9,3X,A9)') '#[00]StepNo', '[01]t','[02]dt',&
              '[03]nSubsteps','[04]dt_min','[05]dt_max','[06]nr_nh','[07]nr_th'
            close(funit)
          endif !exist endif
        endif !MASTER_PE endif
      endif !EUV bands endif
    end do


    first_call = .false.
    !Return to Cool.F90 (where the first call will be made from; so the block of code above is basically all for initialising)
    return
#else
    call Driver_abortFlash("[VETIonize]: Either IONY or IHA_SPEC variables have to be present to use VETIonize.")
#endif
  endif

  if (.not. useEUVIonize)then 
    !Commented here; uncomment if you want to test ray-tracer without ionisation
    !call rt_setOpacity()
    !call RadTrans(blockCount_, blockList_, dt_substep, pass=0)
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

  cumulative_dt = 0
  nr_its = 0
  nr_substeps = 0
  dt_substep_min = huge(dt_substep_min)
  dt_substep_max = tiny(dt_substep_max)
  nr_nh = 0
  nr_th = 0
  dt_substep = dt
  do
    !If there are multiple Ion bands do the below as fluxes and opacities would be overwritten
    !Note: This would work by saving absorbed energy as UEUV_VAR in the higher energy band, and HEUV_VAR in lower
    if(multiple_ionbands .and. current_band .eq. 'EUV_15P2_INFTY') then
      call rt_computeIonH() !Compute the heating + cooling term
      !Get the substep timestep now (important to prevent double counting of radiation pressure etc.)
      call SubstepTimestep()
      if(hy_useHydro) then 
        call rt_ionMomentum(blockCount_,blockList_,dt_substep,time)
#ifdef REIR_VAR
        if(dusttoGasRatio .gt. 0) call rt_ionDustRep(blockCount_,blockList_,dt_substep,time)
#endif
      endif
#ifdef H2_SPEC
      !Ionize H2: this needs to be done before the 13P6 band is solved; see Baczynski et al. 2015
      if(useH2Ionize) call IonizeH2(blockCount_,blockList_,dt_substep,time)
#endif
      !Now do the solution for the 13P6 band
      current_band = 'EUV_13P6_15P2'
      !Set opacity
      call rt_setOpacity()
      if(ion_implicit .or. nr_substeps .eq. 0) then
        !Call RT every substep if ion_implicit; always call the very first time so that one solution is available that can then be reused if ion_implicit is false
        !Get sink source term
        call rt_sinkInject(blockCount_, blockList_, dt_substep) !This has to be called because the current JSTR_VAR is for the 15P2 band
        call Radiation_Step(dt_substep)
      else
        call UpdateUEUV()
      endif
      !Now momentum and energy for this band
      if(hy_useHydro) then 
        call rt_ionMomentum(blockCount_,blockList_,dt_substep,time)
#ifdef REIR_VAR
        if(dusttoGasRatio .gt. 0) call rt_ionDustRep(blockCount_,blockList_,dt_substep,time)
#endif
      endif
      !Reset current band to above; this is important as some calls below need this
      current_band = 'EUV_15P2_INFTY'
    endif

    if(ion_type .eq. 1 .or. ion_type .eq. 2) then
      if(.not. multiple_ionbands) then 
        call rt_computeIonH() !Compute the heating + cooling term
        call SubstepTimestep() !If multiple bands are there substep timestep already computed above
      endif
      !Restrict substep time to match hydro timestep
      !Eq 21 of Kim et al 2017
      dt_substep = MIN(dt_substep,dt-cumulative_dt)
      dt_substep_min = MIN(dt_substep, dt_substep_min)
      dt_substep_max = MAX(dt_substep, dt_substep_max)
      if(myPE .eq. MASTER_PE) then 
        write(*,"(A,ES13.6,3X,A,ES13.6,1X,A,A2,A)") 'Ionization Solver substep dt = ',dt_substep, &
        & 'Hydro timestep = ',dt, '(', substep_reason, ')'
      endif
      if(substep_reason .eq. 'th') nr_th = nr_th + 1
      if(substep_reason .eq. 'nh') nr_nh = nr_nh + 1
    ! No substepping for implicit methods
    else
      dt_substep = dt
    endif

    !Ionize H2: if multiple ionbands are present, H2 ionization is already performed above
#ifdef H2_SPEC
    if(.not. multiple_ionbands) call IonizeH2(blockCount_,blockList_,dt_substep,time)
#endif

    !Call rate equation solver for hydrogen
    call IonizeRateEquation(blockCount_,blockList_,dt_substep,time)

    ! Heating/cooling + Momentum 
    if(hy_useHydro) then
      call rt_ionHeatCool(blockCount_,blockList_,dt_substep,time)
      !Radiation Pressure. If there are multiple EUV bands, this is already added before. 
      if(.not. multiple_ionbands) call rt_ionMomentum(blockCount_,blockList_,dt_substep,time)
      !Same as above for reprocessed energy to the IR by dust
#ifdef REIR_VAR
      if(.not. multiple_ionbands .and. dusttoGasRatio .gt. 0) call rt_ionDustRep(blockCount_,blockList_,dt_substep,time)
#endif
    endif

    !Explicit methods
    if(ion_type .eq. 1 .or. ion_type .eq. 2) then
      cumulative_dt = cumulative_dt + dt_substep
      if(cumulative_dt .ge. dt) exit
    !Implicit methods
    else
      nr_its = nr_its + 1 !Iteration counter
      call CheckforConvergence(converged)
      if(converged) then
        if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,A)") 'Converged in: ', nr_its, ' iterations'
        exit
      else if(nr_its .gt. ion_conv_maxits) then 
        call Driver_abortFlash("[Ionize]: Cycles between RadTrans and Ionisation did not converge.")
      else
        call UpdateOldIon()
        continue
      endif
    endif

    !Get new LyC opacity in the domain (which changes rapidly every substep)
    call rt_setOpacity()
    !If implicit evolution of ion+rad switched on, get RT solution
    if(ion_implicit) then
      !Update the radiation energy solution (ERAD_VAR)
      if(multiple_ionbands) call rt_sinkInject(blockCount_, blockList_, dt_substep) !Required since JSTR now has the 13P6 output
      call Radiation_step(dt_substep)
    else
      !Update the absorbed energy in cell as opacity changed now; ERAD_VAR is fixed in this case unlike above
      call UpdateUEUV()
    endif
    nr_substeps = nr_substeps + 1
  end do

  !Logging
  if(ion_type .eq. 1 .and. log_ion .and. dr_globalMe .eq. MASTER_PE) then
    outfile = trim('Substep_') // trim(current_band) // trim('.log')
    open(funit, file=trim(outfile), position='APPEND')
    write(funit,'(I0,1X,ES16.9,1X,ES16.9,1X,I0,1X,ES16.9,1X,ES16.9,1X,I0,1X,I0)') dr_nstep, &
        dt, dr_simTime, nr_substeps,dt_substep_min,dt_substep_max,nr_nh, nr_th
    close(funit)
  endif

  !Fill guard cells. TODO: Do I really need this? 
  gcmask(:) = .false.
  !Only update IONY, TAUP, and TGAS
#ifdef IHA_SPEC
  gcmask(IHA_SPEC) = .true.
#ifdef H2_SPEC
  if(useH2Ionize) gcmask(H2_SPEC) = .true.
#endif
#elif defined(IONY_MSCALAR)
  gcmask(IONY_MSCALAR) = .true.
#endif
  gcmask(UEUV_VAR) = .true.
#ifdef TEMP_VAR
  gcmask(TEMP_VAR) = .true.
#elif defined(TGAS_VAR) 
  gcmask(TGAS_VAR) = .true.
#endif
  !Fill Guard cells
  call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS,mask=gcmask,makeMaskConsistent = .true.,&
  selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)

  call Timers_stop("Ionize")

CONTAINS 

SUBROUTINE SubstepTimestep
  use RadTrans_data
  use Grid_interface
  use Eos_interface
  use Driver_data, ONLY: dr_globalMe
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use rt_ionisemodule, ONLY: get_recombination_coefficient
  implicit none

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k, error, substep_PE
  real, dimension(:,:,:,:), pointer :: solnData
  real :: rho, nH, nH0, nHplus, neutfrac, erad, ne, dt_local, dt_temp, nHtotal, iony, alpha_rec_cell, temp
  real :: RecombinationRate, IonizationRate, HeatCoolRate, Eint, del(MDIM)
  real :: IonizationRate_H, IonizationRate_H2, dnH_by_dt, dnH2_by_dt, dEint_by_dt, IonizationRate_13P6
  integer, parameter :: communication_count = 1, tag = 1
  real :: frac_15p2_infty_H2, frac_15p2_infty_H, nH2
  real, dimension(2) :: send_data, receive_data
  character(len=MAX_STRING_LENGTH) :: local_substep_reason
  integer, dimension(MPI_STATUS_SIZE) :: status

  dt_local = huge(dt_substep)
  local_substep_reason = 'No'

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID, del)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          rho      = solnData(DENS_VAR,i,j,k)
          ! Neutral fraction
          !Neutral fraction = nH0/(nH0 + nHplus)
#ifdef IHA_SPEC
          ! Total hydrogen nuclei number density = nHplus + nH0
          nH0 = rho*solnData(IHA_SPEC,i,j,k)/hA
          nHplus = rho*solnData(IHP_SPEC,i,j,k)/hpA
          nH     = nH0 + nHplus
          ! ne = mass fraction * density/m_e
          ne       = solnData(ELEC_SPEC,i,j,k) * rho/(elecA)
#elif defined(IONY_MSCALAR)
          nH = rho/(eos_singleSpeciesA*mH)
          !Neutral fraction
          iony     = solnData(IONY_MSCALAR,i,j,k)
          !Neutral H number density
          nH0      = iony * nH
          !Ionized H and Electron number density
          nHplus       = nH*(1-iony)
          ne = nHplus
#endif
          !Local ionising radiation
#ifdef UEUV_VAR
          !Volumetric ionization rate
          IonizationRate    = solnData(UEUV_VAR,i,j,k)/(hnu)
#else
          !NOTE: If ray-tracer not included then user can define something here
          IonizationRate     = 0.0
#endif

#if defined(H2_SPEC)
          nH2 = rho * solnData(H2_SPEC,i,j,k)/h2A
          !If H2 ionization is switched on, compute the fraction of photons absorbed in the hnu>15.2 eV band by H and H2 respectively
          if(current_band .eq. 'EUV_15P2_INFTY' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2_15p2_infty*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH_15p2_infty*nH0/(solnData(TAUH_VAR,i,j,k))
          !If there is only one EUV band >13.6 eV
          else if(current_band .eq. 'EUV' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH*nH0/(solnData(TAUH_VAR,i,j,k))
          else if(current_band .eq. 'EUV_13P6_15P2' .or. .not. useH2Ionize) then
            frac_15p2_infty_H2 = 0.0
            frac_15p2_infty_H = 1.0
          endif

          IonizationRate_H2 = IonizationRate * frac_15p2_infty_H2
          IonizationRate_H = IonizationRate * frac_15p2_infty_H
            
#else       
          !Pure Hydrogen case
          IonizationRate_H2 = 0.0
          IonizationRate_H = IonizationRate
#endif

#ifdef TEMP_VAR
          temp = solnData(TEMP_VAR,i,j,k)
#elif defined(TGAS_VAR)
          temp = solnData(TGAS_VAR,i,j,k)
#endif
          if(alpha_type .eq. 'default') then
            alpha_rec_cell = get_recombination_coefficient(temp,ion_ots)
          else if(alpha_type .eq. 'constant') then
            alpha_rec_cell = alpha_rec_constant
          else
            print *, alpha_type
            call Driver_abortFlash("[rt_ionise]:alpha_type should be 'constant' or 'default'; check!")
          endif

          RecombinationRate = alpha_rec_cell * ne * nHplus
          dnH_by_dt = IonizationRate_H - RecombinationRate
#ifdef H2_SPEC
          dnH2_by_dt = IonizationRate_H2
          !Contribution to H due to H2 ionization: each H2 gives 1 H
          if(useH2Ionize) dnH_by_dt = dnH_by_dt + 2*dnH2_by_dt
#endif
          !Now consider case where there are multiple bands of ionisation for H
          if(multiple_ionbands) then
            IonizationRate_13P6 = solnData(HEUV_VAR,i,j,k)/hnu
            !Add ionisations from the 15.2 band
            dnH_by_dt = dnH_by_dt + IonizationRate_13P6
          endif

          dt_temp = C_tss * nH/(abs(dnH_by_dt))
          dt_local = min(dt_temp,dt_local)
          if(dt_local/dt .lt. 1.e-50) then
            print *, 'dt_nH too small. dt_local, nH0, &
              & nHplus, dnH_by_dt, dnH2_by_dt, IonizationRate_H, RecombinationRate = ', dt_local, &
              & nH0, nHplus, dnH_by_dt, dnH2_by_dt, IonizationRate_H, RecombinationRate
            call Driver_abortFlash("dt_local is too small.")
          endif
          if(dt_temp .eq. dt_local .and. dt_temp .lt. dt) local_substep_reason = 'nh'

          !Heating/cooling condition
#ifdef EINT_VAR
          HeatCoolRate = solnData(PHHE_VAR,i,j,k)
          Eint = solnData(EINT_VAR,i,j,k)
          dEint_by_dt = HeatCoolRate
          dt_temp = C_therm*Eint/abs(dEint_by_dt)
          if(hy_useHydro) dt_local = min(dt_temp,dt_local)
          if(dt_local/dt .lt. 1.e-50) then
            print *, 'dt_therm too small. dt_local, nH, nH0, ne, Eint, &
              & HeatCoolRate,  = ', dt_local, nH, nH0, ne, Eint, HeatCoolRate
            call Driver_abortFlash("dt_local is too small.")
          endif
          if(dt_temp .eq. dt_local .and. dt_temp .lt. dt) local_substep_reason = 'th'
#endif
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

  !Now communicate and get min timestep and processor containing this value
  send_data(1) = dt_local
  send_data(2) = real(myPE)
  call MPI_AllReduce(send_data,receive_data,communication_count,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,error)
  dt_substep = receive_data(1)
  substep_PE = int(receive_data(2))
  !Now communicate the reason for the substep to Master for printing
  if(substep_PE .eq. myPE) call MPI_SEND(local_substep_reason, len(substep_reason), MPI_CHAR, MASTER_PE, tag, MPI_COMM_WORLD, error)
  if(myPE .eq. MASTER_PE) call MPI_Recv(substep_reason, len(substep_reason), MPI_CHAR, substep_PE, tag, MPI_COMM_WORLD, status, error)
  

end SUBROUTINE SubstepTimestep

SUBROUTINE rt_computeIonH

  use rt_ionisemodule
  implicit none 

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k, error
  real, dimension(:,:,:,:), pointer :: solnData
  real :: neutfrac, tgas, iony, dE_dt, rho, nH0, nHplus, ne, nions_pve, nH
  real :: IonizationRate, Ion_Heating_Rate, Recombination_Cooling_Rate, FF_Cooling_Rate
  real :: IonizationRate_H, IonizationRate_H2, nH2, frac_15p2_infty_H, frac_15p2_infty_H2, IonizationRate_13P6

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
          do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

            rho = solnData(DENS_VAR,i,j,k)
#ifdef IHA_SPEC
            iony = (solnData(IHA_SPEC,i,j,k)/hA)/(solnData(IHA_SPEC,i,j,k)/hA+ &
            & solnData(IHP_SPEC,i,j,k)/hpA)
            nH0 = rho*solnData(IHA_SPEC,i,j,k)/hA
            nHplus = rho*solnData(IHP_SPEC,i,j,k)/hpA
            ne = solnData(ELEC_SPEC,i,j,k) * rho/(elecA) 
            !Assuming only Hplus as the positive ion
            !TODO: If using cooling here when KROME is compiled, consider adding other positive ion species here
            nions_pve = nHplus
#elif defined(IONY_MSCALAR)
            iony = solnData(IONY_MSCALAR,i,j,k)
            nH = rho/(eos_singleSpeciesA*mH)
            nH0      = iony * nH
            nHplus       = nH*(1-iony)
            ne = nHplus
            !Unlike KROME, here there is indeed only H as a positive ion
            nions_pve = nHplus
#endif

#ifdef EINT_VAR
            tgas = solnData(TEMP_VAR,i,j,k)
            !1. Heating due to ionisations; units: erg cm^-3)
            !Volumetric ionization rate: No of ionisations per unit volume per unit time
            IonizationRate    = solnData(UEUV_VAR,i,j,k)/(hnu)
#if defined(H2_SPEC)
            nH2 = rho * solnData(H2_SPEC,i,j,k)/h2A
            !If H2 ionization is switched on, compute the fraction of photons absorbed in the hnu>15.2 eV band by H and H2 respectively
            if(current_band .eq. 'EUV_15P2_INFTY' .and. useH2Ionize) then
              frac_15p2_infty_H2 = ion_sigmaH2_15p2_infty*nH2/(solnData(TAUH_VAR,i,j,k))
              frac_15p2_infty_H = ion_sigmaH_15p2_infty*nH0/(solnData(TAUH_VAR,i,j,k))
            !If there is only one EUV band >13.6 eV
            else if(current_band .eq. 'EUV' .and. useH2Ionize) then
              frac_15p2_infty_H2 = ion_sigmaH2*nH2/(solnData(TAUH_VAR,i,j,k))
              frac_15p2_infty_H = ion_sigmaH*nH0/(solnData(TAUH_VAR,i,j,k))
            else if(current_band .eq. 'EUV_13P6_15P2' .or. .not. useH2Ionize) then
              frac_15p2_infty_H2 = 0.0
              frac_15p2_infty_H = 1.0
            endif
            IonizationRate_H2 = IonizationRate * frac_15p2_infty_H2
            IonizationRate_H = IonizationRate * frac_15p2_infty_H
            
#else       
            !Pure Hydrogen case
            IonizationRate_H2 = 0.0
            IonizationRate_H = IonizationRate
#endif

            !Heating due to ionizations : both due to ionizing atomic and molecular hydrogen
            Ion_Heating_Rate = IonizationRate_H * energyPerIonH + IonizationRate_H2 * energyPerIonH2

            !Now consider case where there are multiple bands of ionisation for H
            if(multiple_ionbands) then
              IonizationRate_13P6 = solnData(HEUV_VAR,i,j,k)/hnu
              !Add heating; energyPerIonH_13p6_15p2 would be set in rt_sinkInject
              Ion_Heating_Rate = Ion_Heating_Rate + IonizationRate_13P6 * energyPerIonH_13p6_15p2
            endif

            !2. Free-free cooling rate= Gamma_ff * n_e * n_pions (units: erg cm^-3 s^-1). Gamma_ff is in rt_ionisemodule 
            FF_Cooling_Rate = ff_cool_coeff(tgas) * ne * nions_pve

            !3. Recombination cooling rate= Gamma_rec * n_e * n_Hp (units: erg cm^-3 s^-1). This is recombination of ionised hydrogen
            Recombination_Cooling_Rate = recombination_cool_coeff(tgas) * ne * nHplus

            !Other terms can be added here or called here.

            ! Net heating - cooling rate
            dE_dt = Ion_Heating_Rate - FF_Cooling_Rate - Recombination_Cooling_Rate
            dE_dt = dE_dt/solnData(DENS_VAR,i,j,k) !Consistent with EINT of FLASH
            solnData(PHHE_VAR,i,j,k) = dE_dt
#endif
            end do
        end do
      end do
      call Grid_releaseBlkPtr(blockID,solnData)
    end do
END SUBROUTINE rt_computeIonH

SUBROUTINE CheckforConvergence(converged)
  use Grid_interface
  implicit none
#include "Flash_mpi.h"
  logical, intent(OUT)   :: converged
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k, error
  integer, parameter :: communication_count = 1
  real, dimension(:,:,:,:), pointer :: solnData
  real :: iony, iony_prev, rel_change,rel_change_local, rel_change_global
  real :: abs_change,abs_change_local, abs_change_global

  rel_change = 0.0
  rel_change_local = 0.0
  rel_change_global = 0.0
  abs_change = 0.0
  abs_change_local = 0.0
  abs_change_global = 0.0
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
#ifdef IHA_SPEC
          iony = (solnData(IHA_SPEC,i,j,k)/hA)/(solnData(IHA_SPEC,i,j,k)/hA+ &
            & solnData(IHP_SPEC,i,j,k)/hpA)
#elif defined(IONY_MSCALAR)
          iony = solnData(IONY_MSCALAR,i,j,k)
#endif
          iony_prev = solnData(IONP_VAR,i,j,k)
          abs_change = ABS(iony- iony_prev)
          rel_change = ABS(iony- iony_prev)/(iony_prev)
          rel_change_local = MAX(rel_change,rel_change_local)
          abs_change_local = MAX(abs_change,abs_change_local)

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do
  ! This will give us max change in the entire domain
  call MPI_AllReduce(rel_change_local,rel_change_global,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
  call MPI_AllReduce(abs_change_local,abs_change_global,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)

  if(rel_change_global .lt. ion_conv_rtol .or. abs_change_global .lt. ion_conv_atol) then 
    converged = .true.
  else
    converged = .false.
  endif

END SUBROUTINE CheckforConvergence

SUBROUTINE UpdateOldIon()
  use Grid_interface
  integer :: b, blockID
  real, dimension(:,:,:,:), pointer :: solnData

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkPtr(blockID,solnData)
#ifdef IHA_SPEC
    solnData(IONP_VAR,:,:,:) = solnData(IHA_SPEC,:,:,:)
#elif defined(IONY_MSCALAR)
    solnData(IONP_VAR,:,:,:) = solnData(IONY_MSCALAR,:,:,:)
#endif
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

END SUBROUTINE UpdateOldIon

SUBROUTINE UpdateUEUV()
  use Grid_interface
  use RadTrans_data, ONLY: rt_speedlt
  integer :: b, blockID, i, j, k, absorbed_var
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData
  real :: ueuv_old

  if(current_band .eq. 'EUV_13P6_15P2' .and. multiple_ionbands) then
    absorbed_var = HEUV_VAR
  else
    absorbed_var = UEUV_VAR
  endif

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          if(multiple_ionbands) then
            !First figure out which band
            if(current_band .eq. rt_band1) then
#ifdef B1ER_VAR
              solnData(absorbed_var,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(B1ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B1ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                  & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band2) then
#ifdef B2ER_VAR
              solnData(absorbed_var,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(B2ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B2ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                  & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band3) then
#ifdef B3ER_VAR
              solnData(absorbed_var,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(B3ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B3ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                  & has to be present in the simulation config.")
#endif
            endif
          else
            solnData(absorbed_var,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(ERAD_VAR,i,j,k)
          endif
        end do 
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

END SUBROUTINE UpdateUEUV

end SUBROUTINE rt_ionise
