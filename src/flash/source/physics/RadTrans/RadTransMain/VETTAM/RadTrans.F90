!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/RadTrans
!!
!! NAME
!!  
!!  RadTrans
!!
!!
!! SYNOPSIS
!! 
!!  RadTrans(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt,
!!          real(IN) :: time)
!!  
!! DESCRIPTION
!!
!!	Radiation Hydrodynamic update using the VETTAM RHD Module. Please cite Menon et al. 2022; MNRAS 512, 401 if using this module in your work 
!!
!! AUTHOR
!! 	Shyam Harimohan Menon (2020-2024). Postdoc at Rutgers University & Flatiron Institute
!!  Email: shyam.menon@rutgers.edu ; smenon@flatironinstitute.org
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

#include "Flash.h"
#include "constants.h"
#include "petsc/finclude/petscksp.h"
SUBROUTINE RadTrans(blockCount_,blockList_,dt,time)

  use RadTrans_data
  use Grid_interface
  use Hydro_data, ONLY: hy_useHydro
  use Eos_interface
  use Driver_data, ONLY: dr_globalMe
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
#ifdef RAYTRACE_3DRT  
  use raytrace_data, ONLY : rt_useRayTrace
#endif  

#ifdef SINK_PART_TYPE
  use Particles_sinkData, ONLY : useSinkParticles, localnpf
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
#endif

#ifdef UEUV_VAR
  use rt_ionisedata, ONLY: multiple_ionbands, useEUVIonize
#endif
#ifdef FSHM_VAR
  use rt_lwdata, ONLY: lwdiss_type
#endif

  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampVarMask
  
  
  implicit none

!#include "Flash_mpi.h"

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time
  integer                                    :: b,blockID, nAllParticles, band_count
  integer,dimension(2,MDIM)                  :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer            :: solnData
  logical                                    :: gcmask(NUNK_VARS), gcMaskLogged = .false., resetlog, dust_addsrcterms
  logical, save                              :: doRTDiffuse, directSinkContribution, sink_useRayTrace, firstCall=.true.
  real, save                                 :: tsinceLastCall
  integer, parameter                         :: max_freqbands = 5
  
  if(firstCall) then

    if(rt_freqbands .gt. max_freqbands) then
      call Driver_abortFlash("Only a maximum of five frequency bands implemented.")
    endif

  
#ifdef RAYTRACE_3DRT
    !Determine if Diffuse RT needs to be done
    if(rt_etens .eq. 'vet') then 
      doRTDiffuse = .true.
    else
      doRTDiffuse = .false.
    endif
#else
    !Safety check
    if(rt_etens .eq. 'vet') then 
      call Driver_abortFlash("rt_etens = vet requires the Ray Tracer to be compiled in.")
    endif

#ifdef SINK_PART_TYPE
    !Check if direct sink contribution handled by Mom method or ray-trace (former recommended)
    call RuntimeParameters_get("sink_useRayTrace", sink_useRayTrace)
    if(.not. useSinkParticles .and. sink_useRayTrace) &
      & call Driver_abortFlash("sink_useRayTrace is set to true, but useSinkParticles is set to false.")
#endif
#endif

#if defined(IONY_MSCALAR) || defined(IHA_SPEC)
    !Make this call to initialise some variables in ionise which will be used below
    call rt_ionise(blockCount_,blockList_,dt,time)
#endif

#ifdef FSHM_VAR
!Make this call to initialise some variables in LW dissociation which will be used below
    call rt_lw(blockCount_,blockList_,dt)
#endif

    !Make time since last call a huge number to ensure that RadTrans is called in the first step
    tsinceLastCall = 1.e99
    firstCall = .false.
  endif !firstCall

  
  !Check if RadTrans is to be used
  if (.not. rt_useRadTrans &
#ifdef RAYTRACE_3DRT
.and. .not. sink_useRayTrace & 
#endif
  )return

  !Check if enough time has passed since last call to RadTrans
  tsinceLastCall = tsinceLastCall + dt
  if(tsinceLastCall .lt. rt_useeverydt) then
    if(dr_globalMe .eq. MASTER_PE) &
      & print *, "Skipping VETTAM since tsinceLastCall", tsinceLastCall, "is less than rt_useeverydt", rt_useeverydt 
    return
  endif

#ifdef SINK_PART_TYPE
  if(localnpf .eq. 0 .and. .not. noSink_VETTAM) then
    if(dr_globalMe .eq. MASTER_PE) &
    & print *, "No sinks in domain and noSink_VETTAM is false. Skipping VETTAM" 
    return
  endif
#endif

  call Timers_start("RadTrans")

  if(rt_useRadTrans) then
    !Setup matrix including solver and PC. This is fixed for all bands.
    call Timers_start("vet_petscSetup")
    call Petsc_init(blockCount_,blockList_)
    call Timers_stop("vet_petscSetup")

    !Save old time solution for the IR band. UV is time-independent, so this exercise is not required for the UV bands.
    !TODO: Note this is done before looping through the bands since ERAD and MOH vars are recycled across the bands. The fact that IR is done last is crucial for this.
    do b = 1, blockCount_
          
      blockID = blockList_(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      !Set old time quantities for IR band. Note: This assumes ERAD and MOH contains the IR solution
      solnData(ERAO_VAR,:,:,:) = solnData(ERAD_VAR,:,:,:)
      solnData(MHXO_VAR,:,:,:) = solnData(MOHX_VAR,:,:,:)
#if NDIM>1
      solnData(MHYO_VAR,:,:,:) = solnData(MOHY_VAR,:,:,:)
#if NDIM>2
      solnData(MHZO_VAR,:,:,:) = solnData(MOHZ_VAR,:,:,:)
#endif
#endif

      !Initialise some variables to zero

#ifdef UEUV_VAR
      solnData(UEUV_VAR,:,:,:) = 0.0
#endif

#ifdef REIR_VAR
      solnData(REIR_VAR,:,:,:) = 0.0
#endif

#ifdef PEFL_VAR
      solnData(PEFL_VAR,:,:,:) = 0.0
#endif

#ifdef ULWD_VAR
      solnData(ULWD_VAR,:,:,:) = 0.0
#endif

      call Grid_releaseBlkPtr(blockID,solnData)
    end do

  endif

  !Reset frequency band counter to the number of bands
  band_count = rt_freqbands
  do
    !Choose current frequency band
    SELECT CASE(band_count)
      CASE(1)
        current_band = rt_band1
#ifdef UEUV_VAR
        resetlog = .not. multiple_ionbands .or. current_band .ne. 'EUV_13P6_15P2'
#else
        resetlog = .true.
#endif
        ! reset iteration logs. Don't reset if multiple Ion bands are present since the saved value will be printed below.
        if(resetlog) then
          cumulative_iter_Band1 = 0
          max_iter_band1 = 0
          min_iter_band1 = rt_maxits
        endif
      CASE(2)
        current_band = rt_band2
        if(resetlog) then
          cumulative_iter_Band2 = 0
          max_iter_band2 = 0
          min_iter_band2 = rt_maxits
        endif
      CASE(3)
        current_band = rt_band3
        if(resetlog) then
          cumulative_iter_Band3 = 0
          max_iter_band3 = 0
          min_iter_band3 = rt_maxits
        endif
      CASE(4)
        current_band = rt_band4
        if(resetlog) then
          cumulative_iter_Band4 = 0
          max_iter_band4 = 0
          min_iter_band4 = rt_maxits
        endif
      CASE(5)
        current_band = rt_band5
        if(resetlog) then
          cumulative_iter_Band5 = 0
          max_iter_band5 = 0
          min_iter_band5 = rt_maxits
        endif
      CASE DEFAULT
        if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,I0)") 'Band Count, No of bands: ', band_count, rt_freqbands
        call Driver_abortFlash("[VETTAM]: Something is wrong; band_count is greater than 5, should not be happening!")
    END SELECT

#ifdef UEUV_VAR
    !If there are multiple EUV bands then the lower EUV band is already taken care of (in rt_ionise), so you can skip
    if(current_band .eq. 'EUV_13P6_15P2' .and. multiple_ionbands) then
      call Print_PETSCConvergence()
      !Reset the iterations counter for this band
      if(current_band .eq. rt_band1) then
        cumulative_iter_Band1 = 0
        min_iter_Band1 = rt_maxits
        max_iter_Band1 = 0
      else if(current_band .eq. rt_band2) then
        cumulative_iter_Band2 = 0
        min_iter_Band2 = rt_maxits
        max_iter_Band2 = 0
      else if(current_band .eq. rt_band3) then
        cumulative_iter_Band3 = 0
        min_iter_Band3 = rt_maxits
        max_iter_Band3 = 0
      else if(current_band .eq. rt_band4) then
        cumulative_iter_Band4 = 0
        min_iter_Band4 = rt_maxits
        max_iter_Band4 = 0
      else if(current_band .eq. rt_band5) then
        cumulative_iter_Band5 = 0
        min_iter_Band5 = rt_maxits
        max_iter_Band5 = 0
      endif
      band_count = band_count - 1
      if(band_count .eq. 0) exit
      cycle
    endif
#endif
    
    directSinkContribution = .true.
#ifdef RAYTRACE_3DRT
    ! Allocate memory and/or initialise stuff for RayTrace
    if(doRTDiffuse .or. sink_useRayTrace .or. current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW' &
      ) call rt_prepRayTrace(blockCount_, blockList_, dt)
#endif

    if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') then
#ifdef FSHM_VAR
      !If purely using raytrace for LW band, call rt_lw here since the lines of code below are unnecessary
      if(lwdiss_type .eq. 'raytrace') then
        call rt_lw(blockCount_, blockList_,dt) !This subroutine does the necessary action for the LW band
        band_count = band_count -1 !Update for next iteration of do loop
        if(band_count .eq. 0) then
          exit !Exit if all bands completed
        else
          cycle ! Nothing else to be done for this band here
        endif
      endif
#else
      call Driver_abortFlash("VETTAM/H2Dissociation not compiled in, but LW band is included.")
#endif
    endif

    !Add Sink/direct-RT Contributions
#ifdef SINK_PART_TYPE
    if(useSinkParticles .and. directSinkContribution) then 
      call Timers_start("vet_sink")
      !Call RT for Sink (Note: This likely will not work correctly. Please contact Shyam for more details.)
      if(sink_useRayTrace) then
#ifdef SINK_RT
        call RadTrans_RayTrace_Sink()
        !Add contribution to gas energy/momentum
        call Timers_start("vet_sinkhydro")
        call rt_sinkHydro(blockCount_,blockList_,dt)
        call Timers_stop("vet_sinkhydro")
#else
        call Driver_abortFlash("sink_useRayTrace is set to true, but RayTrace/SinkRT is not compiled in.")
#endif
      else
        !Gaussian source term to radiation energy (recommended)
        call rt_sinkInject(blockCount_,blockList_,dt)
        !Gas energy/momentum contributions handled by rt_hydro_update in this case
      endif
      call Timers_stop("vet_sink")
    endif
#endif

    !Compute opacities for the respective band
    call Timers_start("vet_opacity")
    call rt_setOpacity()
    call Timers_stop("vet_opacity")
    
#ifdef TEMP_VAR
    !Save old time solution of TEMP_VAR for the IR update. This is done here instead of above as the EUV band can change TEMP_VAR.
    if(current_band .eq. 'IR') then
      do b = 1, blockCount_
            
        blockID = blockList_(b)
        call Grid_getBlkPtr(blockID,solnData)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        !Set old time quantities

        solnData(TMPO_VAR,:,:,:) = solnData(TEMP_VAR,:,:,:)
        call Grid_releaseBlkPtr(blockID,solnData)
      end do
    endif
#endif
  
    !Start 
    if(rt_useRadTrans) then
      ! Set Eddington Tensor
      call Timers_start("vet_setet")
      call rt_setEddTensor()
      call Timers_stop("vet_setet")

      !Fill Guard Cells for ETensors and TAUR as both used for wavespeeds
      gcmask(:) = .false.
      gcmask(XXED_VAR) = .true.
#if NDIM>1
      gcmask(XYED_VAR) = .true.
      gcmask(YYED_VAR) = .true.
#if NDIM>2
      gcmask(XZED_VAR) = .true.
      gcmask(YZED_VAR) = .true.
      gcmask(ZZED_VAR) = .true.
#endif
#endif
      gcmask(TAUR_VAR) = .true.
      gcmask(TAUP_VAR) = .true.
      if (.NOT.gcMaskLogged) call Logfile_stampVarMask(gcmask, .false., '[RadTrans]', 'gcMask')
      call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS,mask=gcmask,makeMaskConsistent = .true.,&
                               selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)
      gcMaskLogged = .TRUE.
      if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Implicit time integration with Petsc Started", "[VETTAM] ")
      call Timers_start("vet_implicit")
      call Radiation_step(dt)
      call Timers_stop("vet_implicit")
      !For ionising radiation solve for the ionisation state
#ifdef UEUV_VAR
      if(current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
          & current_band .eq. 'EUV_15P2_INFTY') then
        call Timers_start("vet_ionize")
        call rt_ionise(blockCount_,blockList_,dt,time)
        if(dr_globalMe.eq.MASTER_PE)  call Logfile_stamp("Ionization Done","[VETTAM]")
        call Timers_stop("vet_ionize")
      endif
      
#endif

#ifdef FSHM_VAR
      if(current_band .eq. 'LW' .or. current_band .eq. 'LYMAN_WERNER') then
        !This is the case where lwdiss_type is not "raytrace"; if it is, it never reaches this part of the code
        if(lwdiss_type .eq. 'raytrace') call Driver_abortFlash("lwdiss_type set to raytrace; should not get here")
        call rt_lw(blockCount_, blockList_,dt) !This subroutine does the necessary action for the LW band
      endif
#endif

      if(rt_update_hydro) then
        call Timers_start("vet_picard")
        !Update dust temperature (only applicable for the IR band)
        if(current_band .eq. 'IR') call rt_dustTemperature(dt)
        !Flag to add energy/momentum due to dust absorption
        dust_addsrcterms = .true.
#ifdef UEUV_VAR
        !Do not add if Photoionization is switched on; if it is rt_ionise.F90 takes care of the radiation pressure and UV reprocessing
        if((current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. current_band .eq. &
          & 'EUV_15P2_INFTY') .and. useEUVIonize) &
          & dust_addsrcterms = .false.
#endif
        ! Now add radiation pressure on dust and energy reprocessing to the IR (for the EUV bands, this is already done in rt_ionMomentum and rt_ionHeatCool)
        if(hy_useHydro .and. dust_addsrcterms) then
          call rt_dustTerms(dt)
        endif
        if(dr_globalMe.eq.MASTER_PE)  call Logfile_stamp("Other (dust/radiation pressure) radiation-hydro coupling done","[VETTAM]")
        !Fill guard cells after hydro quantities updated
        gcmask(:) = .false.
#ifdef TEMP_VAR
        gcmask(TEMP_VAR) = .true.
#endif
        gcmask(VELX_VAR) = .true.
#if NDIM>1
        gcmask(VELY_VAR) = .true.
#if NDIM>2
        gcmask(VELZ_VAR) = .true.
#endif
#endif
        
        ! TODO: Not sure if we need both fill Guard cells and Eos. I think doEos=True will do this for me. 
        if (.NOT.gcMaskLogged) call Logfile_stampVarMask(gcmask, .false., '[rt_hydro_update]', 'gcMask')
        call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS,mask=gcmask,makeMaskConsistent = .true.,&
                                selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged,doEos=.true.,eosMode=MODE_DENS_TEMP)
        gcMaskLogged = .TRUE.
        call Timers_stop("vet_picard")
      endif
    endif

    !Print Convergence info of the implicit matrix solver
    call Print_PETSCConvergence()

    band_count = band_count - 1

    !Update band counter and check exit
    if(band_count .eq. 0) exit
  end do

  !Update the time since the last call to zero
  tsinceLastCall = 0.0

  call Timers_stop("RadTrans")

end SUBROUTINE RadTrans

SUBROUTINE Radiation_Step(dtstep)
  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe
  implicit none

  real, intent(in) :: dtstep
  type(PetscInt) :: reason
  integer :: rad_substeps, substep_counter
  real :: dt_substep, min_dt_substep

  rad_substeps = 1
  dt_substep = dtstep
  !TODO: Hardcoded for now
  min_dt_substep = dtstep/64.
  do
    do substep_counter = 1, rad_substeps
      if(rad_substeps .gt. 1) then 
       if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,A,I0)") 'Radiation Step Subcycle: ', substep_counter,'/', rad_substeps
      endif
      call Petsc_step(dt_substep,reason) !Solve the implicit system (in rt_petsc.F90)
      if(reason .eq. -3) exit
    end do
    !Converged after subcycling
    if(reason .ge. 0) then 
      exit
    !Still did not converge, subcycle further
    else
      rad_substeps = rad_substeps*2
      dt_substep = dtstep/rad_substeps
      if(dt_substep .lt. min_dt_substep) call Driver_abortFlash("[VETTAM]: Radiation step iteration subcycling limit reached.")
    endif
  end do

END SUBROUTINE Radiation_Step