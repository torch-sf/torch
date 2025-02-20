!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_sinkInject
!!
!! NAME
!!  
!!  rt_sinkInject
!!
!!
!! SYNOPSIS
!! 
!!  rt_sinkInject(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt)
!!  
!! DESCRIPTION
!!
!!	Inject radiation energy from sinks into the radiation field
!!
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!  time : the current time
!!***
#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"
SUBROUTINE rt_sinkInject(blockCount_, blockList_, dtstep)
#ifdef SINK_PART_TYPE
  use Particles_sinkData, ONLY : particles_global, localnp, localnpf
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getCellCoords, Grid_getMinCellSize
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_globalMe
  use RadTrans_data, ONLY: current_band, rt_freqbands
#ifdef UEUV_VAR
  use rt_ionisedata, ONLY: energyPerIonH, hnu, energyPerIonH2, energyPerIonH_13p6_15p2, multiple_ionbands
#endif
  implicit none

  integer, INTENT(IN) :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(in) :: dtstep

  real, dimension(:), allocatable :: x, y, z
  integer :: sizeX, sizeY, sizeZ
  real, dimension(:,:,:,:), pointer :: solnData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: lum_sink, radius_cluster, jstr, sigma_sink, sink_pos(NDIM)
  real :: xdist, ydist, zdist, distsqr_point, delx
  integer :: b, blockID, i, j, k, p, istat

  !Ionizing Star Case (e.g.: source/Particles/ParticlesMain/Sink/StellarEvolution/PopIII/)
#ifdef EUVRATE_13P6_15P2_PART_PROP
  integer, parameter :: gather_nprops = 6
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: EUVRATE_13P6_15P2_PART_PROP, EUVRATE_15P2_INFTY_PART_PROP, THERMPION_13P6_15P2_PART_PROP, &
                  &  THERMPION_15P2_INFTY_PART_PROP, THERMPION_15P2_INFTY_H2_PART_PROP, LUMLW_PART_PROP /)
  real :: lum_sink_sum, lum_weighted_thermpionH, lum_weighted_thermpionH2
  !New Star cluster case (e.g.:  /scratch/ek9/sm5890/flash_newnew/flash-rsaa/source/Particles/ParticlesMain/Sink/StellarEvolution/StarClusterNew)
#elif LUMEUV_PART_PROP
  integer, parameter :: gather_nprops = 3
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: LUMEUV_PART_PROP, LUMLW_PART_PROP, LUMPE_PART_PROP /)
  !Star cluster case (e.g.:  /scratch/ek9/sm5890/flash_newnew/flash-rsaa/source/Particles/ParticlesMain/Sink/StellarEvolution/StarCluster)
#elif NION_PART_PROP
    integer, parameter :: gather_nprops = 3
    integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: NION_PART_PROP, NPEP_PART_PROP, EPEP_PART_PROP /)
    real :: lum_sink_sum, lum_weighted_thermpion, energyPerIonH
  !Star case (e.g.: /scratch/ek9/sm5890/flash_newnew/flash-rsaa/source/Particles/ParticlesMain/Sink/StellarEvolution)
#else
  integer, parameter :: gather_nprops = 2
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: STELLAR_RADIUS_PART_PROP, LUMINOSITY_PART_PROP /)
#endif
  real, save :: sigma_star
  logical, save :: firstCall = .true.
#if defined(EION_PART_PROP) || defined(THERMPION_13P6_15P2_PART_PROP) || defined(THERMPION_15P2_INFTY_PART_PROP)
  real :: thermpionH, thermpionH2
#endif


  if(firstCall) then 
    call RuntimeParameters_get("sigma_star", sigma_star)
    firstCall = .false.
    if ((dr_globalMe .eq. MASTER_PE)) &
      & print *, 'Gaussian source input sigma_star = ',sigma_star, ' delx'
    firstCall = .false.
  endif !firstCall

  call Timers_start("sink_gaussian")
  !Gather all sink particles even from other procs
  call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

  !Currently just using minimum cell size globally in domain for dx
  call Grid_getMinCellSize(delx)

  !Reset JSTR_VAR to zero
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkPtr(blockID,solnData)
    !Reset JSTR_VAR to zero
    solnData(JSTR_VAR,:,:,:) = 0.0
#ifdef EIH_VAR
    solnData(EIH_VAR,:,:,:) = 0.0
#endif
#ifdef EIH2_VAR
    solnData(EIH2_VAR,:,:,:) = 0.0
#endif
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

#ifdef EION_PART_PROP
  lum_sink_sum = 0.0
  lum_weighted_thermpion = 0.0
  thermpionH = 0.0
#endif

#ifdef THERMPION_13P6_15P2_PART_PROP
  lum_sink_sum = 0.0
  lum_weighted_thermpionH = 0.0
  lum_weighted_thermpionH2 = 0.0
  thermpionH = 0.0
  thermpionH2 = 0.0
#endif

  !Loop over all sink particles in domain
  do p = 1, localnpf
    sink_pos = particles_global(POSX_PART_PROP:POSZ_PART_PROP,p)
  !Ionizing Star Case (e.g.: source/Particles/ParticlesMain/Sink/StellarEvolution/PopIII/)
#ifdef EUVRATE_13P6_15P2_PART_PROP
    !TODO: Here hnu is kept fixed in both the bands = 18eV. Should use 13.6eV and 18eV maybe respectively.
    if(current_band .eq. 'EUV_13P6_15P2') then
      lum_sink = particles_global(EUVRATE_13P6_15P2_PART_PROP,p) * hnu
    else if(current_band .eq. 'EUV_15P2_INFTY') then
      lum_sink = particles_global(EUVRATE_15P2_INFTY_PART_PROP,p) * hnu
    else if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') then
#ifdef LUMLW_PART_PROP
      lum_sink = particles_global(LUMLW_PART_PROP,p)
#else
      lum_sink = particles_global(NPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)
#endif
    else
      print *, "WARNING: Current band ", current_band, " not considered by Particles_sinkStellarEvolution. &
      & Setting quantities to zero."
      lum_sink = 0.0
    endif

    !Thermal energy per ionisation
#ifdef THERMPION_13P6_15P2_PART_PROP
    lum_sink_sum = lum_sink_sum + lum_sink
    if(current_band .eq. 'EUV_13P6_15P2') then
      lum_weighted_thermpionH = lum_weighted_thermpionH + lum_sink * particles_global(THERMPION_13P6_15P2_PART_PROP,p)
      lum_weighted_thermpionH2 = 0.0
      thermpionH = particles_global(THERMPION_13P6_15P2_PART_PROP,p)
      thermpionH2 = 0.0
    else if (current_band .eq.'EUV_15P2_INFTY') then
      lum_weighted_thermpionH = lum_weighted_thermpionH + lum_sink * particles_global(THERMPION_15P2_INFTY_PART_PROP,p)
      lum_weighted_thermpionH2 = lum_weighted_thermpionH + lum_sink * particles_global(THERMPION_15P2_INFTY_H2_PART_PROP,p)
      thermpionH = particles_global(THERMPION_15P2_INFTY_PART_PROP,p)
      thermpionH2 = particles_global(THERMPION_15P2_INFTY_H2_PART_PROP,p)
    endif
#endif

! Star cluster New (or similar) subgrid case (Particles/ParticlesMain/Sink/StellarEvolution/StarClusterNew/)
#elif defined(LUMEUV_PART_PROP)
    if(current_band .eq. 'EUV') then
#ifdef UEUV_VAR
      lum_sink = particles_global(LUMEUV_PART_PROP,p) !Note here it is not the EUVRate but the EUVLuminosity
#else
      call Driver_abortFlash("VETTAM/Photoionization is not compiled in but the EUV band is included.")
#endif
    else if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') then
      lum_sink = particles_global(LUMLW_PART_PROP,p)
    else if(current_band .eq. 'PE') then
      lum_sink = particles_global(LUMPE_PART_PROP,p)
    else
      lum_sink = 0.0
    endif

! Star cluster (or similar) subgrid case (Particles/ParticlesMain/Sink/StellarEvolution/StarCluster/)
#elif defined(NION_PART_PROP) 
    if(current_band .eq. 'EUV') then
#ifdef UEUV_VAR
      lum_sink = particles_global(NION_PART_PROP,p) * hnu
#ifdef EION_PART_PROP
      lum_sink_sum = lum_sink_sum + lum_sink
      lum_weighted_thermpion = lum_weighted_thermpion + lum_sink * particles_global(EION_PART_PROP,p)
      thermpionH = particles_global(EION_PART_PROP,p)
#endif
#else
      call Driver_abortFlash("VETTAM/Photoionization is not compiled in but the EUV band is included.")
#endif
    else if(current_band .eq. 'FUV') then 
      lum_sink = particles_global(NPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)
    else if(current_band .eq. 'IR' .and. rt_freqbands .eq. 1) then
      lum_sink = particles_global(NPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p) !Assume entire FUV luminosity of cluster reprocessed to IR at the subgrid level
    else if(current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW') then
#ifdef LUMLW_PART_PROP
      lum_sink = particles_global(LUMLW_PART_PROP,p)
#else
      lum_sink = particles_global(NPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)*particles_global(EPEP_PART_PROP,p)
#endif
    else
      lum_sink = 0.0
    endif
  ! Star particle subgrid case (Particles/ParticlesMain/Sink/StellarEvolution/)
#else
  lum_sink = particles_global(LUMINOSITY_PART_PROP,p)
#endif


    sigma_sink = sigma_star * delx

    !Loop over all blocks and cells
    do b = 1, blockCount_
      blockID = blockList_(b)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      call Grid_getBlkPtr(blockID,solnData)
      
      sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
      sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
      sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
      allocate(x(sizeX), stat=istat)
      allocate(y(sizeY), stat=istat)
      allocate(z(sizeZ), stat=istat)
      call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., z, sizeZ)
      call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., y, sizeY)
      call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., x, sizeX)

      do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
          do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

            xdist = x(i) - sink_pos(1)
            ydist = y(j) - sink_pos(2)
            zdist = z(k) - sink_pos(3)
            distsqr_point = xdist**2 + ydist**2 + zdist**2

            jstr = (lum_sink/(2.*PI*sigma_sink**2)**1.5) * exp(-(distsqr_point)/(2.*sigma_sink**2)) !Smoothed gaussian source
            solnData(JSTR_VAR,i,j,k) = solnData(JSTR_VAR,i,j,k) + jstr
            !If variable to store heating rate per ionization (of H) is present
#ifdef EIH_VAR
            solnData(EIH_VAR,i,j,k) = solnData(EIH_VAR,i,j,k) + jstr * thermpionH
#endif

#ifdef EIH2_VAR
            solnData(EIH2_VAR,i,j,k) = solnData(EIH2_VAR,i,j,k) + jstr * thermpionH2
#endif

          end do !i
        end do !j
      end do !k
      call Grid_releaseBlkPtr(blockID,solnData)
    end do !b
  enddo !p

#ifdef EION_PART_PROP
  !Compute luminosity-weighted mean thermal energy deposited per ionisation
  if(lum_sink_sum .gt. 0) energyPerIonH = lum_weighted_thermpion/lum_sink_sum
#endif

!Pop III case
#ifdef THERMPION_13P6_15P2_PART_PROP
  if(lum_sink_sum .gt. 0) then
    !Don't overwrite energyPerIonH if multiple_ionbands is on
    if(multiple_ionbands .and. current_band .eq. 'EUV_13P6_15P2') then
      energyPerIonH_13p6_15p2 = lum_weighted_thermpionH/lum_sink_sum
    else
      energyPerIonH = lum_weighted_thermpionH/lum_sink_sum
      energyPerIonH2 = lum_weighted_thermpionH2/lum_sink_sum
    endif
  endif
#endif

!If EIH_VAR (EIH2_VAR) is present, compute the mean energy per ionization of H (H2)
#if defined(EIH_VAR) || defined(EIH2_VAR)
  do b = 1, blockCount_
      blockID = blockList_(b)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      call Grid_getBlkPtr(blockID,solnData)
      do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
          do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
            if(solnData(JSTR_VAR,i,j,k) .gt. 0) then 
#ifdef EIH_VAR
              solnData(EIH_VAR,i,j,k) = solnData(EIH_VAR,i,j,k)/solnData(JSTR_VAR,i,j,k)
#endif
#ifdef EIH2_VAR
              solnData(EIH2_VAR,i,j,k) = solnData(EIH2_VAR,i,j,k)/solnData(JSTR_VAR,i,j,k)
#endif
            else
#ifdef EIH_VAR
              solnData(EIH_VAR,i,j,k) = 0.0
#endif
#ifdef EIH2_VAR
              solnData(EIH2_VAR,i,j,k) = 0.0
#endif
            endif
          end do !i
        end do !j
      end do !k
      call Grid_releaseBlkPtr(blockID,solnData)
  end do
#endif

  call Timers_stop("sink_gaussian")

#else
    integer, INTENT(IN) :: blockCount_
    integer, dimension(blockCount_), intent(IN) :: blockList_
    real, intent(in) :: dtstep
    !If Sinks not compiled in this is just a dummy routine
    return

#endif

END SUBROUTINE rt_sinkInject
