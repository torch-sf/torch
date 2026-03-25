!! Photodissociation due to LW photons of H2 molecules
!! AUTHOR
!!  Shyam Harimohan Menon (2023)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks 
!!  dt : the current timestep
!!  time : the current time
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
SUBROUTINE rt_lw(blockCount_,blockList_,dt)
  use Grid_interface
  use Eos_interface
  use Grid_data, ONLY: gr_smallx
  use Hydro_data, ONLY: hy_useHydro
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use rt_ionisedata, ONLY: h2A, hA
  use rt_lwdata

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt
  integer                                    :: b,blockID, nAllParticles, myPE, i, j, k
  integer,dimension(2,MDIM)                  :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer            :: solnData
  logical                                    :: gcmask(NUNK_VARS), gcMaskLogged = .false.
  logical, save                              :: first_call = .true.

  if(first_call) then
#ifndef RAYTRACE_3DRT
  call Driver_abortFlash("Ray-tracer needs to be compiled for the LW dissociation to work.")
#endif

#ifndef H2_SPEC
  call Driver_abortFlash("[rt_lw]: The Multispecies variable H2_SPEC needs to be defined for this implementation.")
#endif
    !Read runtime params
    call RuntimeParameters_get("useH2Dissociate",useH2Dissociate)
    call RuntimeParameters_get("useHIshield",useHIshield)
    call RuntimeParameters_get("lwdiss_type",lwdiss_type)
    call RuntimeParameters_get("bfive",bfive)
    call RuntimeParameters_get("fpump",fpump)
    call RuntimeParameters_get("energyPerDissociation",energyPerDissociation)
    call RuntimeParameters_get("avgEnergyLW",avgEnergyLW)
    first_call = .false.
    !Return to RadTrans.F90 (where the first call will be made from; so the block of code above is basically all for initialising)
    return
  end if

  if (.not. useH2Dissociate)then
    return
  endif

  call Timers_start("LWDiss")

  SELECT CASE(lwdiss_type)
    CASE("raytrace")
      !Get the LW flux and shielding factor with the RT
      call RayTrace_LW()

    CASE("hybrid") ! This uses the solution from the moment eqns for the flux, and the shielding factor calculated with the sink RT
      !Get the shielding factor with the RT
      call RayTrace_LW()
      !Now combine fshield from RT and net (dust-attenuated) flux from the moment solve to get shielded LW flux
      call Compute_LWFL()

    CASE("noraytrace")
      call Driver_abortFlash("[rt_lw]: The noraytrace Sobolev-based fshield implementation for LW dissociation not yet implemented")
    CASE DEFAULT
      call Driver_abortFlash("[rt_lw]: lwdiss_type is unrecognisable")

  END SELECT

  !Now change the H2 fraction and add heating due to dissociation
  call H2Dissociate()

  !Fill guard cells.
  gcmask(:) = .false.
  !Only update H2_SPEC
  gcmask(H2_SPEC) = .true.
  !Fill Guard cells
  call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS,mask=gcmask,makeMaskConsistent = .true.,&
  selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)

  call Timers_stop("LWDiss")

CONTAINS 

SUBROUTINE H2Dissociate()
  use rt_lwmodule
  implicit none
  real :: lwfl, H2DissRate, H2DissHeat, H2PumpHeat, nH2_old, nH2_new, nH, nH_new, dnH, sum,n,energyPerPump

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          nH2_old = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
          
          !Dissociation
          lwfl = solnData(LWFL_VAR,i,j,k) !Lyman Werner self-shielded flux
          H2DissRate = lwfl/(2.e7 * avgEnergyLW) * 5.8e-11 ! Dissociation rate in per unit second
          nH2_new = nH2_old*exp(-H2DissRate*dt)
          solnData(H2_SPEC,i,j,k) = max(nH2_new * h2A/solnData(DENS_VAR,i,j,k),gr_smallx)

          !Update H abundance
          nH = solnData(DENS_VAR,i,j,k) * solnData(IHA_SPEC,i,j,k)/hA
          dnH = 2*(nH2_old-nH2_new) !Each dissociation produces two H atoms
          nH_new = nH + dnH
          solnData(IHA_SPEC,i,j,k) = nH_new * hA/solnData(DENS_VAR,i,j,k)

          !Renormalise
          !Normalise species fraction to maintain closure relation to machine precision
          ! sum up species fractions
          sum = 0.0
          do n = SPECIES_BEGIN, SPECIES_END
            sum = sum + max(gr_smallx, solnData(n,i,j,k))
          enddo

          ! re-normalise sum of species fractions to 1
          do n = SPECIES_BEGIN, SPECIES_END
            solnData(n,i,j,k) = solnData(n,i,j,k) / sum
          enddo

          if(hy_useHydro) then
            !Heating due to dissociations
            H2DissHeat = H2DissRate * nH2_old * energyPerDissociation * dt
            !Heating due to UV pumping
            
            energyPerPump = get_Epump(solnData(TEMP_VAR,i,j,k),nH2_old,nH)
            H2PumpHeat = fpump * H2DissRate * nH2_old * energyPerPump * dt
#ifdef EINT_VAR
            solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + (H2DissHeat+H2PumpHeat)/solnData(DENS_VAR,i,j,k)
#ifdef ENER_VAR
            solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + (H2DissHeat+H2PumpHeat)/solnData(DENS_VAR,i,j,k)
#endif
#endif
          endif
        end do
      end do
    end do
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

END SUBROUTINE H2Dissociate

SUBROUTINE Compute_LWFL()
  use rt_lwmodule
  implicit none
  real :: fshield, totalflux

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          fshield = solnData(FSHM_VAR,i,j,k)
          !Include shielding due to atomic H if useHIShield included
          if(useHIshield) fshield = fshield * solnData(FSHA_VAR,i,j,k) 
          totalflux = solnData(MOHX_VAR,i,j,k)**2
#if NDIM>1
          totalflux = totalflux + solnData(MOHY_VAR,i,j,k)**2
#if NDIM>2
          totalflux = totalflux + solnData(MOHZ_VAR,i,j,k)**2
#endif
#endif
          totalflux = SQRT(totalflux)

          !Now set the self-shielded flux = totalflux * fshield
          solnData(LWFL_VAR,i,j,k) = totalflux * fshield
          
        end do
      end do
    end do
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do



END SUBROUTINE Compute_LWFL

END SUBROUTINE rt_lw