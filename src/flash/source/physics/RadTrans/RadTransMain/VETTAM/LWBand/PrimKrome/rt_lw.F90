!! Photodissociation due to LW photons of H2 molecules
!! AUTHOR
!!  Shyam Harimohan Menon (2023)
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


  call Timers_stop("LWDiss")

CONTAINS 

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
