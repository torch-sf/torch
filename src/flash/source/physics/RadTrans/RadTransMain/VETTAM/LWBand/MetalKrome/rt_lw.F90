!! Photochemical rates from the LW band: Photodissociation of H2 and CO, and the Photoionization of C
!! AUTHOR
!!  Shyam Harimohan Menon (2024)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
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
#if defined C_SPEC || defined CO_SPEC
  use Multispecies_interface, ONLY : Multispecies_getProperty
#endif

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
#ifdef C_SPEC
    call Multispecies_getProperty(C_SPEC,A,cA)
#endif

#ifdef CO_SPEC
    call Multispecies_getProperty(CO_SPEC,A,coA)
#endif
    first_call = .false.
    !Return to RadTrans.F90 (where the first call will be made from; so the block of code above is basically all for initialising)
    return
  end if

  call Timers_start("LWRates")

  SELECT CASE(lwdiss_type)
    CASE("raytrace")
      !Get the LW flux and shielding factor with the RT
      call RayTrace_LW()

    CASE("hybrid") ! This uses the solution from the moment eqns for the flux, and the shielding factor calculated with the sink RT
      !Get the shielding factors with the RT
      call RayTrace_LW()      
      !Get the total dust-attenuated LW flux which we have from the Petsc solve
      call Compute_LWFL()
    CASE("noraytrace")
      call Driver_abortFlash("[rt_lw]: The noraytrace Sobolev-based fshield implementation for LW dissociation not yet implemented")
    CASE DEFAULT
      call Driver_abortFlash("[rt_lw]: lwdiss_type is unrecognisable")

  END SELECT

  !Now combine fshield from RT and net (dust-attenuated) flux from the moment solve to get shielded LW flux
  call Compute_Rates()

  call Timers_stop("LWRates")

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

          !Get the total dust-attenuated LW flux which we have from the Petsc solve
          totalflux = solnData(MOHX_VAR,i,j,k)**2
#if NDIM>1
          totalflux = totalflux + solnData(MOHY_VAR,i,j,k)**2
#if NDIM>2
          totalflux = totalflux + solnData(MOHZ_VAR,i,j,k)**2
#endif
#endif
          totalflux = SQRT(totalflux)

          !Now set the self-shielded flux = totalflux * fshield
          solnData(LWFL_VAR,i,j,k) = totalflux
          
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do



END SUBROUTINE Compute_LWFL


SUBROUTINE Compute_Rates()
  use rt_lwmodule
  implicit none
  real :: fshield, totalflux
  real, parameter :: H2dissrate_ISRF = 5.7e-11, COdissrate_ISRF = 2.4e-10, Cdissrate_ISRF = 3.5e-10 !Gong et al 2017 (in per s), Table 2
  real, parameter :: J_LW_ISRF = 3.e-5 !Mean intensity of ISRF in LW band (Eq. 10 Kim et al. 2023; erg s^-1 cm^-2 sr^-1)

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          !Now set the self-shielded flux = totalflux * fshield
          solnData(DIH2_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(FSHM_VAR,i,j,k)
          if(useHIshield) solnData(DIH2_VAR,i,j,k) = solnData(DIH2_VAR,i,j,k) * solnData(FSHA_VAR,i,j,k)
          solnData(IC_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(IC_VAR,i,j,k)
          solnData(DICO_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(DICO_VAR,i,j,k)

          !Now convert this to dissociation rates by i) converting to mean intensity, ii)scaling with the ISRF value, and iii) multiplying by the optically-thin rate
          solnData(DIH2_VAR,i,j,k) = solnData(DIH2_VAR,i,j,k) / (4*PI * J_LW_ISRF) * H2dissrate_ISRF
          solnData(IC_VAR,i,j,k) = solnData(IC_VAR,i,j,k) / (4*PI * J_LW_ISRF) * Cdissrate_ISRF
          solnData(DICO_VAR,i,j,k) = solnData(DICO_VAR,i,j,k) / (4*PI * J_LW_ISRF) * COdissrate_ISRF
          
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do



END SUBROUTINE Compute_Rates

END SUBROUTINE rt_lw
