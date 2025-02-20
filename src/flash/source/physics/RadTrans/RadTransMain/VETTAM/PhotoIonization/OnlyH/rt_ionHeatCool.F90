!! Contains and calls the set of routines coupling the heating/cooling of photoionised gas.
!!
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
#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"
SUBROUTINE rt_ionHeatCool(blockCount_,blockList_,dt,time)
  use rt_ionisedata
  use rt_ionisemodule
  use Grid_interface
  use Eos_interface
  use RadTrans_data, ONLY: rt_speedlt

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k, error
  real, dimension(:,:,:,:), pointer :: solnData
  real :: IonizationRate, Ion_Heating, Recombination_Cooling, FF_Cooling, dE, Tneutral, Tion, tgas, iony

  !TODO: Can make this a runtime parameter
  Tneutral = 10.0 !Neutral gas temperature
  Tion     = 1.e4 !Fully ionised gas temperature
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
          do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
            !An ideal gas (or similar) equation of state: Explicit heating/cooling introduced here
#ifdef EINT_VAR
            dE = solnData(PHHE_VAR,i,j,k) * dt
            solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + dE
#ifdef ENER_VAR
            solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + dE
#endif
            !Two-temperature Isothermal approximation (i.e. when compiled with physics/Hydro/HydroMain/split/Bouchut/IonizeIsothermal or similar)
#elif defined(TGAS_VAR)
            iony = solnData(IONY_MSCALAR,i,j,k)
            tgas = (1.-iony)*Tion + iony*Tneutral
            solnData(TGAS_VAR,i,j,k) = tgas
#else
            call Driver_abortFlash("[VET_Ionise]: neither TEMP_VAR or TGAS_VAR is defined. Unsure how to account for ionisation-thermal state coupling.")
#endif

          end do
        end do
      end do
    !Update thermal pressure
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do
end SUBROUTINE rt_ionHeatCool
