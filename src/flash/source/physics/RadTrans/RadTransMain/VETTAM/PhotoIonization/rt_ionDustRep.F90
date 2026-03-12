!! Computes the amount of energy in the EUV radiation field reprocessed to the IR by dust
!!
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
#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"

#ifdef REIR_VAR
SUBROUTINE rt_ionDustRep(blockCount_,blockList_,dt,time)
  use rt_ionisedata
  use rt_ionisemodule
  use Grid_interface
  use Eos_interface
  use RadTrans_data, ONLY: rt_speedlt, rt_band1, rt_band2, rt_band3, current_band

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k
  real, dimension(:,:,:,:), pointer :: solnData
  real :: opac_dust, erad
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          if(multiple_ionbands .and. .not. ion_implicit) then
            if(current_band .eq. rt_band1) then
#ifdef B1ER_VAR
              erad = solnData(B1ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B1ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band2) then
#ifdef B2ER_VAR
              erad = solnData(B2ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B2ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band3) then
#ifdef B3ER_VAR
              erad = solnData(B3ER_VAR,i,j,k)
#else
              call Driver_abortFlash("B3ER missing; When there are multiple H-ionisation bands the corresponding BnEr &
                                & has to be present in the simulation config.")
#endif
            endif
          else
            erad = solnData(ERAD_VAR,i,j,k)
          endif

          !Opacity due to dust alone
          opac_dust = solnData(TAUP_VAR,i,j,k) - solnData(TAUH_VAR,i,j,k)
          solnData(REIR_VAR,i,j,k) = solnData(REIR_VAR,i,j,k) + &
            & opac_dust * rt_speedlt * erad * dt
        end do
      end do
    end do
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

end SUBROUTINE rt_ionDustRep
#endif