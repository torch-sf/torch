!! Contains and calls the set of routines coupling the momentum of EUV radiation.
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
SUBROUTINE rt_ionMomentum(blockCount_,blockList_,dt,time)
  use rt_ionisedata
  use rt_ionisemodule
  use Grid_interface
  use Eos_interface
  use RadTrans_data, ONLY: rt_speedlt, current_band, rt_band1, rt_band2, rt_band3

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b, blockID, i, j, k, error
  real, dimension(:,:,:,:), pointer :: solnData
  real :: opac_total, ekin_old, ekin_new, mohx, mohy, mohz, d_ekin, dMomx, dMomy, dMomz
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          dMomx = 0.0
          dMomy = 0.0
          dMomz = 0.0

          ! This is the opacity due to dust + the opacity due to gas
          opac_total = solnData(TAUR_VAR,i,j,k)
          !Store Kinetic Energy
          ekin_old = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_old = ekin_old+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_old = ekin_old + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif

          if(multiple_ionbands .and. .not. ion_implicit) then
            !MOHX
            if(current_band .eq. rt_band1) then
#ifdef B1FX_VAR
              mohx = solnData(B1FX_VAR,i,j,k)
#else
              call Driver_abortFlash("B1FX missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z & 
                                  & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band2) then
#ifdef B2FX_VAR
              mohx = solnData(B2FX_VAR,i,j,k)
#else
              call Driver_abortFlash("B2FX missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z & 
                                  & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band3) then
#ifdef B3FX_VAR
              mohx = solnData(B3FX_VAR,i,j,k)
#else
              call Driver_abortFlash("B3FX missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z & 
                                  & has to be present in the simulation config.")
#endif
            endif
          else  
            mohx = solnData(MOHX_VAR,i,j,k)
          endif
#if NDIM>1
          if(multiple_ionbands .and. .not. ion_implicit) then
            !MOHY
            if(current_band .eq. rt_band1) then
#ifdef B1FY_VAR
              mohy = solnData(B1FY_VAR,i,j,k)
#else
              call Driver_abortFlash("B1FY missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band2) then
#ifdef B2FY_VAR
              mohy = solnData(B2FY_VAR,i,j,k)
#else
              call Driver_abortFlash("B2FY missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band3) then
#ifdef B3FY_VAR
              mohy = solnData(B3FY_VAR,i,j,k)
#else
              call Driver_abortFlash("B3FY missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            endif
          else  
            mohy = solnData(MOHY_VAR,i,j,k)
          endif
#if NDIM>2
          if(multiple_ionbands .and. .not. ion_implicit) then
            !MOHZ
            if(current_band .eq. rt_band1) then
#ifdef B1FZ_VAR
              mohz = solnData(B1FZ_VAR,i,j,k)
#else
              call Driver_abortFlash("B1FZ missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band2) then
#ifdef B2FZ_VAR
              mohz = solnData(B2FZ_VAR,i,j,k)
#else
              call Driver_abortFlash("B2FZ missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            else if(current_band .eq. rt_band3) then
#ifdef B3FZ_VAR
              mohz = solnData(B3FZ_VAR,i,j,k)
#else
              call Driver_abortFlash("B3FZ missing; When there are multiple H-ionisation bands the corresponding BnFx/y/z &
                                   & has to be present in the simulation config.")
#endif
            endif
          else  
            mohz = solnData(MOHZ_VAR,i,j,k)
          endif
#endif
#endif


          
          !Momentum gain 
          dMomx = dMomx + dt * opac_total * mohx/rt_speedlt
#if NDIM>1
          dMomy = dMomy + dt * opac_total * mohy/rt_speedlt
#if NDIM>2
          dMomz = dMomz + dt * opac_total * mohz/rt_speedlt
#endif
#endif

          !Update velocities
          solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + dMomx/solnData(DENS_VAR,i,j,k)
#if NDIM>1
          solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + dMomy/solnData(DENS_VAR,i,j,k)
#if NDIM>2
          solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + dMomz/solnData(DENS_VAR,i,j,k)  
#endif 
#endif

          !New Kinetic Energy
          ekin_new = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_new = ekin_new+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_new = ekin_new + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif
          ! compute energy injection
          d_ekin = ekin_new - ekin_old

#ifdef ENER_VAR
          ! update the total energy to be consistent with new velocities
          solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + d_ekin
#endif

        end do
      end do
    end do
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do
END SUBROUTINE rt_ionMomentum


