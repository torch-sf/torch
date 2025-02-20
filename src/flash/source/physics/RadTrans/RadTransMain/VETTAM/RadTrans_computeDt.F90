!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RadTrans_computeDt
!!
!! NAME
!!  RadTrans_computeDt
!!
!! SYNOPSIS
!!  RadTrans_computeDt(integer(IN)  :: blockID
!!                 integer(IN)  :: blkLimits(2,MDIM)
!!                 integer(IN)  :: blkLimitsGC(2,MDIM)
!!                 real,pointer :: solnData(:,:,:,:)
!!                 real(OUT)    :: dt_stir
!!                 real(OUT)    :: dt_minloc(5))
!!
!! DESCRIPTION
!!  compute a timestep limiter for VETTAM RHD
!!
!! ARGUMENTS
!!  blockID       --  local block ID
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_stir       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!! SEE ALSO
!!  Driver_computeDt
!!
!!***

subroutine RadTrans_computeDt(blockID,                      &
                        blkLimits,blkLimitsGC,        &
                        solnData,                     &
                        dt_rt, dt_minloc)

use Driver_interface, ONLY: Driver_getMype, Driver_getSimTime
use Grid_interface, ONLY : Grid_getDeltas
use RadTrans_data
use Hydro_data, ONLY : hy_meshMe, hy_cfl

#include "constants.h"
#include "Flash.h"
#ifdef MAGX_VAR
use Hydro_data, ONLY: hy_bref
#endif

implicit none

!! arguments
integer, intent(IN)   :: blockID
integer, intent(IN),dimension(2,MDIM)::blkLimits, blkLimitsGC
real, pointer           :: solnData(:,:,:,:)
real, intent(INOUT)     :: dt_rt
integer, intent(INOUT)  :: dt_minloc(5)


!! User defined 
real, dimension(MDIM)   :: cellsize
integer :: dr_myPE, i,j,k, temploc(5)
real :: cs2, v_max, rad2, dt_ltemp, dt_temp, delx
#ifdef MAGX_VAR
  real :: va2
#endif

!!===================================================================


if ((.not.rt_useRadTrans).or.(.not.rt_compute_Dt)) return

call Driver_getMype(GLOBAL_COMM, dr_myPE)
! initialize the timestep from this block to some obscenely high number
dt_temp = HUGE(0.0)
call Grid_getDeltas(blockID,cellsize)

SELECT CASE(rt_dt_type)

!Light crossing time of smallest cell
CASE(1)
  call Grid_getMinCellSize(delx)
  dt_rt = rt_dtFactor * delx/rt_speedlt
  dt_minloc(1) = 0
  dt_minloc(2) = 0
  dt_minloc(3) = 0
  dt_minloc(4) = 0
  dt_minloc(5) = 0

!Radiation modified CFL condition
CASE(2) 

 do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
   do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

       ! sound speed squared
#ifdef GAMC_VAR
       cs2 = solnData(GAMC_VAR,i,j,k) * solnData(PRES_VAR,i,j,k) / solnData(DENS_VAR,i,j,k)
#else
#ifdef PRES_VAR
       ! isothermal case
       cs2 = solnData(PRES_VAR,i,j,k) / solnData(DENS_VAR,i,j,k)
#else
       ! isothermal case (sound speed = 1)
       cs2 = 1.0
#endif
#endif

!Compute rad contribution
!TODO: Is the cellsize in delta in IAXIS or other? 
#ifdef ERAD_VAR
      rad2 = (4./9.)*solnData(ERAD_VAR,i,j,k)*(1.0 - exp(-1.0 * solnData(TAUR_VAR,i,k,k)*cellsize(IAXIS)))/ &
      solnData(DENS_VAR,i,j,k)
#else 
      rad2 = 0.0
#endif


!Account for magnetic fields, if present.
#ifdef MAGX_VAR
       ! Alfven speed squared
       va2 = ( solnData(MAGX_VAR,i,j,k)**2 + &
               solnData(MAGY_VAR,i,j,k)**2 + &
               solnData(MAGZ_VAR,i,j,k)**2 ) / (hy_bref*solnData(DENS_VAR,i,j,k))
       ! fast MHD wave speed
       v_max = sqrt(cs2 + va2 + rad2)
#else
       v_max = sqrt(cs2 + rad2)
#endif


       ! compute minimum timestep based on MHD waves
       dt_ltemp = cellsize(IAXIS) / (abs(solnData(VELX_VAR,i,j,k)) + v_max)
#if NDIM >= 2
       dt_ltemp = min(dt_ltemp, cellsize(JAXIS) / (abs(solnData(VELY_VAR,i,j,k)) + v_max))
#endif
#if NDIM == 3
       dt_ltemp = min(dt_ltemp, cellsize(KAXIS) / (abs(solnData(VELZ_VAR,i,j,k)) + v_max))
#endif
       if (dt_ltemp < dt_temp) then
          dt_temp = dt_ltemp
          temploc(1) = i
          temploc(2) = j
          temploc(3) = k
          temploc(4) = blockID
          temploc(5) = hy_meshMe
       endif

       if (isnan(dt_ltemp)) call Driver_abortFlash("[Hydro]: Unphysical state in Hydro! Aborting!")

      enddo
   enddo
enddo

! apply CFL and return the minimum timestep
dt_temp = hy_cfl * dt_temp
if (dt_temp < dt_rt) then
   dt_rt = dt_temp
   dt_minloc = temploc
endif

CASE DEFAULT 
  call Driver_abortFlash("[VET] : Check rt_dt_type in flash.par.")
END SELECT

return

end subroutine RadTrans_computeDt
