!!****f* source/physics/sourceTerms/Heat/Heat_computeDt
!!
!! NAME
!!  
!!  Heat_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Heat_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real(IN):  x(:), 
!!                  real(IN): dx(:), 
!!                  real(IN): uxgrid(:),
!!                  real(IN): y(:), 
!!                  real(IN): dy(:), 
!!                  real(IN): uygrid(:), 
!!                  real(IN):  z(:), 
!!                  real(IN): dz(:), 
!!                  real(IN): uzgrid(:), 
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for heating source term solver.
!! 
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  x, y, z         three, directional coordinates
!!  d*              deltas in each {*=x, y z} directions
!!  u*grid          velocity of grid expansion in {*=x, y z} directions
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Heat_computeDt (blockID, &
                              x, dx, uxgrid, &
                              y, dy, uygrid, &
                              z, dz, uzgrid, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )
#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif

  return
end subroutine Heat_computeDt

