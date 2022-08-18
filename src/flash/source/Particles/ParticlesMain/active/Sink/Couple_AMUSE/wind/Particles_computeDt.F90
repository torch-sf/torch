!!****if* source/Particles/ParticlesMain/Particles_computeDt
!!
!! NAME
!!  Particles_computeDt
!!
!! SYNOPSIS
!!  Particles_computeDt( integer(in) :: blockID,
!!                       real(inout) :: dt_part,
!!                       integer(5)(inout) :: dt_minloc)
!!
!! DESCRIPTION
!!   Timestep computation routine for the particle unit.
!!   This routine sets the timestep by requiring that
!!   particles travel no more than some fraction pt_dtFactor
!!   during a single step.  Even for time integration
!!   schemes that can deal with pt_dtFactor > 1, you should
!!   not set pt_dtFactor to be larger than half the number
!!   of guard cells.  Otherwise particles that leave the block
!!   that "owns" them may overshoot the block's immediate
!!   neighbors, causing problems when particle data are transmitted
!!   to the neighbors.
!!
!!
!! ARGUMENTS
!!
!!    blockID:         local block ID
!!    
!!    dt_part:         variable to hold timestep constraint
!!    dt_minloc(5):    array to hold limiting zone info:  index[1-3] are zone
!!                     indices (i,j,k); index[4]=block ID; index[5]=processor 
!!                     number. The zone indices indicate the zone containing
!!                     the particle whose velocity restricted the timestep.
!!
!! PARAMETERS
!!
!!   pt_dtFactor:      REAL [default 0.5] A factor multiplying dx/|v| to limit
!!                       the movement of particles outside a single block. 
!!   pt_small:         REAL [default 1.0E-10] if velocities are larger than this
!!                       value, time stepping may be limited
!!
!!
!! NOTES
!!
!!   Note that this routine does NOT calculate a timestep limitation
!!   based on the theoretical stability limits of the time integration
!!   routines.  It SHOULD, but these limits were deemed to be unlikely
!!   to be restrictive based on the high velocities involved in normal
!!   astronomical simulations.
!!
!!***
!!===============================================================================
  
subroutine Particles_computeDt (blockID, dt_part, dt_minloc)

  !================================================================

  use Particles_data, ONLY:  pt_dtFactor, pt_numLocal, particles, pt_small, pt_meshMe
  use Particles_windData, ONLY : min_wind_dt
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
    Grid_getBlkBoundBox, Grid_getDeltas
  use Particles_interface, only: Particles_sinkComputeDt

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, INTENT(in)    :: blockID
  
  real, INTENT(inout)    :: dt_part
  integer, INTENT(inout) :: dt_minloc(5)


  integer         :: i, lb
  real            :: bsize(MDIM), delta(MDIM), lowerBound(MDIM), boundBox(2,MDIM)
  real            :: dtx, dty, dtz, dtnew, velxabs, velyabs, velzabs

  !===============================================================================

!!! Okay now we set this first to the wind dt, then we check
!!! the sink dt, with that function picking the smaller of the two.
!!! The rest of the regular particle timestepping is commented out
!!! since these are integrated in AMUSE. -JW
 
 dt_part = min(dt_part,min_wind_dt)

!  !! In Flash3, we don't assume a regular grid.  Ask the package for its dimensions
!  call Grid_getBlkPhysicalSize(blockID,bsize)   ! physical size of the block in each direction
!  call Grid_getBlkBoundBox(blockID,boundBox)    ! physical bounding box of the block
!  lowerBound = boundBox(1,:)                        


!  call Grid_getDeltas(blockID,delta)



!  !! loop over all local particles
!  do i = 1, pt_numLocal
  
!     !! find the local block of each particle
!     lb = int(particles(BLK_PART_PROP,i))
     
!     if (lb == blockID) then
        
!        dtx   = HUGE(1.0)
!        dty   = HUGE(1.0)
!        dtz   = HUGE(1.0)
        
!        !! calculate the time needed to move a particle through the block
!        !!      in each direction:  blocksize / particle velocity
        
        
!        velxabs = abs(particles(VELX_PART_PROP,i))
!        if (velxabs > pt_small) then
!           dtx = delta(1) / velxabs
!        endif
        
!        if (NDIM > 1) then
!           velyabs = abs(particles(VELY_PART_PROP,i))
!           if (velyabs > pt_small) then
!              dty = delta(2) / velyabs
!           endif
!        end if

!        if (NDIM > 2) then
!           velzabs = abs(particles(VELZ_PART_PROP,i))
!           if (velzabs > pt_small) then
!              dtz = delta(3) / velzabs
!           endif
           
!        end if
           

!        !! decrease by the restraining factor
        
!        dtnew = pt_dtFactor * min(dtx, dty, dtz)
!        if (dtnew < dt_part) then
!           dt_part = dtnew
!           !! information about where the minimum restriction took place
!           dt_minloc(1) = int((particles(POSX_PART_PROP,i)-lowerBound(1))/delta(1)) + NGUARD+1
!           if (NDIM > 1) &
!                dt_minloc(2) = int((particles(POSY_PART_PROP,i)-lowerBound(2))/delta(2)) + NGUARD+1
!           if (NDIM > 2) &
!                dt_minloc(3) = int((particles(POSZ_PART_PROP,i)-lowerBound(3))/delta(3)) + NGUARD+1
!           dt_minloc(4) = lb
!           dt_minloc(5) = pt_meshMe
!        endif
!     endif
     
!  enddo
!  !
!  !===============================================================================
  call Particles_sinkComputeDt(blockID,dt_part,dt_minloc)
  return
  
end subroutine Particles_computeDt
