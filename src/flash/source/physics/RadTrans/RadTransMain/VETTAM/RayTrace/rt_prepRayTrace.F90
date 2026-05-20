!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/rt_prepRayTrace
!!
!!  NAME 
!!
!!  rt_prepRayTrace
!!
!!  SYNOPSIS
!!
!!  call rt_prepRayTrace( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt)
!!
!!  DESCRIPTION 
!!      This subroutine does some operations to prepare the code for the ray-trace. 
!!      This needs to be called before the sink/diffuse ray trace is called.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!
!!***
#undef DEBUG_GRID_GCMASK
subroutine rt_prepRayTrace(nblk, blklst, dt)

  use raytrace_data, ONLY: rt_useRayTrace

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
       Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_DIRICHLET

  use RadTrans_interface, ONLY: RadTrans_planckInt, RadTrans_sumEnergy
  use Eos_interface, ONLY: Eos_wrapped, Eos_guardCells
  use Opacity_interface, ONLY: Opacity
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use RadTrans_RayTrace_3DRT, ONLY : FirstStep, rad_resetData

#ifdef DEBUG_GRID_GCMASK
  use Logfile_interface, ONLY: Logfile_stampVarMask
#endif

  
  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt

  ! Local variables:
  integer :: lb, i, j, k
  integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)
  real    :: opac
  LOGICAL :: GcMask(NUNK_VARS)
  real, pointer   :: solnData(:,:,:,:)

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,parameter :: gcMaskLogged =.TRUE.
#endif
  !=========================================================================
  if(.not.rt_useRayTrace) return

  call Timers_start("RadTrans_Setup") 
  call FirstStep()
  call RadTrans_update_source_function(nblk, blklst)


#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(rt_gcMask, .FALSE., '[RayTrace]', 'gcNeed')
  end if
#endif

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

  GcMask(:) = .FALSE.
  GcMask(SOUR_VAR) = .TRUE.
  GcMask(OPAC_VAR) = .TRUE.
  GcMask(TAUP_VAR) = .TRUE.

  call rad_resetData()

  !
  ! This is the raytracer for the diffuse radiation field
  ! In this routine, we compute the formal solution of the
  ! radition field which is stored in 'mean'
  !
  !Communicate SOUR, OPAC and TAUP for 2 Guard cells - Important for the ray-trace
  call Grid_fillGuardCells(CENTER, ALLDIR,doEos=.FALSE.,minLayers=2,&
         maskSize=NUNK_VARS,mask=GcMask)

  ! Update the radiation field, it's actually just a formal step.
  ! Important when scattering is present, not relavant otherwise.
#ifdef LAMB_VAR
  call RadTrans_update_radiation_field(nblk, blklst)
#endif
  
  call Timers_stop("RadTrans_Setup") 

  return

end subroutine rt_prepRayTrace
