!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_potentialListOfBlocks
!!
!! NAME 
!!
!!     Gravity_potentialListOfBlocks
!!
!! SYNOPSIS
!!
!!  call Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                     integer(IN) :: blockList(blockCount),
!!                            optional,integer(IN) :: potentialIndex)
!!
!! DESCRIPTION
!!
!!      This routine computes the gravitational potential on all
!!      blocks specified in the list, for the gravity implementations
!!      (i.e., various Poisson implementations), which make use of it
!!      in computing the gravitational acceleration.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.  For some implementation of Gravity,
!!      in particular with Barnes-Hut tee solver, additional combinations
!!      of boundary conditions may be supported.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!   potentialIndex : If present, determines which variable in UNK to use
!!                    for storing the updated potential.  If not present,
!!                    GPOT_VAR is assumed.
!!                    Presence or absense of this optional dummy argument
!!                    also determines whether some side effects are enabled
!!                    or disabled; see discussion of two modes under NOTES
!!                    below.
!!
!! NOTES
!!
!!  Gravity_potentialListOfBlocks can operate in one of two modes:
!!  * automatic mode  - when called without the optional potentialIndex.
!!    Such a call will usually be made once per time step, usually
!!    from the main time advancement loop in Driver_evolveFlash.
!!    Various side effects are enabled in this mode, see SIDE EFFECT below.
!!
!!  * explicit mode  - when called with the optional potentialIndex.
!!    The potential is stored in the variable explicitly given, and
!!    side effects like saving the previous potential in GPOL_VAR
!!    and updating some sink particle state and properties are
!!    suppressed.
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  gravitational potential.  Invokes a solver (of the Poisson equation)
!!  if necessary. On return, if potentialIndex is not present,
!!     GPOT_VAR:  contains potential for the current simulation time.
!!     GPOL_VAR (if defined): contains potential at the previous simulation time.
!!  On return, if potentialIndex is present, the UNK variable given by
!!  potentialIndex contains the newly computed potential.
!!
!!  May affect other variables related to particle properties if particles
!!  are included in the simulation.  In particular,
!!     PDEN_VAR (if defined): may get updated to the current density from
!!                particles if particles have mass.
!!
!!  There are additional side effects if sink particles are used.
!!  These effects happen by calls to Particles_sinkAccelGasOnSinks and
!!  Particles_sinkAccelSinksOnGas, which may update sink particle
!!  properties and additional UNK variables that store
!!  accelerations. The calls are only made in automatic mode.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Gravity implementation.
!!  The following information is subject to change without notice.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!     IMGM_VAR (image mass)
!!     IMGP_VAR (image potential)
!!  For the Multipole implementation:
!!     (none)
!!
!!***

!! NOTE: In the AMUSE coupled version gas->sink and sink->gas gravity
!! is handled using a gravity bridge. So no particle-gas interactions
!! should happen in this routine. They are all commented out. -JW

!!REORDER(4): solnVec

subroutine Gravity_potentialListOfBlocks(blockCount,blockList, potentialIndex)


#include "Flash.h"
#ifdef PFFT_WITH_MULTIGRID
  use gr_hgPfftData, ONLY: gr_hgbcTypes
#endif

  use Gravity_data, ONLY : grav_poisfact, grav_temporal_extrp, grav_boundary, &
       grav_unjunkPden, &
       useGravity, updateGravity, grv_meshComm
  use Cosmology_interface, ONLY : Cosmology_getRedshift, &
       Cosmology_getOldRedshift
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_updateGridVar, &
       Particles_sinkAccelSinksOnGas, Particles_sinkAccelGasOnSinks
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_notifySolnDataUpdate, &
       Grid_solvePoisson
  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN), optional :: potentialIndex


  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  integer       :: ierr

  real          :: redshift, oldRedshift
  real          :: scaleFactor, oldScaleFactor
  real          :: invscale, rescale
  integer       :: lb
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  integer       :: density
  integer       :: newPotVar
  logical       :: saveLastPot

  saveLastPot = (.NOT. present(potentialIndex))
  if (present(potentialIndex)) then
     newPotVar = potentialIndex
  else
     newPotVar = GPOT_VAR
  end if

  call Cosmology_getRedshift(redshift)
  call Cosmology_getOldRedshift(oldRedshift)
  
  scaleFactor = 1./(1.+redshift)
  oldScaleFactor = 1./(1.+oldRedshift)
  
  invscale = 1./scaleFactor**3

! Rescaling factor to try and keep initial guess at potential close to
! final solution (in cosmological simulations).  Source term in Poisson
! equation has 1/a(t)^3 in it; and in linear theory (Omega=1, matter dom.)
! the comoving density increases as a(t), so comoving peculiar potential
! (which is what we are calculating here) should vary as 1/a(t)^2.  For
! noncosmological simulations this has no effect, since oldscale = scale = 1.

  rescale = (oldScaleFactor/scaleFactor)**2

!=========================================================================

  if(.not.useGravity) return
  
  if(.not.updateGravity) return

  call Timers_start("gravity Barrier")
  call MPI_Barrier (grv_meshComm, ierr)
  call Timers_stop("gravity Barrier")

  call Timers_start("gravity")

#ifdef PFFT_WITH_MULTIGRID
  bcTypes = gr_hgbcTypes
#else
  bcTypes = grav_boundary
#endif

  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == ISOLATED)
     bcTypes = GRID_PDE_BND_ISOLATED
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where
  bcValues = 0.
     
  if (grav_temporal_extrp) then
     
     call Driver_abortFlash("shouldn't be here right now")
     !call extrp_initial_guess( igpot, igpol, igpot )
     
  else 
     
     do lb = 1, blockCount
        call Grid_getBlkPtr(blocklist(lb), solnVec)
#ifdef GPOL_VAR
        if (saveLastPot) solnVec(GPOL_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:)
#endif
        solnVec(newPotVar,:,:,:) = solnVec(newPotVar,:,:,:) * rescale

        ! CTSS - We should also be storing the old sink particle accelerations:
        !!! Not in AMUSE coupled version! I'm commenting out. - JW
!#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
!        if (saveLastPot) then   !... but only if we are saving the old potential - kW
!           solnVec(SGXO_VAR,:,:,:) = solnVec(SGAX_VAR,:,:,:)
!           solnVec(SGYO_VAR,:,:,:) = solnVec(SGAY_VAR,:,:,:)
!           solnVec(SGZO_VAR,:,:,:) = solnVec(SGAZ_VAR,:,:,:)
!        end if
!#endif

        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
     enddo
     
#ifdef GPOL_VAR
     if (saveLastPot) call Grid_notifySolnDataUpdate( (/GPOL_VAR/) )
#endif

  endif

! Poisson is solved with the total density of PDEN_VAR + DENS_VAR 
  density=DENS_VAR
! This only gets called if there are active particles.
!!! Not in AMUSE coupled version! I'm commenting out. - JW
!#ifdef PDEN_VAR
!#ifdef MASS_PART_PROP
!  call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR)
!  if (.NOT. grav_unjunkPden) call Grid_notifySolnDataUpdate( (/PDEN_VAR/) )
!  density = PDEN_VAR
!#ifdef DENS_VAR           
!  do lb = 1, blockCount
!     call Grid_getBlkPtr(blocklist(lb), solnVec)
!     solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
!          solnVec(DENS_VAR,:,:,:)
!     call Grid_releaseBlkPtr(blocklist(lb), solnVec)
!  enddo
!#endif
!#endif
!#endif

  invscale=grav_poisfact*invscale
  call Grid_solvePoisson (newPotVar, density, bcTypes, bcValues, &
       invscale)
  call Grid_notifySolnDataUpdate( (/newPotVar/) )

! Un-junk PDEN if it exists and if requested.
!!! Not in AMUSE coupled version! I'm commenting out. - JW
!#ifdef PDEN_VAR
!  if (grav_unjunkPden) then
!     density = PDEN_VAR
!#ifdef DENS_VAR           
!     do lb = 1, blockCount
!        call Grid_getBlkPtr(blocklist(lb), solnVec)
!        solnVec(density,:,:,:) = solnVec(density,:,:,:) - solnVec(DENS_VAR,:,:,:)
!        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
!     enddo
!     if (density .NE. PDEN_VAR) call Grid_notifySolnDataUpdate( (/density/) )
!#endif
!  end if
!#endif

!  if (.NOT. present(potentialIndex)) then
! Compute acceleration of the sink particles caused by gas
! We must call Particles_sinkAccelGasOnSinks() here first, because it
! updates particles_global positions and masses by calling pt_sinkGatherGlobal(),
! required in Particles_sinkAccelSinksOnGas(). But we avoid another call to
! pt_sinkGatherGlobal() by calling Particles_sinkAccelGasOnSinks() here first.
!!! Not in AMUSE coupled version! I'm commenting out. - JW
!     call Particles_sinkAccelGasOnSinks()
! Compute acceleration of the gas caused by sink particles
!     call Particles_sinkAccelSinksOnGas(blockCount,blockList)
! This is both done here to make the computation symmetric (CTSS)
!  end if


#ifdef USEBARS
  call MPI_Barrier (grv_meshComm, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_potentialListOfBlocks
