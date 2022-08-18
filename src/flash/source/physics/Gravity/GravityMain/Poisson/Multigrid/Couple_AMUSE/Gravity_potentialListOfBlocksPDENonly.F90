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

!!REORDER(4): solnVec

subroutine Gravity_potentialListOfBlocksPDENonly(blockCount,blockList, potentialIndex, which_particles)


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
       Particles_updateGridVarMassiveOnly, Particles_updateGridVarSinksOnly

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_notifySolnDataUpdate, &
       Grid_solvePoisson, Grid_getSingleCellCoords, &
       Grid_getDeltas, Grid_getBlkIndexLimits, Grid_mapParticlesToMesh
       
  use Particles_sinkData, ONLY : particles_local, localnp, sink_maxSinks
       
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN), optional :: potentialIndex
  integer, intent(IN), optional :: which_particles


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
  logical, save :: firstcall
  
  integer       :: blkLimits(2,NDIM), blkLimitsGC(2,NDIM)
  real          :: part_mass, del(3), cell_com(3), cell_coords(3)
  integer       :: ii, jj, kk, part_block_id
  integer       :: part_type

  saveLastPot = (.NOT. present(potentialIndex))
  ! This version is for PDEN only, exit with an error if we are
  ! not storing in SGPT_VAR (Sink on Gas PoTential)  here. -JW
  if (present(potentialIndex)) then
     newPotVar = potentialIndex
     !print*, "newPotVar =", newPotVar, " and SGPT_VAR =", SGPT_VAR
  else
     !newPotVar = GPOT_VAR
     call Driver_abortFlash("Gravity_potentialListOfBlocksPDENonly: &
                            You're not using the BGPT_VAR here!")
  end if

  if (present(which_particles)) then
     part_type = which_particles
  else
     part_type = 0
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

  bcTypes = grav_boundary
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
     
  !else 
     
  !   do lb = 1, blockCount
  !      call Grid_getBlkPtr(blocklist(lb), solnVec)
!#ifdef GPOL_VAR
!        if (saveLastPot) solnVec(GPOL_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:)
!#endif
    !    solnVec(newPotVar,:,:,:) = solnVec(newPotVar,:,:,:) * rescale
        
        !print*, "Before mapping max PDEN is ", maxval(solnVec(PDEN_VAR, :, :, :))

        ! CTSS - We should also be storing the old sink particle accelerations:
!#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
!        if (saveLastPot) then   !... but only if we are saving the old potential - kW
!           solnVec(SGXO_VAR,:,:,:) = solnVec(SGAX_VAR,:,:,:)
!           solnVec(SGYO_VAR,:,:,:) = solnVec(SGAY_VAR,:,:,:)
!           solnVec(SGZO_VAR,:,:,:) = solnVec(SGAZ_VAR,:,:,:)
!        end if
!#endif

  !      call Grid_releaseBlkPtr(blocklist(lb), solnVec)
  !   enddo
     
!#ifdef GPOL_VAR
!     if (saveLastPot) call Grid_notifySolnDataUpdate( (/GPOL_VAR/) )
!#endif

  endif

! Poisson is solved with the total density of PDEN_VAR + DENS_VAR 
!  density=DENS_VAR

! This only gets called if there are active particles.
! Note that we have disabled the check that prevents sink masses from being
! mapped to the PDEN_VAR on the grid, since that's what we want here.
! Remember, with great power comes great responsibilty! - JW
!#ifdef PDEN_VAR
!#ifdef MASS_PART_PROP

if (part_type == ACTIVE_PART_TYPE) then
  call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR) ! Just active/massive.
else if (part_type == SINK_PART_TYPE) then
  call Particles_updateGridVarSinksOnly(MASS_PART_PROP, PDEN_VAR) ! Just sinks.
else
  call Driver_abortFlash("[Gravity_potentialListOfBlocksPDENonly]: Improper part type given!")
end if

! Let's try to map it by hand and see what we get:

!  call Grid_mapParticlesToMesh( &
!    particles_local(:,:), NPART_PROPS, localnp, sink_maxSinks, &
!    MASS_PART_PROP, PDEN_VAR, 0)

!  part_block_id = particles_local(BLK_PART_PROP,1)
!  call Grid_getDeltas(part_block_id,del)

!  if (.NOT. grav_unjunkPden) call Grid_notifySolnDataUpdate( (/PDEN_VAR/) )
  density = PDEN_VAR
!#ifdef DENS_VAR
!  part_mass = 0.0
!  do lb = 1, blockCount
!     call Grid_getBlkPtr(blocklist(lb), solnVec)
!     call Grid_getBlkIndexLimits(blocklist(lb), blkLimits, blkLimitsGC, CENTER)
!!     solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
!!          solnVec(DENS_VAR,:,:,:)

!    do ii=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!      do jj=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!        do kk=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
!          if (solnVec(density,ii,jj,kk) .gt. 1e-30) then
          
!            call Grid_getSingleCellCoords([ii,jj,kk],blocklist(lb), CENTER, &
!                                     EXTERIOR, cell_coords)
            
!            cell_com = cell_com + cell_coords*solnVec(density,ii,jj,kk)
            
!            part_mass = part_mass + solnVec(density,ii,jj,kk)
            
!          end if

!        end do
!      end do
!    end do

!     call Grid_releaseBlkPtr(blocklist(lb), solnVec)
!  enddo
  
!  cell_com = cell_com / part_mass
!  part_mass = part_mass*(del(1)**3.0)
  
!  print*, "After mapping total PDEN*cellVol is ", part_mass
!  print*, "After mapping total cell mapped particle center of mass is ", cell_com
!  print*, "Particle location is ", particles_local(POSX_PART_PROP, 1), &
!                                   particles_local(POSY_PART_PROP, 1), &
!                                   particles_local(POSZ_PART_PROP, 1)

!#endif
!#endif
!#endif

  invscale=grav_poisfact*invscale
  call Grid_solvePoisson (newPotVar, density, bcTypes, bcValues, &
       invscale)
  
  !!!! IMPORTANT!!!! Actually, this is wrong, as a look at Driver_evolve
  !!!! shows. - JW
  ! FLASH solves the potential at the present gas density field, then calculates
  ! the evolution. So at the end of an evolve step the stored potential field
  ! DOES NOT MATCH the current density field of the gas. Therefore we FIRST
  ! must calculate the correct potential for the current gas density field.
  
  ! Note this means this kick step should go first (sinks -> gas kick)
  ! because the gas -> sinks kick maps the current acceleration to the sinks
  ! which is in turn calculated from the SAVED potential. -JW
  
  ! So here we calculate the current potential field so that it is updated
  ! properly. - JW

!    call Grid_solvePoisson( GPOT_VAR, DENS_VAR, bcTypes, bcValues, &
!         invscale)
!    !firstcall = .false.
!    call Grid_notifySolnDataUpdate( (/GPOT_VAR/) )

       
  ! Get the difference between the gas potential and the gas+particle potential.
  ! and store that so we have only the potential due to the sink particles.
  ! We later use this to difference for the acceleration sink -> gas.
  
!  do lb=1, blockCount
!      call Grid_getBlkPtr(blocklist(lb), solnVec)
!      solnVec(newPotVar,:,:,:) = solnVec(newPotVar,:,:,:) - &
!          solnVec(GPOT_VAR,:,:,:)
      
!  end do
       
  call Grid_notifySolnDataUpdate( (/newPotVar/) )

! Note this call shouldn't be needed since Particles_UpdateGridVar
! zeros out the value in PDEN before it calculates a new one... -JW
! But we do it anyway just to be sure. - JW

! Un-junk PDEN if it exists and if requested.

!#ifdef PDEN_VAR
!  if (grav_unjunkPden) then
!     density = PDEN_VAR
!#ifdef DENS_VAR           
!     do lb = 1, blockCount
!        call Grid_getBlkPtr(blocklist(lb), solnVec)
!        solnVec(density,:,:,:) = 0.0
!        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
!     enddo
!     call Grid_notifySolnDataUpdate( (/density/) )
!#endif
!  end if
!#endif

  if (.NOT. present(potentialIndex)) then
    ! Compute acceleration of the sink particles caused by gas and vice versa
    ! call Particles_sinkAccelGasOnSinksAndSinksOnGas()
  end if


#ifdef USEBARS
  call MPI_Barrier (grv_meshComm, ierr)
#endif  
  call Timers_stop ("gravity")
  
  return
end subroutine Gravity_potentialListOfBlocksPDENonly
