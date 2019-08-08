!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement.
!!
!! ARGUMENTS
!!
!!  none
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_eosModeNow
  use tree, ONLY : newchild, refine, derefine, stay, nodetype
!!$  use physicaldata, ONLY : force_consistency
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Particles_interface, only: Particles_sinkMarkRefineDerefine

!!! For my modified refinement only in the box/sphere I specify. -JW
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_Data, ONLY : dr_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
  
!!! For my modified refinement only in the box I specify. -JW
  logical, save :: first_call = .true.
  logical, save :: use_zoom
!!! Boundary of box where refinement is allowed.
  real, dimension(7), save :: specs
  integer, save :: l_max_outside, num_specs, zoom_criterion
!!! zoom_type = 1 for rectangle, 2 for sphere.  
  character(len=15), save :: zoom_type

!!! Check and see if we are zooming. If so, pick the proper
!!! type and then get the info for Grid_markRefineSpecialized.

  if (first_call) then
	
    call RuntimeParameters_get("use_zoom", use_zoom)
    call RuntimeParameters_get("zoom_type", zoom_type)
    
    zoom_type = trim(zoom_type)

    if (use_zoom .and. (zoom_type == "RECTANGLE")) then
		
		zoom_criterion = RECTANGLE
		num_specs = 7

		call RuntimeParameters_get("refine_xmin", specs(1))
		call RuntimeParameters_get("refine_xmax", specs(2))
		call RuntimeParameters_get("refine_ymin", specs(3))
		call RuntimeParameters_get("refine_ymax", specs(4))
		call RuntimeParameters_get("refine_zmin", specs(5))
		call RuntimeParameters_get("refine_zmax", specs(6))
		call RuntimeParameters_get("outer_lmax", l_max_outside)
		call RuntimeParameters_get("only_in_box", specs(7))

		if (dr_globalMe .eq. MASTER_PE) &
           print*, "Using rectangle zoom. Boundaries are:", &
           specs(1), specs(2), specs(3), specs(4), specs(5), specs(6) 
		

    else if (use_zoom .and. (zoom_type == "INRADIUS")) then
		
		zoom_criterion = INRADIUS
		num_specs = 4

		call RuntimeParameters_get("x_center", specs(1))
		call RuntimeParameters_get("y_center", specs(2))
		call RuntimeParameters_get("z_center", specs(3))
		call RuntimeParameters_get("sph_radius", specs(4))
		call RuntimeParameters_get("outer_lmax", l_max_outside)

		if (dr_globalMe .eq. MASTER_PE) &
           print*, "Using sphere zoom. Center and radius are:", &
           specs(1), specs(2), specs(3), specs(4)

	else if (use_zoom .and. ((zoom_type /= "RECTANGLE") .and. &
	                         (zoom_type /= "INRADIUS"))) then
		call Driver_abortFlash("Grid_markRefineDerefine: Zoom set to true but zoom_type /= RECTANGLE or INRADIUS.")

	end if

    first_call = .false.

  end if


  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.
!!$  gcMask(NUNK_VARS+1:maskSize) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

!!$  force_consistency = .FALSE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
     gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
     
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! For PARAMESH2, call gr_markRefineDerefine here if it hasn't been called above.
  ! This is necessary to make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  
  call Particles_sinkMarkRefineDerefine()
   
  !!! Now apply refinement only in the box I specify. Note
  !!! this routine DEREFINES everything outside the box/sphere
  !!! and allows normal refinement inside. This results in
  !!! "zooming in." - JW
  
  if (use_zoom) &
  call Grid_markRefineSpecialized(zoom_criterion, num_specs, specs, l_max_outside)

  ! When the flag arrays are passed to Paramesh for processing, only leaf
  ! blocks should be marked. - KW
  where (nodetype(:) .NE. LEAF)
     refine(:)   = .false.
     derefine(:) = .false.
  end where
  
  return
end subroutine Grid_markRefineDerefine

