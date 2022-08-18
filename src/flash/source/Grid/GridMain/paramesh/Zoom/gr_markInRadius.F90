!!****if* source/Grid/GridMain/paramesh/gr_markInRadius
!!
!! NAME
!!  gr_markInRadius
!!
!!  
!! SYNOPSIS 
!!  gr_markInRadius(real(in) :: ic, 
!!                  real(in) :: jc,
!!                  real(in) :: kc, 
!!                  real(in) :: radius, 
!!                  integer(in) :: lref) 
!!  
!! PURPOSE 
!!  Refine all blocks containing points within a circular/spherical region of
!!  given radius about a given point (xc,yc,zc).  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.  
!!  
!! ARGUMENTS 
!!  ic -   Center of the interval/circle/sphere : IAXIS
!!  jc -                                          JAXIS
!!  kc -                                          KAXIS
!!               (Coordinates for nonexistent dimensions are ignored.)
!!  radius -       Radius of the region 
!!  lref  -        If > 0, bring all qualifying blocks to this level of refinement.
!!                 If <= 0, refine qualifying blocks once.
!!  
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!  
!!  
!!***

subroutine gr_markInRadius(ic, jc, kc, radius, lref)

!-------------------------------------------------------------------------------
  use tree, ONLY : refine, derefine, lrefine, bsize, coord, lnblocks, nodetype
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry
#include "constants.h"
#include "Flash.h"
  implicit none

! Arguments

  real, intent(IN)      :: ic, jc, kc, radius
  integer, intent(IN)   :: lref

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: bxl, bxr, byl, byr, bzl, bzr
  real                  :: dist2, xdist2, ydist2, zdist2
  integer               :: b
!  integer  :: num_refined, num_derefined, at_lref
  integer :: nsteps ! Number of steps since last derefinement.
                    ! Useful for preventing derefinement too quickly
                    ! in an area (especially for restarts). -JW
  logical, save :: derefined_last_step = .false.
  logical       :: derefined_this_step
  integer, parameter :: derefine_steps = 5, & ! # of steps between derefinement
                        max_total_nsteps = 50 ! total # of steps to slowly derefine over
  integer,save  :: total_nsteps = 0           ! total simulation steps taken so far
  
!  num_refined = 0
!  num_derefined = 0
!  at_lref = 0
  
  derefined_this_step = .false.
  nsteps = nsteps + 1

  if((gr_geometry == CARTESIAN).or.(gr_geometry == CYLINDRICAL)) then
     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = 0.5*bsize(:,b)
           
           bxl = blockCenter(1) - blockSize(1) - ic
           bxr = blockCenter(1) + blockSize(1) - ic
           if (NDIM > 1) then
              byl = blockCenter(2) - blockSize(2) - jc
              byr = blockCenter(2) + blockSize(2) - jc
           else
              byl = 0.
              byr = 0.
           endif
           if ((NDIM == 3).and.(gr_geometry==CARTESIAN)) then
              bzl = blockCenter(3) - blockSize(3) - kc
              bzr = blockCenter(3) + blockSize(3) - kc
           else
              bzl = 0.
              bzr = 0
           endif
           
! Find minimum distance from (ic,jc,kc) for each dimension.  For each
! coordinate, if both "left" and "right" distances have the same sign,
! then the smaller magnitude is the minimum.  Otherwise (ic,jc,kc) is
! contained within the interval for that dimension, so the minimum is 0.
! Nonexistent dimensions have had all distances set to zero, so they are
! ignored.

           if (bxl*bxr > 0.) then
              xdist2 = min( bxl**2, bxr**2 )
           else
              xdist2 = 0.
           endif
           if (byl*byr > 0.) then
              ydist2 = min( byl**2, byr**2 )
           else
              ydist2 = 0.
           endif
           if (bzl*bzr > 0.) then
              zdist2 = min( bzl**2, bzr**2 )
           else
              zdist2 = 0.
           endif

!!! Modified by me. -JW
! Now compute the minimum distance to (ic,jc,kc) and compare it to the
! specified radius.  If it is less than this radius, then the block contains
! at least part of the interval/circle/sphere and nothing is done.
! If it is outside this radius it is marked for derefinement. - JW

           dist2 = xdist2 + ydist2 + zdist2    ! Currently assumes Cartesian
           ! or 2D axisymmetric (r-z)
           ! or 1D spherical (r)
           if (dist2 > radius**2) then
              
              if (lrefine(b) > lref ) then
! Never refine outside the radius. - JW
                 refine(b)   = .false.
! If we derefined last time, don't do it this time. - JW
! If we haven't taken at least nsteps, don't derefine. - JW
! Note the number of steps should be a parameter for the 
! par file, not hard coded at 5. I'll get to it eventually! - JW
                 if ((derefined_last_step .or. (nsteps .lt. derefine_steps)) & 
                     .and. (total_nsteps .lt. max_total_nsteps)) then
                   derefine(b) = .false.
                 else
                   derefine(b) = .true.
                   derefined_this_step = .true.
                 !num_derefined = num_derefined + 1
                 end if
              else if (lrefine(b) == lref) then
                 refine(b)   = .false.
                 !at_lref = at_lref + 1
              !else if (lref <= 0) then
              !   refine(b) = .true.
              endif
              
           else
              !num_refined = num_refined + 1
           
           endif
           
           ! End of leaf-node block loop
        endif
     end do
  elseif((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL)) then

     do b = 1, lnblocks
        if(nodetype(b) == LEAF) then
           blockCenter(:) = coord(:,b)
           blockSize(:) = bsize(:,b)
           
           bxl = blockCenter(1) - blockSize(1) - ic
           bxr = blockCenter(1) + blockSize(1) - ic
           
           if (bxl*bxr > 0.) then
              dist2 = min( bxl, bxr )
           else
              dist2 = 0.
           endif
           
           if (dist2 > radius) then
              
              if (lrefine(b) > lref ) then
                 refine(b)   = .false.
                 derefine(b) = .true.
              else if (lrefine(b) == lref) then
                 refine(b) = .false.
              !else if (lref <= 0) then
              !   refine(b) = .true.
              endif
              
           endif
           
           
        endif
     end do
  else
     call Driver_abortFlash("MarkRefine: geometry spec is wrong")
     !-------------------------------------------------------------------------------
  end if
  
!  print*, "Number of derefined blocks = ", num_derefined
!  print*, "Number of blocks inside radius not touched = ", num_refined
!  print*, "Number of blocks at lowest refinement already = ", at_lref

! If we derefined something, reset the counter for steps since derefinement. - JW
  if (derefined_last_step) nsteps = 0
! Keep track of the total number of steps we've taken.
  total_nsteps = total_nsteps + 1
  
  return
end subroutine gr_markInRadius
