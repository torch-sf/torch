!!****if* source/Grid/GridMain/paramesh/gr_markInRectangle
!!
!! NAME
!!  gr_markInRectangle
!!
!! SYNOPSIS
!!  gr_markInRectangle(real (in) :: ilb, 
!!                             real (in) :: irb, 
!!                             real (in) :: jlb, 
!!                             real (in) ::jrb, 
!!                             real (in) ::klb, 
!!                             real (in) ::krb, 
!!                        integer (in) :: lref, 
!!                        integer (in) :: contained)
!!
!! PURPOSE
!!  Refine blocks containing any points within a given rectangular 
!!  region having lower left coordinate (xlb,ylb,zlb) and upper right
!!  coordinate (xrb,yrb,zrb).  "Rectangular" is interpreted on a
!!  dimension-by-dimension basis: the region is an interval, rectangle,
!!  or rectangular parallelipiped in 1/2/3D Cartesian geometry; the
!!  rectangular cross-section of a rectangular torus in 2D axisymmetric
!!  (r-z) cylindrical geometry; an annular wedge in 2D polar (r-theta)
!!  cylindrical geometry; or an annulus in 1D spherical (r) geometry.
!!  Either blocks are brought up to a specific level of refinement or
!!  each block is refined once.
!!
!! ARGUMENTS
!!   ilb - lower bounding coordinate along i 
!!   irb - upper bounding coordinate along i 
!!   jlb - lower bounding coordinate along j
!!   jrb - upper bounding coordinate along j 
!!   klb - lower bounding coordinate along k 
!!   krb - upper bounding coordinate along k 
!!   lref -   If > 0, bring all qualifying blocks to this level of refinement.
!!            If <= 0, refine qualifying blocks once.
!!   contained - If /= 0, refine only blocks completely contained within
!!               the rectangle; otherwise refine blocks with any overlap.
!!
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!  
!!
!!
!!***

subroutine gr_markInRectangle(ilb, irb, jlb, jrb, klb, krb, lref, contained)

!-------------------------------------------------------------------------------

  use tree, ONLY: refine, derefine, lrefine, bsize, coord, nodetype, lnblocks
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_data, ONLY : gr_geometry


#include "constants.h"
#include "Flash.h"

  implicit none

! Arguments

  real, intent(IN)    :: ilb, irb, jlb, jrb, klb, krb
  integer, intent(IN) :: lref, contained

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: xl, xr, yl, yr, zl, zr
  integer               :: b
  logical               :: x_in_rect, y_in_rect, z_in_rect

#ifdef DEBUG
  if((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL))&
       call Driver_abortFlash("markRefineInRectangle : wrong geometry")
  if((gr_geometry==CYLINDRICAL).and.(NDIM==3))&
       call Driver_abortFlash("markRefineInRectangle : not valid in 3d for cylindrical")
#endif

  do b = 1, lnblocks
     if (nodetype(b) == LEAF) then
        
        blockCenter = coord(:,b)
        blockSize  = 0.5 * bsize(:,b)
        
        xl = blockCenter(1) - blockSize(1)
        xr = blockCenter(1) + blockSize(1)
        if (NDIM > 1) then
           yl = blockCenter(2) - blockSize(2)
           yr = blockCenter(2) + blockSize(2)
        endif
        if (NDIM == 3) then
           zl = blockCenter(3) - blockSize(3)
           zr = blockCenter(3) + blockSize(3)
        endif
        
        ! For each dimension, determine whether the block overlaps the specified
        ! rectangle.  Nonexistent dimensions are ignored.  This method assumes
        ! Cartesian coordinates (or the cross-section of a rectangular torus in
        ! 2D axisymmetric coordinates, or an annulus in 1D spherical coordinates).
        
        if (contained /= 0) then       ! only refine if completely contained
           
           x_in_rect =   ((xl >= ilb) .and. (xr <= irb))
           
           if (NDIM >= 2) then
              y_in_rect = ((yl >= jlb) .and. (yr <= jrb))
           else
              y_in_rect = .true.
           endif
           
           if (NDIM == 3) then
              z_in_rect = ((zl >= klb) .and. (zr <= krb))
           else
              z_in_rect = .true.
           endif
           
        else                           ! refine if any overlap with rectangle
           
           x_in_rect =   .not. ((xr <= ilb) .or. (xl >= irb))
           
           if (NDIM >= 2) then
              y_in_rect = .not. ((yr <= jlb) .or. (yl >= jrb))
           else
              y_in_rect = .true.
           endif
           
           if (NDIM == 3) then
              z_in_rect = .not. ((zr <= klb) .or. (zl >= krb))
           else
              z_in_rect = .true.
           endif
           
        endif

!!! I'll let FLASH do the regular refinement inside the box.
!!! Note this means this is called after gr_markRefineDerefine. - JW
        
        ! Refine the block if all of the dimensions overlap/are contained.
        
!        if (x_in_rect .and. y_in_rect .and. z_in_rect) then
!           
!           if (lrefine(b) < lref ) then
!              refine(b)   = .true.
!              derefine(b) = .false.
!           else if (lrefine(b) == lref) then
!              derefine(b) = .false.
!           else if (lref <= 0) then
!              refine(b) = .true.
!           endif
!           
!        endif

        if (.not. x_in_rect .and. .not. y_in_rect .and. &
          & .not. z_in_rect) then

!!! If you're not in the box and more refined than lref, derefine.
        
            if (lrefine(b) > lref)  then
            
              derefine(b) = .true.
              refine(b)   = .false.

!!! If you're not in the box and at lref, don't allow refinement.
!!! Note we do allow derefinement however.
              
            else if (lrefine(b) == lref) then
            
!              derefine(b) = .false.
              refine(b)   = .false.
!              stay(b)     = .true.
              
            end if
            
        end if
        
        ! End of leaf-node block loop
        
     endif
  end do
  
  !-------------------------------------------------------------------------------
  
  return
end subroutine gr_markInRectangle
