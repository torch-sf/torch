!!****if* source/Grid/GridMain/paramesh/gr_markOutRectangle
!!
!! NAME
!!  gr_markOutRectangle
!!
!! SYNOPSIS
!!  gr_markOutRectangle(real (in) :: ilb, 
!!                             real (in) :: irb, 
!!                             real (in) :: jlb, 
!!                             real (in) ::jrb, 
!!                             real (in) ::klb, 
!!                             real (in) ::krb, 
!!                        integer (in) :: lref, 
!!                        integer (in) :: contained)
!!
!! PURPOSE
!!  Refine blocks containing any points outside of a given rectangular 
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

subroutine gr_markOutRectangle(ilb, irb, jlb, jrb, klb, krb, lref, contained)

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
  logical               :: x_out_rect, y_out_rect, z_out_rect

#ifdef DEBUG
  if((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL))&
       call Driver_abortFlash("markRefineOutRectangle : wrong geometry")
  if((gr_geometry==CYLINDRICAL).and.(NDIM==3))&
       call Driver_abortFlash("markRefineOutRectangle : not valid in 3d for cylindrical")
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
        
           
           x_out_rect =   ((xl <= ilb) .and. (xr <= ilb)) .or. ((xl >= irb) .and. (xr >= irb))
           
           if (NDIM >= 2) then
           	y_out_rect =   ((yl <= jlb) .and. (yr <= jlb)) .or. ((yl >= jrb) .and. (yr >= jrb))
           endif
           
           if (NDIM == 3) then
           	z_out_rect =   ((zl <= klb) .and. (zr <= klb)) .or. ((zl >= krb) .and. (zr >= krb))
           endif
           

        
        ! Refine the block if all of the dimensions overlap/are contained.
        
        if (x_out_rect .or. y_out_rect .or. z_out_rect) then
	    if (lrefine(b) > lref) then
              print*,"Derefining block in [gr_markOutRectangle]"
              refine(b)   = .false.   
              derefine(b) = .true.
           endif
	    if (lrefine(b) == lref) then
              refine(b)   = .false.   
              derefine(b) = .false.
   	    endif	      
           
        endif
        
        ! End of leaf-node block loop
        
     endif
  end do
  
  !-------------------------------------------------------------------------------
  
  return
end subroutine gr_markOutRectangle
