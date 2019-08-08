!!****if* source/Grid/GridMain/paramesh/gr_markInRectangle
!!
!!  
!! special refinement for the stratified box setup
!! fixed to 5 levels
!! starts at highest supplied refined level and decreases by one in each lower resolution grid
!! if dynamic range of refinements is too low, lower resolution levels are truncated to lowest 
!! refinement level
!! allows for dynamic refinement levels 
!! origin has to be in box center
!! TODO clean up dimensional checks (not needed)
!! TODO clean up x and y checks (not needed)
!!
!!***

subroutine gr_markInZ(limits,lref)

!-------------------------------------------------------------------------------

  use tree, ONLY: refine, derefine, lrefine, bsize, coord, nodetype, lnblocks, lrefine_max, lrefine_min
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_geometry,gr_imin,gr_imax,gr_jmin,gr_jmax

#include "constants.h"
#include "Flash.h"

  implicit none

! Arguments

  real, dimension(7) ,intent(IN)    :: limits
  integer, intent(IN) :: lref
  integer	:: i

! Local data

  real, dimension(MDIM) :: blockCenter, blockSize
  real                  :: xl, xr, yl, yr, zl, zr
  integer               :: b
  logical               :: x_in_rect, y_in_rect, z_in_rect

	real 									:: xlowb1,ylowb1,zlowb1
	real 									:: xhighb1,yhighb1,zhighb1

	real 									:: xlowb2,ylowb2,zlowb2l,zlowb2u
	real 									:: xhighb2,yhighb2,zhighb2l,zhighb2u

	real 									:: xlowb3,ylowb3,zlowb3l,zlowb3u
	real 									:: xhighb3,yhighb3,zhighb3l,zhighb3u

	real 									:: xlowb4,ylowb4,zlowb4l,zlowb4u
	real 									:: xhighb4,yhighb4,zhighb4l,zhighb4u

	real 									:: xlowb5,ylowb5,zlowb5l,zlowb5u
	real 									:: xhighb5,yhighb5,zhighb5l,zhighb5u

	real 									:: xlowb6,ylowb6,zlowb6l,zlowb6u
	real 									:: xhighb6,yhighb6,zhighb6l,zhighb6u

	real 									:: xlowb7,ylowb7,zlowb7l,zlowb7u
	real 									:: xhighb7,yhighb7,zhighb7l,zhighb7u

  integer               :: lref1,lref2,lref3,lref4,lref5,lref6,lref7,refdiff

#ifdef DEBUG
  if((gr_geometry==POLAR).or.(gr_geometry==SPHERICAL))&
       call Driver_abortFlash("markRefineInRectangle : wrong geometry")
  if((gr_geometry==CYLINDRICAL).and.(NDIM==3))&
       call Driver_abortFlash("markRefineInRectangle : not valid in 3d for cylindrical")
#endif


! dim check
	if (NDIM .lt. 3) return

! this is ugly but meh
! 1 is most central region 
! 5 is outer most in z
! origin has to be in box center
!	xlowb1 	 = gr_imin*1.1
!	ylowb1 	 = gr_jmin*1.1
	zlowb1 	 = -limits(1)

!	xhighb1   = gr_imax*1.1
!	yhighb1   = gr_jmax*1.1
	zhighb1  =  limits(1)

	lref1 	 = lref
  if( lref1 .gt. lrefine_max) lref1 = lrefine_max

!	xlowb2 	 = gr_imin*1.1
!	ylowb2 	 = gr_jmin*1.1
!	xhighb2  = gr_imax*1.1
!	yhighb2  = gr_jmax*1.1

! zlowb2u < zhighb2u
	zlowb2u  = limits(1) 
	zhighb2u = limits(2)

! zlowb2l > zhighb2l
	zlowb2l  = -limits(2)
	zhighb2l = -limits(1)

	lref2 	 = lref1 - 1
  if( lref2 .lt. lrefine_min) lref2 = lrefine_min

!	xlowb3   = gr_imin*1.1
!	ylowb3   = gr_jmin*1.1
!	xhighb3  = gr_imax*1.1
!	yhighb3  = gr_jmax*1.1
 
	zlowb3u  =  limits(2)
	zhighb3u =  limits(3)
	zlowb3l  = -limits(3)
	zhighb3l = -limits(2)

	lref3 	 = lref2 - 1
  if( lref3 .lt. lrefine_min) lref3 = lrefine_min

!	xlowb4 	 = gr_imin*1.1
!	ylowb4 	 = gr_jmin*1.1
!	xhighb4  = gr_imax*1.1
!	yhighb4  = gr_jmax*1.1

	zlowb4u  =  limits(3)
	zhighb4u =  limits(4)
	zlowb4l  = -limits(4)
	zhighb4l = -limits(3)

	lref4 	 = lref3 - 1
  if( lref4 .lt. lrefine_min) lref4 = lrefine_min

!	xlowb5 	 = gr_imin*1.1
!	ylowb5 	 = gr_jmin*1.1
!	xhighb5  = gr_imax*1.1
!	yhighb5  = gr_jmax*1.1

	zlowb5u  =  limits(4)
	zhighb5u =  limits(5)
	zlowb5l  = -limits(5)
	zhighb5l = -limits(4)

	lref5 	= lref4 - 1
  if( lref5 .lt. lrefine_min) lref5 = lrefine_min

	zlowb6u  =  limits(5)
	zhighb6u =  limits(6)
	zlowb6l  = -limits(6)
	zhighb6l = -limits(5)

	lref6 	= lref5 - 1
  if( lref6 .lt. lrefine_min) lref6 = lrefine_min


	zlowb7u  =  limits(6)
	zhighb7u =  limits(7)
	zlowb7l  = -limits(7)
	zhighb7l = -limits(6)

	lref7 	= lref6 - 1
  if( lref7 .lt. lrefine_min) lref7 = lrefine_min


  do b = 1, lnblocks
     if (nodetype(b) == LEAF) then
        
        blockCenter = coord(:,b)
        blockSize  = 0.5 * bsize(:,b)
        
!        xl = blockCenter(1) - blockSize(1)
!        xr = blockCenter(1) + blockSize(1)

!        if (NDIM > 1) then
!           yl = blockCenter(2) - blockSize(2)
!           yr = blockCenter(2) + blockSize(2)
!        endif
        
        ! For each dimension, determine whether the block overlaps the specified
        ! rectangle.  Nonexistent dimensions are ignored.  This method assumes
        ! Cartesian coordinates (or the cross-section of a rectangular torus in
        ! 2D axisymmetric coordinates, or an annulus in 1D spherical coordinates).
        
! only refine if completely contained 
! zcut 1, inner contiguos region
!        x_in_rect =   ((xl >= xlowb1) .and. (xr <= xhighb1))
           
!        if (NDIM >= 2) then
!        	y_in_rect = ((yl >= ylowb1) .and. (yr <= yhighb1))
!        else
!        	y_in_rect = .true.
!        endif

! calculate blockcenter at target refinement 
				refdiff = 0

				if(lrefine(b) < lref1) then
					refdiff = lref1-lrefine(b)
				endif

				if(blockCenter(3) .gt. 0.0 ) then
					zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
					zr = zl
				else
					zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
					zr = zl
				endif
! check if it is inside
       	z_in_rect = (((zl >= zlowb1) .and. (zr <= zhighb1)))

!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
	      if (z_in_rect) then
					if (lrefine(b) < lref1 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref1) then
! do nothing
						derefine(b) = .false.
						refine(b)   = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
! don't derefine if higher refinement, allows for additional dynamic levels
					else if (lrefine(b) > lref1) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref1 <= 0) then
						refine(b) = .true.
      		endif
				endif

! zcut 2
!        x_in_rect =   ((xl >= xlowb2) .and. (xr <= xhighb2))
!           
!        if (NDIM >= 2) then
!        	y_in_rect = ((yl >= ylowb2) .and. (yr <= yhighb2))
!        else
!        	y_in_rect = .true.
!        endif
! calculate blockcenter at target refinement 
				refdiff = 0

				if(lrefine(b) < lref2) then
					refdiff = lref2-lrefine(b)
				endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif
! check if it is inside
        	z_in_rect = (((zl >= zlowb2u) .and. (zr <= zhighb2u))) .or. &
											(((zl >= zlowb2l) .and. (zr <= zhighb2l)))


!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
      if (z_in_rect) then
					if (lrefine(b) < lref2 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref2) then
						refine(b)   = .false.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref2) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref2 <= 0) then
						refine(b) = .true.
      		endif
      	endif

! zcut 3
!        x_in_rect =   ((xl >= xlowb3) .and. (xr <= xhighb3))
!           
!        if (NDIM >= 2) then
!        	y_in_rect = ((yl >= ylowb3) .and. (yr <= yhighb3))
!        else
!        	y_in_rect = .true.
!        endif

					refdiff = 0

					if(lrefine(b) < lref3) then
						refdiff = lref3-lrefine(b)
					endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif

! check if it is inside
        	z_in_rect = (((zl >= zlowb3u) .and. (zr <= zhighb3u))).or. &
											(((zl >= zlowb3l) .and. (zr <= zhighb3l)))

!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
	      if (z_in_rect) then
					if (lrefine(b) < lref3 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref3) then
						derefine(b) = .false.
						refine(b)   = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref3) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
	     		else if (lref3 <= 0) then
						refine(b) = .true.
      		endif
      	endif

! zcut 4
!        x_in_rect =   ((xl >= xlowb4) .and. (xr <= xhighb4))
           
!        if (NDIM >= 2) then
!        	y_in_rect = ((yl >= ylowb4) .and. (yr <= yhighb4))
!        else
!        	y_in_rect = .true.
!        endif
           
					refdiff = 0

					if(lrefine(b) < lref4) then
						refdiff = lref4-lrefine(b)
					endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif

! check if it is inside
        	z_in_rect = (((zl >= zlowb4u) .and. (zr <= zhighb4u))).or. &
											(((zl >= zlowb4l) .and. (zr <= zhighb4l)))

!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
	      if (z_in_rect) then
					if (lrefine(b) < lref4 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref4) then
						refine(b)   = .false.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref4) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref4 <= 0) then
						refine(b) = .true.
      		endif
      	endif
! zcut 5
!        x_in_rect =   ((xl >= xlowb5) .and. (xr <= xhighb5))
!           
!        if (NDIM >= 2) then
!        	y_in_rect = ((yl <= ylowb5) .and. (yr >= yhighb5))
!        else
!        	y_in_rect = .true.
!        endif
           
					refdiff = 0

					if(lrefine(b) < lref5) then
						refdiff = lref5-lrefine(b)
					endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif

! check if it is inside
        	z_in_rect = (((zl >= zlowb5u) .and. (zr <= zhighb5u))).or. &
											(((zl >= zlowb5l) .and. (zr <= zhighb5l)))

!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
	      if (z_in_rect) then
					if (lrefine(b) < lref5 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref5) then
						refine(b)   = .false.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref5) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref5 <= 0) then
						refine(b) = .true.
      		endif
      	endif

! calculate blockcenter at target refinement 
				refdiff = 0

				if(lrefine(b) < lref6) then
					refdiff = lref6-lrefine(b)
				endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif
! check if it is inside
        	z_in_rect = (((zl >= zlowb6u) .and. (zr <= zhighb6u))) .or. &
											(((zl >= zlowb6l) .and. (zr <= zhighb6l)))


!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
      if (z_in_rect) then
					if (lrefine(b) < lref6 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref6) then
						refine(b)   = .false.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref6) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref6 <= 0) then
						refine(b) = .true.
      		endif
      	endif

! calculate blockcenter at target refinement 
				refdiff = 0

				if(lrefine(b) < lref7) then
					refdiff = lref7-lrefine(b)
				endif

					if(blockCenter(3) .gt. 0.0 ) then
						zl = blockCenter(3)+(0.5**refdiff-1.0)*blockSize(3)
						zr = zl
					else
						zl = blockCenter(3)+(-(0.5**refdiff)+1.0)*blockSize(3)
						zr = zl
					endif
! check if it is inside
        	z_in_rect = (((zl >= zlowb7u) .and. (zr <= zhighb7u))) .or. &
											(((zl >= zlowb7l) .and. (zr <= zhighb7l)))


!	      if (x_in_rect .and. y_in_rect .and. z_in_rect) then
      if (z_in_rect) then
					if (lrefine(b) < lref7 ) then
						refine(b)   = .true.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
        	else if (lrefine(b) == lref7) then
						refine(b)   = .false.
						derefine(b) = .false.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lrefine(b) > lref7) then
						refine(b)	 = .false.
						derefine(b) = .true.
! this belongs to the local slab, cycle out so no other slab interferes
						cycle
      		else if (lref7 <= 0) then
						refine(b) = .true.
      		endif
      	endif

        ! End of leaf-node block loop        
     endif
  end do
  
  !-------------------------------------------------------------------------------
  
  return
end subroutine gr_markInZ
