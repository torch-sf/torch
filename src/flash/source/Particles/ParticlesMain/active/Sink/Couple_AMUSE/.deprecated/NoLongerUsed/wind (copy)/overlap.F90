

!!! Overlap between a cell and a sphere or right cylinder.
!!! Adapted from the function of the same name written in Fortran 77
!!! by David Clarke, originally included in ZEUS-MP.

!!! This function does a "Monte Carlo" like integration with a
!!! uniform cartesian grid within a single cell.

!!! Joshua Wall
!!! Drexel University
!!! June 2016

!!! inputs:
!!! integer :: ishp [1 for sphere, 2 for right cylinder] 
!!! real    :: rad  [ radius of the sphere/cylinder]
!!! real    :: obj_center [ center location of the sphere/cylinder]
!!! real    :: cell_bot [ bottom/left points of cell in x,y,z]
!!! real    :: cell_top [ top/right points of cell in x,y,z]
!!! integer :: nstep [ number of points in x,y,z to sample in the cell]

!!! output:
!!! real    :: overlap [the fractional volume overlap of the object and cell]

subroutine overlap(ishp, rad, center, cell_bot, cell_top, nsteps, overlap_vol)

  implicit none

  integer, parameter :: dp = kind(1.d0)

  integer, intent(in)     :: ishp, nsteps
  real(dp), intent(in)        :: rad, center(3), cell_bot(3), cell_top(3)
  real(dp), intent(out)       :: overlap_vol

  real(dp)                    :: dx, dy, dz, r, factor
  real(dp), dimension(nsteps) :: xsq, ysq, zsq
  integer                 :: i, j, k
  integer                 :: inside_count

! Sphere sets factor = 1.0, right cylinder sets to 0.0.

  if (ishp .eq. 1) then
    factor = 1.0_dp
  else
    factor = 0.0_dp
  end if

! Get the deltas and positions to sample across the cell.

  dx = (cell_top(1) - cell_bot(1)) / real(nsteps)
  dy = (cell_top(2) - cell_bot(2)) / real(nsteps)
  dz = (cell_top(3) - cell_bot(3)) / real(nsteps)
  
  do i=1, nsteps

    xsq(i) = (cell_bot(1) + (0.5_dp+real(i-1))*dx - center(1))**2.0_dp
    ysq(i) = (cell_bot(2) + (0.5_dp+real(i-1))*dy - center(2))**2.0_dp
    zsq(i) = (cell_bot(3) + (0.5_dp+real(i-1))*dz - center(3))**2.0_dp
    
  end do
  
! Now sample these to see if they are inside the object.

  inside_count = 0

  do i=1, nsteps
    do j=1, nsteps
      do k=1, nsteps
      
        r = sqrt(factor*xsq(i) + ysq(j) + zsq(k))

        if (r .le. rad) inside_count = inside_count + 1
        
      end do
    end do
  end do
  
  ! Fractional volume (unitless).
  overlap_vol = real(inside_count) / real(nsteps**3.0_dp)
  !if (overlap_vol .gt. 0.0) print*, "YO! Overlap is = ", overlap_vol
  ! Actual volume (in physical units).
  ! overlap_vol =  overlap_vol * (cell_top(1) - cell_bot(1))*(cell_top(2) - cell_bot(2))*(cell_top(3) - cell_bot(3))
  
end subroutine overlap
