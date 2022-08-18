module Particles_sort

contains

subroutine sortInd(array, index_array, ascending)

real, dimension(:), intent(inout)    :: array
!integer, intent(in)                  :: sort_axis ! which column to sort rows on

logical, intent(in), optional        :: ascending 

integer, dimension(:), intent(inout) :: index_array ! returned sorted index array.

logical                              :: asc_order
integer                              :: num_rows, row, new_min_row
real                                 :: buffer

if (present(ascending)) then
  asc_order = ascending
else
  asc_order = .true.
end if

num_rows = size(array)

do row=1, num_rows
    index_array(row) = row
end do

if (asc_order) then

do row=1, num_rows

    new_min_row = minloc(array(row:num_rows), dim=1) + row - 1
    
    ! Now pivot the rows to move this one to the proper location
    buffer              = array(row)
    array(row)          = array(new_min_row)
    array(new_min_row)  = buffer
    
    buffer              = index_array(row)
    index_array(row)    = index_array(new_min_row)
    index_array(new_min_row) = buffer
    
    !print*, "min_row =", new_min_row
    
end do

else

do row=1, num_rows

    new_min_row = maxloc(array(row:num_rows), dim=1) + row - 1
    
    ! Now pivot the rows to move this one to the proper location
    buffer              = array(row)
    array(row)          = array(new_min_row)
    array(new_min_row)  = buffer
    
    buffer              = index_array(row)
    index_array(row)    = index_array(new_min_row)
    index_array(new_min_row) = buffer
    
    !print*, "min_row =", new_min_row
    
end do

end if

end subroutine sortInd

end module Particles_sort
