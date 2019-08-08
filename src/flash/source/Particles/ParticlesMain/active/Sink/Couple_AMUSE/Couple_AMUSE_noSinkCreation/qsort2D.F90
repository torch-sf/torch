!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! A routine to quickly sort a matrix by row (low to high) based on
!!! comparisons of a specific column given by sort_column.
!!!
!!! Based on Nick Lomuto's scheme, as implemented by Jon Bentley.
!!! Found in numerous locations online, including the quicksort wiki
!!! article.
!!!
!!! matrix      : The matrix, which will be sorted in place.
!!! sort_column : The column to use to sort.
!!! lo          : Starting index to sort from.
!!! hi          : Ending index to sort to.
!!!
!!! Josh Wall, 11-2015, Drexel University
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module qsort2D

contains

recursive subroutine qsortRows(matrix, sort_column, lo, hi)
implicit none
real, dimension(:,:), intent(inout) :: matrix
integer, intent(in) :: sort_column, lo, hi

integer :: element, p

element = sort_column

if (lo .lt. hi) then

  call partition(matrix, element, lo, hi, p)
  call qsortRows(matrix, element, lo, p-1)
  call qsortRows(matrix, element, p+1, hi)

end if

end subroutine qsortRows

subroutine partition(matrix, element, lo, hi, i)
implicit none
real, dimension(:,:), intent(inout) :: matrix
integer, intent(in) :: element, lo, hi
integer, intent(out) :: i
integer :: j
real    ::  pivot

  pivot = matrix(hi, element)
  i = lo
  
  do j=lo, hi-1
  
    if (matrix(j, element) .le. pivot) then
    
        call swap(matrix(j,:), matrix(i,:))
        i = i + 1
    
    end if
  
  end do
  
  call swap(matrix(i,:),matrix(hi,:))
  
end subroutine partition

subroutine swap(a, b)
implicit none
real, dimension(:), intent(inout) :: a,b
real, dimension(size(a)) :: temp

temp = a
a    = b
b    = temp

end subroutine swap

end module qsort2D
