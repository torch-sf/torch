module gaus

contains

subroutine gaussian(x, gaus, mu_in, var_in)

#include "constants.h"

implicit none

!real, dimension(:), intent(in) :: x
real, intent(in) :: x
real, optional, intent(in) :: mu_in, var_in ! Center of gaus, variance
!real, dimension(:), allocatable, intent(out) :: gaus
real, intent(out) :: gaus

!real, dimension(:), allocatable :: y
real :: mu, var, y
real :: gaus_factor, norm_factor ! 1 / sqrt(2*pi*sig**2), normalization factor

!allocate(gaus(size(x)))
!allocate(y(size(x)))

!if (present(mu_in)) then
  mu = mu_in
!else
!  mu = sum(x)/real(size(x))
!end if

!if (present(var_in)) then
  var = var_in
!else
!  var = sum(x**2.0) / real(size(x)) 
!  var = var - mu**2.0
!end if

y = x - mu

gaus_factor = 1 / sqrt(2.0*PI*var)

gaus = gaus_factor * exp(-(y**2.0 / (2.0*var))) ! Calculate gaussian
!norm_factor = 1 / sum(gaus) ! Find normalization
!gaus = norm_factor*gaus ! Answer

!deallocate(y)

end subroutine gaussian

end module gaus
