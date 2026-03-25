!! Some functions regarding the LW photodissociation of H2.
module rt_lwmodule

#include "Flash.h"
#include "constants.h"

implicit none

contains

!
!===============================================================================
!
function get_fshield_H2(NH2,bfive)

!
! Returns the H2 self-shielding function
! Eq 12 of Wolcott-Green, Haiman and Bryan 2011: note this is slightly different from DB function
  implicit none
  real, intent(in) :: NH2, bfive
  real :: get_fshield_H2

  get_fshield_H2 = 0.965/(1+(NH2/(5.e14*bfive)))**1.1 + &
                   0.035/((1. + (NH2/5.e14))**0.5) * exp(-8.5 * 1.e-4 * (1. + (NH2/5.e14))**0.5)
  return
end function get_fshield_H2

function get_fshield_H(NH,bfive)
!
! Returns the self-shielding due to Lyman-alpha lines on the LW band
! Eq 15 of Wolcott-Green, Haiman and Bryan 2011
  implicit none
  real, intent(in) :: NH, bfive
  real :: get_fshield_H

  get_fshield_H = max( 1/(1+(NH/(2.85e23)))**1.6 * exp(-0.15 * (NH/(2.85e23))), 1e-15 )
  return
end function get_fshield_H

function get_fshield_C(NH2,NC,bfive)
!
! Returns the shielding factor for C using the treatment of Tielens & Hollenbach (1985)
! Eq 9 in Gong, Ostriker & Wolfire 2017
  implicit none
  real, intent(in) :: NH2, NC, bfive
  real :: get_fshield_C

  get_fshield_C = exp(-NC * 1.6e-17) * exp(-NH2 * 2.8e-22)/(1 + (2.8e-22 * NH2))
  return
end function get_fshield_C

!*****************************
  !2D interpolation at (x0,y0) for (x(:), y(:)) in z(:,:)
  !Added by Piyush Sharda in 2024 for CO shielding
function interpolate2DCO(x, y, z, x0, y0)
  implicit none
  real*8 :: interpolate2DCO
  real*8 :: x(:), y(:), z(size(x), size(y))
  real*8 :: x0, y0
  real*8 :: f
  integer :: i, j
  real*8 :: t, u

  ! Find indices i and j such that x(i) <= x0 < x(i+1) and y(j) <= y0 < y(j+1)
  i = 1
  do while (i < size(x) - 1 .and. x0 > x(i + 1))
    i = i + 1
  end do

  j = 1
  do while (j < size(y) - 1 .and. y0 > y(j + 1))
    j = j + 1
  end do

  ! Compute interpolation weights
  t = (x0 - x(i)) / (x(i + 1) - x(i))
  u = (y0 - y(j)) / (y(j + 1) - y(j))

  ! Perform bilinear interpolation
  interpolate2DCO = (1 - t) * (1 - u) * z(i, j) + t * (1 - u) * z(i + 1, j) + &
      (1 - t) * u * z(i, j + 1) + t * u * z(i + 1, j + 1)

end function interpolate2DCO

function get_fshield_CO(NH2,NCO,bfive)
!
! Returns the shielding factor for CO using tabulated data from Visser et al. 2009, compiled by Gong et al. 2017
! Tabulated data procured from: https://github.com/munan/pdr/blob/master/shielding.cpp
  implicit none
  real, intent(in) :: NH2, NCO, bfive
  real :: get_fshield_CO
  real :: x(8), y(6), z(8, 6), clipped_x,clipped_y

  x = (/ 0d0, 13d0, 14d0, 15d0, 16d0, 17d0, 18d0, 19d0 /) !N_CO
  y = (/ 0d0, 19d0, 20d0, 21d0, 22d0, 23d0 /) !N_H2
  z(:, 1) = (/ 1d0, 8.080d-1, 5.250d-1, 2.434d-1, 5.467d-2, 1.362d-2, 3.378d-3, 5.240d-5 /)
  z(:, 2) = (/ 8.176d-1, 6.347d-1, 3.891d-1, 1.787d-1, 4.297d-2, 1.152d-2, 2.922d-3, 4.662d-4 /)
  z(:, 3) = (/ 7.223d-1, 5.624d-1, 3.434d-1, 1.540d-1, 3.515d-2, 9.231d-3, 2.388d-3, 3.899d-4 /)
  z(:, 4) = (/ 3.260d-1, 2.810d-1, 1.953d-1, 8.726d-2, 1.907d-2, 4.768d-3, 1.150d-3, 1.941d-4 /)
  z(:, 5) = (/ 1.108d-2, 1.081d-2, 9.033d-3, 4.441d-3, 1.102d-3, 2.644d-4, 7.329d-5, 1.437d-5 /)
  z(:, 6) = (/ 3.938d-7, 3.938d-7, 3.936d-7, 3.923d-7, 3.901d-7, 3.893d-7, 3.890d-7, 3.875d-7 /)

  !Clip logNCO and logNH2 to the ranges in the data
  clipped_x = max(x(1), min(log10(NCO), x(8)))
  clipped_y = max(y(1), min(log10(NH2), y(6)))
  get_fshield_CO = 1d1**interpolate2DCO(x, y, log10(z), clipped_x, clipped_y)
  return
end function get_fshield_CO

function get_Epump(T, nH2, nH)
!
! Returns the heat deposited per UV pumping event -- depends on temperature
!
  implicit none
  real, intent(in) :: T, nH2, nH
  real :: get_Epump

  get_Epump = 2 * 1.6021764620000066e-12 * get_Cdex(T,nH2,nH)/(get_Cdex(T,nH2,nH) + 2.e-7)
  return


end function get_Epump

function get_Cdex(T, nH2, nH)
!
! Returns the quantity Cdex (Eq 47 of Baczynski 15) which represents the collisional de-excitation rate
!
  implicit none
  real, intent(in) :: T, nH2, nH
  real :: get_Cdex

  get_Cdex = 1.e-12 * (1.4 * exp(-18100/(T+1200)) * nH2 + exp(-1000/T)*nH) * SQRT(T)
  return
end function get_Cdex

end module rt_lwmodule
