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
