!! Some functions regarding the reaction and cooling rate coefficients for ionisation and recombination of Hydrogen.
module rt_ionisemodule

#include "Flash.h"
#include "constants.h"

implicit none

contains

!
!===============================================================================
!
function get_recombination_coefficient(temp,on_the_spot)
!
! Returns the recombination coefficient. If on the spot approximation is switched on, then the Case-B coefficient is returned.
!
  implicit none
  
  real, intent(in) :: temp
  logical, intent(in) :: on_the_spot
  real :: get_recombination_coefficient

  if(on_the_spot) then
    get_recombination_coefficient = recombination_coefficient_B(temp)
  else
    get_recombination_coefficient = recombination_coefficient_A(temp)
  endif

  return

end function get_recombination_coefficient

!
!===============================================================================
!
function get_recombination_ground(temp,on_the_spot)
!
! Returns the recombination coefficient for transitions directly to the ground state (for diffuse emission)
!
  implicit none
  
  real, intent(in) :: temp
  logical, intent(in) :: on_the_spot
  real :: get_recombination_ground

  if(on_the_spot) then
    get_recombination_ground = 0.0
  else
    get_recombination_ground = recombination_coefficient_A(temp) - recombination_coefficient_B(temp)
  endif

  return

end function get_recombination_ground

!
!===============================================================================
!
elemental function recombination_coefficient_A(temp)
!
! Returns the Case-A recombination coefficient as a function of temperature
! Uses the value reported in Abel et al. 1997, which fits data from Ferland et al 1992. This is also used in KROME.
! This returns the value of the coefficient in unit of cm^3s^-1
!
  implicit none
  
  real, intent(in) :: temp
  real :: recombination_coefficient_A

  if(temp .lt. 5.5d3) then 
    recombination_coefficient_A = 3.92d-13*(1./(temp*8.617652504247027e-05))**0.6353d0
  else
    recombination_coefficient_A = exp(-28.61303380689232d0-0.7241125657826851d0*log(temp*8.617652504247027d-05)- &
            0.02026044731984691d0*log(temp*8.617652504247027d-05)**2d0- 0.002380861877349834d0*log(temp*8.617652504247027d-05) &
            **3d0-0.0003212605213188796d0*log(temp*8.617652504247027d-05)**4d0-0.00001421502914054107d0 &
            *log(temp*8.617652504247027d-05)**5d0 + 4.989108920299513e-6*log(temp*8.617652504247027d-05)**6d0+ &
            5.755614137575758e-7*log(temp*8.617652504247027d-05)**7d0-1.856767039775261e-8* &
            log(temp*8.617652504247027d-05)**8d0-3.071135243196595e-9*log(temp*8.617652504247027d-05)**9d0)
  endif

  return

end function recombination_coefficient_A

!===============================================================================
!
elemental function recombination_coefficient_B(temp)
!
! Returns the Case-B recombination coefficient as a function of temperature
! Uses the value reported in Hui & Gnedin (1997)
! This returns the value of the coefficient in unit of cm^3s^-1
!
  implicit none
  
  real, intent(in) :: temp
  real :: recombination_coefficient_B
  
  recombination_coefficient_B = 2.753d-14*(315614/temp)**(1.5)/(1.d0 + (315614/temp/2.74d0)**(0.407))**(2.242)

  return

end function recombination_coefficient_B

!
!===============================================================================
!
elemental function recombination_cool_coeff(temp)
!
! Returns the cooling rate coefficient for recombinations Gamma_rec. The cooling rate per vol is then Gamma_rec n_e n_Hplus
! Uses the value reported adopted in KROME which uses fits used in Cen+1992
!
  implicit none
  
  real, intent(in) :: temp
  real :: recombination_cool_coeff

  recombination_cool_coeff = 8.7d-27 * sqrt(temp) * (temp/1.d3)**(-0.2)/(1.d0 + (temp/1.d6)**0.7)

  return

end function recombination_cool_coeff

!
!===============================================================================
!
elemental function ff_cool_coeff(temp)
!
! Returns the cooling rate coefficient for free-free cooling Gamma_ff. The cooling rate per vol is then Gamma_ff n_e n_plus, where n_plus is all +ve ions
! Uses the value from Osterbrock & Ferland 2006, with a Gaunt factor correction of 1.5 (peak value) to be consistent with KROME
!
  implicit none
  
  real, intent(in) :: temp
  real :: ff_cool_coeff

  ff_cool_coeff = 1.42e-27*sqrt(temp)*1.5d0

  return

end function ff_cool_coeff



end module rt_ionisemodule


