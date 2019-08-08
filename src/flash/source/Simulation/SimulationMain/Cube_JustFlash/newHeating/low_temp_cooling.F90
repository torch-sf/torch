!! Low temperature cooling
!! Juan Ibañez-Mejia 2014

!!****f* 
!! NAME
!!  
!!  low_temp_cooling
!!
!!
!! SYNOPSIS
!! 
!!  call low_temp_cooling(rho, energy, Z, G0, cool)
!! 
!! 
!! DESCRIPTION
!!
!!	Given the local volume density, temperature, metallicity and Interstellar ratiadion field strength
!!  in terms of the Draine's field, the function returns the specific cooling rate. This parameterized
!!  curve accounts for collisional excitation of C+ with H2 and a transition to CO cooling at higher
!!	densities. It does not account for dust cooling at densities highet than 10^4.5 cm^-3.
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine low_temp_cooling(rho, temp, Z, cool)
!	subroutine low_temp_cooling(rho, energy, Z, cool)
!
	implicit none

	! Input: rho    -- mass density    [cgs units]
	! Input: energy -- specific energy [cgs units]
	! Input: Z      -- metallicity     [solar units, 1 ==> solar metallicity]
	!
	! Output: lambda -- cooling rate per unit mass [cgs units]

	real, intent(IN)	:: rho, temp, Z
	real, intent(OUT)	:: cool

	real :: n, ntot, n3, yn, LTE_factor
	real :: CII_fraction, CII_cooling, CO_cooling

	
	! Assumptions made in this routine:
	! 1) Fully molecular hydrogen
	! 2) He:H ratio same as in local ISM
	! 3) gamma = 5/3
	
	! Declare local constants 
	REAL, PARAMETER :: abhe  = 0.1
	REAL, PARAMETER :: gamma = 5.0 / 3.0
	REAL, PARAMETER :: mp    = 1.6726e-24
	REAL, PARAMETER :: kb    = 1.38066e-16
	
	! Lower threshold for cooling.
	if(temp .le. 10.0) then
		cool = 1.0e-35
		return
	end if
	  
	! Compute number density of H nuclei
      n  = rho / ((1d0 + 4d0 * abhe) * mp)
      n3 = n / 1d3

	! Compute total particle number density
      ntot = (0.5d0 + abhe) * n

	! number density of hydrogen nuclei
	  yn = rho / (1.4 * mp)

	! Compute temperature
    !  temp = (rho * energy) / ((gamma - 1d0) * ntot * kb)

	! Cooling from CII excitation
	! - Account for transition from low-density rate to LTE rate at n ~ 2d3 cm^-3
      LTE_factor   = 2d0 / (2d0 + n3)  

	! - Account for conversion from C+ to CO at high densities
	!   N.B. At metallicities << 1, this should be modified so that conversion 
	!   occurs at higher density -- see Glover & Clark (2012, MNRAS, 426, 327)
      CII_fraction = 1d0 / (1d0 + n3**2)

	! CII cooling rate -- assumes collisions with HI dominate, neglects minor 
	! temperature dependence of CII-HI collisional de-excitation rate
      CII_cooling = 1.4d-27 * exp(-9.2d1 / temp) * Z * CII_fraction * LTE_factor * yn**2

	! CO cooling rate
	!  - Assumes all CO in LTE
	!  - Neglects opacity, freeze-out
	!  - Approximately reproduces rates shown in Neufeld et al (1995, ApJS, 100, 132),
	!    at least for temperatures T = 10-20 K; may become inaccurate at much higher T
      CO_cooling = 1.4d-30 * temp**3 * exp(-6d0 / temp) * Z * (1d0 - CII_fraction) * yn
	
	! Total cooling rate [erg/(s cm^3)]
      cool = CII_cooling + CO_cooling

	! convert to specific cooling rate [erg/(s g)].
      cool = cool / rho

!	if(temp .le. 100.0) then
!		PRINT*,'cooling rate is:', cool
!	endif

      return
end
