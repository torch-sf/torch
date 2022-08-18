!! New Parameterized heating.
!! Juan Ibañez-Mejia 2014

!!****f* 
!! NAME
!!  
!!  radloss
!!
!!
!! SYNOPSIS
!! 
!!  old:
!!  call low_temp_heating(rho, Z, G0, heat)
!! 
!!	Now:
!!	Parameterized_heating(zz, rho, Z, G0, pe_h, cr_h, heat)
!! 
!! DESCRIPTION
!!
!!		Parameterized heating from Simon Glover non-equilibrium chemistry simulations. Accounts for
!!		Dust shielding and molecular cooling with different species dominating at different densities.
!!		An empiric column density function is used to account for the dust shielding.
!!
!! ARGUMENTS
!!
!!		Input: zz		-- distance to the midplane [cm]
!!		Input: rho    	-- mass density [g cm-3]
!!		Input: energy 	-- specific energy [cgs units]
!!		Input: Z      	-- metallicity [solar units, 1 ==> solar metallicity]
!!		Input: G0     	-- interstellar radiation field strength [units of Draine field]
!!		  			   		Draine's field is 1.7 Habings = 1d-24
!!		Input: stratify -- Stratify the photoelectric heating and the cosmic ray heating ?
!!		Input: peh(2)	-- Photoelectric Heating stratification scale height.
!!		Input: crh(2)	-- Cosmic rays stratification scale height. 		   
!!		
!!		Output: lambda 	-- cooling rate per unit mass [cgs units]
!!
!!
!!***

subroutine Parameterized_heating(zz, rho, Z, G0, stratify, pe_h, cr_h, heat)

	implicit none

	real, intent(IN)	:: zz, rho
	real, intent(IN)	:: Z, G0
	real, intent(IN)	:: pe_h, cr_h 
	logical, intent(IN)	:: stratify
	real, intent(OUT)	:: heat

	real :: n, n3
	real :: cosmic, photoelectric
	real :: AVeff

	! Assumptions made in this routine:
	! 1) Fully molecular hydrogen
	! 2) He:H ratio same as in local ISM
	! 3) gamma = 5/3
	
	! Declare local constants 
	REAL, PARAMETER :: abhe  = 0.1
	REAL, PARAMETER :: gamma = 5.0 / 3.0
	REAL, PARAMETER :: mp    = 1.6726e-24
	REAL, PARAMETER :: kb    = 1.38066e-16
	

	! Compute number density of H nuclei
	  n  = rho / ((1d0 + 4d0 * abhe) * mp)
	  n3 = n / 1d3

	! Compute approximate visual extinction
	  AVeff = (Z * n3)**0.4


	! Photoelectric heating term
	  photoelectric = G0 * exp(-2.5d0 * AVeff) * Z * n
	
	! Cosmic ray heating rate -- assumes canonical cosmic ray ion. rate
	  cosmic = 6.4d-28 * n
	
	if(stratify) then
		! Stratify the photoelectric heating & the cosmic rays
		photoelectric = photoelectric * exp(-abs(zz) / pe_h)
		cosmic = cosmic * exp(-abs(zz)/cr_h)
	else
		! dont stratify the heating
		photoelectric = photoelectric
		cosmic = cosmic
	endif

	! Total heating rate [erg/(s cm^3)]
	  heat = photoelectric + cosmic
	
	! convert to specific heating rate [erg/(s g)].
	  heat = heat / rho
	
	return
end