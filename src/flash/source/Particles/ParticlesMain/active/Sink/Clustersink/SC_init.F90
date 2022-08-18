subroutine SC_init(restart)

  use SC_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

	logical, intent(IN) :: restart

  	call RuntimeParameters_get('SF_delay',sc_SF_delay)
  	call RuntimeParameters_get('SFE',     sc_SFE)

	! Exponent for the IMF of the broken power law.
	call RuntimeParameters_get('IMFalpha1', sc_IMFalpha1)
	call RuntimeParameters_get('IMFalpha2', sc_IMFalpha2)

	! Mass bounds for a SN explosion.
	call RuntimeParameters_get('IMF_SN_low', sc_IMFSNCut1)
	call RuntimeParameters_get('IMF_SN_high', sc_IMFSNCut2)

	! Mass cuts for the IMF.
	call RuntimeParameters_get('IMF_lowMass',  sc_IMFMassLim1)
	call RuntimeParameters_get('IMF_intMass',  sc_IMFMassLim2)
	call RuntimeParameters_get('IMF_highMass', sc_IMFMassLim3)

	! Single star mass cut
	call RuntimeParameters_get('Single_star_cut', sc_singleStarMassCut)

	! Models for the stellar evolution
	! I'll probably drop the different models as I will keep only 1.
	call RuntimeParameters_get('atmos_model', sc_atmos_model)
	call RuntimeParameters_get('wind_model', sc_wind_model)
	call RuntimeParameters_get('iophot_switch', sc_iophot_switch)
	call RuntimeParameters_get('SNRate_switch', sc_SNRate_switch)
	call RuntimeParameters_get('Wind_switch', sc_wind_switch)
	call RuntimeParameters_get('cluster_Metallicity', sc_cluster_metal)
	call RuntimeParameters_get('cluster_spectra', sc_cluster_spectra)
	
	! Calculate the IMF k_prime and k constants from: continuity and completeness(total cluster mass).
	k_prime0 = sc_IMFMassLim2**( sc_IMFalpha1 - sc_IMFalpha2 )
	
	k0 = 1./( k_prime0/(2-sc_IMFalpha1) * ( sc_IMFMassLim2**(2-sc_IMFalpha1) - sc_IMFMassLim1**(2-sc_IMFalpha1) ) &
	&	   	+ 	    1/ (2-sc_IMFalpha2) * ( sc_IMFMassLim3**(2-sc_IMFalpha2) - sc_IMFMassLim1**(2-sc_IMFalpha2) ))

	
	! Calculate IMF probability distribution function constants kp1 and kp2 assuming continuity of the
	! PDF and normalized to 1
	kp2 = 1 / ( &
	& 		sc_IMFMassLim2**(sc_IMFalpha1 - sc_IMFalpha2) / (1 - sc_IMFalpha1)* 		  &
	& 		(sc_IMFMassLim2**(1 - sc_IMFalpha1) - sc_IMFMassLim1**(1 - sc_IMFalpha1))+	  &
	& 		(sc_IMFMassLim3**(1 - sc_IMFalpha2) - sc_IMFMassLim2**(1 - sc_IMFalpha2))/ 	  &
	& 		(1 - sc_IMFalpha2)	 														  &
	& 	)
	
	kp1 = kp2 * sc_IMFMassLim2**(sc_IMFalpha1 - sc_IMFalpha2)

	! Division between the probability distributions functions for the different power law slopes of the IMF
	X01 = kp1/(1-sc_IMFalpha1)*(sc_IMFMassLim2**(1-sc_IMFalpha1)-sc_IMFMassLim1**(1-sc_IMFalpha1))
	
	
! JCIM. sample should always be stochastic.
	! Method used to sample the IMF: 0-stochastic, 1-complete, 2-single star.
	!allocate(IMFsample(nsinks))
	!IMFsample(:) = 0 
	
! JCIM. I also dont need this.	
	! Number of sampled stellar populations in a particle.
	!allocate(SinkSample_status(nsinks))
	!Sinksample_status(:) = 0
	
! JCIM. I think I dont need this anymore. This information is saved in the sink particle.	
	!allocate(starMass(nsinks,MaxStochStars))
	!allocate(numStars(nsinks))
	!allocate(numSNStars(nsinks))
	!allocate(clusterMass(nsinks))
	!allocate(criticalMass(nsinks))
	!allocate(criticalup(nsinks))
	!
	!starMass(:,:)   = 0
	!numStars(:) 	= 0
	!numSNStars(:)	= 0
	!clusterMass(:) 	= 0
	!criticalMass(:) = IMF_HighMass + 1

! JCIM. Do I need this ?
	! SN table, holds all the massive stars.
	!allocate(SNTable(nsinks, 100, 2))   	! 100 -> just guessing a number for the upper limit of SNs.
	!SNTable(:,:,:)	= 0	
	
	!allocate(freqRad(numSBWaves))
	!allocate(ioRadOrig(numSBWaves))
	
! JCIM. should uncomment for testing.
	
	!! Arrays used to test the results
	!allocate(ioRadNum(nsinks))
	!allocate(Winde(nsinks))
	!allocate(SNEnergy (nsinks))
	!allocate(photHeatE(nsinks))
	!allocate(avgphotE (nsinks))
	!allocate(iophotE (nsinks))
	!allocate(numSN(nsinks))
	!
	!ioRadNum(:) = 0
	!Winde(:) = 0
	!SNEnergy (:) = 0
	!photHeatE(:) = 0
	!avgphotE (:) = 0
	!iophotE (:) = 0
	!numSN(:) = 0
	
	! binned imf.
	!massiveIMFOnly = use_massiveIMFOnly
	!IMFbins = numIMFbins

! JCIM. I dont need this.
	! This saved all the sampled stars.
	!allocate(StarMassArray(nsinks,IMFbins, 2)) !should be initialized to the number of bins desired (bins num should be set in param).
	!StarMassArray(:,:,:)  = 0
	
	
	!! Calculate the IMF k_prime and k constants from: continuity and completeness(total cluster mass).
	!k_prime0 = IMFMassLim(2)**( IMFxponent(1) - IMFxponent(2) )
	!
	!k0 = 1./( k_prime0/(2-IMFxponent(1)) * ( IMFMassLim(2)**(2-IMFxponent(1)) - IMFMassLim(1)**(2-IMFxponent(1)) ) &
	!&	   	+ 	    1/ (2-IMFxponent(2)) * ( IMFMassLim(3)**(2-IMFxponent(2)) - IMFMassLim(1)**(2-IMFxponent(2)) ))
    !
	!
	!! Calculate IMF probability distribution function constants kp1 and kp2 assuming continuity of the
	!! PDF and normalized to 1
	!kp2 = 1 / ( &
	!& 		IMFMassLim(2)**(IMFxponent(1) - IMFxponent(2)) / (1 - IMFxponent(1))* &
	!& 		(IMFMassLim(2)**(1 - IMFxponent(1)) - IMFMassLim(1)**(1 - IMFxponent(1)))+	  &
	!& 		(IMFMassLim(3)**(1 - IMFxponent(2)) - IMFMassLim(2)**(1 - IMFxponent(2)))/ 	  &
	!& 		(1 - IMFxponent(2))	 														  &
	!& 	)
	!
	!kp1 = kp2 * IMFMassLim(2)**(IMFxponent(1) - IMFxponent(2))
    !
	!! Division between the probability distributions functions for the different power law slopes of the IMF
	!X01 = kp1/(1-IMFxponent(1))*(IMFMassLim(2)**(1-IMFxponent(1))-IMFMassLim(1)**(1-IMFxponent(1)))
	

	! Allocate the stellar evolution tables.
	! call get_SB_tables()

end subroutine SC_init
