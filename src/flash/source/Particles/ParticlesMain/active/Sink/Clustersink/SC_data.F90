!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkData
!!
!! NAME
!!
!!    Particles_sinkData
!!
!! SYNOPSIS
!!
!!    Particles_sinkData()
!!
!! DESCRIPTION
!!
!!    Module to hold local variables and data types for sink particle unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

module SC_data

  implicit none

#include "constants.h"
#include "Particles.h"
#include "Flash.h"

	real, save :: sc_SF_delay, sc_SFE
	
		!integer,save	:: accretionSwitch = 0

		!integer, save	:: numinterv
		!integer, save	:: numinterv2

		! Parameters of the IMF, exponents, mass cut for SN expl, low, cut and high mass limits.
		real, save	:: sc_IMFalpha1, sc_IMFalpha2
		real, save	:: sc_IMFSNCut1, sc_IMFSNCut2, 
		real, save 	:: sc_IMFMassLim1, sc_IMFMassLim2, sc_IMFMassLim3


		! Arrays containing the masses of the stars per sink particle, total cluster mass and total number of stars.]
		real,allocatable,dimension(:,:,:),save		:: SNTable
		real,allocatable,dimension(:,:),save		:: starMass
		real,allocatable,dimension(:,:,:),save		:: StarMassArray
		real,allocatable,dimension(:),save			:: clusterMass
		integer,allocatable,dimension(:),save		:: numStars, numSNStars

		! Arrays storing the mass of the last star that exploded as a SN.
		real,allocatable,dimension(:),save			:: criticalMass
		real,allocatable,dimension(:),save			:: criticalup

		! Conditions for the maximum mass in the cloud and maximum number of stars in the arary 
		! to be sampled stochastically
		real,save :: MaxStochMass
		integer,save :: MaxStochStars

		! Maximum mass in sink to be treated as a single star. Above this threshold value cluster sink
		real, save :: StarSinkMassCut

		! Switches to implement different models in SB 99.
		integer,save	:: atmos_model, wind_model
		integer,save	:: iophot_switch, SNRate_switch, StellarWind_switch, SynSpectra_switch
		integer,save	:: Metallicity

		! normalization constants to sample the IMF and
		! probability distribution constants for the stochastic sample.
		real,save	:: k_prime0, k0
		real,save	:: kp1, kp2, X01	

		! map StarBurst 99 into FLASH4 bins.
		integer, parameter	:: numSBWaves = 118
		integer, parameter	:: SB99bins = 3000

		! Starburst 99 ionizing Radiation.
		real,allocatable,dimension(:),save :: freqRad
		real,allocatable,dimension(:),save :: ioRadOrig


		! Arrays to test the output
		integer,	save	:: numFeedbackPoints
		real,allocatable,dimension(:),save	:: ioRadNum, WindE, SNEnergy, avgphotE, photHeatE, numSN, iophotE
	!	real,allocatable,dimension(:,:,:),save	:: ioRad, WindE, SNEnergy, avgphotE, photHeatE

		! Binned IMF.
		logical, save :: massiveIMFOnly  
		integer, save :: IMFbins	
	
	

end module
