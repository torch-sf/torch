!!****if* source/Particles/ParticlesMain/Particles_data
!!
!! NAME
!!    Particles_data
!!
!! SYNOPSIS
!!    Particles_data()
!!
!! DESCRIPTION
!!    Module to hold local variables and data types for Particles unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!                                    in the setup.
!!    pt_maxPerProc  INTEGER [1000]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!                                   don't move farther than one block in each step
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation (EstiMidpoint):
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!    pt_numParticlesWanted  INTEGER [100]  for withDensity -- requested number of particles per block
!!
!!***
!*******************************************************************************


module Particles_data
!===============================================================================

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#include "GridParticles.h"
!-------------------------------------------------------------------------------

!! The particles data structure is an real array of size (NPART_PROPS,MAX_PARTICLES_PER_PROCESSOR)
!! It is allocated at runtime, as MAX_PARTICLES_PER_PROCESSOR is a runtime parameter = pt_maxPerProc
   real, save, allocatable, dimension(:,:), target :: particles


  integer, save   :: pt_numLocal, pt_numLost
  logical, save   :: pt_keepLostParticles

! Run-time parameters, also described in Particles_init and above under PARAMETERS

  logical, save   :: useParticles = .true.   !if including Particles in simulation, default is true
  real, save      :: pt_dtFactor              ! a multiplying factor for time step limits
  real, save      :: pt_small                 ! a small number for velocity comparisons
  real, save      :: pt_dtChangeTolerance     ! percentage change threshold that
                                              ! controls when some schemes fall back to Euler

! Run-time parameters, local copies defined from some other Unit

  real, save      :: pt_xmin, pt_xmax, pt_ymin, pt_ymax, pt_zmin, pt_zmax ! from Grid/common
!                                  physical domain lower and upper bound in x/y/z
  integer, save    :: pt_geometry                            ! from Grid/common, an integer and string
  character(len=MAX_STRING_LENGTH), save :: pt_str_geometry  ! designation of the geometry configuration

! Run-time parameters, for grid initialization
  integer, save      :: pt_maxPerProc  ! for Lattice


! Guard Cell masks, initialized in Particles_init
  logical, save      :: pt_gcMaskForAdvance(NUNK_VARS) !This one is for Particles_advance
  logical, save      :: pt_gcMaskForWrite(NUNK_VARS+NDIM*NFACE_VARS) !This one for Particles_updateAttributes
  integer,save       :: pt_gcMaskSizeForAdvance=NUNK_VARS
  integer,save       :: pt_gcMaskSizeForWrite=NUNK_VARS+NDIM*NFACE_VARS
! Local variable to control timestepping
  logical, save      :: pt_restart

  integer,save         :: pt_globalMe, pt_globalNumProcs, pt_globalComm
  integer,save         :: pt_meshMe, pt_meshNumProcs, pt_meshComm
  logical, save        :: pt_posInitialized, pt_velInitialized, pt_resetTag
  integer,save         :: pt_logLevel
  integer,save         :: pt_numAtOnce
  integer,dimension(PT_MAX_ATTRIBUTES),save :: pt_attributes
  integer,dimension(PART_ATTR_DS_SIZE,PT_MAX_ATTRIBUTES), save :: pt_meshVar
  integer, save        :: pt_numAttributes
  integer, save        :: pt_velNumAttrib, thermNumAttrib
  integer, dimension(PART_ATTR_DS_SIZE,MDIM),save :: pt_velAttrib,pt_velPredAttrib
  integer, dimension(PART_ATTR_DS_SIZE,2) :: therm_Attrib
  integer, dimension(MDIM),save :: pt_posAttrib,pt_posPredAttrib
  integer, dimension(PART_TYPE_INFO_SIZE,NPART_TYPES),save :: pt_typeInfo
  


!! Paramters needed only with Lattice initialization
  integer, save      :: pt_numX, pt_numY, pt_numZ
  real, save         :: pt_initialXMin, pt_initialXMax, &
       pt_initialYMin, pt_initialYMax, pt_initialZMin, pt_initialZMax, &
       pt_initialRadius

!! Parameter needed only for With Density initialization

  integer, save      :: pt_pRand, pt_numParticlesWanted 
  real, save         :: pt_totalMass, pt_totalVolume, pt_averageDensity
  integer, dimension(GRPT_ALL),save      :: pt_indexList
  integer, save      :: pt_indexCount

!! Parameters needed only with NParPerZone initialisation
  integer, save      :: pt_parPerSide
  real, save, allocatable, dimension(:,:,:) :: pt_tracerPosX
  real, save, allocatable, dimension(:,:,:) :: pt_tracerPosY
  real, save, allocatable, dimension(:,:,:) :: pt_tracerPosZ

  logical, save :: pt_reduceGcellFills

!! Parameter needed for generating the unique consecutive tags in parallel
  integer, save :: pt_startTagNumber=0
end module Particles_data
