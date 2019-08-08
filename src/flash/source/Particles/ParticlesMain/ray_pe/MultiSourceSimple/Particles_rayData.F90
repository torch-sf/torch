!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   data for ray tracing
!!

module Particles_rayData
!===============================================================================

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#include "GridParticles.h"
!-------------------------------------------------------------------------------

!! Parameters needed only with HEALPix initialization
  integer*4, save  :: ph_initHPlevel
!! HEALPIx tag arrays
  integer*4, save, dimension(0:1023) :: pix2x, pix2y
  !short int
  integer*2, save, dimension(128)    :: x2pix, y2pix
  !short int, with one integer offset (weird HEALPix thing)
  integer*2, save, dimension(0:127)  :: x2pix1, y2pix1

!! ray traversal stuff	
  real,save    :: ph_sampling
  real,save    :: ph_locSampling
  real,save    :: ph_periodicBoxL
  logical,save :: ph_inBlockSplit
  logical,save :: ph_rotRays
  logical,save :: ph_radPressure
  real, save :: speedoflight

! for mpi and sorting
  integer,save :: ph_maxNRays, ph_localRays, ph_radOutput
! default value for outflow boundaries
  integer,save :: ph_transProp
  integer,parameter :: ph_numIntProp = 6, ph_numRealProp = 11, ph_numProp = 16
  
! Radiation bin / photon type. Ion >= 13.6 eV, 5.6 < PE < 13.6 eV - JW
  integer, parameter :: ion_photon = 1, pe_photon = 2

! rot matrix, columns 1,2,3
  real,save    :: aa, ab, ac
  real,save    :: ba, bb, bc
  real,save    :: ca, cb, cc

! index into particle array, radiation properties, 10, real
  integer, save:: inion, ieion, irad, ihnum

! position and direction, 6, real
  integer, save:: iposx, iposy, iposz
  integer, save:: ivelx, ively, ivelz

! photon properties, 5, real
! sigh is low energy atomic H ionisation
! sigh2 is high energy H_2 ionisation
! sighi is high energy H ionisation
  integer, save:: isigh

! ray properties, 5, integer
  integer, save:: ihlev, istpd, isid, itype

! position in data structure, 2, integer
  integer, save :: iblk, iproc

  logical, save :: xperiodic, yperiodic, zperiodic, ph_periodic
  logical, save :: xhydroper, yhydroper, zhydroper, useRadTransfer
! for mixed BC
  integer, save :: ph_BCcase 

! array for local raytracing
  real,save,allocatable,dimension(:,:) :: raysRealProp

! array for local raytracing
  integer,save,allocatable,dimension(:,:) :: raysIntProp

! indices into transport buffers 
  integer, save :: itnion, ithnum, itinfo, itrad, itid, itblk, ittype

  integer, save :: iteion, itvelx, itvely, itvelz, itposx, itposy, itposz, itsigh, itstpd      

! global box sizes
  real,save    :: glDX, glDY, glDZ

! H2 dissociation parameters
!  real,save    :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10
!  real,save    :: e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e0, emax

! could use more deegeeets
  real, parameter :: evToErg = 1.60217657e-12
  real, parameter :: MBTocm2 = 1e-18

! Terminate the FUV bin rays if they are less than background FUV at local
! shielding level.
  logical, save :: early_term_FUV, ph_EUVonDust
! ray properties
!  integer, save	:: ph_numDest

! delay ionising radiation by some time to clear out mol. H
!  real,save :: ph_delayIon

! holds target processor ID
!  integer,save,allocatable,dimension(:) :: rayDestBufTarget

! arrays for ray MPI
!  real,save,allocatable,dimension(:,:) :: rayDestBuf, raySourceBuf

end module Particles_rayData
