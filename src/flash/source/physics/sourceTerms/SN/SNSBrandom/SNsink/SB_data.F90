!!****if* source/physics/sourceTerms/Heat/SN+Heat+Cool/Heat_data
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!
!!
!!
!! DESCRIPTION
!!
!!
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

!!TODO READ IN RHO AMBIENT AND PRESSURE AMBIENT FOR HEATING ROUTING

module SB_data

#include "Flash.h"

   use SN_data, only: he_nsndt

  implicit none

  !this should not be initialised like that, but from nsnmax-nsnmin
  !this should be at least sb_nsnmax-sb_nsnmin, 100 should be enough 
  !make number of maximum SN allocatable
  real, save, DIMENSION(100) :: sb_edge

	!sb creation rate, vertical stratification 
  real,    save	:: sb_tsb, sb_hstarb, sb_MaxV
  logical, save	:: useSB,  sb_trackV, sb_useSBrandom, sb_accTracer, sb_useSNsink
  logical, save	:: sn_tracerAccFrac, sn_sinkBulkMotion

! sb for life
  real,    save	:: sb_life, sn_outflowFrac, sn_cloudMassPerStar, sn_SFdelay
  integer, save	:: sb_nsnmax, sb_nsnmin, sb_nSN, sb_SBmax
!  integer,save	:: tracerMemID

! number of simultaneous sb sn is 1000
  integer, parameter :: sbMax = 1000

! could also allocate this in SB_init as number of simultaneous SB 
!  real,allocatable, dimension(:,:) :: sbPOS
  real, save, dimension(1:4,1:1000) :: sbPOS 

! buffer for MPI
!  real,allocatable, dimension(:,:) :: sbPOSbuff
  real, save, dimension(1:4,1:1000) :: sbPOSbuff 

end module SB_data
