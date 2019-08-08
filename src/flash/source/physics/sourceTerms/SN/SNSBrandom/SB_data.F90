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

  implicit none

  !this should not be initialised like that, but from nsnmax-nsnmin
	!this should be at least sb_nsnmax-sb_nsnmin, 100 should be enough 
  ! make number of maximum SN allocatable
  real, DIMENSION(100), save	:: sb_edge

	!sb creation rate, vertical stratification 
  real,save	:: sb_tsb, sb_hstarb, sb_MaxV
  logical,save	:: useSB,  sb_trackV, sb_useSBrandom

! sb for life, yo
  real,save	:: sb_life
  integer,save	:: sb_nsnmax, sb_nsnmin, sb_nSN, sb_SBmax

!  integer,parameter :: sbMax = 1000

! number of simultaneous sb sn is 1000
! could also allocate this in SB_init as number of simultaneous SB 

! allocatable is not persistant enough in scope
  real, save, allocatable, dimension(:,:) :: sbPOS
!  real, save,  dimension(1:4,1:1000) :: sbPOS

! buffer for MPI
  real, save, allocatable, dimension(:,:) :: sbPOSbuff
!  real, save,  dimension(1:4,1:1000) :: sbPOSbuff

end module SB_data
