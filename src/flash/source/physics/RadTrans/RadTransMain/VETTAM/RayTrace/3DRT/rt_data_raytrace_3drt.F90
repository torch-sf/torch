!!****if* source/physics/RadTrans/RadTransMain/HybridChar3DRTFLD/RayTrace/3DRT/rt_data_raytrace_3drt
!!
!!  NAME 
!!    rt_data
!!
!!  SYNOPSIS
!!    use rt_data
!!
!!  DESCRIPTION 
!!    Stores data for HybridChar3DRTFLD
!!
!!***

#include "Flash.h"

module rt_data_raytrace_3drt
  implicit none
  
  ! The direction (phi,theta) in cartesian coordinates
  ! for the characteristic
  !
  real,    save       :: domainSizeX,domainSizeY,domainSizeZ

  ! Values for the use of solid angle groups. The number of
  ! angle groups defines for how many angles the face values
  ! are computed before they are communicated with an
  ! mpi_allgather command. Less communication is faster,
  ! but is also very memory consuming. One has to keep this
  ! in mind and choose the number of angle groups very carefully.
  !
  integer, save       :: nrOfAnglesPerGroup
  integer, save       :: nrOfAngles ! = rt_nPhi*rt_nTheta

  !
  ! The solid angle element. Should be constant and depends
  ! on the method of angular discretization.
  !
  real,    save       :: dOmega ! = 4pi/nrOfAngles
  !
  ! The following values are used, if periodic boundary conditions are
  ! invoked. We basically copy the outgoing radiation of the opposite 
  ! boundary and use it as incoming intensities. This can be done several
  ! times (e.g. shallow polar rays in the planar atmosphere problem).
  ! TODO: This is very buggy yet, also doesn't work together with angle groups.
  ! TODO: Use independen boundary conditions for radiation 
  !
  integer, save       :: nrOfBoundIter
  logical, save       :: x_periodic = .false. 
  logical, save       :: y_periodic = .false. 
  logical, save       :: z_periodic = .false.

end module rt_data_raytrace_3drt
