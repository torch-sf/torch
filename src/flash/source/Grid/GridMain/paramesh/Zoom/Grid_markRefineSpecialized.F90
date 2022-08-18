!!***if* source/Grid/GridMain/paramesh/Grid_markRefineSpecialized
!!
!! NAME
!!  Grid_markRefineSpecialized
!!
!! SYNOPSIS
!!  Grid_markRefineSpecialized(integer(IN) :: criterion,
!!                              integer(IN) :: size,
!!                              real(IN)    :: specs(size),
!!                              integer(IN) :: lref )
!!
!! DESCRIPTION
!!  The routine provides an interface to a collection of routines
!!  that define very specialized refinement criteria. The currently
!!  supported options are:
!!
!!   THRESHOLD   : when a specific variable is below or above a
!!                 threshold
!!   ELLIPSOID   : The blocks that fall within the specified ellipsoid
!!   RECTANGLE   : The blocks that fall within the specified rectangle
!!   INRADIUS    : The blocks that fall within the specified radius
!!   WITHRADIUS  : The blocks that fall on the specified radius
!!
!! ARGUMENTS
!!  criterion - the creterion on which to refine
!!  size      - size of the specs data structure
!!  specs     - the data structure containing information specific to 
!!              the creterion
!!              For THRESHOLD 
!!                 specs(1) = real(variable_name), for example
!!                           if variable is density, then
!!                           specs(1)=real(DENS_VAR)
!!                 specs(2) = the threshold value
!!                 specs(3) = if < 0 refine if variable < threshold
!!                            if > 0 refine if variable > threshold
!!              For ELLIPSOID
!!                 specs(1:3) = center of the ellipsoid
!!                 specs(4:6) = the semimajor axes of the ellipsoid
!!              For INRADIUS
!!                 specs(1:3) = center of the circle/sphere
!!                 specs(4)   = the radius
!!              For WITHRADIUS
!!                 specs(1:3) = center of the circle/sphere
!!                 specs(4)   = the radius
!!              For RECTANGLE
!!                 specs(1:6) = bounding coordinates of rectangle
!!                 specs(7)   = if 0 refine block with any overlap
!!                              if /= refine only blocks fully
!!                              contained in the rectangle
!!
!!
!!  lref      - If > 0, bring selected blocks to this level of refinement.
!!              If <= 0, refine qualifying blocks once.
!!
!! NOTES
!! 
!!  This collection of routines has not been tested well and can be
!!  used as a guideline for a user's implementation.
!!
!!  Non-Cartesian geometries may not be supported in the default
!!  implementations of geometric criteria; the level of support depends
!!  on the routine that implements a given criterion.
!!
!!***

subroutine Grid_markRefineSpecialized (criterion,size,specs,lref)

  implicit none
#include "constants.h"
! Arguments

  integer, intent(IN) :: criterion
  integer, intent(IN) :: size
  real,dimension(size),intent(IN) :: specs
  integer, intent(IN) ::  lref
  integer :: contained, var, icmp

  select case (criterion)
     case(THRESHOLD)
        var = int(specs(1))
        icmp = int(specs(3))
        call gr_markVarThreshold(var,specs(2),icmp,lref)
     case(ELLIPSOID)
        call gr_markEllipsoid(specs(1),specs(2),specs(3),specs(4),&
             specs(5),specs(6),lref)
     case(INRADIUS)
        call gr_markInRadius(specs(1),specs(2),specs(3),specs(4),lref)
     case(WITHRADIUS)
        call gr_markWithRadius(specs(1),specs(2),specs(3),specs(4),lref)
     case(RECTANGLE)
        contained=int(specs(7))
        call gr_markInRectangle(specs(1),specs(2),specs(3),specs(4),&
             specs(5),specs(6),lref,contained)
  end select
  return
end subroutine Grid_markRefineSpecialized
