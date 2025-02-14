!!***if* source/Grid/GridMain/paramesh/Grid_markDerefineSpecialized
!!
!! NAME
!!  Grid_markDerefineSpecialized
!!
!! SYNOPSIS
!!  Grid_markDerefineSpecialized(integer(IN) :: criterion,
!!                              integer(IN) :: size,
!!                              real(IN)    :: specs(size),
!!                              integer(IN) :: lref )
!!
!! DESCRIPTION
!!  The routine provides an interface to a collection of routines
!!  that define very specialized Derefinement criteria. The currently
!!  supported options are:
!!   
!!   RECTANGLE   : The blocks that fall within the specified rectangle
!!
!! ARGUMENTS
!!  criterion - the creterion on which to Derefine
!!  size      - size of the specs data structure
!!  specs     - the data structure containing information specific to 
!!              the creterion
!!              For RECTANGLE
!!                 specs(1:6) = bounding coordinates of rectangle
!!                 specs(7)   = if 0 Derefine block with any overlap
!!                              if /= Derefine only blocks fully
!!                              contained in the rectangle
!!
!!
!!  lref      - If > 0, bring selected blocks to this level of Derefinement.
!!              If <= 0, Derefine qualifying blocks once.
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

subroutine Grid_markDerefineSpecialized (criterion,size,specs,lref)

  implicit none
#include "constants.h"
! Arguments

  integer, intent(IN) :: criterion
  integer, intent(IN) :: size
  real,dimension(size),intent(IN) :: specs
  integer, intent(IN) ::  lref
  integer :: contained, var, icmp

  select case (criterion)
     case(RECTANGLE)
        contained=int(specs(7))
        call gr_markOutRectangle(specs(1),specs(2),specs(3),specs(4),&
             specs(5),specs(6),lref,contained)
  end select
  return
end subroutine Grid_markDerefineSpecialized
