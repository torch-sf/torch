! get the distance between two points,
! accounting for periodic boundary conditions
!
! a simple convenience wrapper for gi_distanceVector
!
! not guaranteed to work with 1-D or 2-D.
! assumes Cartesian coordinates

function gi_distance (xx, yy, zz, xloc, yloc, zloc)

  use gi_interface, ONLY : gi_distanceVector

  implicit none

  real, intent(IN) :: xx, yy, zz, xloc, yloc, zloc
  real :: dx, dy, dz, gi_distance

  call gi_distanceVector(xx, yy, zz, xloc, yloc, zloc, dx, dy, dz)

  gi_distance = sqrt(dx**2 + dy**2 + dz**2)

end function gi_distance
