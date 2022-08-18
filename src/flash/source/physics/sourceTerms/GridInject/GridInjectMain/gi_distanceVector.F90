! get the signed components of distance vector from
! (xloc,yloc,zloc) to (xx,yy,zz)
! accounting for periodic boundary conditions
!
! not guaranteed to work with 1-D or 2-D.
! assumes Cartesian coordinates

subroutine gi_distanceVector (xx, yy, zz, xloc, yloc, zloc, dx, dy, dz)

  use Grid_data, ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax, &
                        gr_domainBC

  implicit none
#include "constants.h"

  real, intent(IN) :: xx, yy, zz, xloc, yloc, zloc
  real, intent(OUT) :: dx, dy, dz
  real :: dplus, dminus
  integer, dimension(LOW:HIGH,MDIM) :: boundary

  call Grid_getDomainBC(boundary)

  dx = xx-xloc

  if (boundary(LOW,IAXIS)==PERIODIC .and. boundary(HIGH,IAXIS)==PERIODIC) then
    dplus  = xx-(xloc+gr_imax-gr_imin)
    dminus = xx-(xloc-gr_imax+gr_imin)
    if (abs(dplus) < abs(dx)) then
      dx = dplus
    end if
    if (abs(dminus) < abs(dx)) then
      dx = dminus
    end if
  end if

  dy = yy-yloc

  if (boundary(LOW,JAXIS)==PERIODIC .and. boundary(HIGH,JAXIS)==PERIODIC) then
    dplus  = yy-(yloc+gr_jmax-gr_jmin)
    dminus = yy-(yloc-gr_jmax+gr_jmin)
    if (abs(dplus) < abs(dy)) then
      dy = dplus
    end if
    if (abs(dminus) < abs(dy)) then
      dy = dminus
    end if
  end if

  dz = zz-zloc

  if (boundary(LOW,KAXIS)==PERIODIC .and. boundary(HIGH,KAXIS)==PERIODIC) then
    dplus  = zz-(zloc+gr_kmax-gr_kmin)
    dminus = zz-(zloc-gr_kmax+gr_kmin)
    if (abs(dplus) < abs(dz)) then
      dz = dplus
    end if
    if (abs(dminus) < abs(dz)) then
      dz = dminus
    end if
  end if

  return
end subroutine gi_distanceVector
