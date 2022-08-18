!!****if* source/Grid/GridMain/Grid_coordTransfm
!!
!! NAME
!!
!!  Grid_coordTransfm
!!
!! SYNOPSIS
!!
!!  call Grid_coordTransfm(real(IN)    :: x,
!!                         real(IN)    :: y,
!!                         real(IN)    :: z,
!!                         real(OUT)   :: xout,
!!                         real(OUT)   :: yout,
!!                         real(OUT)   :: zout,
!!                OPTIONAL,integer(IN) :: geometryIn,
!!                OPTIONAL,integer(IN) :: geometryOut,
!!                OPTIONAL,integer(IN) :: ndim,
!!                OPTIONAL,real(IN)    :: velI,
!!                OPTIONAL,real(IN)    :: velJ,
!!                OPTIONAL,real(IN)    :: velK,
!!                OPTIONAL,real(OUT)   :: velIOut,
!!                OPTIONAL,real(OUT)   :: velJOut,
!!                OPTIONAL,real(OUT)   :: velKOut)
!!
!!
!! DESCRIPTION
!!
!!  Convert between Cartesian and other coordinates.
!!
!!  Supports the FLASH geometries CARTESIAN, CYLINDRICAL, SPHERICAL, and POLAR.
!!
!! ARGUMENTS
!!
!!  x - First input coordinate.
!!  y - Second input coordinate.
!!  z - Third input coordinate.
!!  xout - First Cartesian coordinate.
!!  yout - Second Cartesian coordinate.
!!  zout - Third Cartesian coordinate.
!!  geometryIn - What is the input geometry? (default is the FLASH geometry established by
!!                                            configuration and runtime parameter)
!!  geometryOut - What is the output geometry? (default is the FLASH geometry established by
!!                                            configuration and runtime parameter)
!!  ndim - dimensionality (default N_DIM)
!!  velI,velJ,velK - input velocity components
!!  velIOut,velJOut,velKOut - output velocity components
!!
!!
!!***

#include "constants.h"


subroutine Grid_coordTransfm(x,y,z, xout,yout,zout, geometryIn,geometryOut, ndimArg, &
     velI,velJ,velK,velIOut,velJOut,velKOut)
  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_data, ONLY: gr_meshMe, gr_geometry
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: xout,yout,zout
  integer,OPTIONAL,intent(IN) :: geometryIn
  integer,OPTIONAL,intent(IN) :: geometryOut
  integer,OPTIONAL,intent(IN) :: ndimArg
  real,OPTIONAL,intent(IN) :: velI,velJ,velK
  real,OPTIONAL,intent(OUT) :: velIOut,velJOut,velKOut

  logical,save :: first_call = .TRUE.
  integer :: geoIn, geoOut, ndim

  if (present(geometryIn)) then
     geoIn = geometryIn
  else
     geoIn = -1
  end if
  if (geoIn == -1) geoIn = gr_geometry
     
  if (present(geometryOut)) then
     geoOut = geometryOut
  else
     geoOut = -1
  end if
  if (geoOut == -1) geoOut = gr_geometry
     
  if (present(ndimArg)) then
     ndim = ndimArg
  else
     ndim = N_DIM
  end if

#ifdef DEBUG_GRID
  if (first_call) then
     print*,'geoIn=',geoIn
     print*,'geoOut=',geoOut
  end if
#endif

  if (present(velIOut).AND.present(velI)) velIOut = velI
  if (present(velJOut).AND.present(velJ)) velJOut = velJ
  if (present(velKOut).AND.present(velK)) velKOut = velK
  if (geoOut==geoIn) then
     if (gr_meshMe==MASTER_PE) then
        if (first_call) print*,'Input geometry ok.'
     end if
     xout = x
     yout = y
     zout = z
  else if (geoIn==CARTESIAN) then
     if (geoOut==SPHERICAL) then
        call gr_coordTransfmSph(x,y,z, xout,yout,zout, ndim)
     else if (geoOut==CYLINDRICAL) then
        call gr_coordTransfmCyl(x,y,z, xout,yout,zout, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
     else
        call Driver_abortFlash("[Grid_coordTransfm] invalid input geometry")
     end if
  else if (geoOut==CARTESIAN) then
     if (geoIn==CYLINDRICAL) then
        call gr_coordTransfmFromCyl(x,y,z, xout,yout,zout, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
     else
        call Driver_abortFlash("[Grid_coordTransfm] invalid output geometry")
     end if
  else
     call Driver_abortFlash("[Grid_coordTransfm] invalid geometry")
  end if

  first_call = .FALSE.


contains
!!****if* source/Grid/GridMain/gr_coordTransfmSph
!!
!! NAME
!!
!!  gr_coordTransfmSph
!!
!! SYNOPSIS
!!
!!  call gr_coordTransfmSph(real(IN)    :: x,
!!                        real(IN)    :: y,
!!                        real(IN)    :: z,
!!                        real(OUT)   :: r,
!!                        real(OUT)   :: theta,
!!                        real(OUT)   :: phi,
!!                        integer(IN) :: ndim)
!!
!!
!! DESCRIPTION
!!
!!  Convert from Cartesian to spherical coordinates.
!!
!! ARGUMENTS
!!
!!  x - First Cartesian coordinate.
!!  y - Second Cartesian coordinate.
!!  z - Third Cartesian coordinate.
!!  r - First spherical coordinate.
!!  theta - Second spherical coordinate.
!!  phi - Third spherical coordinate.
!!  ndim - dimensionality
!!
!!
!!***
subroutine gr_coordTransfmSph(x,y,z, r,theta,phi, ndim)
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: r,theta,phi
  integer,intent(IN) :: ndim

  real :: rsq

  Theta = PI*0.5
  phi = 0.0

  rsq = x**2
  if (ndim > 1) rsq = rsq + y**2
  if (ndim > 2) rsq = rsq + z**2
  r = sqrt(rsq)

  if (r == 0.0) return

  if (ndim > 1) then
     theta = acos(z/r)
  end if

  if (ndim > 2) then
     phi = atan2(y,x)
  end if

end subroutine gr_coordTransfmSph

subroutine gr_coordTransfmCyl(x,y,z, r,zout,phi, ndim, velI,velJ,velK,velRout,velZout,velPhiOut)
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: r,zout,phi
  integer,intent(IN) :: ndim
  real,intent(IN)   ,OPTIONAL :: velI,velJ,velK
  real,intent(INOUT),OPTIONAL :: velRout
  real,intent(INOUT),OPTIONAL :: velZout,velPhiOut

  real :: rsq, velSq, velIJ

  zout = 0.0
  phi = 0.0

  rsq = x**2
  if (ndim > 1) rsq = rsq + y**2
  r = sqrt(rsq)

  if (ndim > 1) then
     velSq = 0.0
     if (present(velI)) velSq = velI**2
     if (present(velJ)) velSq = velSq + velJ**2
     if (present(velK)) velSq = velSq + velK**2
     zout = z
     if (velSq.NE.0.0) then
!!$        velIJ = sqrt(velSq)
        if (r .NE. 0.0) then
           velRout                          = (velI*x + velJ*y) / r
           if(present(velPhiOut)) velPhiout = (velJ*x - velI*y) / r
        else
           velRout                          = 0.0
           if(present(velPhiOut)) velPhiout = 0.0
        end if
        if(present(velZOut  ) .AND. present(velK) ) velZout   = velK
     end if
  end if

  if (r == 0.0) return

  if (ndim > 2) then
     phi = atan2(y,x)
  end if

end subroutine gr_coordTransfmCyl

subroutine gr_coordTransfmFromCyl(r,zcyl,phi, x,y,z, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
  implicit none
  real,intent(OUT) :: x,y,z
  real,intent(IN)  :: r,zcyl,phi
  integer,intent(IN) :: ndim
  real,OPTIONAL,intent(IN) :: velI,velJ,velK
  real,OPTIONAL,intent(INOUT) :: velIOut,velJOut,velKOut

  real :: xphi
  real :: rsq, velSq, velIJ

  xphi = phi
  if (ndim < 3) xphi = 0.0

  z = zcyl              !this is set from "Y", the second coordinate in FLASH cylindrical coordinates
  y = r * sin(xphi)
  x = r * cos(xphi)

  if (present(velKOut)) then
     velKOut = velJ
  end if
  if (present(velJOut)) then
     velJOut = velI * sin(xphi)  +  velK * cos(xphi)
  end if
  if (present(velIOut)) then
     velIout = velI * cos(xphi)  -  velK * sin(xphi)
  end if


end subroutine gr_coordTransfmFromCyl

end subroutine Grid_coordTransfm

!#ifdef DEBUG
!program test
!  implicit none
!
!  real :: x,y,z,r,theta,phi
!  do while(.true.)
!     print*,'Enter x,y,z:'
!     read*,x,y,z
!     call gr_coordTransfmSph(x,y,z,r,theta,phi,3)
!     print*,'R    =',r
!     print*,'theta=',theta,' (',theta*180/PI,' degrees)'
!     print*,'phi  =',phi  ,' (',phi  *180/PI,' degrees)'
!  end do
!
!end program test
!#endif
