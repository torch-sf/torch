
module RadTrans_HybridCharModule

#include "Flash.h"
#include "constants.h"

implicit none
#include "interpolate_bezier.h"

INTERFACE equal
MODULE PROCEDURE equal_rr,equal_ri,equal_ir
END INTERFACE

INTERFACE obeybox
MODULE PROCEDURE obeybox_real, obeybox_int
END INTERFACE

  real, parameter :: infinity  = 1.0e+100
  real, parameter :: precision = 1.e-9 !1.e+5*EPSILON(infinity)
  integer, dimension(NDIM), parameter :: KD = 1 !(/1,K2D,K3D/)
  integer, dimension(NDIM), parameter :: NB &
#if K3D
    = (/NXB,NYB,NZB/)
#elif K2D
    = (/NXB,NYB/)
#else
    = (/NXB/)
#endif
!
! io for radiation logfiles
!
  integer, parameter :: io      = 100
  integer, parameter :: io_conv = 101
  integer, parameter :: io_sour = 102
!
  type threeDReal
     real :: x, y, z
  end type threeDReal
  type threeDInt
     integer :: i, j, k
  end type threeDInt
!
! The current maximal residual of the temperature iteration 
!
  real, save                :: T_resid_max
  real, save                :: T_resid_max_global
!
! The accuracy for the dust temperature iteration   
!
  real, parameter           :: T_change_max = 1.d-4
!
! Lambda iteration 
!
  real, parameter           :: lambdaIterAccuracy = 1.d-4
  real, save                :: dSourceMax,dMeanMax
  integer(8), save          :: nrOf2ndOrders,nrOf3rdOrders
!
contains
!
!===============================================================================
!
elemental function sb_law(temp)
!
! Stefan-Boltzmann Law for a grey LTE source function
!
  implicit none
  !
  real, intent(in) :: temp
  real :: sb_law
  !
  ! grey LTE source function
  ! S = B = F/pi = sigmaSB/pi * T^4
  !  
  sb_law = 1.8049444d-05 * temp**4  
  !
  return  
!
end function sb_law
!
!===============================================================================
!
elemental function bplanck(temp,nu)
   implicit none
   real, intent(in) :: temp, nu
   real :: bplanck
   if(temp.eq.0.d0) then
      bplanck = 0.d0
      return
   endif
   bplanck = 1.47455d-47 * nu * nu * nu /                    &
             (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
   return
end function bplanck
!
!===============================================================================
!
elemental logical function isinf(x)
   implicit none
   real, intent(in) :: x

   isinf = x.gt.infinity
end function isinf
!
!===============================================================================
!
logical elemental function equal_rr(a,b)
   implicit none
   real, intent(in) :: a,b

   equal_rr = (abs(a).lt.precision.and.abs(b).lt.precision) &
              .OR. (abs(a-b) .le. (precision * max(abs(a),abs(b))))

end function equal_rr

logical elemental function equal_ri(a,b)
   implicit none
   real, intent(in) :: a
   integer, intent(in) :: b
   equal_ri = equal_rr(a,REAL(b))
end function equal_ri

logical elemental function equal_ir(a,b)
   implicit none
   integer, intent(in) :: a
   real, intent(in) :: b
   equal_ir = equal_rr(REAL(a),b)
end function equal_ir
!
!===============================================================================
!-------------------------------------------------------------------------
!             SUBROUTINE FOR OLSON & KUNASZ (1987) QUADRATURE
!-------------------------------------------------------------------------
#ifndef CHECKS_RT
elemental &
#endif
subroutine qdr_olsonkunasz(src_prev,src_curr,src_next, &
                          dtau_prev,dtau_next,order, &
                          u,v,w,res)
  implicit none
  real, intent(in)  :: src_prev,src_curr,src_next
  real, intent(in)  :: dtau_prev,dtau_next
  real, intent(out) :: u,v,w,res
  integer, intent(in) :: order
  real :: xp,e0,e1,e2
  real :: q_max,q

#ifdef CHECKS_RT
  if(order.eq.2.and.dtau_prev.LE.0..OR.&
     order.eq.3.and.(dtau_prev.LE.0..OR.dtau_next.LE.0.)) then
    print *, "order:     ", order
    print *, "dtau_prev: ", dtau_prev
    print *, "dtau_next: ", dtau_next
    call Driver_abortFlash("qdr_olsonkunasz: dtau_prev or dtau_next zero or below.")
  end if
#endif

  ! Precalculate e^{-\Delta\tau}
  !
  xp = exp(-dtau_prev)
  !
  !
  !
  q_max = 0.5d0*(src_prev+src_curr)*dtau_prev
  !
  ! Calculate
  !
  select case(order)
  case(1)
     !
     ! First order integration
     !
     u  = 0.5d0*(1.d0-xp)
     v  = 0.5d0*(1.d0-xp)
     w  = 0.d0
  case(2)
     !
     ! Second order integration
     !
     u  = (1.d0-(1.d0+dtau_prev)*xp)/dtau_prev
     v  = (dtau_prev-1.d0+xp)/dtau_prev
     w  = 0.d0
  case(3)
     !
     ! Third order integration
     !
     e0 = 1.d0-xp
     e1 = dtau_prev-e0
     e2 = dtau_prev*dtau_prev-2.d0*e1
     u  = e0+(e2-(2.d0*dtau_prev+dtau_next)*e1)/(dtau_prev*(dtau_prev+dtau_next))
     v  = ((dtau_prev+dtau_next)*e1-e2)/(dtau_prev*dtau_next)
     w  = (e2-dtau_prev*e1)/(dtau_next*(dtau_prev+dtau_next))
#ifdef CHECKS_RT
  case default
     write(*,*) "qdr_olsonkunasz: ERROR: Do not know order = ",order
     stop
#endif
  end select
  !
  ! Now return the new intensity
  !
  q = u*src_prev + v*src_curr + w*src_next
  !
#ifdef CHECKS_RT
  if(isNaN(q)) then
     write(6,*) 'WARNING: NaN error in qdr_olsonkunasz'
     write(6,*) 'xp',xp
     write(6,*) 'dtau_prev',dtau_prev
     write(6,*) 'dtau_next',dtau_next
     write(6,*) 'u,v,w',u,v,w
     write(6,*) 'order',order
     write(6,*) 'u', 0.5d0*(1.d0-xp)
     write(6,*) 'v', 0.5d0*(1.d0-xp)
  endif
#endif
  ! Prevent overshoots or undershoots
    if(dtau_prev<1.e-6) then
      q = q_max
    end if
  SELECT CASE(order)
  CASE(2)
    !
    ! Prevent undershoots
    !
    if(dtau_prev<1.e-6) then
      res = q_max
      !v = 0.5d0*dtau_prev !(q_max for src_curr=1, src_prev=0)
      ! limits:
      ! for tau->0: v->0
      ! for tau->oo: v->1
      ! make sure these are obeyed:
      !v = MAX(MIN(v, 1.), 0.)
    else
      res = q
    endif
  CASE(3)
    res = max(q,0.d0)
  CASE DEFAULT
    res = q
  END SELECT
end subroutine qdr_olsonkunasz


#ifndef CHECKS_RT
elemental &
#endif
subroutine qdr_bezier(src_prev,src_curr,src_next, &
                      dtau_prev,dtau_next, &
                      a,b, &
                      u,v,w,res)
  implicit none
  real, intent(in)  :: src_prev,src_curr,src_next
  real, intent(in)  :: dtau_prev,dtau_next
  real, intent(in)  :: a, b
  real, intent(out) :: u,v,w,res
  real :: g0,Sc
  real :: q

  ! Prevent overshoots or undershoots
  if(dtau_prev<1.e-6) then
    q = 0.5d0*(src_prev+src_curr)*dtau_prev
  else
    ! Third order integration on bezier interpolation
    ! with support for integration limits: 0 <= a <= b <= 1
    ! See Hayek et al. (2010) for a version with a=0, b=1.
    ! Maybe a monotone cubic interpolation would even be better:
    ! https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
    ! see Fritsch and Carlson (1979): http://epubs.siam.org/doi/10.1137/0717021
    g0 = get_g0(a,b,dtau_prev)
    Sc = get_Sc(src_curr,src_next,src_prev,dtau_prev,dtau_next)
    IF((src_curr-src_prev)*(src_next-src_curr).lt.0.) THEN
      u = Psiu_1(a, b, g0, dtau_prev)
      v = Psi0_1(a, b, g0, dtau_prev)
      w = 0.
    ELSE IF(Sc.lt.min(src_prev,src_curr).or.&
            Sc.gt.max(src_prev,src_curr)) THEN
      u = Psiu_2(a, b, g0, dtau_prev)
      v = Psi0_2(a, b, g0, dtau_prev)
      w = 0.
    ELSE
      u = Psiu_0(a, b, g0, dtau_prev, dtau_next)
      v = Psi0_0(a, b, g0, dtau_prev, dtau_next)
      w = Psid_0(a, b, g0, dtau_prev, dtau_next)
    END IF
    q = u*src_prev + v*src_curr + w*src_next
  end if
  !
#ifdef CHECKS_RT
  if(isNaN(q)) then
     write(6,*) 'WARNING: NaN error in qdr_bezier'
     write(6,*) 'dtau_prev',dtau_prev
     write(6,*) 'dtau_next',dtau_next
     write(6,*) 'u,v,w',u,v,w
  endif
#endif
  res = max(q,0.d0)
end subroutine qdr_bezier
!
!===============================================================================
!
elemental logical function valid(a1,a2,b1,b2,a,b)
   implicit none
   real,intent(in) :: a1,a2,b1,b2,a,b

   valid = ((a1.le.a .or. equal(a1,a)) .and. (a.le.a2 .or. equal(a,a2)) .and.&
            (b1.le.b .or. equal(b1,b)) .and. (b.le.b2 .or. equal(b,b2)))
end function valid
!

pure subroutine obeyintervall(a,v)
  implicit none
  real, dimension(2), intent(in) :: a
  real, intent(inout) :: v
  if(v.lt.a(1)) then
!    if(equal(v,a(1))) &
      v = a(1)
  else if(v.gt.a(2)) then
!    if(equal(a(2),v)) &
      v = a(2)
  end if
end subroutine obeyintervall

elemental subroutine obeybox_real(b0,b1,x)
  implicit none
  real, intent(in) :: b0,b1
  real, intent(inout) :: x
  if(x.lt.b0) then
    if(equal(x,b0)) &
      x = b0
  else if(x.gt.b1) then
    if(equal(b1,x)) &
      x = b1
  end if
end subroutine obeybox_real

elemental subroutine obeybox_int(b0,b1,x)
  implicit none
  integer, intent(in) :: b0,b1
  real, intent(inout) :: x
  call obeybox_real(real(b0),real(b1),x)
end subroutine obeybox_int

#ifndef CHECKS_RT
pure &
#endif
subroutine getAngles(ipix,randvec_theta,randvec_phi,randrot,theta,phi,dirFact)
  !
  use raytrace_data, ONLY: rt_dirX, rt_dirY, rt_dirZ, rt_healpix_nside, &
       rt_nPhi, rt_nTheta
  
  use rt_data_raytrace_3drt, ONLY : nrOfAngles

  use healpix

  implicit none
  !
  
  integer, intent(in) :: ipix
  real, intent(in)    :: randvec_theta,randvec_phi,randrot
  real, intent(inout) :: theta,phi
  real, intent(inout), dimension(NDIM) :: dirFact
  real                :: dir, mu
  integer             :: ipix1, iPhi, iTheta
  logical             :: flipFlop
  real                :: dirX,dirY,dirZ
  flipFlop = .true.
  !
  if(flipFlop) then
     if(mod(ipix,2).eq.0) then
        ipix1 = ipix / 2
     else
        ipix1 = nrOfAngles - ipix/2
     endif
  else
     ipix1 = ipix
  endif
  !
  select case(NDIM)
     !
  case(3)
     !
     if(rt_healpix_nSide.ge.0) then
        !
        ! use healpix tesselation
        !
        call pix2ang_ring(rt_healpix_nSide,ipix1-1,theta,phi)

        !Randomly rotate angle
        call rotate_angle(theta,phi,randvec_theta,randvec_phi,randrot)
        theta = MOD(theta,PI)
        phi = MOD(phi,2*PI)

        ! Get Cartesian direction vector
        !
        dirX = sin(theta) * cos(phi)
        dirY = sin(theta) * sin(phi)
        dirZ = cos(theta)
           !
     elseif(rt_healpix_nSide==-1) then
        !
        ! use direction vector from parameter context
        !
        dirX = rt_dirX ! call get_parm_from_context(global_parm_context, 'rt_dirX', dirX)
        dirY = rt_dirY ! call get_parm_from_context(global_parm_context, 'rt_dirY', dirY)
        dirZ = rt_dirZ ! call get_parm_from_context(global_parm_context, 'rt_dirZ', dirZ)
        !
        dir  = sqrt(dirX**2+dirY**2+dirZ**2)
        !
        dirX = dirX / dir 
        dirY = dirY / dir 
        dirZ = dirZ / dir
        !
        theta = acos(dirZ)
        !
        if(dirX.gt.precision) then
           phi = atan(dirY/dirX)
        elseif(abs(dirX).lt.precision) then
           phi = sign(1.0,dirY) * 0.5 * PI
        elseif(dirX.lt.-1.0*precision) then
           phi = atan(dirY/dirX) + PI
        elseif(dirX.lt.-1.0*precision.and.dirY.lt.-1.0*precision) then
           phi = atan(dirY/dirX) - PI
        endif
        !
     else
        !
        ! use old discretization 
        !
        iPhi    = 1+int((ipix1-1)/rt_nTheta)
        iTheta  = ipix1-(iPhi-1)*rt_nTheta
        phi = 2.d0*PI/dble(2.d0*rt_nPhi) + 2.d0*PI/dble(rt_nPhi) * (iPhi-1)
        !phi = 2.d0*PI/dble(rt_nPhi-1) * (iPhi-1)
        if(rt_nTheta==1) then
           theta = 0.5 * PI ! we are in the xy-plane
           mu    = 0.0
        elseif(rt_nTheta==2) then
           ! along +z and -z direction
           mu    = 2.0 / dble(rt_nTheta-1) * (ipix1-1) - 1.0
           theta = acos(mu)
        else
           mu     = 1.0 / dble(rt_nTheta) + 2.0 / dble(rt_nTheta) * (iTheta-1) - 1.0
           !mu     = 2.0 / dble(rt_nTheta-1) * (iTheta-1) - 1.0
           theta = acos(mu)
        endif
        !
        ! Get Cartesian direction vector
        !
        dirX = sin(theta) * cos(phi)
        dirY = sin(theta) * sin(phi)
        dirZ = cos(theta)
        !
     endif
     !
  case(2)
     !
     ! 2d case 
     !
     dirZ  = 0.0
     theta = 0.5 * PI ! we are in the xy-plane
     !
     if(nrOfAngles>1) then
        phi  = 2.0*PI/nrOfAngles * (ipix1-1)
     elseif(nrOfAngles==1) then
        dirX = rt_dirX ! call get_parm_from_context(global_parm_context, 'rt_dirX', dirX)
        dirY = rt_dirY ! call get_parm_from_context(global_parm_context, 'rt_dirY', dirY)
        dir  = sqrt(dirX**2+dirY**2)
        dirX = dirX / dir 
        dirY = dirY / dir
        phi  = atan2(dirY,dirX)
     endif
  case(1)
     !
     dirZ  = 0.0
     dirY  = 0.0
     theta = 0.5 * PI
     !
     if(nrOfAngles==2) then
        phi = PI * (ipix-1)
        dirX = cos(phi)
     else
        dirX = rt_dirX ! call get_parm_from_context(global_parm_context, 'rt_dirX', dirX)
        dirX = sign(1.0,dirX)
     endif
     !
  end select
  !
  ! check for precision
  !
  if(abs(dirX) < precision)  dirX = 0.0
  if(abs(dirY) < precision)  dirY = 0.0
  if(abs(dirZ) < precision)  dirZ = 0.0
  if(abs(phi) < precision)   phi = 0.0
  if(abs(theta) < precision) theta = 0.0

  dirFact = (/ dirX, dirY, dirZ /)
  ! 
end subroutine getAngles

! Rotate vector theta,phi by angle_rot about axis theta_rotax, phi_rotax
pure subroutine rotate_angle(theta,phi,theta_rotax,phi_rotax,angle_rot)
  implicit none 
  real, intent(inout) :: theta, phi
  real, intent(in) :: theta_rotax,phi_rotax,angle_rot
  real :: ux,uy,uz,sinp,cosp,cosp1,x,y,z,xdash,ydash,zdash

  !Cartesian components of original vector
  x = sin(theta)*cos(phi)
  y = sin(theta)*sin(phi)
  z = cos(theta)

  !Cartesian components of rotation axis
  ux = sin(theta_rotax)*cos(phi_rotax)
  uy = sin(theta_rotax)*sin(phi_rotax)
  uz = cos(theta_rotax)

  ! Some definitions
  sinp = sin(angle_rot)
  cosp = cos(angle_rot)
  cosp1 = 1-cos(angle_rot)

  !Apply rotation matrix
  xdash = x*(cosp+ux**2*cosp1) + y*(ux*uy*cosp1-uz*sinp)+z*(ux*uz*cosp1+uy*sinp)
  ydash = x*(uy*ux*cosp1+uz*sinp) + y*(cosp+uy**2*cosp1)+z*(uy*uz*cosp1-ux*sinp)
  zdash = x*(uz*ux*cosp1-uy*sinp) + y*(uz*uy*cosp1+ux*sinp) + z*(cosp+uz**2*cosp1)

  !This would be in range -pi to pi
  phi = ATAN2(ydash,xdash)
  !Convert to to pi range
  phi = phi+PI
  theta = ACOS(zdash/(sqrt(xdash**2+ydash**2+zdash**2)))


end subroutine rotate_angle

  pure function check_incoords(bndmin, bndmax, idx)
!-------------------------------------------------------------------------------
    implicit none
    logical :: check_incoords
    real, dimension(NDIM), intent(in) :: bndmin, bndmax, idx
!-------------------------------------------------------------------------------
    check_incoords = .NOT.(ANY(bndmin.GT.idx) .OR. ANY(bndmax.LT.idx))
    !check_incoords = ALL(bndmin.LE.idx) .AND. ALL(bndmax.GE.idx)
  end function check_incoords

!
!===============================================================================
!

end module RadTrans_HybridCharModule
