!
! The Healpix Module
!
! We extracted just the few lines of code for the pix2ang_ring subroutine
! from the pix_tools module of the original healpix package.
! We only need the subroutine that gives us the angular coordinates
! for a pixel of pixel index ipix.
!
! 15/10/14 Lars Buntemeyer
!
module Healpix

  implicit none

  contains

  !=======================================================================
  !     pix2ang_ring
  !
  !     renders theta and phi coordinates of the nominal pixel center
  !     for the pixel number ipix (RING scheme)
  !     given the map resolution parameter nside
  !
  !     NOTE: Extracted from pix_tools.F90
  !     We got rid of all the healpix dependent datatypes, works still
  !     very accurately and works for our purposes.
  !
  !=======================================================================

  elemental subroutine pix2ang_ring(nside, ipix, theta, phi)

    implicit none

    INTEGER, INTENT(IN)  :: nside
    INTEGER, INTENT(IN)  :: ipix
    REAL,    INTENT(OUT) :: theta, phi

    INTEGER ::  nl2, nl4, iring, iphi
    INTEGER ::  npix, ncap, ip
    REAL ::  fodd, dnside
    INTEGER, PARAMETER :: ns_max4=8192     ! 2^13
    INTEGER, PARAMETER :: ns_max=ns_max4 ! largest nside available
    REAL, parameter :: half = 0.5d0
    REAL, parameter :: PI = 3.14159265359d0
    REAL, parameter :: HALFPI = 1.57079632679d0
    !real(kind=dp), parameter :: one  = 1.000000000000000_dp
    !real(kind=dp), parameter :: three = 3.00000000000000_dp
    REAL, parameter :: threehalf = 1.5d0
    character(len=*), parameter :: code = "pix2ang_ring"
    !-----------------------------------------------------------------------
!    if (nside > ns_max4) then
!!       print*,code,"> nside out of range"
!       stop
!    endif
    npix = nside2npix(nside)       ! total number of points
!    if (ipix <0 .or. ipix>npix-1) then
!!       print*,code,"> ipix out of range"
!       stop
!    endif

    nl2  = 2*nside
    ncap = nl2*(nside-1) ! points in each polar cap, =0 for nside =1
    dnside = real(nside)

    if (ipix < ncap) then ! North Polar cap -------------

!       iring = nint( sqrt( (ipix+1) * half ), kind=MKD) ! counted from North pole
       iring = (cheap_isqrt(2*ipix+2) + 1)/2
       iphi  = ipix - 2*iring*(iring - 1)

!       theta = ACOS( one - (iring/dnside)**2 / three )
       theta = 2.d0 * asin(iring / (sqrt(6.d0)*dnside))
       phi   = (real(iphi) + half) * HALFPI/iring

    elseif (ipix < npix-ncap) then ! Equatorial region ------

       ip    = ipix - ncap
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = iand(ip, nl4-1)

       fodd  = half * ( iand(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (threehalf*dnside) )
       phi   = (real(iphi) + fodd) * HALFPI / dnside

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix
!       iring = nint( sqrt( ip * half ), kind=MKD)     ! counted from South pole
       iring = (cheap_isqrt(2*ip) + 1) / 2
       iphi  = 2*iring*(iring + 1) - ip

!       theta = ACOS( (iring/dnside)**2 / three  - one)
       theta = PI - 2.d0 * asin(iring / (sqrt(6.d0)*dnside))
       phi   = (real(iphi) + half) * HALFPI/iring

    endif

    return

    contains


    pure function nside2npix(nside) result(npix_result)
      !=======================================================================
      ! given nside, returns npix such that npix = 12*nside^2
      !  nside should be a power of 2 smaller than ns_max
      !  if not, -1 is returned
      ! EH, Feb-2000
      ! 2009-03-04: returns i8b result, faster
      !=======================================================================
      INTEGER             :: npix_result
      INTEGER, INTENT(IN) :: nside
      INTEGER :: npix, ns
      CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
      !=======================================================================
      if (nside < 1 .or. nside > ns_max .or. iand(nside-1,nside) /= 0) then
         !print*,code,": Nside=",nside," is not a power of 2."
         ! Since we want to make everything pure in these routines, we can not
         ! print. If an invalid nside was chosen, default to the next possible
         ! smaller power of 2
         ns = 2**(FLOOR(LOG(REAL(nside))/LOG(2.)))
      else
         ns = nside
      endif
      npix = (12*ns)*ns
      npix_result = npix
      return
    end function nside2npix

    elemental function cheap_isqrt(lin) result (lout)
      integer, intent(in) :: lin
      integer :: lout, diff
      real :: dout, din
      lout = floor(sqrt(dble(lin))) ! round-off error may offset result
      diff = lin - lout*lout ! test Eq (1)
      if (diff <0)      lout = lout - 1
      if (diff >2*lout) lout = lout + 1
      return
    end function cheap_isqrt

  end subroutine pix2ang_ring

end module Healpix
