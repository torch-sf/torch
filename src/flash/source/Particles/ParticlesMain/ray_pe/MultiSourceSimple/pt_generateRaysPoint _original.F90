!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   loops over all source particles and generates initial rays
!!   
!!
!! Input: 
!!   dt: current simulation timestep

subroutine pt_generateRaysPoint(dt)

! use Particles_data,ONLY : particles, pt_typeInfo
  use Particles_rayData
  use HEALPixModule, ONLY :  mk_pix2xy, mk_xy2pix, mk_xy2pix1, pix2vec_nest
  use mtmod
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  real,intent(in) :: dt
  integer         :: p, nside, i
  integer(kind=8) :: n_photons, j
  real            :: ionizingN, dissN, heatN, ionizingNH
  real            :: xrot, yrot, zrot, wrot
  real            :: u1, u2, u3, magrot

! direction vector of heal pix (center position)
  real, dimension(NDIM) :: dir
  real, dimension(NPART_PROPS) :: sourcedata

! no local rays
  ph_localRays = 0
  nside = 2**(ph_initHPlevel-1)

! get number of pixels, misleadingly named n_photons
! for level 4 this is 768, -1 because we count from 0
  n_photons = 12*nside*nside-1

! generate new rotation matrix
  if(ph_rotRays) then

! random values
    u1 = grnd()
    u2 = grnd()
    u3 = grnd()

! magical uniform sampling of SO(3), see ...
    wrot = sqrt(1d0-u1)*sin(2d0*PI*u2)
    xrot = sqrt(1d0-u1)*cos(2d0*PI*u2)
    yrot = sqrt(u1)*sin(2d0*PI*u3)
    zrot = sqrt(u1)*cos(2d0*PI*u3)

!normalize quaternion, should be one but root gives errars
    magrot = sqrt(xrot*xrot+yrot*yrot+zrot*zrot+wrot*wrot)
    xrot   = xrot/magrot
    yrot   = yrot/magrot
    zrot   = zrot/magrot
    wrot   = wrot/magrot

! quaternion v' = q v q*
! = (w^2-|q|^2) V + 2 (v.q) Q + 2 w (qXv)
! precalculate some stuff 
    aa = (xrot*xrot + wrot*wrot - yrot*yrot - zrot*zrot)
    ab = 2d0*(yrot*xrot - wrot*zrot)
    ac = 2d0*(zrot*xrot + wrot*yrot)

    ba = 2d0*(xrot*yrot + wrot*zrot)
    bb = (yrot*yrot + wrot*wrot - xrot*xrot - zrot*zrot) 
    bc = 2d0*(zrot*yrot - wrot*xrot) 

    ca = 2d0*(xrot*zrot - wrot*yrot)
    cb = 2d0*(zrot*yrot + wrot*xrot) 
    cc = (zrot*zrot + wrot*wrot - yrot*yrot - xrot*xrot)

!alternative calculation
! u1:
!dir(1) *(xrot*xrot + wrot*wrot - yrot*yrot - zrot*zrot) +
!dir(2) *2d0*(yrot*xrot - wrot*zrot) +
!dir(3) *2d0*(zrot*xrot + wrot*yrot)

! u2:
!dir(1) *2d0*(xrot*yrot + wrot*zrot) +
!dir(2) *(yrot*yrot + wrot*wrot - xrot*xrot - zrot*zrot) +
!dir(3) *2d0*(zrot*yrot - wrot*xrot)

! u3:
!dir(1) *2d0*(xrot*zrot - wrot*yrot) +
!dir(2) *2d0*(zrot*yrot + wrot*xrot) +
!dir(3) *(zrot*zrot + wrot*wrot - yrot*yrot - zrot*zrot)

! u1 = 2d0*(dir(1)*xrot + dir(2)*yrot + dir(3)*zrot)*xrot   			+ &
! (wrot*wrot - xrot*xrot - yrot*yrot - zrot*zrot)*dir(1) + &
! 2d0*wrot*(yrot*dir(3) - zrot*dir(2))

! u2 = 2d0*(dir(1)*xrot + dir(2)*yrot + dir(3)*zrot)*yrot   			+ &
! (wrot*wrot - xrot*xrot - yrot*yrot - zrot*zrot)*dir(2) + &
! 2d0*wrot*(zrot*dir(1) - xrot*dir(3))

! u3 = 2d0*(dir(1)*xrot + dir(2)*yrot + dir(3)*zrot)*zrot   			+ &
! (wrot*wrot - xrot*xrot - yrot*yrot - zrot*zrot)*dir(3)	+ &
! 2d0*wrot*(xrot*dir(2) - yrot*dir(1))
  endif

#ifdef ACTIVE
  call pt_sourceUpdateActive()
  
  p = ph_localRays
! loop over sources
  i = 0
! this should take the particle as argument to fully decouple from number of particle properties
  call pt_assignRayActive(i, sourcedata)
  
  do while (i .ge. 0)
    ionizingN   = sourcedata(NION_PART_PROP)

! ionizing photons at surface of star per HEALPix pixel in this timestep
    ionizingN   = dt*ionizingN/(n_photons+1)

! generate rays
    do j = 0, n_photons
      p = p + 1

! direction of photon package rays
      call pix2vec_nest(int(nside,4), j, dir)

      if (p > ph_maxNRays) then
        print *,'PARAMETER ph_maxNRays is set to ', ph_maxNRays,p
        print *,'  To avoid this crash, redimension bigger in your flash.par'
        print *,'  or buy more ram'
        call Driver_abortFlash &
        ("pt_generateRays:  Exceeded max # of rays/processor!")
      endif

!   same block as source initially
      raysIntProp(iblk, p)   = sourcedata(BLK_PART_PROP)
      raysIntProp(isid, p)   = sourcedata(TAG_PART_PROP)
      raysIntProp(iproc, p)  = sourcedata(PROC_PART_PROP)
 
      if(ph_rotRays) then

        u1 = dir(1)*aa + dir(2)*ab + dir(3)*ac
        u2 = dir(1)*ba + dir(2)*bb + dir(3)*bc
        u3 = dir(1)*ca + dir(2)*cb + dir(3)*cc

        dir(1) = u1
        dir(2) = u2
        dir(3) = u3

! renormalize, because there is a slight error from rotation
        magrot = sqrt(sum(dir*dir))
        dir = dir/magrot
      else
! fix numerical precision issue, only happens with no rotation and 
! rays along cardinal directions, i.e. exactly 1 for x, y or z coordinate
        if(dir(1) .eq. 1.0 .or. dir(1) .eq. -1.0 ) then
          dir(2) = 0.0
          dir(3) = 0.0
        endif

       if(dir(2) .eq. 1.0 .or. dir(2) .eq. -1.0 ) then
          dir(1) = 0.0
          dir(3) = 0.0
        endif

       if(dir(3) .eq. 1.0 .or. dir(3) .eq. -1.0 ) then
          dir(2) = 0.0
          dir(1) = 0.0
        endif
      endif

!   column, distance from source, healpix level and number
      raysRealprop(inion,p)  = ionizingN
      raysRealProp(ieion,p)  = sourcedata(EION_PART_PROP)
      raysRealProp(isigh,p)  = sourcedata(SIGH_PART_PROP)
      raysRealProp(irad,p)   = 0d0
      raysIntProp(ihlev, p)  = nside

! int*4 might not be enough, save as real?
      raysRealProp(ihnum, p) = j
      raysIntProp(istpd, p)  = -1

!   for rays this is the sources position
      raysRealprop(iposx,p)  = sourcedata(POSX_PART_PROP)
!   for rays this is the direction of the ray
      raysRealprop(ivelx,p)  = dir(1)

      raysRealprop(iposy,p)  = sourcedata(POSY_PART_PROP)
      raysRealprop(ively,p)  = dir(2)

      raysRealprop(iposz,p)  = sourcedata(POSZ_PART_PROP)
      raysRealprop(ivelz,p)  = dir(3)
    enddo
    call pt_assignRayActive(i, sourcedata)
  enddo
  ph_localRays = p
#endif

! check for other source types if applicable
  return
end subroutine pt_generateRaysPoint
