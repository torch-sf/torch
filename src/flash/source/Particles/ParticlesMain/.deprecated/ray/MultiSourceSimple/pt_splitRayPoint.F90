!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014

!! Description:
!!   splits a ray into 4 child rays using the HEALPix library
!!   raytracing is continued using one of the child rays
!!
!! Input: 
!!   pid         : parent ray index
!!   low         : lower bounds of block  
!!   up          : upper bounds of block
!!   radius      : ray distance from source 
!!   IsChild     : flag if returned ray is child of split ray 
!!   Nnew        : number of generated rays
!!   EndPhotonID : index of last created ray
!!   xpoint      : intersection point of ray and zone
!!   lowZ        : lower zone boundary
!!   upZ         : upper zone boundary
!!   indZ        : index of zone inside its block
!!   zonesize    : size of the zone
!!
!! TODO replace p counter with ph_localRays
!#define DEBUG
#define PROC 1

subroutine pt_splitRayPoint(pid, low, up, radius, IsChild, Nnew, xpoint, lowZ, upZ, indZ, zonesize)

  use HEALPixModule
  use Particles_rayData
  use tree, only : bnd_box
  use Driver_interface, ONLY : Driver_abortFlash

! global boundary box
  use Grid_data, ONLY :  gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax
  use pt_rayAsyncComm, only : ph_meshMe, ph_sendRay, ph_progressComm

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in)  :: pid
  real,  intent(in)    :: radius, zonesize
  logical, intent(out) :: IsChild
  integer, intent(out) :: Nnew
  real,dimension(NDIM), intent(inout) :: xpoint
  real,dimension(NDIM), intent(in)    :: lowZ, upZ, low, up
  integer,dimension(NDIM), intent(in) :: indZ

  integer,dimension(NDIM) :: newindZ
  integer    :: i, j
  integer*4  :: blocko, hpl, x, y, z, newBID, blockproc
  integer*8  :: hpid, p
  real,dimension(NDIM) :: s, st
  real,dimension(NDIM) :: norm, diff
  real :: normmag, rad
  real :: blockBounds(LOW:HIGH,1:MDIM) !could be updated during splitting
  real,dimension(3)  :: newlow, newup

! for neighbours
  integer :: numNegh, face, sourceID, info, tmp, lvl, EndPhotonID

! for 3d 4 high res neighbours and 3 properties
  integer, dimension(3,4)   :: negh
  integer, dimension(MDIM,4):: neghCornerID
  integer,dimension(2) :: sind
  real,dimension(2) :: o
! parent ray properties
  real :: step, ionizingN, ionizingE
  real :: sigmaH

  logical :: isGone, leave
  real, dimension(ph_transProp) :: oneRay

  integer :: countray

  s(1) = raysRealProp(iposx,pid)
  s(2) = raysRealProp(iposy,pid)
  s(3) = raysRealProp(iposz,pid)

  Nnew = 0
  IsChild = .false.
  blocko  = raysIntProp(iblk,pid)

! increase HEALPix refinement level, this is offset l0 p1 -> l1 p4, room for first three higher ref pixel
  hpid      = raysRealProp(ihnum,pid)*4
  hpl       = raysIntProp(ihlev,pid)*2
  sourceID  = raysIntProp(isid,pid)

! work on parent ray
  ionizingN = raysRealProp(inion,pid)/4d0
  ionizingE = raysRealProp(ieion,pid)

  sigmaH    = raysRealProp(isigh,  pid)

! get last used particle index
  p = ph_localRays 
  EndPhotonID = ph_localRays

!============================================
!==== create child rays
!============================================
! particle memory is already allocated just make active
  do i = 0, 3
! get new direction, check which block it belongs to
    call pix2vec_nest(hpl, hpid+i, norm)
    if(ph_rotRays) then

! quaternion multiplication, use precalculated values
      st(1) = norm(1)*aa + norm(2)*ab + norm(3)*ac
      st(2) = norm(1)*ba + norm(2)*bb + norm(3)*bc
      st(3) = norm(1)*ca + norm(2)*cb + norm(3)*cc

! assign rotated direction vector
      norm = st

! renormalize, because there is a slight error from rotation
      normmag = sqrt(sum(norm*norm))
      norm = norm/normmag
    endif

! updated absolute position in domain
    diff = (s + radius*norm)

! check against global domain and periodic BCs
    if(ph_periodic) then
      if(.not. xperiodic) then
! check x global domain boundaries
        if(  (diff(1) .lt. gr_imin) .or. (diff(1) .gt. gr_imax) ) then
          cycle
        endif
      endif

      if(.not. yperiodic) then
! check y global domain boundaries
        if(   (diff(2)  .lt. gr_jmin)  .or. (diff(2) .gt. gr_jmax) ) then
          cycle 
        endif
      endif

      if(.not. zperiodic) then
! check z global domain boundaries
        if(   (diff(3)  .lt. gr_kmin)  .or. (diff(3) .gt. gr_kmax) ) then
          cycle 
        endif
      endif
    else
! check all boundaries at once
      if((diff(1) .lt. gr_imin) .or. (diff(1) .gt. gr_imax) .or. (diff(2)  .lt. gr_jmin) &
        & .or. (diff(2) .gt. gr_jmax) .or. (diff(3) .lt. gr_kmin) .or. (diff(3) .gt. gr_kmax) ) then
        cycle
      endif
    endif

! check if split ray leaves zone, let's hope it's not exactly on an edge (very unlikely)
    newIndZ = floor( (diff - low)/zonesize)

!============================================
!==== split ray outside current block
!============================================
! compound logicals are slow ...
    if( (newIndZ(1) .gt. NXB-1) .or. (newIndZ(2) .gt. NYB-1) .or. (newIndZ(3) .gt. NZB-1) .or. & 
      (newIndZ(1) .lt. 0)     .or. (newIndZ(2) .lt. 0)     .or. (newIndZ(3) .lt. 0) ) then

! check which block to take, binary search tree would be better but meh
! only one zone is skipped
      if ( newIndZ(1) .lt. 0 ) then
        x = 1
! change newIndZ to local values 
      elseif(newIndZ(1) .gt. NXB-1) then
        x = 3
      else
        x = 2
      endif

      if ( newIndZ(2) .lt. 0 ) then
        y = 1 !face 3
      elseif( newIndZ(2) .gt. NYB-1) then
        y = 3 !face 4
      else
        y = 2
      endif

      if ( newIndZ(3) .lt. 0 ) then
        z = 1 !face 5
      elseif( newIndZ(3) .gt. NZB-1 ) then
        z = 3 !face 6
      else
        z = 2
      endif

! in total 26 neighbouring blocks
! 6 faces
      if(x .eq. 1 .and. y .eq. 2 .and. z .eq. 2) then
! -x
        face = 1
      else if (x .eq. 3 .and. y .eq. 2 .and. z .eq. 2) then
! +x
        face = 2
      else if (x .eq. 2 .and. y .eq. 1 .and. z .eq. 2) then
! -y
        face = 3
      else if (x .eq. 2 .and. y .eq. 3 .and. z .eq. 2) then
! +y
        face = 4
      else if (x .eq. 2 .and. y .eq. 2 .and. z .eq. 1) then
! -z
        face = 5
      else if (x .eq. 2 .and. y .eq. 2 .and. z .eq. 3) then
! +z
        face = 6
      endif

! diagonals and MORE, share edge with current block
      if(x .eq. 2 .and. y .ne. 2 .and. z .ne. 2) then
! x constant, ring of 4 possible blocks
        face = 7
      else if (x .ne. 2 .and. y .eq. 2 .and. z .ne. 2) then
! y constant, ring of 4 possible blocks
        face = 8
      else if (x .ne. 2 .and. y .ne. 2 .and. z .eq. 2) then
! z constant, ring of 4 possible blocks
        face = 9
      endif

! here cases where a vertex is shared
! if split criterion is fine enough these cases never happen
! should only give back one neighbour in ptfindNegh, so don't treat explicitly
!    if      (x .ne. 2 .and. y .ne. 2 .and. z .ne. 2) then
!     face = 10
!    endif

! fix source position to old value 
      st = s

      call gr_ptFindNegh(blocko, (/x,y,z/), negh, neghCornerID, numNegh)

      if(numNegh .gt. 1) then
! get stepsize from parent, not most elegant way, but no communication needed
        blockBounds = bnd_box(:,:,blocko)

! higher rest boundary extents
        step    = (blockBounds(HIGH,IAXIS) - blockBounds(LOW,IAXIS))/2d0
        newlow  = blockBounds(LOW,IAXIS:KAXIS)
        newup   = blockBounds(HIGH,IAXIS:KAXIS)

        select case (face)
! use diff position to find out to which child block the ray belongs
! case 1-6: 4 possible child blocks
          case (1)
! /2,4,6,8/, /ymin,zmin,ymin,zmax,ymax,zmin,ymax,zmax/
            o(1) = diff(2) - newlow(2)
            o(2) = diff(3) - newlow(3)

            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1

            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif

            newBID    = negh(1,j)
            blockproc = negh(2,j)

          case (2)
            o(1) = diff(2) - newlow(2)
            o(2) = diff(3) - newlow(3)
            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1

            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif

! /1,3,5,7/, /ymin,zmin,ymin,zmax,ymax,zmin,ymax,zmax/
            newBID    = negh(1,j)
            blockproc = negh(2,j)

          case (3)
            o(1) = diff(1) - newlow(1)
            o(2) = diff(3) - newlow(3)

            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1

            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif

            newBID    = negh(1,j)
            blockproc = negh(2,j)

! /3,4,7,8/, /xmin,zmin,xmax,zmin,xmin,zmin,xmax,zmax/
          case (4)
            o(1) = diff(1) - newlow(1)
            o(2) = diff(3) - newlow(3)
            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1
 
            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif
 
            newBID    = negh(1,j)
            blockproc = negh(2,j)

! /1,2,5,6/, /xmin,zmin,xmax,zmin,xmin,zmin,xmax,zmax/
          case (5)
            o(1) = diff(1) - newlow(1)
            o(2) = diff(2) - newlow(2)
            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1
  
            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif

            newBID    = negh(1,j)
            blockproc = negh(2,j)

! /5,6,7,8/, /xmin,ymin,xmax,ymin,xmin,zmax,xmax,zmax/
          case (6)
            o(1) = diff(1) - newlow(1)
            o(2) = diff(2) - newlow(2)
            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1
            if(sind(2) .eq. 2) sind(2) = 1

            if(sind(2) .gt. 0) then
              j = sum(sind) + 2
            else
              j = sum(sind) + 1
            endif
   
            newBID    = negh(1,j)
            blockproc = negh(2,j)

! cases 7-9: pick upper or lower neighbour 
          case (7)
            j = (diff(1) - newlow(1))/step
            if (j .eq. 2) j = 1
            j = j + 1 
            newBID    = negh(1,j)
            blockproc = negh(2,j)

          case (8)
            j = (diff(2) - newlow(2))/step
            if (j .eq. 2) j = 1
            j = j + 1 
            newBID    = negh(1,j)
            blockproc = negh(2,j)

          case (9)
            j = (diff(3) - newlow(3))/step
            if (j .eq. 2) j = 1
            j = j + 1 
            newBID = negh(1,j)
            blockproc = negh(2,j)
! case 10 only one possible child block per corner
! the numnegh should be one and can be treated in next if block
          end select
        endif

! same or lower res
        if(numnegh .eq. 1) then
          newBID    = negh(1,1) !Neigh(1,face,blocko)
          blockproc = negh(2,1) !Neigh(2,face,blocko)

! move back ray
        else if(numnegh .eq. 0) then
          newBid = -20
        endif
 
        if(newBID .le. -20) then
! continue with next photon
          cycle
        endif

! periodic boundary conditions
! check if ray was split over boundary
        if(ph_periodic) then
          if(xperiodic) then
! check faces, and flip source position if necessary
            if( x .eq. 1 ) then
              if(low(1) .eq. gr_imin ) then
                st(1)    = s(1) + glDX
              endif
            endif

            if( x .eq. 3 ) then
              if(up(1) .eq. gr_imax )  then
                st(1) = s(1) - glDX
              endif
            endif
          endif

          if(yperiodic) then
            if( y .eq. 1 ) then
              if(low(2) .eq. gr_jmin ) then
                st(2) = s(2) + glDY
              endif
            endif

          if( y .eq. 3 ) then
            if(up(2) .eq. gr_jmax )  then
              st(2) = s(2) - glDY
            endif
          endif
        endif

        if(zperiodic) then
          if( z .eq. 1 ) then
            if(low(3) .eq. gr_kmin ) then
              st(3) = s(3) + glDZ
            endif
          endif

          if( z .eq. 3 ) then
            if(up(3) .eq. gr_kmax )  then
              st(3) = s(3) - glDZ
            endif
          endif
        endif
      endif

      p = p + 1
      if (p > ph_maxNRays) then
        call Driver_abortFlash &
          ("pt_splitRay:  Exceeded max # of rays/processor!")
      endif

! this is always added, even if ray is marked for transport
      ph_localRays = ph_localRays + 1
      Nnew        = Nnew + 1

! unique Particle tag, not sure if needed, for single core good enough
! never going to look for specific ray anyway, also why not use HEALPix ID?
!   rays(itag,p)  = p

! Photon's current block number
      raysIntProp(iblk,p)  = newBID
      raysIntProp(isid,p)  = sourceID

      raysRealProp(ivelx,p) = norm(1)
      raysRealProp(ively,p) = norm(2)
      raysRealProp(ivelz,p) = norm(3)

      raysRealProp(ihnum,p) = hpid+i
      raysIntProp(ihlev,p) = hpl

! move to current edge of block
! assumes moved ray is close to old blockboundary, no explicit check is performed 
! to make sure a zone is skipped, as it is very unlikely
                       !call moveRayEdgeBlock(blocko, norm, s, diff, rad)

! copy parent data
      raysRealProp(irad,p)   = radius
      raysRealProp(inion,p)  = ionizingN
      raysRealProp(ieion,p)  = ionizingE

      raysRealProp(isigh,p)  = sigmaH

      raysRealProp(iposx,p) = st(1)
      raysRealProp(iposy,p) = st(2)
      raysRealProp(iposz,p) = st(3)
      raysIntProp(iproc,p) = ph_meshMe

! check if local and mark for transport if not 
      if (blockproc .ne. ph_meshMe) then
! mark as frozen until transported, not needed
!    raysIntProp(istpd,p) = 1
! might be right block might be not
        raysIntProp(iblk,p)  = newBID
        raysIntProp(iproc,p) = blockproc

! buffer reals 
        oneRay(itnion)  = raysRealProp(inion,  p)
        oneRay(ithnum)  = raysRealProp(ihnum,  p)
        oneRay(itrad)   = raysRealProp(irad,   p)
 
        oneRay(itvelx)  = raysRealProp(ivelx,  p)
        oneRay(itvely)  = raysRealProp(ively,  p)
        oneRay(itvelz)  = raysRealProp(ivelz,  p)

        oneRay(itposx)  = raysRealProp(iposx,  p)
        oneRay(itposy)  = raysRealProp(iposy,  p)
        oneRay(itposz)  = raysRealProp(iposz,  p)

        oneRay(itsigh)  = raysRealProp(isigh,  p)  

        oneRay(iteion)  = raysRealProp(ieion,  p)  
        oneRay(itstpd)  = raysRealProp(istpd,  p)
 
        oneRay(itinfo)  = raysIntProp(ihlev,p)
        oneRay(itblk)   = raysIntProp(iblk,p)
        oneRay(itid)    = raysIntProp(isid,p)

! clean up freed up slot
        ph_localRays = ph_localRays - 1

! empty particle slot as it is now copied to destBuf
        raysIntProp(:,p) = -1
        raysRealProp(:,p) = -1

! roll back the local crude memory array pointer
        p = p - 1
        Nnew = Nnew - 1

        countray = ph_localRays
        call ph_sendRay(oneRay, blockproc, isGone)

! try again, and again, again and again, possibly again
! this only happens if send buffer is full and has to be cleared for new ray
        do while (.not. isGone)
          call ph_progressComm(.false.,leave,.false.)
          call ph_sendRay(oneRay, blockproc, isGone)
        end do

! received rays adjust counter
        if(countray .ne. ph_localRays) then
#ifdef DEBUG
          print*,'recv in splitRayPoint',countray,ph_localRays
#endif
! change pointer to accurately point at free slot
          p = ph_localRays
        endif

! photon leaves local domain
! also part of debug, as later sanity check would fail if leave is true
        cycle
      else
        raysIntProp(istpd,p) = -1
      endif
    else

!============================================
!==== split ray inside current block
!============================================

! same zone or different zone but same block
! check a bit redundant
      if( all(newIndZ .eq. IndZ) ) then
! overwrite parent
        if (.not. IsChild) then
          raysRealProp(ivelx,pid) = norm(1)
          raysRealProp(ively,pid) = norm(2)
          raysRealProp(ivelz,pid) = norm(3)
          raysRealProp(inion,pid) = ionizingN
          raysRealProp(ieion,pid) = ionizingE

! do new ray current cell intersection
! split was done at cell entry so move ray to egde pointing TO source
! this should treat artefacts from splitting, for ionizing photons 
! and improve the column based result a bit 
! call moveRayEdge(indZ, blocko, norm, s, diff, rad)
          call moveRayEdge2(indZ, blocko, norm, s, diff, rad)

          raysRealProp(irad,pid)    = rad
          raysRealProp(isigh,pid)   = sigmaH

! implicit type conversion here to double
          raysRealProp(ihnum,pid) = hpid+i
          raysIntProp(ihlev,pid) = hpl

! update current ray position for blockRayTrace
          xpoint = diff

          IsChild = .true.
! continue with next photon
          cycle
        else
! parent was already overwritten just append the other in-block photons
          p = p + 1
          if (p > ph_maxNRays) then
            call Driver_abortFlash &
            ("pt_splitRay:  Exceeded max # of rays/processor!")
          endif

          ph_localRays = ph_localRays + 1
          Nnew = Nnew + 1

! unique Particle tag, not sure if needed, for single core good enough
!     rays(itag,p)  = p

! Photon's current block number
          raysIntProp(iblk,p)   = blocko
          raysIntProp(isid,p)   = sourceID

          raysRealProp(ivelx,p) = norm(1)
          raysRealProp(ively,p) = norm(2)
          raysRealProp(ivelz,p) = norm(3)

          raysRealProp(ihnum,p) = hpid+i
          raysIntProp(ihlev,p)  = hpl

! do new ray current cell intersection
! split was done at cell entry so move ray to egde pointing TO source
! this should treat artefacts from splitting, for ionizing photons 
! and improve the column based result a bit 
          call moveRayEdge2(indZ, blocko, norm, s, diff, rad)

! copy parent data
          raysRealProp(irad,p)  = rad
          raysRealProp(inion,p) = ionizingN
          raysRealProp(ieion,p) = ionizingE
          raysRealProp(isigh,p)  = sigmaH

          raysRealProp(iposx,p) = s(1)
          raysRealProp(iposy,p) = s(2)
          raysRealProp(iposz,p) = s(3)

          raysIntProp(iproc,p) = ph_meshMe
          raysIntProp(istpd,p) = -1
! no need to check if it is the same processor, as it is the same block
! continue with next photon

          cycle
        endif
      else
! not the same zone
! still overwrite parent
        if(.not. IsChild) then
          raysRealProp(ivelx,pid)  = norm(1)
          raysRealProp(ively,pid)  = norm(2)
          raysRealProp(ivelz,pid)  = norm(3)
          raysRealProp(inion,pid)  = ionizingN
          raysRealProp(ieion,pid)  = ionizingE

      !call moveRayEdge(newindZ, blocko, norm, s, diff, rad)
          call moveRayEdge2(newindZ, blocko, norm, s, diff, rad)

          raysRealProp(irad,pid)   = rad
          raysRealProp(isigh,pid)  = sigmaH

! implicit type conversion here
          raysRealProp(ihnum,pid)  = hpid+i
          raysIntProp(ihlev,pid)   = hpl

! update current ray position for blockRayTrace
          xpoint = diff

          IsChild = .true.
! continue with next photon

          cycle
        else
! parent was already overwritten just append the other in-block photons
          p = p + 1
          if (p > ph_maxNRays) then
            call Driver_abortFlash &
              ("pt_splitRay:  Exceeded max # of rays/processor!")
          endif

          ph_localRays = ph_localRays + 1
          Nnew = Nnew + 1

! unique Particle tag, not sure if needed, for single core good enough
!     rays(itag,p)  = p

! Photon's current block number
          raysIntProp(iblk,p)  = blocko
          raysIntProp(isid,p)  = sourceID

          raysRealProp(ivelx,p) = norm(1)
          raysRealProp(ively,p) = norm(2)
          raysRealProp(ivelz,p) = norm(3)

          raysRealProp(ihnum,p) = hpid+i
          raysIntProp(ihlev,p)  = hpl

!call moveRayEdge(newindZ, blocko, norm, s, diff, rad)
          call moveRayEdge2(newindZ, blocko, norm, s, diff, rad)
! copy parent data
          raysRealProp(irad,p)   = rad
          raysRealProp(inion,p)  = ionizingN
          raysRealProp(ieion,p)  = ionizingE

          raysRealProp(isigh,p)  = sigmaH

          raysRealProp(iposx,p) = s(1)
          raysRealProp(iposy,p) = s(2)
          raysRealProp(iposz,p) = s(3)

          raysIntProp(iproc,p) = ph_meshMe
          raysIntProp(istpd,p) = -1
! no need to check if it is the same processor, as it is the same block
          cycle
! continue with next photon
        endif ! end isChild
      endif ! end same/diff zone
    endif ! end same/diff block
  enddo 

!============================================
!==== no child ray inside current block, pick next one
!============================================

! don't know if this is needed, if the code breaks, look here first 
  if(.not. IsChild) then

! fill new free slot
    if(ph_localRays .gt. 1 .and. pid .ne. ph_localRays) then 
! move last ray data to current slot
      raysintProp (:,pid) = raysIntProp (:,ph_localRays)
      raysRealProp(:,pid) = raysRealProp(:,ph_localRays)

! empty last particle slot
      raysintProp (:,ph_localRays) = -1
      raysRealProp(:,ph_localRays) = -1
! update exit point
      xpoint = raysRealProp(iposx:iposz,pid)+raysRealProp(ivelx:ivelz,pid)*raysRealProp(irad,pid)
    else 
! empty particle slot as it is now copied to destBuf
      raysIntProp(:,pid) = -1
      raysRealProp(:,pid) = -1
    endif

    ph_localRays = ph_localRays - 1
  endif

  return
contains

! move to upper or lower boundary
  subroutine moveRayEdge(indZ, blockID, dir, pos, xpoint,radius)
  
    use Particles_rayData
    use tree, ONLY : bnd_box, lrefine

! to calculate zone sizes 
    use Grid_data, ONLY : gr_delta

    implicit none

#include "Flash.h"
#include "constants.h"

    integer,dimension(NDIM), intent(in) :: indZ
    integer, intent(in) :: blockID
    real, dimension(NDIM), intent(in) :: dir,pos
    real, dimension(NDIM), intent(inout) :: xpoint
    real, intent(out) :: radius
    real :: blockBounds(LOW:HIGH,1:MDIM) !could be updated during splitting
    real :: xH0, xHp, tnear, tfar
    real, dimension(NDIM) :: lowZ, upZ, div
    real, dimension(NDIM) :: zonesize

    real, dimension(6) :: tall, tallabs
    integer :: tmin

    blockBounds = bnd_box(:, :, blockID)

! assumed cubic zones
    zonesize = gr_delta(1:MDIM,lrefine(blockID))

! find zone boundaries
    lowZ = indZ*zonesize     + blockBounds(LOW,:)  
    upZ  = (indZ+1)*zonesize + blockBounds(LOW,:)

    div = 1. /dir

! check first boundary to encounter
    tnear = -1e99
    tfar  =  1e99

    tall(1) = (lowZ(1) - pos(1)) * div(1)
    tall(2) = (upZ(1)  - pos(1)) * div(1)

    tnear = max(tnear, min(tall(1), tall(2)))
    tfar  = min(tfar,  max(tall(1), tall(2)))

    tall(3) = (lowZ(2) - pos(2)) * div(2)
    tall(4) = (upZ(2)  - pos(2)) * div(2)

    tnear = max(tnear, min(tall(3), tall(4)))
    tfar  = min(tfar,  max(tall(3), tall(4)))

    tall(5) = (lowZ(3) - pos(3)) * div(3)
    tall(6) = (upZ(3)  - pos(3)) * div(3)

    tnear = max(tnear, min(tall(5), tall(6)))
    tfar  = min(tfar,  max(tall(5), tall(6)))

! update xpoint and radius and so on
! current intersection of block and ray
      !xpoint = tnear*dir 
    radius = tfar !sqrt(sum(xpoint*xpoint))
    xpoint = radius*dir + pos 

  end subroutine moveRayEdge

  subroutine moveRayEdgeBlock(blockID, dir, pos, xpoint,radius)

    use Particles_rayData
    use tree, ONLY : bnd_box

    implicit none

#include "Flash.h"
#include "constants.h"

    integer, intent(in) :: blockID
    real, dimension(NDIM), intent(in) :: dir, pos
    real, dimension(NDIM), intent(inout) :: xpoint
    real, intent(out) :: radius
    real :: blockBounds(LOW:HIGH,1:MDIM) !could be updated during splitting
    real :: xH0, xH2, xHp, tnear, tfar
    real, dimension(NDIM) :: lowZ, upZ, div
    real, dimension(NDIM) :: zonesize

    real, dimension(6) :: tall
    integer :: tmin

    blockBounds = bnd_box(:, :, blockID)

! find block boundaries
    lowZ = blockBounds(LOW,:)  
    upZ  = blockBounds(HIGH,:)
    div = 1. /dir

! check first boundary to encounter
    tnear = -1e99
    tfar  =  1e99 

    tall(1) = (lowZ(1) - pos(1)) * div(1)
    tall(2) = (upZ(1)  - pos(1)) * div(1)

    tnear = max(tnear, min(tall(1), tall(2)))
    tfar  = min(tfar,  max(tall(1), tall(2)))

    tall(3) = (lowZ(2) - pos(2)) * div(2)
    tall(4) = (upZ(2)  - pos(2)) * div(2)

    tnear = max(tnear, min(tall(3), tall(4)))
    tfar  = min(tfar,  max(tall(3), tall(4)))

    tall(5) = (lowZ(3) - pos(3)) * div(3)
    tall(6) = (upZ(3)  - pos(3)) * div(3)

    tnear = max(tnear, min(tall(5), tall(6)))
    tfar  = min(tfar,  max(tall(5), tall(6)))

!      if(tfar .lt. tnear) print*, 'miss'

    radius = tfar !sqrt(sum(xpoint*xpoint))
    xpoint = radius*dir + pos
  end subroutine moveRayEdgeBlock
!==========================
! move to the closest edge
  subroutine moveRayEdge2(indZ, blockID, dir, pos,xpoint,radius)

    use Particles_rayData
    use tree, ONLY : bnd_box, lrefine

! to calculate zone sizes 
    use Grid_data, ONLY : gr_delta

    implicit none

#include "Flash.h"
#include "constants.h"

    integer,dimension(NDIM), intent(in) :: indZ
    integer, intent(in) :: blockID
    real, dimension(NDIM), intent(in) :: dir,pos
    real, dimension(NDIM), intent(inout) :: xpoint
    real, intent(out) :: radius
    real :: blockBounds(LOW:HIGH,1:MDIM) !could be updated during splitting
    real :: xH0, xH2, xHp, tnear, tfar
    real, dimension(NDIM) :: lowZ, upZ, div
    real, dimension(NDIM) :: zonesize

    real, dimension(6) :: tall, tallabs
    integer :: tmin

    blockBounds = bnd_box(:, :, blockID)

! assumed cubic zones
    zonesize = gr_delta(1:MDIM,lrefine(blockID))

! find zone boundaries
    lowZ = indZ*zonesize     + blockBounds(LOW,:)  
    upZ  = (indZ+1)*zonesize + blockBounds(LOW,:)

    div = 1. /dir

! calculate the number of zones the ray traverses
    tall(1:3) = (lowZ - xpoint) * div 
    tall(4:6) = (upZ  - xpoint) * div

! find the absolute minimum value of t and use that
! could be moved ahead or could be moved back (probs back)
    tmin  = minloc(abs(tall),1)
    tnear = tall(tmin)

    if(tmin .lt. 4) then
      radius = (lowZ(tmin) - pos(tmin))*div(tmin)
    else
      tmin = tmin -3 
      radius = (upZ(tmin) - pos(tmin))*div(tmin)
    endif
    xpoint = radius*dir + pos


  end subroutine moveRayEdge2
end subroutine pt_splitRayPoint
