!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014

!! Description:
!!   raytracing inside a block
!!   calculates intersection of zones and the ray 
!!   calls splitting routine to adapt number of rays to mesh 
!!   calls rate calculation for later radiation transport
!!   hardcoded for 3d
!!
!! Input: 
!!   p            : ray index
!!   blocko       : block index 
!!   blockBounds  : bounding box of the block
!!   dr           : ray segment length
!!   xpoint       : intersection point 
!!   face         : the face of the block the rays leaves
!!   dt           : current simulation timestep
!!   facedir      : direction the ray leave the block

!!  TODO remove final debug code

!! Note: After reading the paper appendix and calculating some of this,
!!       I'm adding my own comments to alleviate confusion. - JW
#define DEBUG
recursive subroutine pt_blockRayTracePoint(p,blocko,blockBounds,dr,xpoint,face,dt,facedir)

  use Particles_rayData
  use tree, ONLY : bnd_box, lrefine
  
  use Driver_interface, ONLY : Driver_abortFlash

! to calculate zone sizes 
  use Grid_data, ONLY : gr_delta, gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax
! access to state variables, for testing only, should be accessed through actual RT routine
  use Grid_interface, ONLY : Grid_getPointData, Grid_putPointData
! for debug
#ifdef DEBUG
  use pt_rayAsyncComm, only : ph_meshMe
#endif
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: p
  real, intent(inout) :: blockBounds(LOW:HIGH,1:MDIM) !could be updated during splitting
  integer, intent(inout) :: blocko !blockid could be updated during splitting, photons could be created
  real, dimension(NDIM), intent(inout) :: xpoint ! intersection point used for finding neighbors
  integer, dimension(NDIM), intent(out) :: facedir ! intersection point used for finding neighbors
  real, intent(in) :: dt
! not used facdir instead
  integer, intent(out) :: face

! not sure this is used anywhere
  real, intent(out) :: dr
  real, dimension(NDIM) :: n, s, ns !normal vector, source pos., distance to source from zone center, normal zc to s
  real, dimension(NDIM) :: zonesize

! variables needed for box-ray intersection
  real, dimension(NDIM) :: div, low, up, lowB, upB, tmin, tmax
  real, dimension(NDIM) :: ttmax, st, o, ray, nsq, oldray
  real :: radius, tnear, step, newRadius

! zone related stuff, sind = step index, tsind = total step index
  integer, dimension(NDIM) :: sind, tsind

  real :: AngleFaceQuotient, Eion, Nion, Vpix
  real :: sH

  real :: solidA!, Apix
  logical :: IsChild, inside, stopp
  integer :: facenear, Nnew, boundCheck

! for radiation pressure
  real :: dirx, diry, dirz

! bin of the ray - JW
  integer :: ph_type

  radius = raysRealProp(irad,p)

! xpoint is only needed for splitting as ray restarts in middle of the block
! not used here until overwritten by ray, legacy code for good reason
  s = xpoint

! direction saved in velocity property to reduce number of added properties
! Note these are the direction cosines. - JW
  n(1)  = raysRealProp(ivelx,p)
  n(2)  = raysRealProp(ively,p)
  n(3)  = raysRealProp(ivelz,p)

! current ray position
! source position (its IS ACTUALLY source position - JW)
  st(1) = raysRealProp(iposx,p)
  st(2) = raysRealProp(iposy,p)
  st(3) = raysRealProp(iposz,p)

! stupid root, numerical errors add up in radius so use xpoint
! 99% sure xpoint is the current location of the ray. - JW
! I don't think this is necessarily the entrance point of the ray
! into the block. -JW
  oldray = (xpoint-st)

  dr  = 0d0 ! reset just to be sure
  nsq = n*n

! intersection points

! Not actually intersection points, instead these are the
! lower and upper block boundaries. We'll use these to find the
! closest *exit* boundary, which is used to calculate the exit point. - JW
  lowB(1)  = blockBounds(LOW ,IAXIS)
  lowB(2)  = blockBounds(LOW ,JAXIS)
  lowB(3)  = blockBounds(LOW ,KAXIS)

  upB(1)   = blockBounds(HIGH,IAXIS)
  upB(2)   = blockBounds(HIGH,JAXIS)
  upB(3)   = blockBounds(HIGH,KAXIS)

! assumed cubic zones
  zonesize = gr_delta(1:MDIM,lrefine(blocko))
  Vpix = zonesize(1)*zonesize(2)*zonesize(3) ! Not actually a pixel, its a grid zone -JW

! Added a line to error out if the grid zones aren't cubes. - JW
#ifdef DEBUG
  if (zonesize(1) .ne. zonesize(2) .or. &
      zonesize(3) .ne. zonesize(2) .or. &
      zonesize(1) .ne. zonesize(3) ) &
      call Driver_abortFlash("FERVENT assumes cubic zones!")
#endif

! could also be saved directly as particle property
! 1 / n is the inverse of the direction cosines. - JW
  div   = 1. /  n

! calculate the number of zones the ray traverses
! This actually is the distance from the source to each boundary (upper
! and lower) in each direction, divided by the direction cosine. - JW
! Note also this is always a positive number, since its divided by n. - JW
  tmin  = (lowB - st) * div
  tmax  = (upB  - st) * div

! This finds the sides most distance from the source position in each
! direction (i.e. the sides the ray exits.) - JW
! ttmin  = min(tmin, tmax)
  ttmax = max(tmin, tmax)

! min and maxval are F95 intrinsic functions!, face gives the face the ray cuts
! And here we find the closest side to the source that the ray crosses on
! exit. Note this is the actual side the ray uses to exit the block, all
! other lines that correspond to the block boundaries are crosses outside
! the actual block. - JW
! We can then use this one exit side to create the entire exit point from
! the direction cosines. - JW

  tmin  = minloc(ttmax, 1) ! Index in the array of the closest exit side (x,y or z) - JW

! type conversion here
! And the actual value of this boundary. - JW
  facenear = tmin(1)
  tnear = ttmax(facenear) ! Note again, this is the nearest *exit* face

  face  = facenear

! check how many and which zones are cut
! zone size, assumes cubic zones
! This cubic zone assumption means the code quietly fails if zones aren't cubic! - JW
  step = zonesize(IAXIS)

! current intersection of block and ray
! Now using the direction cosines construct the entire exit point for the
! ray (the intersection of the block and ray *exit*.) - JW
  ray = tnear*n + st ! This is the distance from the source to the exit point. - JW
  newRadius = tnear

! corner of the block is origin of block coordinate system
  o(1)  = blockBounds(LOW, IAXIS)
  o(2)  = blockBounds(LOW, JAXIS)
  o(3)  = blockBounds(LOW, KAXIS)

! gives steps to traverse by dividing end position of ray with stepsize, gives last zone position in block
! this has to decide to numerical precision where the ray enters and leaves the block
! the initial and final zone IDs are hardcoded to avoid 'meandering rays'

! Right, this one is much clearer. The distance in each direction through
! the block is the (ray - o) and the number of cells is this divided by
! the number of cells in a block. For some reason I can't grasp why this is
! indexed like in C, starting at zero and going to NXB-1... - JW

! Note also the closest exit face is assumed to be completely traversed, so
! if facenear = 1 (x boundary) sind(1) = NXB-1, etc. But this is not
! always true if the ray only cuts the corner of a block, so I'm not sure
! about it... it could be okay as long as checks are made to see if the
! ray has left the block in this direction - JW
  select case (facenear)
    case (1)
      if(n(1) .ge. 0e0) then ! if the ray has an x component - JW
        sind(1)  = NXB-1
        sind(2)  = floor( (ray(2) - o(2))/step)
        sind(3)  = floor( (ray(3) - o(3))/step)
      else ! Put in zero by hand because of truncation error I'm assuming - JW
        sind(1)  = 0 !floor( (ray(1) - o(1) +1d0 )/step )
        sind(2)  = floor( (ray(2) - o(2) )/step)
        sind(3)  = floor( (ray(3) - o(3) )/step)
      endif
    case (2)
      if(n(2) .ge. 0e0) then 
        sind(1)  = floor( (ray(1) - o(1) )/step)
        sind(2)  = NYB-1
        sind(3)  = floor( (ray(3) - o(3) )/step)
      else
        sind(1)  = floor( (ray(1) - o(1) )/step)
        sind(2)  = 0
        sind(3)  = floor( (ray(3) - o(3) )/step)
      endif
    case (3)
      if(n(3) .ge. 0e0) then
        sind(1)  = floor( (ray(1) - o(1) )/step )
        sind(2)  = floor( (ray(2) - o(2) )/step )
        sind(3)  = NZB-1 
      else
        sind(1)  = floor( (ray(1) - o(1) )/step )
        sind(2)  = floor( (ray(2) - o(2) )/step )
        sind(3)  = 0
      endif
  end select

! fix overshoot, sometimes it is too accurate 16.0000.. = 16 which is one zone too much
! this is always one less than the actual index of the zone, as it is used to calculate zoneboundaries
! |1|2|3|4|5|6|7|8| ... zone indexing 
! 0 1 2 3 4 5 6 7   ... zone boundary indexing
  if(sind(1) .ge. NXB ) sind(1) = NXB-1
  if(sind(2) .ge. NYB ) sind(2) = NYB-1
  if(sind(3) .ge. NZB ) sind(3) = NZB-1

! So then here I think we are finding which
! cell the ray is currently in. Note
! int truncates towards zero always. - JW

! NOTE!!! I'm currently getting an error in pt_solveZone that
! this number is less one even though its passed as tsind+1 in the
! call to pt_solveZone! - JW

! This being the current position in the block, it should never be
! negative since the block origin is the lower left corner always. - JW
  tsind = int((xpoint - o)/step)

! As above, if we are = NXB then int(NXB)=NXB, so subtract 1. - JW
  if(tsind(1) .ge. NXB ) tsind(1) = NXB-1
  if(tsind(2) .ge. NYB ) tsind(2) = NYB-1
  if(tsind(3) .ge. NZB ) tsind(3) = NZB-1

! just face cases, no corner or edge
! 2 is center, 3 is upper, 1 is lower face
! Now we are doing things with the block faces again... - JW
  select case (facenear)
    case (1)
      facedir = (/3,2,2/)
      face = 2
      if(n(facenear) .lt. 0) then 
        facedir = (/1,2,2/)
        face = 1 
!      finalface = 1
      endif
    case (2)
      facedir = (/2,3,2/)
      face = 4
!        finalface = 4
      if(n(facenear) .lt. 0) then 
        facedir = (/2,1,2/)
        face = 3
!          finalface = 3
      endif
    case (3)
      facedir = (/2,2,3/)
      face = 6
!        finalface = 6
      if(n(facenear) .lt. 0) then 
        facedir = (/2,2,1/)
        face = 5 
!          finalface = 5
      endif
  end select

! point where ray leaves block
!!! NOTE: We just changed the meaning of xpoint (WHY?! just use ray!) - JW
  xpoint = ray

! find order for the traversal of the zones, gives a hook for gpu ray tracer, 
! as all ray segments in the zones can be calculated at the same time

! These are the cells the ray is leaving and entering (depends on ray
! orientation whether its entering low or up). - JW

  o(1)  = blockBounds(LOW ,IAXIS)
  o(2)  = blockBounds(LOW ,JAXIS)
  o(3)  = blockBounds(LOW ,KAXIS)

  low(1) = o(1)+zonesize(1)*tsind(1)
  low(2) = o(2)+zonesize(2)*tsind(2)
  low(3) = o(3)+zonesize(3)*tsind(3)

  up(1)  = o(1)+zonesize(1)*(tsind(1)+1)
  up(2)  = o(2)+zonesize(2)*(tsind(2)+1)
  up(3)  = o(3)+zonesize(3)*(tsind(3)+1)

! to reduce checks later, find octant the rays moves into
! On the cell size, not the block. - JW
  if (n(1) .ge. 0.0 .and. n(2) .ge. 0.0 .and. n(3) .ge. 0.0) boundCheck = 1 !x_u, y_u, z_u
  if (n(1) .ge. 0.0 .and. n(2) .lt. 0.0 .and. n(3) .ge. 0.0) boundCheck = 2 !x_u, y_l, z_u
  if (n(1) .ge. 0.0 .and. n(2) .ge. 0.0 .and. n(3) .lt. 0.0) boundCheck = 3 !x_u, y_u, z_l
  if (n(1) .lt. 0.0 .and. n(2) .ge. 0.0 .and. n(3) .ge. 0.0) boundCheck = 4 !x_l, y_u, z_u
  if (n(1) .ge. 0.0 .and. n(2) .lt. 0.0 .and. n(3) .lt. 0.0) boundCheck = 5 !x_u, y_l, z_l
  if (n(1) .lt. 0.0 .and. n(2) .ge. 0.0 .and. n(3) .lt. 0.0) boundCheck = 6 !x_l, y_u, z_l 
  if (n(1) .lt. 0.0 .and. n(2) .lt. 0.0 .and. n(3) .ge. 0.0) boundCheck = 7 !x_l, y_l, z_u 
  if (n(1) .lt. 0.0 .and. n(2) .lt. 0.0 .and. n(3) .lt. 0.0) boundCheck = 8 !x_l, y_l, z_l 


! DEBUG Lets check before all the case calls. - JW
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*, "Before case select tsind wrong!"
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! here we decide if in the block loop only upper or lower bounds are checked
! lots of duplicate code here, 8 cases for all occtants

! NOTE: Now we are recasting tmax to be for a single cell, not the whole block. - JW
  select case(boundCheck)
    case (1)! x, y, z up
      tmax(1) = (up(1)  - st(1)) * div(1)
      tmax(2) = (up(2)  - st(2)) * div(2)
      tmax(3) = (up(3)  - st(3)) * div(3)

! sloop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)
! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)
! omegaaaa, this is estimated from equator of healpix sphere
! solidA = area of sphere / number of healpix elements - JW
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA
! this is Azone/omega

! compare cell area to solid angle (not sure why the r^2 is
! left until after the check to see if we are allowing ray
! splitting, I'd think none of these extra calcs happen
! unless this check is on... also we could push this
! behind a preprocessor check for faster code). - JW
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)
! check for splitting before continuing
        if (ph_inBlockSplit) then
! get pixel size at current radius
! now we multiply solid angle by r^2... - JW
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))
! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif
! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            endif
          endif ! AngleFaceQuotient
        endif ! inBlockSplit

! choose x, y or z
! Here we decide which face of the cell the ray leaves (this is no longer a
! reference for a block face). - JW
        facenear = minloc(tmax,1)
        tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
        Eion  = raysRealProp(ieion,p)
        Nion  = raysRealProp(inion,p)
        sH    = raysRealProp(isigh, p)
        dr = tnear-radius

! add dr to particle property
        raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
        stopp = .false.

        dirx = raysRealProp(ivelx,p)
        diry = raysRealProp(ively,p)
        dirz = raysRealProp(ivelz,p)
        ph_type = raysIntProp(itype,p)

        call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                          stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
        if(stopp) then
! mark as dead
          blocko = -20
          exit
        endif

! set new values
        raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
        if(tsind(facenear) .eq. sind(facenear)) then
          inside = .false.
          exit
        endif
! If not the last step, hop to the next cell
! in the direction of the ray. - JW
        if( n(facenear) .lt. 0)then
          tsind(facenear) = tsind(facenear) - 1
        else
          tsind(facenear) = tsind(facenear) + 1
        endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
        if( facenear .eq. 1 ) then 
! low(1)  = o(1)+zonesize(1)*tsind(1)
          up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
          tmax(1) = (up(1)  - st(1)) * div(1)
        else if( facenear .eq. 2 ) then 
  ! low(2)  = o(2)+zonesize(2)*tsind(2)
          up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
          tmax(2) = (up(2)  - st(2)) * div(2)
        else
  ! low(3)  = o(3)+zonesize(3)*tsind(3)
          up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
          tmax(3) = (up(3)  - st(3)) * div(3)
        endif
      enddo
    case(2)
!x_u, y_l, z_u
      tmax(1) = (up(1)  - st(1)) * div(1)
      tmax(2) = (low(2) - st(2)) * div(2)
      tmax(3) = (up(3)  - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)

! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
        if (ph_inBlockSplit) then
! get pixel size at current radius
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif

! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)

! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)

              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            endif
          endif ! AngleFaceQuotient
        endif ! inBlockSplit

! choose x, y or z
        facenear = minloc(tmax,1)
        tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
        Eion  = raysRealProp(ieion,p)
        Nion  = raysRealProp(inion,p)
        sH    = raysRealProp(isigh, p)

        dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
        raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
        stopp = .false.

        dirx = raysRealProp(ivelx,p)
        diry = raysRealProp(ively,p)
        dirz = raysRealProp(ivelz,p)
        ph_type = raysIntProp(itype,p)

        call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                          stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
        if(stopp) then
! mark as dead
          blocko = -20
          exit
        endif

! set new values
        raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
        if(tsind(facenear) .eq. sind(facenear)) then
          inside = .false.
          exit
        endif

        if( n(facenear) .lt. 0)then
          tsind(facenear) = tsind(facenear) - 1
        else
          tsind(facenear) = tsind(facenear) + 1
        endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
        if( facenear .eq. 1 ) then 
! low(1)  = o(1)+zonesize(1)*tsind(1)
          up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
          tmax(1) = (up(1)  - st(1)) * div(1)
        else if( facenear .eq. 2 ) then 
          low(2)  = o(2)+zonesize(2)*tsind(2)
   !  up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
          tmax(2) = (low(2)  - st(2)) * div(2)
        else
  ! low(3)  = o(3)+zonesize(3)*tsind(3)
          up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
          tmax(3) = (up(3)  - st(3)) * div(3)
        endif
      enddo
    case(3)
!x_u, y_u, z_l
      tmax(1) = (up(1) - st(1)) * div(1)
      tmax(2) = (up(2) - st(2)) * div(2)
      tmax(3) = (low(3) - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)

! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
        if (ph_inBlockSplit) then
! get pixel size at current radius
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then

!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif

! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            endif
          endif ! AngleFaceQuotient
        endif ! inBlockSplit

! choose x, y or z
        facenear = minloc(tmax,1)
        tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
        Eion  = raysRealProp(ieion,p)
        Nion  = raysRealProp(inion,p)
        sH    = raysRealProp(isigh, p)
        dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
        raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
        stopp = .false.

        dirx = raysRealProp(ivelx,p)
        diry = raysRealProp(ively,p)
        dirz = raysRealProp(ivelz,p)
        ph_type = raysIntProp(itype,p)

        call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                          stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
        if(stopp) then
! mark as dead
          blocko = -20
          exit
        endif

! set new values
        raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
        if(tsind(facenear) .eq. sind(facenear)) then
          inside = .false.
          exit
        endif

        if( n(facenear) .lt. 0)then
          tsind(facenear) = tsind(facenear) - 1
        else
          tsind(facenear) = tsind(facenear) + 1
        endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
        if( facenear .eq. 1 ) then 
! low(1)  = o(1)+zonesize(1)*tsind(1)
          up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
          tmax(1) = (up(1)  - st(1)) * div(1)
        else if( facenear .eq. 2 ) then 
    !  low(2)  = o(2)+zonesize(2)*tsind(2)
          up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
          tmax(2) = (up(2)  - st(2)) * div(2)
        else
          low(3)  = o(3)+zonesize(3)*tsind(3)
  !   up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
          tmax(3) = (low(3)  - st(3)) * div(3)
        endif
      enddo
    case(4)
!x_l, y_u, z_u
      tmax(1) = (low(1)  - st(1)) * div(1)
      tmax(2) = (up(2)  - st(2)) * div(2)
      tmax(3) = (up(3) - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)
! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
        if (ph_inBlockSplit) then

! get pixel size at current radius
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))
! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif

! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)

! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            endif
          endif ! AngleFaceQuotient
        endif ! inBlockSplit

! choose x, y or z
        facenear = minloc(tmax,1)
        tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
        Eion  = raysRealProp(ieion,p)
        Nion  = raysRealProp(inion,p)
        sH    = raysRealProp(isigh, p)

        dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
        raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
        stopp = .false.

        dirx = raysRealProp(ivelx,p)
        diry = raysRealProp(ively,p)
        dirz = raysRealProp(ivelz,p)
        ph_type = raysIntProp(itype,p)

        call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                          stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
        if(stopp) then
! mark as dead
          blocko = -20
          exit
        endif

! set new values
        raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
        if(tsind(facenear) .eq. sind(facenear)) then
          inside = .false.
          exit
        endif

        if( n(facenear) .lt. 0)then
          tsind(facenear) = tsind(facenear) - 1
        else
          tsind(facenear) = tsind(facenear) + 1
        endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
        if( facenear .eq. 1 ) then 
          low(1)  = o(1)+zonesize(1)*tsind(1)
     !up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
          tmax(1) = (low(1)  - st(1)) * div(1)
        else if( facenear .eq. 2 ) then 
    !  low(2)  = o(2)+zonesize(2)*tsind(2)
          up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
          tmax(2) = (up(2)  - st(2)) * div(2)
        else
  ! low(3)  = o(3)+zonesize(3)*tsind(3)
          up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
          tmax(3) = (up(3)  - st(3)) * div(3)
        endif
      enddo
    case(5)
!x_u, y_l, z_l
      tmax(1) = (up(1)  - st(1)) * div(1)
      tmax(2) = (low(2) - st(2)) * div(2)
      tmax(3) = (low(3) - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)

! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
        if (ph_inBlockSplit) then
! get pixel size at current radius
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif

! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            endif
          endif ! AngleFaceQuotient
        endif ! inBlockSplit

! choose x, y or z
        facenear = minloc(tmax,1)
        tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
        Eion  = raysRealProp(ieion,p)
        Nion  = raysRealProp(inion,p)
        sH    = raysRealProp(isigh, p)

        dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
        raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
        stopp = .false.

        dirx = raysRealProp(ivelx,p)
        diry = raysRealProp(ively,p)
        dirz = raysRealProp(ivelz,p)
        ph_type = raysIntProp(itype,p)

        call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                          stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
        if(stopp) then
! mark as dead
          blocko = -20
          exit
        endif

! set new values
        raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
        if(tsind(facenear) .eq. sind(facenear)) then
          inside = .false.
          exit
        endif

        if( n(facenear) .lt. 0)then
          tsind(facenear) = tsind(facenear) - 1
        else
          tsind(facenear) = tsind(facenear) + 1
        endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
        if( facenear .eq. 1 ) then 
! low(1)  = o(1)+zonesize(1)*tsind(1)
          up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
          tmax(1) = (up(1)  - st(1)) * div(1)
        else if( facenear .eq. 2 ) then 
          low(2)  = o(2)+zonesize(2)*tsind(2)
   !  up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
          tmax(2) = (low(2)  - st(2)) * div(2)
        else
          low(3)  = o(3)+zonesize(3)*tsind(3)
!     up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
          tmax(3) = (low(3)  - st(3)) * div(3)
        endif
      enddo
    case(6)
!x_l, y_u, z_l
      tmax(1) = (low(1) - st(1)) * div(1)
      tmax(2) = (up(2)  - st(2)) * div(2)
      tmax(3) = (low(3) - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
      inside = .true.
      do while(inside)

! update the traversed radius
        radius = raysRealProp(irad,p)
        solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
        solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
        solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
        if (ph_inBlockSplit) then

! get pixel size at current radius
          AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
          if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
            call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
            if(.not. IsChild) then
! no split ray found 
              if (raysIntProp(iblk,p) .lt. 0) then
                return
              endif
! update blockID
              blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
              blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
              return
            else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
              call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
            return
          endif
        endif ! AngleFaceQuotient
      endif ! inBlockSplit

! choose x, y or z
      facenear = minloc(tmax,1)
      tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
      Eion  = raysRealProp(ieion,p)
      Nion  = raysRealProp(inion,p)
      sH    = raysRealProp(isigh, p)

      dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
      raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

      !oldray = ray
      !Vpix = zonesize(1)*zonesize(2)*zonesize(3)

! call rate calculation
      stopp = .false.

      dirx = raysRealProp(ivelx,p)
      diry = raysRealProp(ively,p)
      dirz = raysRealProp(ivelz,p)
      ph_type = raysIntProp(itype,p)

      call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                        stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
      if(stopp) then
! mark as dead
        blocko = -20
        exit
      endif

! set new values
      raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
      if(tsind(facenear) .eq. sind(facenear)) then
        inside = .false.
        exit
      endif

      if( n(facenear) .lt. 0)then
        tsind(facenear) = tsind(facenear) - 1
      else
        tsind(facenear) = tsind(facenear) + 1
      endif
! expand zone boundaries in direction of leaving ray, saves greatly on computation
      if( facenear .eq. 1 ) then 
        low(1)  = o(1)+zonesize(1)*tsind(1)
     !up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
        tmax(1) = (low(1)  - st(1)) * div(1)
      else if( facenear .eq. 2 ) then 
    !  low(2)  = o(2)+zonesize(2)*tsind(2)
        up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
        tmax(2) = (up(2)  - st(2)) * div(2)
      else
        low(3)  = o(3)+zonesize(3)*tsind(3)
     !up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
        tmax(3) = (low(3)  - st(3)) * div(3)
      endif
    enddo
  case(7)
!x_l, y_l, z_u 
    tmax(1) = (low(1) - st(1)) * div(1)
    tmax(2) = (low(2) - st(2)) * div(2)
    tmax(3) = (up(3)  - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
    inside = .true.
    do while(inside)
! update the traversed radius
      radius = raysRealProp(irad,p)
      solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
      solidA = 4d0*PI/(12d0*solidA*solidA)
! for solving save solid angle
!      Apix   = solidA

! this is Azone/omega
      solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
      if (ph_inBlockSplit) then
! get pixel size at current radius
        AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)
        if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
          call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
          if(.not. IsChild) then
! no split ray found 
            if (raysIntProp(iblk,p) .lt. 0) then
              return
            endif
! update blockID
            blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
            blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
            call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
            return
          else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
            call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
            return
          endif
        endif ! AngleFaceQuotient
      endif ! inBlockSplit

! choose x, y or z
      facenear = minloc(tmax,1)
      tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
      Eion  = raysRealProp(ieion,p)
      Nion  = raysRealProp(inion,p)
      sH    = raysRealProp(isigh, p)

      dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
      raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
      stopp = .false.

      dirx = raysRealProp(ivelx,p)
      diry = raysRealProp(ively,p)
      dirz = raysRealProp(ivelz,p)
      ph_type = raysIntProp(itype,p)

      call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                        stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
      if(stopp) then
! mark as dead
        blocko = -20
        exit
      endif

! set new values
      raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
      if(tsind(facenear) .eq. sind(facenear)) then
        inside = .false.
        exit
      endif

      if( n(facenear) .lt. 0)then
        tsind(facenear) = tsind(facenear) - 1
      else
        tsind(facenear) = tsind(facenear) + 1
      endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
      if( facenear .eq. 1 ) then 
        low(1)  = o(1)+zonesize(1)*tsind(1)
     !up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
        tmax(1) = (low(1)  - st(1)) * div(1)
      else if( facenear .eq. 2 ) then 
        low(2)  = o(2)+zonesize(2)*tsind(2)
   !  up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
        tmax(2) = (low(2)  - st(2)) * div(2)
      else
    ! low(3)  = o(3)+zonesize(3)*tsind(3)
        up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
        tmax(3) = (up(3)  - st(3)) * div(3)
      endif
    enddo
  case(8)
!x_l, y_l, z_l
    tmax(1) = (low(1) - st(1)) * div(1)
    tmax(2) = (low(2) - st(2)) * div(2)
    tmax(3) = (low(3) - st(3)) * div(3)

! loop over zones inside block until boundary is reached -> GPU?
    inside = .true.
    do while(inside)

! update the traversed radius
      radius = raysRealProp(irad,p)
      solidA = raysIntProp(ihlev,p)

! omegaaaa, this is estimated from equator of healpix sphere
      solidA = 4d0*PI/(12d0*solidA*solidA)

! this is Azone/omega
      solidA = zonesize(IAXIS)*zonesize(IAXIS)/(solidA)

! check for splitting before continuing
      if (ph_inBlockSplit) then
! get pixel size at current radius
        AngleFaceQuotient = (solidA)/max(1.0d0,radius*radius)

        if(AngleFaceQuotient .lt. ph_locsampling) then
!============================================
!==== split ray
!============================================
          call pt_splitRayPoint(p, lowB, upB, radius, IsChild, Nnew, xpoint, low, up, tsind, zonesize(1))

! child flag from splitting, if ray was replaced by split one: 1 -> 2,3,4,5 raytracing is continued with 2 
          if(.not. IsChild) then
! no split ray found 
            if (raysIntProp(iblk,p) .lt. 0) then
              return
            endif
! update blockID
            blocko = raysIntProp(iblk,p)
! new block could be different refinement, update to be sure
            blockBounds = bnd_box(:,:,blocko)
! resume block ray tracing
            call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
            return
          else
! resume block ray tracing in same block, i.e. split ray stays in this block
! xpoint updated in pt_splitRay
            call pt_blockRayTracePoint(p, blocko, blockBounds, dr, xpoint, face, dt, facedir)
            return
          endif
        endif ! AngleFaceQuotient
      endif ! inBlockSplit

! choose x, y or z
      facenear = minloc(tmax,1)
      tnear    = tmax(facenear)

!============================================
!==== prepare rates calculation
!============================================
  
      Eion  = raysRealProp(ieion,p)
      Nion  = raysRealProp(inion,p)
      sH    = raysRealProp(isigh, p)

      dr = tnear-radius!sqrt(sum(ns*ns))

! add dr to particle property
      raysRealProp(irad,p) = raysRealProp(irad,p)  + dr

! call rate calculation
      stopp = .false.

      dirx = raysRealProp(ivelx,p)
      diry = raysRealProp(ively,p)
      dirz = raysRealProp(ivelz,p)
      ph_type = raysIntProp(itype,p)

      call pt_solveZone(dt, blocko, tsind+1, dr, Eion, sH,Nion, Vpix, &
                        stopp, dirx, diry, dirz, ph_type, zonesize(IAXIS))
! check if photons of the ray have been absorbed, or stopped, annihilated a planet etc.
      if(stopp) then
! mark as dead
        blocko = -20
        exit
      endif

! set new values
      raysRealProp(inion,p)  = Nion

! DEBUG 
#ifdef DEBUG
if(any(tsind .lt. -1) .or. any(tsind .gt. NXB)) then
  print*,'proc ID', ph_meshMe
  print*,'start point', s,p
  print*,'initpos',sind
  print*,'exit point', xpoint
  print*,'outside bounds',tsind,facenear
  print*,'current pos',raysRealProp(iposx:iposz,p)+n*radius
  print*,'lower blockbounds', bnd_box(LOW,:,blocko)
  print*,'upper blockbounds', bnd_box(HIGH,:,blocko)
  print*,'direction vector',raysRealProp(ivelx:ivelz,p)
  print*,'source position',raysRealProp(iposx:iposz,p)
  print*,'radius',raysRealProp(irad,p)
 stop
endif
#endif

! check if it was the last step
      if(tsind(facenear) .eq. sind(facenear)) then
        inside = .false.
        exit
      endif

      if( n(facenear) .lt. 0)then
        tsind(facenear) = tsind(facenear) - 1
      else
        tsind(facenear) = tsind(facenear) + 1
      endif

! expand zone boundaries in direction of leaving ray, saves greatly on computation
      if( facenear .eq. 1 ) then 
        low(1)  = o(1)+zonesize(1)*tsind(1)
     !up(1)   = o(1)+zonesize(1)*(tsind(1)+1)
! update intersection
        tmax(1) = (low(1)  - st(1)) * div(1)
      else if( facenear .eq. 2 ) then 
        low(2)  = o(2)+zonesize(2)*tsind(2)
   !  up(2)   = o(2)+zonesize(2)*(tsind(2)+1)
! update intersection
        tmax(2) = (low(2)  - st(2)) * div(2)
      else
        low(3)  = o(3)+zonesize(3)*tsind(3)
     !up(3)   = o(3)+zonesize(3)*(tsind(3)+1)
! update intersection
        tmax(3) = (low(3)  - st(3)) * div(3)
      endif
    enddo 
  end select

! for numerical imprecision the radius could be recalculated here, but is extra root
 !s = xpoint - st
  raysRealProp(irad,p) = newRadius !sqrt(sum(s*s))

! check for periodicity and adjust accordingly 
! this is after the block traversal, assumes even number of blocks
! periodic boundary condition check
! if ray traversed 1.5 times in any periodic direction stop it
! 1.499 leaves some slop to not traverse neighbour block if distance is exact
  if(ph_periodic) then !global switch
! maybe nested cases are better? 
    if(xperiodic) then
      if( abs(raysRealProp(ivelx,p)*raysRealProp(irad,p)) .gt. ph_periodicBoxL*glDX ) then
        blocko = -20
        return
      endif
    endif

! maybe nested cases are better? 
    if(yperiodic) then
      if( abs(raysRealProp(ively,p)*raysRealProp(irad,p)) .gt. ph_periodicBoxL*glDY ) then
        blocko = -20
        return
      endif
    endif

! maybe nested cases are better? 
    if(zperiodic) then
      if( abs(raysRealProp(ivelz,p)*raysRealProp(irad,p)) .gt. ph_periodicBoxL*glDZ ) then
        blocko = -20
        return
      endif
    endif
  endif
  return
end subroutine pt_blockraytracePoint
