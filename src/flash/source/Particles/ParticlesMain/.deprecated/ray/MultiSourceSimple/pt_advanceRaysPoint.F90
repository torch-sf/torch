!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014
!!
!!
!! Description:
!!   main driver of raytracing algorithm
!!   changes the particles data structure
!!   calls MPI communication for ray particles
!!   loops over all particles and calls raytracing routine for blocks
!!
!! Input: 
!!   dt: current simulation timestep
!!

!! TODO figure out where call to merge-sort to be most effective
!#define DEBUG
subroutine pt_advanceRaysPoint(dt)

  use Particles_rayData!, ONLY : ph_sortMerge, istpd, irad, raysRealProp, &
  use Driver_interface, ONLY : Driver_abortFlash
                                 
  use pt_rayAsyncComm, ONLY : ph_meshMe, ph_progressComm,  ph_CommCheckInterval, ph_sendRay

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax

  use tree, ONLY : bnd_box

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#include "GridParticles.h"

 include "Flash_mpi.h" 

  real, intent(in) :: dt

  integer   :: i, j
  logical   :: activeRays, padvance, leave, firstTrace, noMoreRays
  real      :: dr, radius, timestep, step
  integer   :: p_begin

! as new photons can be created p_count changes
  integer   :: p_countFinal, face

! paramesh block properties variables
  integer   :: blocko, newBID, blockproc, p_count
  real      :: blockBounds(LOW:HIGH,1:MDIM)
  real      :: posx, posy, posz
  real, dimension(NDIM) :: diff, newlow, newup, norm

! for neighbours
  integer :: numNegh, info, lvl, tmp

! for 3d 4 high res neighbours and 3 properties
  integer, dimension(3,4)   :: negh
  integer, dimension(MDIM,4):: neghCornerID
  integer, dimension(MDIM)  :: facedir
  integer,dimension(2)      :: sind
  real, dimension(2)        :: o

! for async comm.
  integer :: raysUntilComm

! for transport
  logical   :: isGone
! 1 ray buffer
  real, dimension(ph_transProp) :: oneRay

! indices for photon particles in data structure
  p_begin = 1
  p_count = ph_localRays

  p_countFinal = ph_localRays

! raytracing is ON
  activeRays = .true.
  i = p_begin

! init number of rays to traverse until we check communication
  raysUntilComm = ph_CommCheckInterval

! active while there are any rays on any subdomain, so global flag
  do while (activeRays)
! only raytrace if there are rays locally, otherwise go to MPI and wait for something to do *thumb twiddle*
    if(ph_localRays .eq. 0 ) then
! locally done
      activeRays = .false.
! don't go to photon loop, but go to MPI
      padvance   = .false.
! this ensures that MPI is still executed though
    else
! reset
      newBID = -100
      padvance   = .true.
      firstTrace = .true.
    endif

! local rays are traversed with splitting, so generation of new rays
    do while(padvance)
      if(raysUntilComm .eq. 0) then

! returns if one only core is used
! the boolean flags if only full send buffers are to be communicated
! if doCounter is .true. then we get a MPI crash
! best guess at the moment is mismatch in allreduce calls
        call ph_progressComm(.false., noMoreRays, .false.)

! reset counter to comm
        raysUntilComm = ph_CommCheckInterval
      endif
! traverse inside current cell block, (hopefully there will be no prison riots)
      blocko = raysIntProp(iblk,i)

! sanity check if photon is inside domain, could also be done with physical boundaries of the local domain
! paramesh defines -20 as outside domain flag, should be in constants.h but its not, actually it is now in 4.0 hooray
#ifdef DEBUG
      if( blocko .le. -20) then
        padvance = .false.
        call Driver_abortFlash("pt_advanceRaysPoint: iterated to an empty slot!")
        exit
      endif
#endif


! get boundaries of block for ray/face intersection
      blockBounds = bnd_box(:,:,blocko)

! get ray position
! for repeat calls we know the position already from pt_blockRayTrace
      if(firstTrace) then
! recalc here?
! distance from source
        radius = raysRealProp(irad,i)
        diff(1) = raysRealProp(iposx,i) + raysRealProp(ivelx,i)*radius
        diff(2) = raysRealProp(iposy,i) + raysRealProp(ively,i)*radius
        diff(3) = raysRealProp(iposz,i) + raysRealProp(ivelz,i)*radius

        firstTrace = .false.
      endif
!============================================
!==== raytrace inside the current block
!============================================
! trace and split rays inside block, see pt_blockRayTracePoint for details
      call pt_blockRayTracePoint(i, blocko, blockBounds, dr, diff, face,  dt, facedir)

! check if ray is done 
      if(blocko .le. -20) then
        padvance = .false.

! fill new free slot
        if(ph_localRays .gt. 1 .and. i .ne. ph_localRays) then
! move last ray data to current slot
          raysintProp (:,i) = raysIntProp (:,ph_localRays)
          raysRealProp(:,i) = raysRealProp(:,ph_localRays)

! empty last particle slot
          raysintProp (:,ph_localRays) = -1
          raysRealProp(:,ph_localRays) = -1

          ph_localRays = ph_localRays - 1
        else
        ! last local ray
          raysRealProp(:,i) = -1
          raysIntProp(:,i) = -1
          ph_localRays = ph_localRays - 1
        endif

! skip rest of routine
        exit
      endif

! this checks for a rare case, as ray splitting is done in pt_blockRayTrace,
! if no child rays are created and current ray is last in index, the current slot will be 
! empty, so bail out
      if( raysIntProp(iblk,i) .lt. 0) then
        padvance = .false.
        exit
      endif

!============================================
!==== find Neighbour block, check AMR
!============================================
      call gr_ptFindNegh(blocko,facedir,negh,neghCornerID,numNegh)

! get current bounding box
      blockBounds = bnd_box(:,:,blocko)

! higher refined region
      if(numNegh .gt. 1) then 
        newBid = -1
! size of block
        step    = (blockBounds(HIGH,IAXIS) - blockBounds(LOW,IAXIS))/2d0
        newlow  = blockBounds(LOW,IAXIS:KAXIS)
        newup   = blockBounds(HIGH,IAXIS:KAXIS)

! each face sees different child block order
        select case (face)
! use position of ray to find where it belongs  
          case (1)
! index in tree structure /2,4,6,8/, physical coordinates /ymin,zmin,ymin,zmax,ymax,zmin,ymax,zmax/
! in negh structure:
! z,z,y
! |----|----|
! | 3  |  4 |
! |----|----| --> j selects right face
! | 1  |  2 |
! |----|----|y,x,x

! ordering is always low low, low max, max max, max low, clockwise in coordinates, so case i and i+1 are the same
            o(1) = diff(2) - newlow(2)
            o(2) = diff(3) - newlow(3)
            sind = floor(o/step)
            if(sind(1) .eq. 2) sind(1) = 1 ! y
            if(sind(2) .eq. 2) sind(2) = 1 ! z

! transform to 1..4 index,
            if(sind(2) .gt. 0) then ! check z
               j = sum(sind) + 2
            else
               j = sum(sind) + 1
            endif

            newBID    = negh(1,j)
            blockproc = negh(2,j)

          case (2)
            o(1) = diff(2) - newlow(2) ! y
            o(2) = diff(3) - newlow(3) ! z
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

! /1,2,3,4/, /xmin,ymin,xmax,ymin,xmin,zmax,xmax,zmax/
          end select
        endif 

! same or lower resolution
        if(numnegh .eq. 1) then
          newBID    = negh(1,1) !Neigh(1,face,blocko)
          blockproc = negh(2,1) !Neigh(2,face,blocko)

        else if(numnegh .eq. 0) then
           newBid = -20
        endif

        if(xhydroper) then
! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
          if(blockBounds(LOW,IAXIS) .eq. gr_imin   .and. face .eq. 1 ) then
            newBid = -20
          endif

          if(blockBounds(HIGH,IAXIS) .eq. gr_imax  .and. face .eq. 2 ) then
            newBid = -20
          endif
        endif

        if(yhydroper) then
! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
          if(blockBounds(LOW,JAXIS) .eq. gr_jmin   .and. face .eq. 3 ) then
            newBid = -20
          endif

          if(blockBounds(HIGH,JAXIS) .eq. gr_jmax  .and. face .eq. 4 ) then
            newBid = -20
          endif
        endif

        if(zhydroper) then
! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
          if(blockBounds(LOW,KAXIS) .eq. gr_kmin   .and. face .eq. 5 ) then
            newBid = -20
          endif

          if(blockBounds(HIGH,KAXIS) .eq. gr_kmax  .and. face .eq. 6 ) then
            newBid = -20
          endif
        endif
!      endif

! it's dead, jim, ray left global domain
        if( newBID .le. -20 ) then
          padvance = .false.
! fill new free slot
          if(ph_localRays .gt. 1 .and. i .ne. ph_localRays) then
! move last ray data to current slot
            raysintProp (:,i) = raysIntProp (:,ph_localRays)
            raysRealProp(:,i) = raysRealProp(:,ph_localRays)
! empty last particle slot
            raysintProp (:,ph_localRays) = -1
            raysRealProp(:,ph_localRays) = -1
            ph_localRays = ph_localRays - 1
        else
        ! last local ray
          raysRealProp(:,i) = -1
          raysIntProp(:,i) = -1
          ph_localRays = ph_localRays - 1
        endif
! skip processor check
        exit
      endif

! check if it is on same core
      if(blockproc .ne. ph_meshMe) then
        padvance = .false.

! mark as stopped until transport
        raysIntProp(istpd,i) = 1
! might be right block might not be, should be though if tree up to date
        raysIntProp(iblk,i)  = newBID
        raysIntProp(iproc,i) = blockproc

! buffer reals , yummy inions
        oneRay(itnion)  = raysRealProp(inion, i)
        oneRay(ithnum)  = raysRealProp(ihnum, i)
        oneRay(itrad)   = raysRealProp(irad,  i)

        oneRay(itnion)  = raysRealProp(inion,  i)
        oneRay(itvelx)  = raysRealProp(ivelx,  i)
        oneRay(itvely)  = raysRealProp(ively,  i)
        oneRay(itvelz)  = raysRealProp(ivelz,  i)

        oneRay(itsigh)  = raysRealProp(isigh,  i)  
        oneRay(iteion)  = raysRealProp(ieion,  i)  
        oneRay(itstpd)  = raysRealProp(istpd,  i)

! for periodic boundaries
        if(ph_periodic) then
          select case(ph_BCcase)
            case (1)! x
! no change
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposy) = raysRealProp(iposy,i)
              oneRay(itposz) = raysRealProp(iposz,i)

              if(blockBounds(LOW,IAXIS)  .eq. gr_imin .and. face .eq. 1 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) + glDX
              endif

              if(blockBounds(HIGH,IAXIS) .eq. gr_imax .and. face .eq. 2 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) - glDX
              endif
   

            case (2)! y
! no change
              oneRay(itposy) = raysRealProp(iposy,i)
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposz) = raysRealProp(iposz,i)

              if(blockBounds(LOW,JAXIS)  .eq. gr_jmin .and. face .eq. 3 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) + glDY
              endif

              if(blockBounds(HIGH,JAXIS) .eq. gr_jmax .and. face .eq. 4 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) - glDY
              endif

            case (3)! z
! no change
              oneRay(itposz) = raysRealProp(iposz,i)
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposy) = raysRealProp(iposy,i)

              if(blockBounds(LOW,KAXIS)  .eq. gr_kmin .and. face .eq. 5 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) + glDZ
              endif

              if(blockBounds(HIGH,KAXIS) .eq. gr_kmax .and. face .eq. 6 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) - glDZ
              endif

            case (4)! x,y
! no change
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposy) = raysRealProp(iposy,i)
              oneRay(itposz) = raysRealProp(iposz,i)

              if(blockBounds(LOW,IAXIS)  .eq. gr_imin .and. face .eq. 1 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) + glDX
              endif

              if(blockBounds(HIGH,IAXIS) .eq. gr_imax .and. face .eq. 2 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) - glDX
              endif

              if(blockBounds(LOW,JAXIS)  .eq. gr_jmin .and. face .eq. 3 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) + glDY
              endif

              if(blockBounds(HIGH,JAXIS) .eq. gr_jmax .and. face .eq. 4 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) - glDY
              endif

            case (5)! x,z
! no change
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposz) = raysRealProp(iposz,i)
              oneRay(itposy) = raysRealProp(iposy,i)

              if(blockBounds(LOW,IAXIS)  .eq. gr_imin .and. face .eq. 1 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) + glDX
              endif

              if(blockBounds(HIGH,IAXIS) .eq. gr_imax .and. face .eq. 2 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) - glDX
              endif

              if(blockBounds(LOW,KAXIS)  .eq. gr_kmin .and. face .eq. 5 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) + glDZ
              endif

              if(blockBounds(HIGH,KAXIS) .eq. gr_kmax .and. face .eq. 6 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) - glDZ
              endif

            case (6)! y,z
! no change
              oneRay(itposy) = raysRealProp(iposy,i)
              oneRay(itposz) = raysRealProp(iposz,i)
              oneRay(itposx) = raysRealProp(iposx,i)

              if(blockBounds(LOW,JAXIS)  .eq. gr_jmin .and. face .eq. 3 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) + glDY
              endif

              if(blockBounds(HIGH,JAXIS) .eq. gr_jmax .and. face .eq. 4 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) - glDY
              endif

              if(blockBounds(LOW,KAXIS)  .eq. gr_kmin .and. face .eq. 5 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) + glDZ
              endif

              if(blockBounds(HIGH,KAXIS) .eq. gr_kmax .and. face .eq. 6 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) - glDZ
              endif

            case (7)! x,y,z
! no change
              oneRay(itposx) = raysRealProp(iposx,i)
              oneRay(itposy) = raysRealProp(iposy,i)
              oneRay(itposz) = raysRealProp(iposz,i)

              if(blockBounds(LOW,IAXIS)  .eq. gr_imin .and. face .eq. 1 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) + glDX
              endif

              if(blockBounds(HIGH,IAXIS) .eq. gr_imax .and. face .eq. 2 ) then
                  oneRay(itposx) = raysRealProp(iposx,i) - glDX
              endif

              if(blockBounds(LOW,JAXIS)  .eq. gr_jmin .and. face .eq. 3 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) + glDY
              endif

              if(blockBounds(HIGH,JAXIS) .eq. gr_jmax .and. face .eq. 4 ) then
                  oneRay(itposy) = raysRealProp(iposy,i) - glDY
              endif

              if(blockBounds(LOW,KAXIS)  .eq. gr_kmin .and. face .eq. 5 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) + glDZ
              endif

              if(blockBounds(HIGH,KAXIS) .eq. gr_kmax .and. face .eq. 6 ) then
                  oneRay(itposz) = raysRealProp(iposz,i) - glDZ
              endif
            end select
        else
          oneRay(itposx)  = raysRealProp(iposx,  i)
          oneRay(itposy)  = raysRealProp(iposy,  i)
          oneRay(itposz)  = raysRealProp(iposz,  i)
        endif

        oneRay(itinfo)  = raysIntProp(ihlev,i)
        oneRay(itblk)   = raysIntProp(iblk,i)
        oneRay(itid)    = raysIntProp(isid,i)

        call ph_sendRay(oneRay, blockproc, isGone)

! try again, and again, again and again, possibly again
! this only happens if send buffer is full and has to be cleared for new ray
        do while (.not. isGone)
          call ph_progressComm(.false.,leave,.false.)
          call ph_sendRay(oneRay, blockproc, isGone)
        end do

! fill new free slot
        if(ph_localRays .gt. 1 .and. i .ne. ph_localRays) then 
! move last ray data to current slot
          raysintProp (:,i) = raysIntProp (:,ph_localRays)
          raysRealProp(:,i) = raysRealProp(:,ph_localRays)

! empty last particle slot
          raysintProp (:,ph_localRays) = -1
          raysRealProp(:,ph_localRays) = -1
        else 
! empty particle slot as it is now copied to destBuf
          raysIntProp(:,i) = -1
          raysRealProp(:,i) = -1
        endif

        ph_localRays = ph_localRays - 1

        exit
      endif

! not a refined region, you could say quite unsophisticated, also on local domain
! new block id
      raysIntProp(iblk,i) = newBID

! check for periodicity and adjust accordingly 
! this is after the block traversal, assumes even number of blocks
! periodic boundary condition check
! if ray traversed 1.5 times in any periodic direction stop it
! 1.499 leaves some slop to not traverse neighbour block if distance is exact
!      if(ph_periodic) then !global switch
      blockBounds = bnd_box(:,:,newBID)        
      norm = raysRealProp(ivelx:ivelz,i)

! maybe nested cases are better? 
! x 
      if(xperiodic) then
! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
        if(blockBounds(LOW,IAXIS) .eq. gr_imin) then
! ray points into domain and entered face is global boundary
          if( (norm(1) .gt. 0.0) .and. (facedir(1) .eq. 3)) then
            raysRealProp(iposx,i) = raysRealProp(iposx,i) - glDX
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
          endif
        endif

        if((blockBounds(HIGH,IAXIS) .eq. gr_imax)) then
          if((norm(1) .lt. 0.0) .and. (facedir(1) .eq. 1)  ) then
            raysRealProp(iposx,i) =  raysRealProp(iposx,i) + glDX
! recalc diff
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
! overwrite critical position
          endif
        endif
      endif
! y
      if(yperiodic) then

! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
        if(blockBounds(LOW,JAXIS) .eq. gr_jmin) then
! ray points into domain and entered face is global boundary
          if( (norm(2) .gt. 0.0) .and. (facedir(2) .eq. 3)) then
            raysRealProp(iposy,i) = raysRealProp(iposy,i) - glDY
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
          endif
        endif

        if((blockBounds(HIGH,JAXIS) .eq. gr_jmax)) then
          if((norm(2) .lt. 0.0) .and. (facedir(2) .eq. 1)  ) then
            raysRealProp(iposy,i) = raysRealProp(iposy,i) + glDY
! recalc diff
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
          endif
        endif
      endif

! z
      if(zperiodic) then
! check if it passed periodic boundary 
! rays from other cores are flipped in pt_reconstructRays
        if(blockBounds(LOW,KAXIS) .eq. gr_kmin) then
! ray points into domain and entered face is global boundary
          if( (norm(3) .gt. 0.0) .and. (facedir(3) .eq. 3)) then
            raysRealProp(iposz,i) = raysRealProp(iposz,i) - glDZ
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
          endif
        endif

        if((blockBounds(HIGH,KAXIS) .eq. gr_kmax)) then
          if((norm(3) .lt. 0.0) .and. (facedir(3) .eq. 1)  ) then
            raysRealProp(iposz,i) = raysRealProp(iposz,i) + glDZ
! recalc diff
            diff = raysRealProp(ivelx:ivelz,i)*raysRealProp(irad,i)+raysRealProp(iposx:iposz,i)
          endif
        endif
      endif
!      endif
    enddo ! photon traversal loop

! go to next particle
    i = i + 1

! counter until it's time to check for communication
! counts terminated and domain leaving photons
    raysUntilComm = raysUntilComm - 1

!============================================
!==== communication check
!============================================
! in principle could use any field instead of nion
    if(raysIntProp(iblk,i) .lt. 0) then
! we might end up here with rays left locally, by refilling and moving the last ray to 
! previous positions, if so go back to main loop
      if(ph_localRays .gt. 0 ) then 
        activeRays = .true.
! start from the top
        i = 1
      else 
! this stays false if no photons are going to be received in communication step
        activeRays = .false.
! absolute index, should also work with stopped rays
        p_countFinal = ph_localRays

! check globally for active rays, if there are some, some might be received by this domain
! so keep inside this loop by setting activeRays to true
        call ph_progressComm(.true.,leave,.true.)
! if we received rays make this routine active again and start iterating over new arrivals
        activeRays = .not. leave

        if(activeRays) then
! start from the top
          i = 1
        else
! rewind i index to point to same slot next iteration, this keeps the loop active and waiting
          i = i - 1
! rewind comm counter otherwise each empty iteration decrements
          raysUntilComm = raysUntilComm + 1
        endif
! use counter as flag
      endif ! leftover check
    endif ! 
  enddo
  return
end subroutine pt_advanceRaysPoint
