
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! subroutine Particle_yieldInjection
!!!
!!! Authors: Eric Andersson
!!!          American Museum of Natural History
!!!          Fall 2023
!!!
!!! A routine for metal injection of SN. Adapted from energy injection.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine Particles_yieldInjection(itracer, injectYieldIn, injectMassIn, xloc, yloc, zloc)

!#define DEBUG2

  use Grid_data, ONLY: gr_meshComm, gr_meshMe
  
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getMinCellSize, &
      Grid_notifySolnDataUpdate, Grid_getBlkIDFromPos
      
  use Eos_interface, ONLY : Eos_wrapped
  
  use pt_enerInjInterface, only : overlap, sphere_and_cell_frac
  
  
  use tree, ONLY: nodetype, coord, bsize, lnblocks, refine, derefine, stay
  
#include "Flash.h"
#include "constants.h"
  
  implicit none
  
#include "Flash_mpi.h"
  
  integer, parameter :: dp = kind(1.d0)
  integer :: ierr
  
  integer, intent(in) :: itracer ! Tracer field indices starts counting at 1.
  real(dp),intent(in) :: injectMassIn, xloc, yloc, zloc
  real(dp),intent(in) :: injectYieldIn
  
  logical :: iHaveInjectBlk
  logical :: snap_to_grid
  
  real(dp) :: loc(3)
  real(dp) :: cell_top(3), cell_bot(3)
  real(dp) :: dVol
  real(dp) :: x, y, z
  
  real(dp),save :: injectRadius, delta(3), SNdelta(3), min_delta
  
  real(dp) :: overlap_frac,  sumOverlap, oldDens, dDens, newDens
  real(dp) :: xcoll, ycoll, zcoll, d2coll
  
  integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), blockID, procID
  integer :: i, j, k, n
  
  real(dp), pointer, dimension(:,:,:,:) :: solndata
  
  ! Indices: blockID, i, j, k, [dimension]
  real(dp), allocatable, dimension(:,:,:,:) :: injectDataOverlap
  integer, allocatable, dimension(:) :: localInjectBlocks
  
  integer :: injBlkNum
  real(dp) :: blkCtr(3), blkSize(3)  ! code requires MDIM=3
  
  real(dp) :: injectYield, injectMass
  real(dp) :: oldTracerField, newTracerField, dTracerField

#ifndef TRACER_FIELDS
  print*, "Function Particles_yieldInjection cannot be called without tracer fields. Something went wrong!"
#endif

  iHaveInjectBlk = .false.
  sumOverlap = 0.0_dp

  snap_to_grid = .false. ! For testing / debugging.

  call Grid_getBlkIDFromPos([xloc, yloc, zloc] ,blockID ,procID, gr_meshComm)

  ! Note only the proc that owns this block has the proper cell size... other procs
  ! may have a block with this number but its the wrong block (and maybe a parent, 
  ! wrong size, etc).

  if (gr_meshMe .eq. procID) then
    call Grid_getDeltas(blockID,SNdelta)
  end if

  ! Send the proper cell size to all procs.
  call MPI_Bcast(SNdelta, 3, MPI_DOUBLE_PRECISION, procID, gr_meshComm, ierr)
  call Grid_getMinCellSize(min_delta)

  ! Place the star in the center of a cell?
  if (snap_to_grid) then
    loc = floor([xloc, yloc, zloc]) + 0.5_dp*SNdelta
  else
    loc = [xloc, yloc, zloc]    
  end if

  injectRadius = 3.0_dp*minval(SNdelta)

#ifdef DEBUG_ENERGY
  if (gr_meshMe == 0) &
    write(*, '(A,ES13.3e3)') "SN injection radius / dx = ", injectRadius/SNdelta(1)
#endif

#ifdef DEBUG2
  write(*, '(A,ES13.3e3, I4)') " injection radius = ", injectRadius, gr_meshMe
  print*, "loc =", loc, gr_meshMe
#endif

  ! count # of blocks which are at least partially within injectRadius, check that
  ! they are maximally refined
  injBlkNum = 0
  do blockID = 1, lnblocks
    if(nodetype(blockID) == LEAF) then
      ! exact collision detection for sphere and rectangular prism
      ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
      call Grid_getBlkCenterCoords(blockID,blkCtr)
      call Grid_getBlkPhysicalSize(blockID,blkSize)
      ! point within block that is closest to SN location
      xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
      ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
      zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))
      if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadius**2) then
        injBlkNum = injBlkNum + 1
        iHaveInjectBlk = .true.
      end if
    end if
  end do

  ! build array of block IDs for all blocks to be injected. Get their center
  ! distances from the injection star. Then for each cell in each block, calculate
  ! its overlap with the injection sphere and store the value.
  if (iHaveInjectBlk) then

#ifdef DEBUG2
    print *, "Found", injBlkNum, "injection blocks on proc ", gr_meshMe
#endif

    allocate(localInjectBlocks(injBlkNum))
    allocate(injectDataOverlap(injBlkNum,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI, &
                                                            GRID_KLO:GRID_KHI))

    localInjectBlocks = 0
    injectDataOverlap = 0.0d0

#ifdef DEBUG2
    print *, "Allocations done"
#endif

    n = 1
    do blockID = 1, lnblocks
      if(nodetype(blockID) == LEAF) then
        ! exact collision detection for sphere and rectangular prism
        ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
        call Grid_getBlkCenterCoords(blockID,blkCtr)
        call Grid_getBlkPhysicalSize(blockID,blkSize)
        ! point within block that is closest to SN location
        xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
        ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
        zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))
        if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadius**2) then
          localInjectBlocks(n) = blockID
          n = n + 1
        end if
      end if
    end do

#ifdef DEBUG2
    print *, "Found injection blocks:", localInjectBlocks, "on proc", gr_MeshMe
#endif
 
    do n = 1, injBlkNum
      blockID = localInjectBlocks(n)
      call Grid_getDeltas(blockID,delta)
      call Grid_getBlkPtr(blockID, solndata)
      do k = GRID_KLO, GRID_KHI
        do j = GRID_JLO, GRID_JHI
          do i = GRID_ILO, GRID_IHI
            
            ! since we have checked that all cells are refined, use
            ! mindelta
            x = (i - NGUARD - NXB/2.0 - 0.5)*delta(1) + coord(1,blockID)
            y = (j - NGUARD - NYB/2.0 - 0.5)*delta(2) + coord(2,blockID)
            z = (k - NGUARD - NZB/2.0 - 0.5)*delta(3) + coord(3,blockID)
            
            ! exact collision detection for sphere and rectangular prism
            ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
            ! point within cell that is closest to SN location
            xcoll = max(x-0.5*delta(1),min(loc(1),x+0.5*delta(1)))
            ycoll = max(y-0.5*delta(2),min(loc(2),y+0.5*delta(2)))
            zcoll = max(z-0.5*delta(3),min(loc(3),z+0.5*delta(3)))
            d2coll = (xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2
            
            ! is cell outside injection sphere?
            if (d2coll > injectRadius**2) then
              cycle
            end if
            
            ! get overlapping volume of inject sphere and this cell,
            ! modified by a tapered center-weighting within overlap(..)
            cell_bot = [ sign(abs(x) - 0.5*delta(1), x), &
                         sign(abs(y) - 0.5*delta(2), y), &
                         sign(abs(z) - 0.5*delta(3), z) ]
            cell_top = [ sign(abs(x) + 0.5*delta(1), x), &
                         sign(abs(y) + 0.5*delta(2), y), &
                         sign(abs(z) + 0.5*delta(3), z) ]
            call overlap(1, injectRadius, loc, cell_bot, &
                         cell_top, 10, overlap_frac)
            
            if (overlap_frac .gt. 0.0d0) then
              sumOverlap = sumOverlap + overlap_frac
              injectDataOverlap(n,i,j,k) = overlap_frac
            end if
          end do
        end do
      end do
      
      call Grid_releaseBlkPtr(blockID, solndata)
    
    end do ! injBlkNum

#ifdef DEBUG2
    print *, "Calculated overlaps"
#endif

  end if ! iHaveInjectBlk

#ifdef DEBUG2
  print*, "Proc ", gr_meshMe, " about to call MPI with sumOverlap = ", sumOverlap
#endif

  call MPI_ALLREDUCE(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

  injectMass   = injectMassIn 
  injectYield = injectYieldIn

  if (iHaveInjectBlk) then
    do n = 1, injBlkNum
      blockID = localInjectBlocks(n)
      call Grid_getDeltas(blockID,delta)
      dVol = product(delta)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      call Grid_getBlkPtr(blockID, solndata)
      
        do k = GRID_KLO, GRID_KHI
          do j = GRID_JLO, GRID_JHI
            do i = GRID_ILO, GRID_IHI
              
              ! Round off errors in calculations lead to small
              ! changes in thermal energy that have big
              ! consequences in EOS calls. So just skip.
              if (injectDataOverlap(n,i,j,k) .le. 0.0_dp) cycle
                
              dDens   = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
              oldDens = solndata(DENS_VAR,i,j,k)
              newDens = solndata(DENS_VAR,i,j,k) + dDens
            
              ! Compute the tracer field density of material which is added to cell.
              oldTracerField = solndata(MASS_SCALARS_BEGIN+(itracer-1),i,j,k)*oldDens ! Tracer field mass per volume
              if(injectYield.LT.0.0) then
                 ! Old field should remain untouched. Apply old field to injected mass
                 dTracerField = solndata(MASS_SCALARS_BEGIN+(itracer-1),i,j,k)*dDens
              else
                  ! Injected mass will update field, if yields=0 -> field is depleted.
                 dTracerField = injectDataOverlap(n,i,j,k)/sumOverlap*injectYield/dVol
              endif
              newTracerField = oldTracerField + dTracerField
              
              ! Update scalar field to new metallicity
              solndata(MASS_SCALARS_BEGIN+(itracer-1), i, j, k) = newTracerField/newDens
            
            end do
          end do
        end do
      
      call Grid_releaseBlkPtr(blockID, solndata)
      call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
    end do ! injBlkNum
  end if ! iHaveInjectBlk

  call Grid_notifySolnDataUpdate()

  if (iHaveInjectBlk) then
    deallocate(localInjectBlocks)
    deallocate(injectDataOverlap)

#ifdef DEBUG2
    print *, "Deallocating done for proc", gr_meshMe
#endif

  endif
  
  call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)

#ifdef DEBUG2
  print *, "Exiting Particles_yieldInjection for proc", gr_meshMe
#endif

end subroutine Particles_yieldInjection
