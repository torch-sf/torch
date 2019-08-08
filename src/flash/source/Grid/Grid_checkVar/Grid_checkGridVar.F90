!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! subroutine Grid_checkTemps(maxTemp, x, y, z, ion_frac)
!!
!! Checks the temperature on all the blocks.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Grid_checkGridVar(gridInd, limit, gridVal, x, y, z, ion_frac)

use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
			   Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_fillGuardCells


implicit none

#include "constants.h"
#include "Flash.h"

integer, intent(in)		:: gridInd ! Index of grid variable
character(len=3), intent(in)	:: limit ! max or min
real, intent(out) 		:: gridVal, x, y, z, ion_frac
logical                         :: overLimit
real                            :: xx, yy, zz, maxStopVal
real				:: f1, f2, f3
real, allocatable               :: xCoord(:), yCoord(:), zCoord(:)
real, pointer                   :: solnData(:,:,:,:)
integer           		:: thisBlock, blockID, xSizeCoord, & 
				   ySizeCoord, zSizeCoord, &
				   blkLimits(2,3), blkLimitsGC(2,3), &
				   i, j, k, &
				   blockList(MAXBLOCKS), blockCount
		     
! Not yet used. Maybe later.
!  if (.not. present(maxStopVal_in)) then
!    maxStopVal = -1.0
!  else
!    maxStopVal = maxStopVal_in
!  end if

  if (limit .eq. 'min') then
    gridVal  = 1d99
  else if (limit .eq. 'max') then
    gridVal  = -1.0d0
  else
    print*, "[Grid_checkGridVar]: Limit not recognized. limit =", limit
  end if
  
  ion_frac = -1.0d0
  f1 = -1.0d0
  f2 = -1.0d0
  f3 = -1.0d0

  call Grid_fillGuardCells(CENTER, ALLDIR)  
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)

  do thisBlock = 1, blockCount
    blockID = blockList(thisBlock)

! Get a pointer to solution data 
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

! Get a pointer to solution data 
    xSizeCoord = blkLimitsGC(HIGH,IAXIS)
    ySizeCoord = blkLimitsGC(HIGH,JAXIS)
    zSizeCoord = blkLimitsGC(HIGH,KAXIS)

! allocate space for dimensions
    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))

    call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)     

! loop over all zones in block
! NOTE: if you loop over guard cells, be sure you filled the cells first!

    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
      zz = zCoord(k)
      do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
          xx = xCoord(i)

	  overLimit = .false.
	  
	  if (gridInd == VELX_VAR) then

	      f1 = max(f1,abs(solnData(gridInd, i, j, k)))
	      f2 = max(f2,abs(solnData(VELY_VAR, i, j, k)))
	      f3 = max(f3,abs(solnData(VELZ_VAR, i, j, k)))

	  else
 
	    if ((solnData(gridInd, i, j, k) .gt. gridVal) .and. (limit == 'max')) then
	      overLimit = .true.
	    else if ((solnData(gridInd, i, j, k) .lt. gridVal) .and. (limit == 'min')) then
	      overLimit = .true.
	    else
	      overLimit = .false.
	    end if
	  
	    if (overLimit) then
	        gridVal  = solnData(gridInd, i, j, k)
#ifdef IHP_SPEC
                ion_frac = solnData(IHP_SPEC, i, j, k)
#endif
	        x = xx
	        y = yy
	        z = zz
	    end if
	    
	  end if  
	end do
      end do
    end do

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  call Grid_releaseBlkPtr(blockID, solnData)   
   
  end do
  
  if (gridInd == VELX_VAR) then
      if (f1 .gt. gridVal) then
        write(*,'(3ES13.3)') f1, f2, f3
        if (f2 .ne. 0.d0) write(*,'(3ES13.3)') f1/f2, f1/f3, f2/f3
	call flush(6)
        if (f2 .ne. 0.d0) then 
	  gridVal = f1/f2
	  if ( f1/f2 .gt. 1.01d0) then
	    print*, f1/f2
	  end if
	end if
      end if
  end if

  return
  
end subroutine Grid_checkGridVar
