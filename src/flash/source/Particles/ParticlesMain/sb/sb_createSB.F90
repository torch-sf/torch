!!
!! NAME
!!
!!   sb_createSB
!!
!! SYNOPSIS
!!
!!  inserts a passive particles as SB
!!
!! ARGUMENTS
!! x,y,z : position of SB
!! time  : current simulation time   
!!
!!
!! DESCRIPTION
!!
!!
!!***

subroutine sb_createSB(x,y,z,time,dt,new)

	use Particles_Data, only: pt_maxPerProc, pt_numLocal, particles, pt_meshMe, pt_typeInfo
	use tree, only: bnd_box
	use gr_ptData, ONLY: gr_ptBlkList, gr_ptBlkCount
	use SB_data
  use mtmodSN

	use Grid_interface, ONLY: Grid_getListOfBlocks, Grid_getPointData

	use SN_data, only: he_SNmapToGrid
	use Grid_data, ONLY : gr_delta
  use Driver_interface, ONLY : Driver_abortFlash

	use tree, ONLY : lrefine, bnd_box

	implicit none

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

	real,  intent(in)			  :: x,y,z, time, dt
	logical,  intent(inout) :: new
	real  :: dxhalf, xx, yy, zz, spin, vSB
	real, dimension(MDIM)  :: newlow, pos, vel
	integer, dimension(MDIM)  :: ind
	logical  :: local
	logical	 :: getGuardCells = .true.
	integer	 :: j, startID, nSB, k, i, l, m, p,nsn

	real, pointer, dimension(:,:,:,:)	:: solnData
	integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
	real, dimension(MDIM) :: del
	real,dimension(NPART_PROPS) :: tmpP
	integer :: xSizeCoord, ySizeCoord, zSizeCoord, zoneindex, blockID
	real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
	real :: blockBounds(LOW:HIGH,1:MDIM)
	real, allocatable,dimension(:) :: dx, dy, dz
	! spacing
  real, dimension(NDIM)	:: gridspace
  integer :: blockList(MAXBLOCKS)
  integer :: tagcounter

! before anything is done find out if SB is local 
	local = .false.

! as random generator is in lockstep all domains have same position info
	pos = (/x,y,z/)
  call Grid_getListOfBlocks(LEAF,blockList,gr_ptBlkCount)
  call gr_findBlock(blockList,gr_ptBlkCount,pos,blockID)

	if(blockID .gt. -1) then
    gridspace = gr_delta(1:MDIM,lrefine(blockID))
    newlow = bnd_box(LOW,:,blockID)
		local = .true.
	endif

! find nearest zone position
	if(local) then
	  ind(1) = floor((x-newlow(1)) / gridspace(1))
!	  pos(1) = gridspace(1)*(zoneindex + 0.5d0)

	  ind(2) = floor((y-newlow(2)) / gridspace(2))
!	  pos(2) = gridspace(2)*(zoneindex + 0.5d0)

	  ind(3) = floor((z-newlow(3)) / gridspace(3))
!	  pos(3) = gridspace(3)*(zoneindex + 0.5d0)

! look inside block for velocity data, don't do CIC or averaging
!		  call Grid_getBlkPtr(blockID, solnData)
!		  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC) !indices for the interior zones and all zones including guard zones
		call Grid_getPointData(blockID, CENTER, VELX_VAR, INTERIOR, ind, vel(1))
		call Grid_getPointData(blockID, CENTER, VELY_VAR, INTERIOR, ind, vel(2))
		call Grid_getPointData(blockID, CENTER, VELZ_VAR, INTERIOR, ind, vel(3))

	endif

	if(.not. local) then
! spin the random number generator to guarantee lockstep
! initial value to guarantee that loop is entered
!		nsn = 0.1
! we don't want any multiple SN from the same SB right away
!		do while(dt/nsn .gt. 0)
!		  call random_number(spin)
			spin = grndSN()

! compare random value to n^2 distribution probability bins, high spin values 
! are needed for many SN in the SB
! .and. comparison looks if it is contained in the probability bin
! assigns number of SN per SB corresponding to probability bin
! loop over all bins
		  do m = 1, sb_nsnmax-sb_nsnmin+1
		    if ( (spin .ge. sb_edge(m) ) .and. ( spin < sb_edge(m+1) ) ) then 
			    nsn = sb_nsnmin+m-1  ! m goes from 1 to nsmax-nsnmin+1
!          exit
		    endif
		  enddo
!		  nsn = sb_life/nsn
!		enddo
! this was done for lockstep don't do anything locally
		new = .false.

! have to call this even for a 0 count to get the MPI right
! update tag offset, this does a couple of MPI calls >(

    call pt_findTagOffset(0,tagcounter)

! out
		return
	endif

! SB is local
	new = .true.

! update particle number
	pt_numLocal = pt_numLocal + 1

	if (pt_numLocal > pt_maxPerProc) then
            call Driver_abortFlash &
 	            ("sb_createSB:  Exceeded max # of particles/processor!")
	endif

!	nsn = 0.1
!	do while(dt/nsn .gt. 0) 
! random number init in SN_init
!	  call random_number(spin)
		spin = grndSN()
! populate SB
!  compare random value to n^2 distribution probability bins, high spin values are needed for many SN in the SB
!  .and. comparison looks if it is contained in the probability bin
!  assigns number of SN per SB corresponding to probability bin
	  do m = 1, sb_nsnmax-sb_nsnmin+1
		  if ( (spin .ge. sb_edge(m) ) .and. ( spin < sb_edge(m+1) ) ) then 
			  nsn = sb_nsnmin+m-1  ! m goes from 1 to nsmax-nsnmin+1
		  endif
	  enddo
! SN rate
	  tmpP(NSN_PART_PROP) = nsn
! to exit while loop, not used further down
!	  nsn =  sb_life/nsn
!	enddo

! creation time, now actually SB ID, look up creation time in output or reconstruct
! as a SN explodes right away this is the same as creation time
	tmpP(TCRT_PART_PROP) = time

! time between SN in SB
	tmpP(TSN_PART_PROP)  = sb_life/(tmpP(NSN_PART_PROP)+1) ! +1 because all SN should go off and we need enough intervals before maximum time is reached

! velocity is initial blk velocity of the surrounding gas
! Limit velocity to an upper value 10 km/s
	if(sb_trackV) then
		vSB = sqrt(sum(vel*vel))

		if(vSB .gt. sb_MaxV ) then
! scale to maximum velocity
			vel = vel*sb_MaxV/vSB
		endif

		tmpP(VELX_PART_PROP) = vel(1)
		tmpP(VELY_PART_PROP) = vel(2)
		tmpP(VELZ_PART_PROP) = vel(3)
	else
		tmpP(VELX_PART_PROP) = 0
		tmpP(VELY_PART_PROP) = 0
		tmpP(VELZ_PART_PROP) = 0
	endif

! position is current SB location
	tmpP(POSX_PART_PROP) = x
	tmpP(POSY_PART_PROP) = y
	tmpP(POSZ_PART_PROP) = z

! find current maximum ID and add 1
  call pt_findTagOffset(1,tagcounter)

! tag and blk id
	tmpP(TAG_PART_PROP)  = tagcounter+1
	tmpP(SBID_PART_PROP) = sb_nSN + 1
	tmpP(BLK_PART_PROP)  = blockID
	tmpP(PROC_PART_PROP) = real(pt_meshMe)
#ifdef TYPE_PART_PROP
	tmpP(TYPE_PART_PROP) = SB_PART_TYPE 
#endif

! update typeInfo
! if already local sb existent, shuffle memory otherwise just append
!//////////////
! as there is a sorting step in Particles_advance befor the time integration, this is nice but too complicated,
! also very hard to debug with additional particle types, so always just append
!//////////////
#ifdef SB_PART_TYPE
	if( pt_typeInfo(PART_LOCAL, SB_PART_TYPE) .gt. 0) then

! just append, will be sorted later on
  	p = pt_numLocal
		particles(:,p) = tmpP
! should be one 
		pt_typeInfo(PART_LOCAL,SB_PART_TYPE) = pt_typeInfo(PART_LOCAL,SB_PART_TYPE) + 1
	else

		p = pt_numLocal
		particles(:,p) = tmpP
		pt_typeInfo(PART_LOCAL,SB_PART_TYPE) = pt_typeInfo(PART_LOCAL,SB_PART_TYPE) + 1
		pt_typeInfo(PART_TYPE_BEGIN,SB_PART_TYPE) = pt_numLocal
	endif 
#endif

	call WriteSBcreate(int(tmpP(SBID_PART_PROP)), int(tmpP(NSN_PART_PROP)), time, tmpP(TSN_PART_PROP),tmpP(POSX_PART_PROP), tmpP(POSY_PART_PROP), &
		tmpP(POSZ_PART_PROP), tmpP(VELX_PART_PROP), tmpP(VELY_PART_PROP), tmpP(VELZ_PART_PROP))

	return

contains

  subroutine WriteSBcreate(SBid, totSN, time, SNinvT, x, y, z, vx, vy, vz)

    use SN_data, ONLY :  he_meshMe

    implicit none

#include "constants.h"
#include "Flash.h"

    real, intent(IN)   :: time, x,y,z, vx, vy, vz, SNinvT
    integer, intent(IN):: SBid,totSN	 ! 1 is I 2 is II and 3 is SN in SB

    integer, parameter :: funit_evol = 15
    character(len=80)  :: outfile = "SBcreate.dat"
    integer            :: i

!   as all processors calculate the SN in lockstep, only the Master process has to output data.

    open(funit_evol, file=trim(outfile), position='APPEND')

    write(funit_evol,'(2(1X,I16),8(1X,ES16.9))') &
			 SBid,  &
	     totSN,	&
	     time,	&
	     SNinvT,&
	     x,  		&
	     y,			&
	     z,			&
	     vx,		&
	     vy,		&
			 vz

    close(funit_evol)

    return
  end subroutine WriteSBcreate
end subroutine sb_createSB
