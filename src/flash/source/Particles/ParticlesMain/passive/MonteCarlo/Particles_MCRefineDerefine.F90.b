!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!
!!
!! SYNOPSIS
!!
!!
!!  
!! DESCRIPTION 
!! goes through to blocks to refine and derefine and 
!! moves any MC tracer particle to a position inside
!! the new block
!! the actual change in blockID is done outside this routine
!!
!! ARGUMENTS
!! 
!! NOTES
!! TODO THIS ASSUMES THAT ALL PASSIVE PARTICLES ARE ADVECTED WITH THE MC METHOD
!!
!!
!!***

subroutine Particles_MCRefineDerefine()

  use tree, ONLY : refine, derefine, stay, lnblocks, lrefine, bnd_box, &
									 coord
  use Particles_data, ONLY : particles, pt_typeInfo
	use pt_advanceMonteCarlo_data
	use Grid_data, ONLY : gr_delta
  use Particles_data, ONLY: pt_xmin, pt_xmax, pt_ymin, &
														pt_ymax, pt_zmin, pt_zmax
  use mtmod

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,iref
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask

	integer :: lb, startID, stopID, partPos, j,i
	integer :: child, rnd 
	integer :: idx, idy, idz
	integer :: p1, p2, p3
  real, dimension(NDIM)	     :: gridspace
	real :: blockBounds(LOW:HIGH,1:MDIM)

! for testing
	real :: lowBx,lowBy,lowBz

! indices for photon particles in data structure
! find start integer for passive particles
! first check if they are defined
	if(PASSIVE_PART_TYPE .gt. 0) then
    startID = pt_typeInfo(PART_TYPE_BEGIN,PASSIVE_PART_TYPE)
    stopID  = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) + startID - 1
! check if there are passive particles that are advected with the MC method
		if(pt_typeInfo(PART_ADVMETHOD,PASSIVE_PART_TYPE) .ne. MONTECARLO) then
			return
		endif
	endif

	j = startID
! loop over all blocks and check if they are refined or derefined 
  print*,'blockloop',lnblocks
	print*,startID,stopID

	idx = 2

	if(stopID .gt. 0) then 
  	do lb = 1,lnblocks

		  if( refine(lb) .or. derefine(lb)) then
!			  print*,'change in ref',lb
!			  print*,'look from',j,stopID
! loop over tracer particle block property
			    do i = j, stopID

					  if(particles(BLK_PART_PROP,i) .eq. lb) then
!						  print*,'found particle in refinement changing block',lb

! get not updated blk data

							if(refine(lb)) then
								print*,'refine',particles(TAG_PART_PROP,i)
! throw dice to see which child block gets the tracer particle
								child = floor(grnd()*8) + 1
						    gridspace = gr_delta(1:MDIM,lrefine(lb))

								idx = particles(VELX_PART_PROP,i)
								idy = particles(VELY_PART_PROP,i)
								idz = particles(VELZ_PART_PROP,i)

								print*,'old ids',idx,idy,idz
								print*,particles(POSX_PART_PROP,i),gridspace(1)*(idx-0.5)+bnd_box(LOW,IAXIS,lb)
								print*,particles(POSY_PART_PROP,i),gridspace(2)*(idy-0.5)+bnd_box(LOW,JAXIS,lb)
								print*,particles(POSZ_PART_PROP,i),gridspace(3)*(idz-0.5)+bnd_box(LOW,KAXIS,lb)
						    gridspace = gr_delta(1:MDIM,lrefine(lb))/4d0

								lowbx = bnd_box(LOW,IAXIS,lb)
								lowby = bnd_box(LOW,JAXIS,lb)
								lowbz = bnd_box(LOW,KAXIS,lb)

! this maps index to new indices, always lower one
! from 1->1,2;  2->3,4; 3->5,6; 4->7,8 for le 4
! from 5->1,2;  6->3,4; 7->5,6; 8->7,8 for ge 5
							  if     (IDx .le. 4 .and. IDy .le. 4 .and. IDz .le. 4 ) then
! x 1..4 to 1...8, y 1..4 to 1...8, z 1..4 to 1...8
									idx = idx*2-1
									idy = idy*2-1
									idz = idz*2-1
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .ge. 5 .and. IDy .le. 4 .and. IDz .le. 4 ) then
! x 5..8 to 1...8, y 1..4 to 1...8, z 1..4 to 1...8
! map x 	
									idx = (idx-4)*2-1
									idy = idy*2-1
									idz = idz*2-1
								  lowbx = bnd_box(LOW,IAXIS,lb)+(bnd_box(HIGH,IAXIS,lb)-bnd_box(LOW,IAXIS,lb))/2d0
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .ge. 5 .and. IDy .ge. 5 .and. IDz .le. 4 ) then
! x 5..8 to 1...8, y 5..8 to 1...8, z 1..4 to 1...8
									idx = (idx-4)*2-1
									idy = (idy-4)*2-1
									idz = idz*2-1
								  lowbx = bnd_box(LOW,IAXIS,lb)+(bnd_box(HIGH,IAXIS,lb)-bnd_box(LOW,IAXIS,lb))/2d0
								  lowby = bnd_box(LOW,JAXIS,lb)+(bnd_box(HIGH,JAXIS,lb)-bnd_box(LOW,JAXIS,lb))/2d0
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .le. 4 .and. IDy .ge. 5 .and. IDz .le. 4 ) then
! x 1..4 to 1...8, y 5..8 to 1...8, z 1..4 to 1...8
									idx = idx*2-1
									idy = (idy-4)*2-1
									idz = idz*2-1
								  lowby = bnd_box(LOW,JAXIS,lb)+(bnd_box(HIGH,JAXIS,lb)-bnd_box(LOW,JAXIS,lb))/2d0
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .le. 4 .and. IDy .le. 4 .and. IDz .ge. 5 ) then
! x 1..4 to 1...8, y 1..4 to 1...8, z 5..8 to 1...8
									idx = idx*2-1
									idy = idy*2-1
									idz = (idz-4)*2-1
								  lowbz = bnd_box(LOW,KAXIS,lb)+(bnd_box(HIGH,KAXIS,lb)-bnd_box(LOW,KAXIS,lb))/2d0
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .ge. 5 .and. IDy .le. 4 .and. IDz .ge. 5 ) then
! x 5..8 to 1...8, y 1..4 to 1...8, z 5..8 to 1...8
									idx = (idx-4)*2-1
									idy = idy*2-1
									idz = (idz-4)*2-1							
								  lowbx = bnd_box(LOW,IAXIS,lb)+(bnd_box(HIGH,IAXIS,lb)-bnd_box(LOW,IAXIS,lb))/2d0
								  lowbz = bnd_box(LOW,KAXIS,lb)+(bnd_box(HIGH,KAXIS,lb)-bnd_box(LOW,KAXIS,lb))/2d0
									print*,'mapped ids',idx,idy,idz
							  else if(IDx .ge. 5 .and. IDy .ge. 5 .and. IDz .ge. 5 ) then
! x 5..8 to 1...8, y 5..8 to 1...8, z 5..8 to 1...8
									idx = (idx-4)*2-1
									idy = (idy-4)*2-1
									idz = (idz-4)*2-1
								  lowbx = bnd_box(LOW,IAXIS,lb)+(bnd_box(HIGH,IAXIS,lb)-bnd_box(LOW,IAXIS,lb))/2d0
								  lowby = bnd_box(LOW,JAXIS,lb)+(bnd_box(HIGH,JAXIS,lb)-bnd_box(LOW,JAXIS,lb))/2d0
								  lowbz = bnd_box(LOW,KAXIS,lb)+(bnd_box(HIGH,KAXIS,lb)-bnd_box(LOW,KAXIS,lb))/2d0

									print*,'mapped ids',idx,idy,idz
							  else if(IDx .le. 4 .and. IDy .ge. 5 .and. IDz .ge. 5 ) then
! x 1..4 to 1...8, y 5..8 to 1...8, z 5..8 to 1...8
									idx = idx*2-1
									idy = (idy-4)*2-1
									idz = (idz-4)*2-1
								  lowby = bnd_box(LOW,JAXIS,lb)+(bnd_box(HIGH,JAXIS,lb)-bnd_box(LOW,JAXIS,lb))/2d0
								  lowbz = bnd_box(LOW,KAXIS,lb)+(bnd_box(HIGH,KAXIS,lb)-bnd_box(LOW,KAXIS,lb))/2d0

									print*,'mapped ids',idx,idy,idz
							  endif

! this change physical position so the particle will be at a center of a refined zone
					      select case (child)
									case(1)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

									case(2)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

										idx = idx + 1
									case(3)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

										idy = idy + 1
									case(4)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)

										idx = idx + 1
										idy = idy + 1
									case(5)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

										idz = idz + 1
									case(6)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

										idx = idx + 1
										idz = idz + 1
									case(7)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

										idy = idy + 1
										idz = idz + 1
									case(8)
										particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)

										idx = idx + 1
										idy = idy + 1
										idz = idz + 1
								end select			

								print*,particles(POSX_PART_PROP,i)/(gridspace(1)*2d0)
								print*,particles(POSY_PART_PROP,i)/(gridspace(2)*2d0)
								print*,particles(POSZ_PART_PROP,i)/(gridspace(3)*2d0)

								print*,particles(POSX_PART_PROP,i),gridspace(1)*2d0*(idx-0.5)+lowbx
								print*,particles(POSY_PART_PROP,i),gridspace(2)*2d0*(idy-0.5)+lowby
								print*,particles(POSZ_PART_PROP,i),gridspace(3)*2d0*(idz-0.5)+lowbz

								particles(VELX_PART_PROP,i) = idx
								particles(VELY_PART_PROP,i) = idy
								particles(VELZ_PART_PROP,i) = idz
							endif

							if(derefine(lb)) then
							  gridspace = gr_delta(1:MDIM,lrefine(lb))

								print*,particles(POSX_PART_PROP,i)/gr_delta(1,lrefine(lb))
								print*,particles(POSY_PART_PROP,i)/gr_delta(2,lrefine(lb))
								print*,particles(POSZ_PART_PROP,i)/gr_delta(3,lrefine(lb))
								print*,gridspace

							  gridspace = gr_delta(1:MDIM,lrefine(lb))/2d0

								idx = particles(VELX_PART_PROP,i)
								idy = particles(VELY_PART_PROP,i)
								idz = particles(VELZ_PART_PROP,i)

!								childnodes  = child(1,1:8,parent(1,lb))
! this might not work if parent is on different domain
								print*,'derefine',particles(TAG_PART_PROP,i),derefine(lb),refine(lb),stay(lb)
!						   	print*,(pt_xmax-pt_xmin)/bnd_box(LOW,IAXIS,lb),mod((pt_xmax-pt_xmin)/bnd_box(LOW,IAXIS,lb),2d0)
!						   	print*,(pt_ymax-pt_ymin)/bnd_box(LOW,JAXIS,lb),mod((pt_ymax-pt_ymin)/bnd_box(LOW,JAXIS,lb),2d0)
!						   	print*,(pt_zmax-pt_zmin)/bnd_box(LOW,KAXIS,lb),mod((pt_zmax-pt_zmin)/bnd_box(LOW,KAXIS,lb),2d0)

									p1 = mod(floor((coord(1,lb)-pt_xmin)/(bnd_box(HIGH,IAXIS,lb)-bnd_box(LOW,IAXIS,lb)))+1,2)
									p2 = mod(floor((coord(2,lb)-pt_ymin)/(bnd_box(HIGH,JAXIS,lb)-bnd_box(LOW,JAXIS,lb)))+1,2)
									p3 = mod(floor((coord(3,lb)-pt_zmin)/(bnd_box(HIGH,KAXIS,lb)-bnd_box(LOW,KAXIS,lb)))+1,2)

									if(p1 .eq. 1) then
										print*,'odd block x'
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
! left of center
										if(mod(idx,2) .eq. 0) then 
											particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										endif
! right of center
										if(mod(idx,2) .eq. 1) then 
											particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										endif

										idx = ceiling(idx/2d0)
									endif

									if(p1 .eq. 0) then
										print*,'even block x'

										if(mod(idx,2) .eq. 0) then
											particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)
										endif

										if(mod(idx,2) .eq. 1) then 
											particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
										endif

										idx = ceiling(idx/2d0)+4
									endif

									if(p2 .eq. 1) then
										print*,'odd block y'
! 1+,2- -> 1; 3+,4- -> 2; 5+,6- -> 3; 7+,8- -> 4
										if(mod(idy,2) .eq. 0) then 
											particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										endif

										if(mod(idy,2) .eq. 1) then 
											particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										endif
										idy = ceiling(idy/2d0)
									endif

									if(p2 .eq. 0) then
										print*,'even block y'
										if(mod(idy,2) .eq. 0) then 
											particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)
										endif

										if(mod(idy,2) .eq. 1) then 
											particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
										endif
										idy = ceiling(idy/2d0)+4
									endif

									if(p3 .eq. 1) then
										print*,'odd block z'
										if(mod(idz,2) .eq. 0) then
											particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)
										endif

										if(mod(idz,2) .eq. 1) then 
											particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
										endif

										idz = ceiling(idz/2d0)
									endif

									if(p3 .eq. 0) then
										print*,'even block z'
										if(mod(idz,2) .eq. 0) then 
											particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)
										endif

										if(mod(idz,2) .eq. 1) then 
											particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
										endif
										idz = ceiling(idz/2d0)+4
									endif

									particles(VELX_PART_PROP,i) = idx
									particles(VELY_PART_PROP,i) = idy
									particles(VELZ_PART_PROP,i) = idz

									print*,particles(POSX_PART_PROP,i)/(gridspace(1)*4d0)
									print*,particles(POSY_PART_PROP,i)/(gridspace(2)*4d0)
									print*,particles(POSZ_PART_PROP,i)/(gridspace(3)*4d0)

!				 					if(abs(particles(POSX_PART_PROP,i)-(gridspace(1)*(idx-0.5)+bnd_box(LOW,IAXIS,newBID))) .gt. 1e6 .or. & 
!				 					   abs(particles(POSY_PART_PROP,i)-(gridspace(2)*(idy-0.5)+bnd_box(LOW,JAXIS,newBID))) .gt. 1e6 .or. &
!										 abs(particles(POSZ_PART_PROP,i)-(gridspace(3)*(idz-0.5)+bnd_box(LOW,KAXIS,newBID))) .gt. 1e6) then

!									  print*,'position wrong'
!										print*,particles(POSX_PART_PROP,i)-(gridspace(1)*(idx-0.5)+bnd_box(LOW,IAXIS,newBID))
!										print*,particles(POSY_PART_PROP,i)-(gridspace(2)*(idy-0.5)+bnd_box(LOW,JAXIS,newBID))
!										print*,particles(POSZ_PART_PROP,i)-(gridspace(3)*(idz-0.5)+bnd_box(LOW,KAXIS,newBID))						
!										stop
!									endif

!								stop

! construct position in parent block, do not count on tree information, as that might not be local
! take global domain boundaries and divide by refinement level

!								particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!								particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!								particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)


! this maps index to new indices
! from 1,2->1; 3,4->2; 5,6->3; 7,8->4 for le 4
! from 1,2->5; 3,4->6; 5,6->7; 7,8->8 for ge 5
!							  if     (IDx .le. 4 .and. IDy .le. 4 .and. IDz .le. 4 ) then
!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) - gridspace(1)*(idx-0.5d0)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) - gridspace(2)*(idy-0.5d0)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) - gridspace(3)*(idz-0.5d0)

! 1 should be + gridspace(ID)/2d0
! 2 should be - gridspace(ID)/2d0
! this might not be the savest way
!									idx = ceiling(idx/2d0)
!									idy = ceiling(idy/2d0)
!									idz = ceiling(idz/2d0)

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)*idx
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)*idy
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)*idz

!							  else if(IDx .ge. 5 .and. IDy .le. 4 .and. IDz .le. 4 ) then
!									idx = ceiling(idx/2d0)+4
!									idy = ceiling(idy/2d0)
!									idz = ceiling(idz/2d0)

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .ge. 5 .and. IDy .ge. 5 .and. IDz .le. 4 ) then
!									idx = ceiling(idx/2d0)+4
!									idy = ceiling(idy/2d0)+4
!									idz = ceiling(idz/2d0)

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .le. 4 .and. IDy .ge. 5 .and. IDz .le. 4 ) then
!									idx = ceiling(idx/2d0)
!									idy = ceiling(idy/2d0)+4
!									idz = ceiling(idz/2d0)

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .le. 4 .and. IDy .le. 4 .and. IDz .ge. 5 ) then
!									idx = ceiling(idx/2d0)
!									idy = ceiling(idy/2d0)
!									idz = ceiling(idz/2d0)+4

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .ge. 5 .and. IDy .le. 4 .and. IDz .ge. 5 ) then
!									idx = ceiling(idx/2d0)+4
!									idy = ceiling(idy/2d0)
!									idz = ceiling(idz/2d0)+4

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .ge. 5 .and. IDy .ge. 5 .and. IDz .ge. 5 ) then
!									idx = ceiling(idx/2d0)+4
!									idy = ceiling(idy/2d0)+4
!									idz = ceiling(idz/2d0)+4

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  else if(IDx .le. 4 .and. IDy .ge. 5 .and. IDz .ge. 5 ) then
!									idx = ceiling(idx/2d0)
!									idy = ceiling(idy/2d0)+4
!									idz = ceiling(idz/2d0)+4

!									particles(VELX_PART_PROP,i) = idx
!									particles(VELY_PART_PROP,i) = idy
!									particles(VELZ_PART_PROP,i) = idz

!									particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + gridspace(1)
!									particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + gridspace(2)
!									particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + gridspace(3)
!							  endif

!							  select case (idx)
!									case(1)										
!									case(2)
!									case(3)
!									case(4)
!									case(5)
!									case(6)
!									case(7)
!									case(8)
!								end select
							endif
					  endif

! exit if we find bigger block number
            if(particles(BLK_PART_PROP,i) .gt. lb) then
!						  print*,'outside block',i,particles(BLK_PART_PROP,i)
  					  j = i - 1 
!							print*,'new index ',j
						  exit
					  endif
! the array is sorted in blockIDs, so move pointer to current blockID
			    enddo
			  endif
! check if we have any particles left to check
				if (lb .gt. particles(BLK_PART_PROP,stopid)) then
					print*,'no particles in ref blocks, exit'
					exit
				endif
		  enddo
	  endif
  return
end subroutine Particles_MCRefineDerefine
