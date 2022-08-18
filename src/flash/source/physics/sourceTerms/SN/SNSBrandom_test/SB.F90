!!****f* source/physics/sourceTerms/Heat/Heat
!!
!! NAME
!!  
!!  SB
!!
!!
!! SYNOPSIS
!! 
!!  
!! DESCRIPTION
!!
!! Supplies SN.F90 with the positions and explosion times to execute SN in SB
!! two parts 
!! 1. creation of SB
!! check if it is time for an SB creation
!! create SB by inserting passive particle
!! 2. check local list of SB particles and check if it's time to do something
!!
!! TODO method for removing SB	
!! TODO improve MPI
!! TODO ONLY WORKS FOR 1 SB per dt
!!
!! COMMENTS
!! assumes SB starts with SN explosion right away
!!
!! ARGUMENTS
!!
!!
!!***
!DEBUG_SN

subroutine SB (expFlag, numExp, dt, time, SBdt, nSB)

	use SB_data
    use mtmodSN

	use SN_data, Only: he_imin, he_imax, he_jmin, he_jmax, he_kmin, he_kmax, &
		 he_stratifySN, he_SNminstep, he_newDt, &
		 he_radialSN, he_erstar2, he_SNmapToGrid

  use Particles_data, ONLY: particles, pt_typeInfo, useParticles, pt_numLocal, pt_meshMe
  use Timers_interface, ONLY : Timers_start, Timers_stop
	use Grid_data, ONLY : gr_delta
	use tree, ONLY : lrefine

	implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Particles.h"
#include "Flash_mpi.h"

! arguments
	integer,intent(INout) :: numExp
	integer,intent(out)		:: nSB
	logical,intent(inout)	:: expFlag
	real,intent(IN) 			:: dt,time
	real,intent(inout) 		:: SBdt

	! SN module scratch variables, only local
	integer :: nsnb, nsnbdt, i, localnSB, j,error,nSBbuff,p, zoneindex, sbid
	logical :: new
	real    :: Tdt, tstarb, tnextb, xc, yc, zc, rc
	real    :: posx,posy,posz
	! spacing
  real, dimension(NDIM)	:: gridspace

! not included

  if (.not. useSB .or. .not. useParticles) return

!=====================================
! generate new SB
!=====================================
	new = .false.

  if(sb_useSBrandom) then
! for general purpose, look back one time step, extend time interval to check backwards for overlap
	  Tdt = time - dt

	  nsnb   = time/sb_tsb ! total number of SB since start of simulation
	  tstarb = nsnb*sb_tsb ! current explosion time for SB
	  nsnbdt = dt/sb_tsb

!	  call random_seed(GET=he_seed(1:he_seedsize))

! next SB creation time
	  tnextb = tstarb + sb_tsb

!  calculate new step to land before explosion add some slop
	  SBdt = 0.5*(tnextb - 0.9*he_SNminstep - time)

!  explosion is inside minimum timestep already
!	if(SBdt .le. 0d0) then
	  if(SBdt .le. he_SNminstep) then
		  SBdt = he_SNminstep
	  endif

! number of SB generated this timestep for output
	  nSB = 0

	  if (Tdt < tstarb) then
		  expFlag = .TRUE.

! create random location for SB
		  do i = numExp + 1, numExp + 1  + nsnbdt

			if(he_radialSN) then
!			  call random_number(xc)
!			  call random_number(yc)	
				xc = grndSN()
				yc = grndSN()
			  	rc = sqrt(xc**2d0 + yc**2d0)
			  	xc = xc/rc
			  	yc = yc/rc
			  	rc = he_erstar2*log(rc) ! use inverse transform, could also be 1-r0h, just to exclude 0
			  	xc = xc*rc
			  	yc = yc*rc
			else
				
! TODO - JCIM: Create the SB in density peaks.
! else if (he_SBpeaks) then
! ...
! else
				
!  determine x location of random SN
!			  call random_number(xc)
				xc = grndSN()
			  	xc = xc*(he_imax-he_imin)
			  	xc = xc + he_imin
!  determine y location of random SN
!			  call random_number(yc)
				yc = grndSN()
			  	yc = yc*(he_jmax-he_jmin)
			  	yc = yc + he_jmin
			endif

!		  call random_number(zc)
				zc = grndSN()
!  MKRJ - exponentially decling supernova rate in the vertical direction
!        the scale height, hstar1, is an input parameter
		  	if(he_stratifySN) then
			  	zc = zc*2.0 - 1.0
			  	zc = ( abs(zc)/zc ) * log(abs(zc))
			  	zc = zc*sb_hstarb
		  	else
				zc = zc*(he_kmax-he_kmin)
				zc = zc+he_kmin
		  	endif

! create sb passive particle, also checks if SB is local
		  	call sb_createSB(xc,yc,zc,time,dt,new)

! position has to be added on every process for SN execution to work
! update total number SB
! runs in lockstep, so update is done on all procs
! only local process actually generates SB in sb_createSB
			sb_nSN = sb_nSN + 1
		enddo
	endif

!=====================================
! check local existing SB for any SN
!=====================================
! loop over all sb and determine the particle IDs
! absolute reference
#ifdef SB_PART_TYPE

! reset buffer
	sbPOSbuff(1:4,1:sb_nSN) = 0d0

	localnSB = pt_typeInfo(PART_TYPE_BEGIN,SB_PART_TYPE) + pt_typeInfo(PART_LOCAL,SB_PART_TYPE) - 1
	do i = pt_typeInfo(PART_TYPE_BEGIN,SB_PART_TYPE), localnSB

	  p = i 

    if(new .and. (i .eq. localnSB)) then
! particles are not one chunk in memory, last one is appended 
	    p  = pt_numLocal
	  endif 

!	  if(particles(TYPE_PART_PROP,p) .ne. SB_PART_TYPE)  then
!	    !print*, 'not a SB particle, abort',pt_meshMe
!	    print*,particles(TYPE_PART_PROP,p) 
!	    print*,particles(TYPE_PART_PROP,1:20) 
!	    print*,particles(BLK_PART_PROP,1:20) 
!	    print*,pt_typeInfo(PART_TYPE_BEGIN,:)
!	    print*,pt_typeInfo(PART_LOCAL,:)
!	    stop
!	  endif
! total number of SB since creation of SB
! time or Tdt
! probably should be ceil(), not sure though
		nsnb   = (time - particles(TCRT_PART_PROP, p)) / particles(TSN_PART_PROP, p)

! SB is already done
		if(nsnb .gt. particles(NSN_PART_PROP, p) - 1 ) then 
! keep all SB for now
! otherwise remove here
      cycle
    endif

! next SN time for new timestep test
		tstarb = nsnb*particles(TSN_PART_PROP, p) + particles(TCRT_PART_PROP, p)
		tnextb = tstarb + particles(TSN_PART_PROP, p)

! calculate new step to land before explosion add some slop
		tnextb = 0.5*(tnextb - 0.9*he_SNminstep - time)

! compare to SB creation derived timestep
		if(tnextb .le. SBdt) then
			SBdt = tnextb
		endif

		if(SBdt .lt. he_SNminstep) then
			SBdt = he_SNminstep
		endif

! test for multiple SN in one timestep, actually should never happen if adaptive time step is on
		nsnbdt = dt/particles(TSN_PART_PROP, p)

		if ( Tdt < tstarb ) then
			expFlag = .TRUE.

! add to global SB explosion array
! TODO add one more slots for number of expl. from same SB, currently only one SN per SB processed
			if(nsnbdt .gt. 0) then
			  print*,'more than one SN in same SB',nsnbdt, dt, particles(TSN_PART_PROP,p),p,localnSB,new
			  print*,'particle props',particles(POSX_PART_PROP:POSZ_PART_PROP,p)
			  print*,p,pt_typeInfo(PART_TYPE_BEGIN,SB_PART_TYPE), localnSB
			  !print*,p,pt_typeInfo(PART_TYPE_BEGIN,PASSIVE_PART_TYPE),pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE)
			  print*,'particle type IDs'
			  print*,particles(TYPE_PART_PROP,1:pt_numLocal)
			  stop 
			endif

			if(he_SNmapToGrid) then
! guaranteed to be local
!				call gr_findBlock(blockList,blockCount,(/x0heat, y0heat, z0heat/),blockID)
! blockID should be particle property and up to date

				gridspace = gr_delta(1:MDIM,lrefine(int(particles(BLK_PART_PROP,p))))

! blockbounds not needed
				zoneindex = floor(particles(POSX_PART_PROP, p) / gridspace(1))
				posx = gridspace(1)*(zoneindex + 0.5d0)

			  	zoneindex = floor(particles(POSY_PART_PROP, p) / gridspace(2))
			  	posy = gridspace(2)*(zoneindex + 0.5d0)

			  	zoneindex = floor(particles(POSZ_PART_PROP, p) / gridspace(3))
			  	posz = gridspace(3)*(zoneindex + 0.5d0)

! write into slot for SB, all slots are unique for the specific SB tag, for low number of SB ok
! i.e. less than thousands, otherwise SB cleanup would help, or more clever MPI
				sbid = int(particles(SBID_PART_PROP, p))
				sbPOSbuff(1,sbid) = posx
				sbPOSbuff(2,sbid) = posy
				sbPOSbuff(3,sbid) = posz
				sbPOSbuff(4,sbid) = particles(SBID_PART_PROP, p)
			else
				sbid = int(particles(SBID_PART_PROP, p))
! write into slot for SB, all slots are unique for the specific SB tag, for low number of SB ok
! i.e. less than thousands, otherwise SB cleanup would help, or more clever MPI
				sbPOSbuff(1,sbid) = particles(POSX_PART_PROP, p)
				sbPOSbuff(2,sbid) = particles(POSY_PART_PROP, p)
				sbPOSbuff(3,sbid) = particles(POSZ_PART_PROP, p)
				sbPOSbuff(4,sbid) = particles(SBID_PART_PROP, p)
			endif

			numExp = numExp + nsnbdt + 1

! new SB
			if (nsnb .lt. 1) then 
			  call WriteSBFeedback(int(particles(SBID_PART_PROP, p)), nsnb+1, time, particles(POSX_PART_PROP, p),particles(POSY_PART_PROP, p), particles(POSZ_PART_PROP, p), &
			  particles(VELX_PART_PROP, p), particles(VELY_PART_PROP, p), particles(VELZ_PART_PROP, p))
			else
! subsequent SN
			  call WriteSBFeedback(int(particles(SBID_PART_PROP, p)), nsnb+1, time, particles(POSX_PART_PROP, p),particles(POSY_PART_PROP, p), particles(POSZ_PART_PROP, p), &
				particles(VELX_PART_PROP, p), particles(VELY_PART_PROP, p), particles(VELZ_PART_PROP, p))
			endif

! number of SN from SB this timestep
			nSB = nSB + 1 
		endif
	enddo
#endif
endif


! MPI step to find all SB on all processes for SN execution
! unfortunately done with allreduce which is bad
	call Timers_start("SB_MPI")

! tell all process that there is a SB SN
! only non zero for the process that the SN goes off on
	nSBbuff = nSB

! nSBbuffers are all added up and saved to nSB
	call MPI_allreduce (nSBbuff, nSB, 1, MPI_INTEGER, MPI_Sum, &
				 MPI_Comm_World, error)

!	sb_nSN = sb_nSN + nSB
!	nsb is 1 everywhere now

! use reduce/gather to send data to main process and bcast to send it back to everyone
! use MPI_reduce on the the 4 1d arrays
! and broadcast them to all afterwards
	call MPI_allReduce (sbPOSbuff(1:4,1:sb_nSN), sbPOS(1:4,1:sb_nSN), 4*sb_nSN, MPI_Double_Precision, MPI_Sum, &
				MPI_Comm_World, error)

! reset buffer
	sbPOSbuff(1:4,1:sb_nSN) = 0d0

! update number of Explosion index 
	numExp = nSB
	if(numExp .gt. 0) then
		expFlag = .true. 
	endif

! could use scatter/gather instead 
	!call MPI_allReduce (sbPOSbuff(1:3,1:sb_nSN+1), sbPOS(1:3,1:sb_nSN+1), 3*sb_nSN, MPI_Double_Precision, MPI_Sum)

  call MPI_Barrier (MPI_Comm_World, error)
  call Timers_stop("SB_MPI")

	return

contains

  subroutine WriteSBFeedback(SBid, nSN, time, x, y, z, vx, vy, vz)

    use SN_data, ONLY :  he_meshMe
!  use Driver_data, ONLY : dr_globalMe

    implicit none

#include "constants.h"
#include "Flash.h"

    real, intent(IN)   :: time, x, y, z, vx, vy, vz
    integer, intent(IN):: nSN, SBid ! 1 is I 2 is II and 3 is SN in SB

    integer, parameter :: funit_evol = 15
    character(len=80)  :: outfile = "SBfeedback.dat"
    integer            :: i


!   as all processors calculate the SN in lockstep, only the Master process has to output data.

    open(funit_evol, file=trim(outfile), position='APPEND')

    write(funit_evol,'(2(1X,I16),7(1X,ES16.9))') &
		SBid,&
	     nSN, 	&
	     time,	&
	     x,  		&
	     y,			&
	     z,			&
	     vx,		&
	     vy,		&
			 vz

    close(funit_evol)

    return
  end subroutine WriteSBFeedback
end subroutine SB
