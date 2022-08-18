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
!! create SB bz inserting passive particle
!! 2. check local list of SB particles and check if it's time to do something
!!
!! TODO method for removing SB
!!
!! COMMENTS
!! assumes SB starts with SN explosion right away
!!
!! ARGUMENTS
!!
!!
!!***

subroutine SB (expFlag, numExp, x, y, z, dt, time, SBdt, nSB)

	use SB_data
	use SN_data, Only: he_imin, he_imax, he_jmin, he_jmax, he_kmin, he_kmax, &
		 he_stratifySN, he_SNminstep, he_newDt, he_nsndt, &
		 he_seed,he_seedsize

  use Particles_data, ONLY: particles, pt_typeInfo, useParticles, pt_numLocal, pt_meshMe
  use Timers_interface, ONLY : Timers_start, Timers_stop


!	use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
!													   Grid_releaseBlkPtr, Grid_getDeltas

!	use Timers_interface, ONLY : Timers_start, Timers_stop
!	use Eos_interface, ONLY : Eos_wrapped, Eos_getAbarZbar

	implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Particles.h"
#include "Flash_mpi.h"

! arguments
	real, DIMENSION(he_nsndt),intent(inout)	:: x, y, z
	integer,intent(INout)			:: numExp
	integer,intent(out)			:: nSB
	logical,intent(inout)			:: expFlag
	real,intent(IN) 			:: dt,time
	real,intent(inout) 			:: SBdt

	! SN module scratch variables, only local
	integer :: nsnb, nsnbdt, i, localnSB, j, pID,error,nSBbuff
	real    :: Tdt, tstarb, tnextb, xc, yc, zc

  if (.not. useSB .or. .not. useParticles) return

!=====================================
! generate new SB
!=====================================

! for general purpose, look back one time step, extend time interval to check backwards for overlap
	Tdt = time - dt

! for field supernovae of Type I
	nsnb   = time/sb_tsb ! total number of SB since start of simulation
	tstarb = nsnb*sb_tsb ! current explosion time for SB
	nsnbdt = dt/sb_tsb

	call random_seed(GET=he_seed(1:he_seedsize))

! next SB creation time
	tnextb = tstarb + sb_tsb

!  calculate new step to land before explosion add some slop
	SBdt = 0.5*(tnextb - 0.9*he_SNminstep - time)

!  explosion is inside minimum timestep already
	if(SBdt .le. 0d0) then
		SBdt = he_SNminstep
	endif

! number of SB generated this timestep for output
	nSB = 0
	!print*,'SB timestep', SBdt

	if (Tdt < tstarb) then
		expFlag = .TRUE.

! create random location for SB
		do i = numExp + 1, numExp + 1  + nsnbdt
!  determine x location of random SN
			call random_number(xc)
			xc = xc*(he_imax-he_imin)
			xc = xc + he_imin
!  determine y location of random SN
			call random_number(yc)
			yc = yc*(he_jmax-he_jmin)
			yc = yc + he_jmin

			call random_number(zc)
!  MKRJ - exponentially decling supernova rate in the vertical direction
!        the scale height, hstar1, is an input parameter
			if(he_stratifySN) then
				zc = zc*2.0 - 1.0
				zc = -( abs(zc)/zc ) * log(abs(zc))
				zc = zc*sb_hstarb
			else
				zc = zc*(he_kmax-he_kmin)
				zc = zc+he_kmin
			endif

! create sb passive particle
			call sb_createSB(xc,yc,zc,time,pID)

! position has to be added on every process for SN execution to work
! update total number SB
			sb_nSN = sb_nSN + 1
!			!he_nSN = he_nSN + 1
!			nSB = nSB + 1 
		enddo

! assumes SB starts with SN explosion right away
	!	numExp = numExp + nsnbdt + 1
	endif

!=====================================
! check local exisiting SB for any SN
!=====================================
! loop over all sb and determine the particle IDs
! absolute reference
	localnSB = pt_typeInfo(PART_TYPE_BEGIN,sb_typeID) + pt_typeInfo(PART_LOCAL,sb_typeID) - 1
	do i = pt_typeInfo(PART_TYPE_BEGIN,sb_typeID), localnSB

! for field supernovae of Type I
! total number of SB since creation of SB
! time or Tdt
		nsnb   = (time - particles(TCRT_PART_PROP, i)) / particles(TSN_PART_PROP, i)

! SB is already done
		if(nsnb .gt. particles(NSN_PART_PROP, i) - 1 ) cycle

! maybe remove SB here

! next SN time for new timestep test
		tstarb = nsnb*particles(TSN_PART_PROP, i) + particles(TCRT_PART_PROP, i)
		tnextb = tstarb + particles(TSN_PART_PROP, i)
		!	print*,'SN in SB', tstarb, nsnb
! calculate new step to land before explosion add some slop
		tnextb = 0.5*(tnextb - 0.9*he_SNminstep - time)

! compare to SB creation derived timestep
		if(tnextb .le. SBdt) then
			SBdt = tnextb
		endif

		if(SBdt .le. 0d0) then
			SBdt = he_SNminstep
		endif

! test for multiple SN in one timestep, actually should never happen if adaptive time step is on
		nsnbdt = dt/particles(TSN_PART_PROP, i)
		!print*, Tdt, tstarb,nsnb, particles(TAG_PART_PROP, i),particles(TCRT_PART_PROP, i),particles(TSN_PART_PROP, i)
		!print*, tstarb + particles(TCRT_PART_PROP, i), Tdt - (tstarb + particles(TCRT_PART_PROP, i))

		if ( Tdt < tstarb ) then
			expFlag = .TRUE.
! add to global SB  explosion array
! TODO add one more slot for number of expl. from same SB, currently only on SN per SB processed
			if(nsnbdt .gt. 0) then
			 print*, 'more than one SN in same SB'
			 stop 
			endif
			do j = numExp + 1, numExp + 1 + nsnbdt
				sbPOSbuff(1,particles(TAG_PART_PROP, i)) = particles(POSX_PART_PROP, i)
				sbPOSbuff(2,particles(TAG_PART_PROP, i)) = particles(POSY_PART_PROP, i)
				sbPOSbuff(3,particles(TAG_PART_PROP, i)) = particles(POSZ_PART_PROP, i)
			enddo

			numExp = numExp + nsnbdt + 1

! subsequent SN
			call WriteSBFeedback(int(particles(TAG_PART_PROP, i)), nsnb, time, particles(POSX_PART_PROP, i),particles(POSY_PART_PROP, i), particles(POSZ_PART_PROP, i), &
				 particles(VELX_PART_PROP, i), particles(VELY_PART_PROP, i), particles(VELZ_PART_PROP, i))

! update total number SN
			!he_nSN = he_nSN + 1

! number of SN from SB this timestep
			nSB = nSB + 1 
		endif
	enddo

! MPI step to find all SB on all processes for SN execution
! unfortunately done with allreduce which is bad
	call Timers_start("SB_MPI")
	call MPI_allReduce (sbPOSbuff(1:3,1:sb_nSN), sbPOS(1:3,1:sb_nSN), 3*sb_nSN, MPI_Double_Precision, MPI_Sum, &
				MPI_Comm_World, error)
! reset buffer
	sbPOSbuff(1:3,1:sb_nSN) = 0d0 
! tell all process that there is a SB SN
	nSBbuff = nSB
	call MPI_allreduce (nSBbuff, nSB, 1, MPI_INTEGER, MPI_Sum, &
				 MPI_Comm_World, error)

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
