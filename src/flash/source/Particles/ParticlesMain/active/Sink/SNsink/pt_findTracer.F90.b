!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkCreateParticle
!!
!! NAME
!!
!!  pt_sinkCreateParticle
!!
!! SYNOPSIS
!!
!!  call pt_sinkCreateParticle(real(in) :: x,y,z,pt
!!                             real(in) :: block_no, MyPE)
!!
!! DESCRIPTION
!!
!! accretes tracers that are in the sink particle block
!! 
!!
!!
!! NOTES
!!
!!
!!***
!TODO check if there has to be an extra sorting step after this

#define DEBUG_SN
subroutine pt_findTracer(block_no, ix, iy, iz, dx, blx, bly, blz, TrAcc, ratio)

  use Particles_data, only : particles, pt_TypeInfo, pt_numLocal
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  real, intent(IN)    :: block_no, ix, iy, iz, dx
! lower blockbounds
  real, intent(IN)    :: blx, bly, blz
! ratio is the quotient of threshold mass to mass in the zone
  real, intent(IN)    :: ratio
  real, intent(INOUT) :: TrAcc
  integer :: i, localnTracer, nAcc
  real :: ox,oy,oz 
! buffer for particle IDs that might be removed
	integer, dimension(NXB*NYB*NZB) :: idBuff	
	logical :: foundBlock = .false.

  localnTracer = pt_typeInfo(PART_TYPE_BEGIN,PASSIVE_PART_TYPE) + pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) - 1

! loop over all local tracers and see if they are contained in zone marked for accretion
! should be sorted so just look the ones containd in the block, but that requires some coding effort

	nAcc = 0
  do i = pt_typeInfo(PART_TYPE_BEGIN,PASSIVE_PART_TYPE), localnTracer
! if sorted correctly then there will be no particles with the required block number after all have been found
		if(foundBlock .and. block_no .ne.  particles(BLK_PART_PROP,i) ) then
			exit	
		endif

    if(block_no .eq.  particles(BLK_PART_PROP,i)) then
		  foundBlock = .true.

! corner of the block is origin of block coordinate system
      ox  = particles(POSX_PART_PROP,i) - blx
      oy  = particles(POSY_PART_PROP,i) - bly
      oz  = particles(POSZ_PART_PROP,i) - blz

! transform to local coordinate system, not sure about this one
      ox  = floor(ox/dx)
      oy  = floor(oy/dx)
      oz  = floor(oz/dx)

#ifdef DEBUG_SN
		  print*,'accreting from:',ix,iy,iz, ox,oy,oz
#endif

      if(ix .eq. ox .and. iy .eq. oy .and. iz .eq. oz) then
			  nAcc   = nAcc + 1
			  idBuff(nAcc) = i
! tracer should be accreted and removed from particle map
! number of accreted tracer particles or something else?
!        TrAcc = TrAcc + 1d0
! good bye tracer
!        pt_numLocal = pt_numLocal - 1
!        pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) - 1   
!        particles(BLK_PART_PROP,i) = -1 
!        particles(TAG_PART_PROP,i) = -1
      endif
   endif
  enddo

! remove particles according to the passed fraction that should be removed 
#ifdef DEBUG_SN
			if(nAcc .gt. 0) then
			  print*,'number of possible tracers to accrete, ratio to accrete',nAcc, ratio
			endif
#endif

	nAcc = ceiling(nAcc*ratio)

#ifdef DEBUG_SN
		if(nAcc .gt. 0) then
		  print*,'number of tracers to accrete',nAcc
		endif
#endif

	do i = 1, nAcc
		TrAcc = TrAcc + 1d0
! good bye tracer
		pt_numLocal = pt_numLocal - 1
		pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) - 1   
		particles(BLK_PART_PROP,i) = -1 
		particles(TAG_PART_PROP,i) = -1
	enddo

  return

end subroutine pt_findTracer
