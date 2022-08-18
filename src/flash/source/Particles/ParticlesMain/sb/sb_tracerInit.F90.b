!!
!! NAME
!!
!!   sb_createSBtracer
!!
!! SYNOPSIS
!!
!!  inserts a passive particles as SB and a sphere of tracer particles around it
!!
!! ARGUMENTS
!! x,y,z : position of SB
!! time  : current simulation time   
!! rad	 : shell radius
!! tag   : this is either the particle tag or the he_nSN number
!! SNtype: 1 for type I, 2 for type 2, 3 for SB SN, 4 for Sink SN 
!!
!!
!! DESCRIPTION
!! generates a passive particle field on a sphere on the SN shell
!! passive tracer particles are then advected with the underlying gas velocity field   
!!
!!***

subroutine sb_tracerInit(nExp, cleanUp)

	use Particles_Data, only: pt_maxPerProc, pt_numLocal, particles, pt_meshMe, pt_typeInfo
	use SB_data

	implicit none

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

  integer,  intent(in)	:: nExp
  logical,  intent(in)	:: cleanUP
  integer 		:: startID, nP, p, maxP, offset
  integer		:: nTracer

print*,'cleanup?',cleanUp
  p       = pt_numLocal
  maxP    = nExp*sb_TracerPerSN
print*,'maxp',maxp
  Ntracer = 0
  nP      = pt_typeInfo(PART_LOCAL, sb_TracerType)
print*,'nP',nP
print*,'numlocal',pt_numLocal
  
! check for space
  if (pt_numLocal+maxP > pt_maxPerProc) then
!                           call abort_flash &
    print*,"  sb_createSB:  Exceeded max # of particles/processor! :(" 
    stop
  endif
! move stuff around to free slots
  if(.not. cleanUp) then 
    if(nP .gt. 0) then

      startID = pt_typeInfo(PART_TYPE_BEGIN, sb_TracerType)
! shuffle memory

! scratch for new SB
! move all particles n slots down (untested!)
      particles(:,(startID+nP+maxP):(pt_numlocal+maxP)) = particles(:,startID+nP:pt_numlocal)
! save index of last freed slot 
      tracerMemID = startID+nP+maxP-1
! update tracer particle data structure in sb_createSBtracer.F90
!    pt_typeInfo(PART_LOCAL, sb_TracerType) = pt_typeInfo(PART_LOCAL, sb_TracerType) + sb_NTracers
!  else
	
!      pt_typeInfo(PART_LOCAL, sb_TracerType) = pt_typeInfo(PART_LOCAL, sb_TracerType) + sb_NTracers
! is just last occupied slot
!      pt_typeInfo(PART_TYPE_BEGIN,sb_TracerType)  = pt_numLocal+1
    else
! prepare for appending new tracers
      pt_typeInfo(PART_TYPE_BEGIN, sb_TracerType) = pt_numlocal+1
    endif
  endif

! close gap caused by unused tracers
  if(cleanUP) then
    startID = pt_typeInfo(PART_TYPE_BEGIN, sb_TracerType)
    print*,'type begin',startID
    print*,'nP,tracerMemid',np,tracerMemID
! number of tracer particles, nP changed between calls
    nTracer = tracerMemID - (startID + nP)

! if negative tracer were all appended  
    if(ntracer .lt. 0) then 
      ntracer = startID+nP-1  
    else 

      if(nTracer .lt. maxP) then
	print*,'closing gap', ntracer,maxP
! last used position by tracer particles
        offset = startID + nP
! check beforehand if start and end slot are part of the gap
        if(particles(TAG_PART_PROP,tracerMemID) .gt. 0 .or. particles(TAG_PART_PROP,tracerMemID-nTracer) .gt. 0) then
	  print*,'gap not actually gap'
	  stop
        endif 
! move particles back  
        particles(:,tracerMemID-nTracer:pt_numLocal-nTracer) = particles(:,tracerMemID+1:pt_numlocal)
      endif
    endif
  endif

  return
end subroutine sb_tracerInit
