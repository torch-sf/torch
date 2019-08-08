!!
!! NAME
!!
!!   sb_createSNtracer
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

subroutine sb_createSNtracer(x, y, z, vx, vy, vz, temp, dens, time, tag, SNtype, blockID)!,rad)

	use Particles_Data, only: pt_maxPerProc, pt_numLocal, particles, pt_meshMe, pt_typeInfo
	use SB_data

	implicit none

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

  real,  intent(in)		:: x,y,z,time,vx,vy,vz,temp,dens!,rad
  integer,  intent(in)		:: tag, SNtype,blockID
  integer :: startID, nP, i, p

! check for space
  if(pt_typeInfo(PART_LOCAL, PASSIVE_PART_TYPE) .gt. 0) then

    startID = pt_numLocal+1
    do i = 0, sb_tracerPerZone-1
! TODO offset here 
      p = startID + i 
!   same block as source initially
      particles(BLK_PART_PROP,p)  = blockID
      particles(PROC_PART_PROP,p) = real(pt_meshMe)
! SB or sink tag, otherwise 0
      particles(TAG_PART_PROP,p)  = tag
      particles(TYPE_PART_PROP,p) = PASSIVE_PART_TYPE

      particles(POSX_PART_PROP,p) = x
      particles(POSY_PART_PROP,p) = y 
      particles(POSZ_PART_PROP,p) = z 

! set from zone 
      particles(VELX_PART_PROP,p) = vx
      particles(VELY_PART_PROP,p) = vy
      particles(VELZ_PART_PROP,p) = vz
! type of origin SN
      particles(TSN_PART_PROP,p)  = SNType
! explosion time 
      particles(TCRT_PART_PROP,p) = time
! free for whatever, temperature 
      particles(NSN_PART_PROP,p)  = temp 
! need one more field for density, could be mass if sinks are in etc. check with some preprocessor code maybe
      particles(MASS_PART_PROP,p) = dens
    enddo

!    pt_typeInfo(PART_LOCAL, sb_TracerType)      = pt_typeInfo(PART_LOCAL, sb_TracerType) + sb_tracerPerZone
!    pt_typeInfo(PART_TYPE_BEGIN, sb_TracerType) = pt_typeInfo(PA
     pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) + sb_tracerPerZone

! is just last occupied slot
     pt_numLocal = pt_numLocal + sb_tracerPerZone
! is just last occupied slot

  else

! just append 
    startID = pt_numLocal+1
    do i = 0, sb_tracerPerZone-1
! TODO offset here 
      p = startID + i 

!      if(particles(TAG_PART_PROP,p) .gt. 0) then
!        print*,'  Particle slot not actually free, but should be, in sb_createSNtracer, append'
!       stop
!     endif
!   same block as source initially
      particles(BLK_PART_PROP,p)  = blockID
      particles(PROC_PART_PROP,p) = real(pt_meshMe)
! SB or sink tag, otherwise 0
      particles(TAG_PART_PROP,p)  = tag
      particles(TYPE_PART_PROP,p) = PASSIVE_PART_TYPE

      particles(POSX_PART_PROP,p) = x
      particles(POSY_PART_PROP,p) = y 
      particles(POSZ_PART_PROP,p) = z 

! set from zone 
      particles(VELX_PART_PROP,p) = vx
      particles(VELY_PART_PROP,p) = vy
      particles(VELZ_PART_PROP,p) = vz
! type of origin SN
      particles(TSN_PART_PROP,p)  = SNType
! explosion time 
      particles(TCRT_PART_PROP,p) = time
! free for whatever, temperature 
      particles(NSN_PART_PROP,p)  = temp 
! need one more field for density, could be mass if sinks are in etc. check with some preprocessor code maybe
      particles(MASS_PART_PROP,p) = dens

    enddo
     pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) = pt_typeInfo(PART_LOCAL,PASSIVE_PART_TYPE) + sb_tracerPerZone

! is just last occupied slot
     pt_numLocal = pt_numLocal + sb_tracerPerZone
! is just last occupied slot
     pt_typeInfo(PART_TYPE_BEGIN,PASSIVE_PART_TYPE) = startID
  endif

  return
end subroutine sb_createSNtracer
