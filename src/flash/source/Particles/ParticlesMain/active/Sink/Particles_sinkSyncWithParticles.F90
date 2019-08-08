!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkSyncWithParticles
!!
!! NAME
!!
!!  Particles_sinkSyncWithParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSyncWithParticles(logical(in) :: sink_to_part)
!!
!! DESCRIPTION
!!
!!  Synchronizes global particle array with sink particle array.
!!
!! ARGUMENTS
!!
!!   sink_to_part -  logical flag indicating whether to sync with global
!!                   particle array or not
!!
!! NOTES
!!
!!   written by John Bachan, 2012
!!   modified to include off-domain support, Christoph Federrath, 2015
!!
!!***

subroutine Particles_sinkSyncWithParticles(sink_to_part)

  use Particles_data, only: pt_numLocal, particles, pt_typeInfo, pt_maxPerProc
  use Particles_sinkData, only: localnp, particles_local, particles_global, &
    pt_sinkParticleProps, sink_maxSinks, sink_offDomainSupport
  use Grid_interface, only: Grid_sortParticles
  use Driver_interface, only: Driver_abortFlash

  implicit none

#include "Particles.h"
#include "Flash.h"

  logical, intent(in) :: sink_to_part
  integer :: i
  integer :: s0, sn
  integer :: particlesPerBlk(MAXBLOCKS,NPART_TYPES)
! I'm correcting this once again to also include BLK_PART_PROP in the
! first call. Note without this fix restarts that include both active
! and sink particles will ERASE THE SINK PARTICLES. - JW
#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP,TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
  ! Now update the pt_typeInfo data structure
  call pt_updateTypeDS(particlesPerBlk)
  
  s0 = pt_typeInfo(PART_TYPE_BEGIN,SINK_PART_TYPE)
  sn = pt_typeInfo(PART_LOCAL,SINK_PART_TYPE)
  
  if (sink_to_part) then

    ! erase existing sink particles from particles
    do i = s0, pt_numLocal-sn
      particles(:,i) = particles(:,i+sn)
    end do
    do i = 1, NPART_TYPES
      if(pt_typeInfo(PART_TYPE_BEGIN,i) > s0) then
        pt_typeInfo(PART_TYPE_BEGIN,i) = pt_typeInfo(PART_TYPE_BEGIN,i) - sn
      end if
    end do
    pt_typeInfo(PART_LOCAL,SINK_PART_TYPE) = 0
    pt_numLocal = pt_numLocal - sn
    ! append new sink particles
    pt_typeInfo(PART_TYPE_BEGIN,SINK_PART_TYPE) = pt_numLocal+1
    pt_typeInfo(PART_LOCAL,SINK_PART_TYPE) = localnp
    particles(:,pt_numLocal+1:pt_numLocal+localnp) = particles_local(:,1:localnp)
    pt_numLocal = pt_numLocal + localnp

  else ! sink_to_part = .false.

    if (.not. allocated(particles_local)) then
      allocate(particles_local(pt_sinkParticleProps, sink_maxSinks))
    end if
    if (.not. allocated(particles_global)) then
      allocate(particles_global(pt_sinkParticleProps, sink_maxSinks))
    end if
    
    localnp = sn
    particles_local(:,1:localnp) = particles(:,s0:s0+sn-1)
    
    ! if we run off-domain sinks, then erase existing sink particles from particles
    ! to allow them to remain active; they are only seen by the sink routines, not the
    ! FLASH particle structures anymore. We only pull sinks into the FLASH particles for I/O
    if (sink_offDomainSupport) then
      do i = s0, pt_numLocal-sn
        particles(:,i) = particles(:,i+sn)
      end do
      do i = 1, NPART_TYPES
        if(pt_typeInfo(PART_TYPE_BEGIN,i) > s0) then
          pt_typeInfo(PART_TYPE_BEGIN,i) = pt_typeInfo(PART_TYPE_BEGIN,i) - sn
        end if
      end do
      pt_typeInfo(PART_LOCAL,SINK_PART_TYPE) = 0
      pt_numLocal = pt_numLocal - sn
    end if

  end if

end subroutine Particles_sinkSyncWithParticles
