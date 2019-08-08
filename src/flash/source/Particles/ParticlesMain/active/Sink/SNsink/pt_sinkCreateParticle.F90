!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkCreateParticle
!!
!! NAME
!!
!!  pt_sinkCreateParticle
!!
!! SYNOPSIS
!!
!!  pno = pt_sinkCreateParticle(real   (in) :: x,
!!                              real   (in) :: y,
!!                              real   (in) :: z,
!!                              real   (in) :: pt,
!!                              integer(in) :: block_no,
!!                              integer(in) :: MyPE)
!!
!! DESCRIPTION
!!
!!   Create new sink particle on local CPU at position x, y, z, time pt,
!!   associated with block_no and MyPE.
!!
!! ARGUMENTS
!!
!!   x - x position of new particle
!!
!!   y - y position of new particle
!!
!!   z - z position of new particle
!!
!!   pt - pt creation time of new particle
!!
!!   block_no - block in which new particle is created
!!
!!   MyPE - Processor ID on which new particle is created
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!
!!	 CB added SNsink Module to FLASH 3/4
!!	 JC ported and modified the SNsink Module to FLASH4.2.2 2015
!!
!!***

function pt_sinkCreateParticle(x, y, z, pt, block_no, MyPE)

  use Particles_sinkData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
!	use SB_Data

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  real, intent(IN)    :: x, y, z, pt
  integer, intent(IN) :: block_no, MyPE

  integer             :: pt_sinkCreateParticle, pno, ire, tag
  logical, parameter  :: debug = .false.


  pt_sinkCreateParticle = 0

  localnp = localnp + 1
  if (localnp .GT. sink_maxSinks) &
       call Driver_abortFlash('pt_sinkCreateParticle: Sink particle number exceeds sink_maxSinks. Increase.')

  pno = localnp

  n_empty = n_empty - 1

  do ire = 1, pt_sinkParticleProps
     particles_local(ire, pno) = 0.
  enddo

  local_tag_number = local_tag_number + 1
  
  particles_local(iptag,pno) = get_tag(local_tag_number, MyPE)
  particles_local(ipcpu,pno) = MyPE
  particles_local(ipblk,pno) = block_no
  particles_local(ipx,pno) = x
  particles_local(ipy,pno) = y
  particles_local(ipz,pno) = z
  particles_local(ipt,pno) = pt
  particles_local(ipvx,pno) = 0.0
  particles_local(ipvy,pno) = 0.0
  particles_local(ipvz,pno) = 0.0
  particles_local(iplx,pno) = 0.0
  particles_local(iply,pno) = 0.0
  particles_local(iplz,pno) = 0.0
  particles_local(ipm,pno) = 0.0
  particles_local(ipmdot,pno) = 0.0
  particles_local(ACCX_PART_PROP,pno) = 0.0
  particles_local(ACCY_PART_PROP,pno) = 0.0
  particles_local(ACCZ_PART_PROP,pno) = 0.0
  particles_local(OACX_PART_PROP,pno) = 0.0
  particles_local(OACY_PART_PROP,pno) = 0.0
  particles_local(OACZ_PART_PROP,pno) = 0.0
  particles_local(ipdtold,pno) = 0.0


!////////SB and or subgridmodel
!/////////////////////////
! call magical_subgrid_model()

! Current:  mass of sink at creation divided by CloudMassPerStar in solar masses and returns number if SN (happens after accretion in Particle_sinkCreateAccrete.F90)
! TODO: Link the subgrid Model to the sink particle creation to give the number of SN.
  particles_local(NSN_PART_PROP,pno)  = 0d0

! creation time, already exists in sink particle reused for tracer accretion 
  particles_local(TCRT_PART_PROP,pno) = 0d0

! time between SN in SB
! Currently: sb_life/( particles_local(NSN_PART_PROP,pno)+1) ! +1 because all SN should go off and we need enough intervals before maximum time is reached
! TODO: Get the SN explosion times from the subgrid model. (either a av SN explosion times or store the actual SN times in an array somewhere.)
  particles_local(TSN_PART_PROP,pno)  = 0d0 

! type for bookeeping, Remember that I am a Sink.
!  particles_local(TYPE_PART_PROP,pno) =  SINK_PART_TYPE
!/////////////////

  tag = int(particles_local(iptag, pno))

  pt_sinkCreateParticle = localnp

  if (debug) then
     print*, "creating sink particle!"
     print*, "x pos=", x
     print*, "y pos=", y
     print*, "z pos=", z
     print*, "creation time=", pt
     print*, "cpu=", MyPE
     print*, "tag=",  int(particles_local(iptag, pno))
  end if

  return

end function pt_sinkCreateParticle
