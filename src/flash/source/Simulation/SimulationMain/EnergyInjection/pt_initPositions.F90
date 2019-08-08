subroutine pt_initPositions (blockID, success)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords, Grid_getPointData
  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc,pt_meshMe
! only for ion time calculation
  use Grid_data, ONLY : gr_delta
  use tree, ONLY : lrefine
  use pt_sourceUtil

  use rt_data, ONLY : eth0, ah0, ev2erg

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, INTENT(in) :: blockID
  logical, intent(out) :: success
  integer :: i, p
  logical :: IsInBlock
  real :: xpos, ypos, zpos, bxLower, byLower, bzLower, bxUpper, byUpper, bzUpper
  real :: xvel, yvel, zvel, blockSize(3), blockCenter(3)

! for source turn on time 
  real :: dens, tmp, fsh, Teff, Nion, Eion, phEnergy
  real :: sigh, ionN, ionE
  real, dimension(NDIM)	:: zonesize
  integer :: j,k,l


! Looking a bit stubby, cause I want sinks instead. - JW

!-------------------------------------------------------------------------------

! Particle slot number (incremented and saved between calls)

!  p = pt_numLocal

!!-------------------------------------------------------------------------------

!! Get locations of block faces.

!  call Grid_getBlkPhysicalSize(blockID, blockSize)
!  call Grid_getBlkCenterCoords(blockID, blockCenter)

!  bxLower    = blockCenter(1) - 0.5*blockSize(1)
!  bxUpper    = blockCenter(1) + 0.5*blockSize(1)
!  if (NDIM >= 2) then
!     byLower = blockCenter(2) - 0.5*blockSize(2)
!     byUpper = blockCenter(2) + 0.5*blockSize(2)
!  endif
!  if (NDIM == 3) then
!     bzLower = blockCenter(3) - 0.5*blockSize(3)
!     bzUpper = blockCenter(3) + 0.5*blockSize(3)
!  endif

!! Loop over both particles and compute their positions.
!  xpos = sim_p1x
!  ypos = sim_p1y
!  zpos = sim_p1z

!  do i = 1, sim_nPtot       ! 1 or 2 particles in orbit

!! Check to see if the particle lies within this block.
!     IsInBlock = (xpos >= bxLower) .and. (xpos < bxUpper)
!     if (NDIM >= 2) &
!          IsInBlock = IsInBlock .and. ((ypos >= byLower) .and. (ypos < byUpper))
!     if (NDIM == 3) &
!          IsInBlock = IsInBlock .and. ((zpos >= bzLower) .and. (zpos < bzUpper))

!! If yes, and adequate particle buffer space is available, initialize it.

!     if (IsInBlock) then
!        p = p + 1
!        if (p > pt_maxPerProc) then
!           call Driver_abortFlash &
!                ("InitParticlePositions:  Exceeded max # of particles/processor!")
!        endif

!        print*,'##found block for sources'
       
!        Eion = sim_Eph*ev2erg
!        sigh = ah0
!        Nion = sim_Nph

!        print*,'##total number of ionising photons [/cm/s]', Nion
!        print*,'##heating energy per ionizing photon [eV]', Eion/1.60217657e-12
!        print*,'##lyman limit cross section', ah0
     
!        particles(NION_PART_PROP,p) = Nion
!        particles(EION_PART_PROP,p) = Eion
!        particles(SIGH_PART_PROP,p) = sigh
  
!        !sim_EionPhot = Eion
!        !sim_sigmaH   = sigh
  
!! Particle current block number.
!        particles(BLK_PART_PROP,p) = real(blockID)
!        particles(PROC_PART_PROP,p) = real(pt_meshMe)

!! Particle mass.
!#ifdef MASS_PART_PROP
!        if (MASS_PART_PROP > 0) particles(MASS_PART_PROP,p) = 1.
!#endif 
!! Particle position and velocity.
!        particles(POSX_PART_PROP,p) = xpos
!        particles(POSY_PART_PROP,p) = ypos
!        particles(POSZ_PART_PROP,p) = zpos

!        particles(VELX_PART_PROP,p) = 0d0
!        particles(VELY_PART_PROP,p) = 0d0
!        particles(VELZ_PART_PROP,p) = 0d0
!     endif 
!  enddo

!! Set the particle database local number of particles.

!  pt_numLocal = p

  success = .true.
  return
end subroutine pt_initPositions
