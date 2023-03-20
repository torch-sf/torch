!******************************************************************************
!
!  Routine:     AddDiskMassAcceleration()
!
!  Description: Adds acceleration due to a galactic potential disk+halo
!  is called from Gravity_accelOneRow
!
!  could use an interface that is then implemented by a specific galactic potential
!  would switching between different models easy, BUT lazynessssss s
!
!  this is a dummy implementation with nothing in it
!  the proper function should be in the Simulation folder, with the needed implementation

subroutine Simulation_addStaticAcceleration (pos, sweepDir, blockID, numCells, grav)
!==============================================================================

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords
!  use Gravity_data, only : onlyStaticGrav  ! onlyStaticGrav for BHTree only -AT 2019 May 09

  use Simulation_data
  
! this is simulation specific, finest imported constants
!  use Simulation_data, only: sim_ka1, sim_kb1, sim_ka2, sim_kb2, sim_Md1, &
!	 				    sim_Md2, sim_haloc, sim_r200, sim_M200, sim_fc

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)	              :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer :: j,k,ii

  real,allocatable,dimension(:)    :: xCenter, yCenter, zCenter
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  logical :: gcell = .true.
  real    :: zHeight, rhoNFW, abszHeight

  if (.not. sim_withStaticGrav) then
    return
  endif

  if(sweepdir .ne. SWEEP_Z) then
    return
  endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

!  sizeX=blkLimitsGC(HIGH,IAXIS)
!  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)

!  allocate(xCenter(sizeX))
!  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))

  call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)

! put here the needed parameters of the static potentials as fortran parameters,
! for read in from .par file use Simulation init/data formalism

    if(sweepdir .eq. SWEEP_Z) then
       do ii = 1, numCells
          zHeight = zCenter(ii)
          abszHeight = abs(zHeight)			
      ! TODO 2019 May 09 - AT need to reimplement sim_NFWh from Juan's stratbox
			!if (abszHeight  .lt. sim_NFWh) then
	      ! gravitational acceleration due to stellar mass distribution
          grav(ii) = grav(ii) - (sim_aParm1*zHeight/dsqrt(zHeight*zHeight + sim_aParm3*sim_aParm3)) - &
                                & sim_aParm2 * zHeight + sim_aParm4 * zHeight * abszHeight !- &
				!& sim_aParm5 * zHeight*zHeight*zHeight
			!else
			!	rhoNFW   = sim_rho_s / ( abszHeight / sim_rs * ( 1d0 + abszHeight / sim_rs )**2d0 )
			!	grav(ii) = grav(ii) -4./3d0 * sim_GravConst * PI * rhoNFW * zheight
			!endif
       enddo
    endif


!  deallocate(xCenter)
!  deallocate(yCenter)
  deallocate(zCenter)

  return

end subroutine Simulation_addStaticAcceleration
