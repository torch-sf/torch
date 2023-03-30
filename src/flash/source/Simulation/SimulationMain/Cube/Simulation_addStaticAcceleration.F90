!!****if* source/Simulation/SimulationMain/Cube/Simulation_addStaticAcceleration
!!
!! NAME
!!
!!   Simulation_addStaticAcceleration
!!
!! SYNOPSIS
!!
!!   Simulation_addStaticAcceleration(pos sweepDir blockID numCells grav)
!!
!! DESCRIPTION
!!
!!   Modifies the grav array by adding other gravitational acceleration
!!   distributions to the existing array. Currently only supports
!!   Z-directional accelerations but can be expanded to add other
!!   components for the X and Y sweeps.
!!
!! ARGUMENTS
!!
!!   pos      - Row indices transverse to the sweep direction
!!   sweepDir - The sweep direction:  allowed values are
!!              SWEEP_X, SWEEP_Y, SWEEP_Z which are defined
!!              in constants.h
!!   blockID  - block ID number
!!   numCells - Number of cells to update in grav array
!!   grav     - array of gravitational acceleration data for each cell.
!!
!! NOTES
!!
!!   o This routinue was originally designed for the StratBox Torch problem by
!!     Juan Ibanez-Mejia. It was subsequently refactored by Aaron Tran in 2019,
!!     and then lightly modified and documented by Sean C. Lewis in 2023 for
!!     use in the VorAMR Torch extension.
!!
!!*** 

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
	      ! gravitational acceleration due to stellar mass distribution (Hill et al. 2012)
          ! grav(ii) = grav(ii) - (sim_aParm1*zHeight/dsqrt(zHeight*zHeight + sim_aParm3*sim_aParm3)) - &
          !                       & sim_aParm2 * zHeight + sim_aParm4 * zHeight * abszHeight !- &
                                 !& sim_aParm5 * zHeight*zHeight*zHeight
          ! polynomial fit for background potential from VorAMR input data - SCL 
          grav(ii) = grav(ii) + sim_aParm1*zHeight**3 + sim_aParm2*zHeight**2 + sim_aParm3*zHeight + sim_aParm4
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
