!!****if* source/Simulation/SimulationMain/EnergyInjection/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId)
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the EnergyInjection test problem.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!!
!!***

subroutine Simulation_initBlock(blockID)
  ! get the needed unit scope data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits,Grid_getCellCoords, &
                              Grid_putPointData
  use Multispecies_interface, ONLY : Multispecies_getSum
  use Eos_interface, ONLY: Eos_wrapped
  use Simulation_data
#ifdef VARY_ATM_FRAC
  use rt_data, ONLY : rt_vary_atomic_frac
#endif
  implicit none

! get all the constants
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

! define arguments and indicate whether they are input or output
  integer, intent(in) :: blockID

! declare all local variables.
  integer :: i, j, k, n
  real :: xx, yy,  zz

  ! arrays to hold coordinate information for the block
  real,allocatable, dimension(:) ::xCoord,yCoord,zCoord

  ! array to get integer indices defining the beginning and the end
  ! of a block.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

! the number of grid points along each dimension
  integer :: sizeX, sizeY, sizeZ

  integer, dimension(MDIM) :: axis
  logical,parameter :: gcell = .true.

  ! these variables store the calculated initial values of physical
  ! variables a grid point at a time.
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, &
       enerZone, ekinZone, tempZone
  real :: sim_mu
  real, dimension(NSPECIES) :: massFrac_box

  massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)    = sim_init_Hp 
  massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)    = (1.0 - sim_init_Hp)

  call Multispecies_getSum(GAMMA, sim_gamma, massFrac_box)

  ! get the integer endpoints of the block in all dimensions
  ! the array blkLimits returns the interior end points
  ! whereas array blkLimitsGC returns endpoints including guardcells
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

! get the size along each dimension for allocation and then allocate
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCoord(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCoord, sizeX)
  call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCoord, sizeY)
  call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCoord, sizeZ)

  !-----------------------------------------------------------------------------
  ! loop over all of the zones in the current block and set the variables.
  !-----------------------------------------------------------------------------
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     zz = zCoord(k) ! coordinates of the cell center in the z-direction

     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        yy = yCoord(j) ! center coordinates in the y-direction

        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx = xCoord(i)

           axis(IAXIS) = i   ! Get the position of the cell in the block
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! Compute the gas energy and set the gamma-values
           ! needed for the equation of  state.
           ! boring temperature

           tempZone = sim_amTemp
           rhoZone  = sim_amNumDens*sim_protonmass
           
           sim_mu = 1.3d0
#ifdef VARY_ATM_FRAC
           if (rt_vary_atomic_frac) then
             if (tempZone > 8e3) then
               sim_mu = 0.61d0
             else if (tempZone > 100.) then
               sim_mu = 1.3d0
             else
               sim_mu = 2.3d0
             end if
           end if
#endif

           presZone  = rhoZone * sim_gasconstant * sim_amTemp / (sim_mu*sim_protonmass)

           ekinZone = 0.0

           velxZone = 0d0
           velyZone = 0d0
           velzZone = 0d0

!          Gamma equation of state.
           enerZone = presZone / (sim_gamma-1.)
           enerZone = enerZone / rhoZone
           enerZone = enerZone + ekinZone

           do n=1,NSPECIES
             call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN+n-1, EXTERIOR, axis, massFrac_box(n))
           enddo

           ! store the variables in the current zone via the Grid_putPointData method
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, enerZone)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tempZone)

           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, TDUS_VAR, EXTERIOR, axis, sim_tdust)
#ifdef VARY_ATM_FRAC
           call Grid_putPointData(blockId, CENTER, ATMU_VAR, EXTERIOR, axis, sim_mu)
#endif
#ifdef MAGX_VAR
           call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, sim_magx)
           call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, sim_magy)
           call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, sim_magz)
#endif
#if NFACE_VARS > 0
           if (sim_killdivb) then
              call Grid_putPointData(blockId, FACEX, MAGX_VAR, EXTERIOR, axis, sim_magx)
              call Grid_putPointData(blockId, FACEY, MAGY_VAR, EXTERIOR, axis, sim_magy)
              if (NDIM == 3) call Grid_putPointData(blockId, FACEZ, MAGZ_VAR, EXTERIOR, axis, sim_magz)
           endif
#endif
        end do
     end do
  end do

! set EOS
  call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
