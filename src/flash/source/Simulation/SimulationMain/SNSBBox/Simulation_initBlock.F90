!!****if* source/Simulation/SimulationMain/SNSBBox/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockId)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!
!! ARGUMENTS
!!
!!  blockId -          The number of the block to initialize
!!
!!
!!***

!!REORDER(4):solnData, face[xy]Data

subroutine Simulation_initBlock(blockId)

  use Simulation_data

  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr
#ifdef IHP_SPEC
  use Multispecies_interface, ONLY : Multispecies_getSum
#endif

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockId
  !!$ ---------------------------------

  integer :: i,j,k,n,istat,sizeX,sizeY,sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: ekinZone, eBZone, eintZone
  real, allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData
  logical :: getGuardCells = .true.

  real, dimension(NSPECIES) :: massFrac_box  ! Multispecies

  ! Overwrites user-input sim_gamma
#ifdef IHP_SPEC
  massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)    = sim_init_Hp
  massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)    = (1.0 - sim_init_Hp)
  call Multispecies_getSum(GAMMA, sim_gamma, massFrac_box)
#endif

  ! compute the maximum length of a vector in each coordinate direction
  ! (including guardcells)
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  ! necessary for NDIM < 3
  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords &
                      (KAXIS, blockId, CENTER,getGuardCells, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords &
                      (JAXIS, blockId, CENTER,getGuardCells, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, getGuardCells, xCoord, sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

          ! temperature is set by Eos_wrapped(...)
          solnData(DENS_VAR,i,j,k)= sim_rho
          solnData(PRES_VAR,i,j,k)= sim_pres

          solnData(VELX_VAR,i,j,k)= 0.0
          solnData(VELY_VAR,i,j,k)= 0.0
          solnData(VELZ_VAR,i,j,k)= 0.0

#ifdef MAGX_VAR
          solnData(MAGX_VAR,i,j,k)= sim_Bx0
          solnData(MAGY_VAR,i,j,k)= sim_By0
          solnData(MAGZ_VAR,i,j,k)= sim_Bz0

          solnData(MAGP_VAR,i,j,k)= .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                   solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
          solnData(DIVB_VAR,i,j,k)= 0.
          eBZone = 0.5 * dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
#else
          eBZone = 0
#endif

          ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                       solnData(VELX_VAR:VELZ_VAR,i,j,k))

         ! specific internal energy
          eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.)/solnData(DENS_VAR,i,j,k)

          solnData(ENER_VAR,i,j,k)=max(eintZone + ekinZone + eBZone, sim_smallP)
          solnData(EINT_VAR,i,j,k)=max(eintZone, sim_smallP)
          solnData(GAMC_VAR,i,j,k)=sim_gamma
          solnData(GAME_VAR,i,j,k)=sim_gamma
#ifdef TDUS_VAR
          solnData(TDUS_VAR,i,j,k) = sim_tdust
#endif
#ifdef IHP_SPEC
! if ionization is excluded from setup, this loop will simply not execute,
! so pre-processor ifdef is not strictly needed.
! but, it is useful to explicitly mark where code hooks into ray-tracing unit
          do n=1,NSPECIES
            solnData(SPECIES_BEGIN+n-1,i,j,k) = massFrac_box(n)
          end do
          !solnData(IHA_SPEC,i,j,k) = massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)
          !solnData(IHP_SPEC,i,j,k) = massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)
#endif
        enddo
     enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID)

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
