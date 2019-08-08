!!****if* source/Simulation/SimulationMain/Cube_Velocity/Simulation_initBlock
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
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data

  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr,&
                             Grid_releaseBlkPtr, Grid_getCellCoords
  use Multispecies_interface, ONLY : Multispecies_getSum

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Flash_mpi.h"

  integer,intent(IN) ::  blockId

  ! blk sizes, coordinates, pointers to cell-centered UNK variables
  ! loop vars to set solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: i, j, k, n, sizeX, sizeY, sizeZ, istat
  logical :: gcell = .true.
  real    :: xx, yy, zz
  real,allocatable,dimension(:) :: xCoord, yCoord, zCoord
  real, dimension(:,:,:,:),pointer :: solnData
  real, dimension(NSPECIES) :: massFrac_box  ! Multispecies
  real :: pres, rho, eint, ek, eB, integratedGz, expV  ! stratbox setup

  ! variables used to init cube (velocity data only), w/ linear interpolation
  integer  :: ix, iy, iz, ixp1, iyp1, izp1
  real     :: dx, dy, dz, p, q, r
  real     :: vx, vy, vz

  massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)    = sim_init_Hp
  massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)    = (1.0 - sim_init_Hp)
  call Multispecies_getSum(GAMMA, sim_gamma, massFrac_box)  ! weighted avg

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)

  ! Grid_putPointData and Grid_getBlkPtr/Grid_releaseBlkPtr work equally well
  ! http://flash.uchicago.edu/pipermail/flash-users/2016-December/002142.html
  call Grid_getBlkPtr(blockId,solnData)

  ! ------------------------------------------------------
  ! init velocities ONLY; code from setup-cube-USM problem
  ! ------------------------------------------------------

  dx = (sim_xMax - sim_xMin) / sim_nCD(IAXIS)
  dy = (sim_yMax - sim_yMin) / sim_nCD(JAXIS)
  dz = (sim_zMax - sim_zMin) / sim_nCD(KAXIS)
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz = zCoord(k)
    iz = (zz - (sim_zMin + 0.5*dz)) / dz
    if (iz < 1) then
      iz = 1
      izp1 = 1
      r = 1.0
    else if (iz >= sim_nCD(KAXIS)) then
      iz = sim_nCD(KAXIS)
      izp1 = sim_nCD(KAXIS)
      r = 1.0
    else
      izp1 = iz + 1
      r = (zz - (sim_zMin+dz*(iz+0.5))) / dz
    endif
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
      yy = yCoord(j)
      iy = (yy - (sim_yMin + 0.5*dy)) / dy
      if (iy < 1) then
        iy = 1
        iyp1 = 1
        q = 1.0
      else if (iy >= sim_nCD(JAXIS)) then
        iy = sim_nCD(JAXIS)
        iyp1 = sim_nCD(JAXIS)
        q = 1.0
      else
        iyp1 = iy + 1
        q = (yy - (sim_yMin+dy*(iy+0.5))) / dy
      endif
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
        xx = xCoord(i)
        ix = (xx - (sim_xMin + 0.5*dx)) / dx
        if (ix < 1) then
          ix = 1
          ixp1 = 1
          p = 1.0
        else if (ix >= sim_nCD(IAXIS)) then
          ix = sim_nCD(IAXIS)
          ixp1 = sim_nCD(IAXIS)
          p = 1.0
        else
          ixp1 = ix + 1
          p = (xx - (sim_xMin+dx*(ix+0.5))) / dx
        endif

        ! linear interpolation to cell position
        vx   = (1.-p)*(1.-q)*(1.-r)*sim_velxArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_velxArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_velxArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_velxArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_velxArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_velxArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_velxArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_velxArr(ixp1, iyp1, izp1)
        vy   = (1.-p)*(1.-q)*(1.-r)*sim_velyArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_velyArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_velyArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_velyArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_velyArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_velyArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_velyArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_velyArr(ixp1, iyp1, izp1)
        vz   = (1.-p)*(1.-q)*(1.-r)*sim_velzArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_velzArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_velzArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_velzArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_velzArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_velzArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_velzArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_velzArr(ixp1, iyp1, izp1)

        solnData(VELX_VAR,i,j,k)=vx
        solnData(VELY_VAR,i,j,k)=vy
        solnData(VELZ_VAR,i,j,k)=vz

      enddo
    enddo
  enddo

  ! -------------------------------------
  ! init stratified box density, pressure
  ! -------------------------------------

  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
    xx = xCoord(i)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
       yy = yCoord(j)
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
          zz = zCoord(k)

          if (sim_useStrat) then

            ! Hydrostatic balance of NFW potential and stratified isothermal gas
            ! Ibanez-Mejia et al. 2016, eqns 5-7
            integratedGz = -sim_aParm1 * dsqrt(zz**2 + sim_aParm3**2) &
                - 0.5e0*sim_aParm2 * zz**2 + 1./3. * sim_aParm4 * zz**3 &
                + sim_aParm1 * sim_aParm3
            expV = dexp(integratedGz * sim_rho / sim_p)
          else
            expV = 1  ! disable stratification
          endif

          pres = max(expV*sim_p, sim_pIGM, sim_smallp)
          rho = max(expV*sim_rho, sim_rhoIGM, sim_smlrho)

          solnData(PRES_VAR,i,j,k)= pres
          solnData(DENS_VAR,i,j,k)= rho
          !solnData(VELX_VAR,i,j,k) = 0.
          !solnData(VELY_VAR,i,j,k) = 0.
          !solnData(VELZ_VAR,i,j,k) = 0.

#ifdef MAGX_VAR
          ! B field stratified in z-direction
          solnData(MAGX_VAR,i,j,k)= sim_magx*sqrt(expV)
          solnData(MAGY_VAR,i,j,k)= sim_magy*sqrt(expV)
          solnData(MAGZ_VAR,i,j,k)= sim_magz*sqrt(expV)
          ! AT20190221 does user need to set MAGP, DIVB?..
          solnData(MAGP_VAR,i,j,k)= .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                   solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
          solnData(DIVB_VAR,i,j,k)= 0.
          eB = 0.5 * dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                 solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
#else
          eB = 0
#endif
          ek = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                   solnData(VELX_VAR:VELZ_VAR,i,j,k))
          eint = pres/(sim_gamma-1.)/rho  ! _specific_ internal energy

          solnData(ENER_VAR,i,j,k) = max(eint + ek + eB, sim_smallP)
          solnData(EINT_VAR,i,j,k) = max(eint, sim_smallP)
          solnData(GAMC_VAR,i,j,k) = sim_gamma
          solnData(GAME_VAR,i,j,k) = sim_gamma
          solnData(TDUS_VAR,i,j,k) = sim_tdust
          do n=1,NSPECIES
            solnData(SPECIES_BEGIN+n-1,i,j,k) = massFrac_box(n)
          end do
          !solnData(IHA_SPEC,i,j,k) = massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)
          !solnData(IHP_SPEC,i,j,k) = massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)
        enddo
     enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID)

  call Grid_releaseBlkPtr(blockID, solnData)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return

end subroutine Simulation_initBlock
