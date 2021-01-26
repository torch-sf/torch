!!****if* source/Simulation/SimulationMain/Cube/Simulation_initBlock
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
!!  a specified block.  This version sets up cube read from file.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine Simulation_initBlock (blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords
  use Eos_interface, ONLY : Eos_wrapped 

 implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
  
  integer,intent(IN) ::  blockId

  
  integer  :: i, j, k, istat, n
  integer  :: ii, jj, kk, ix, iy, iz, ixp1, iyp1, izp1
  real     :: vx, vy, vz, pres, rho, e, ek, gpot
  real     :: xx, yy, zz, dx, dy, dz, p, q, r

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  real, dimension(:,:,:,:),pointer :: solnData, facexData, faceyData, facezData
  integer,dimension(MDIM) :: startingPos
  real          :: del(MDIM)
  real, dimension(NSPECIES) :: massFrac_box

  logical :: gcell = .true.

     
!! ---------------------------------------------------------------------------

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  call Grid_getDeltas(blockId,del)

  call Grid_getBlkPtr(blockId,solnData)

#if NFACE_VARS > 0  
  ! For B-field assignment - SCL 10/2020
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif
  
#ifdef IHP_SPEC
  massFrac_box(IHP_SPEC-SPECIES_BEGIN+1)    = sim_init_Hp 
  massFrac_box(IHA_SPEC-SPECIES_BEGIN+1)    = (1.0 - sim_init_Hp)
#endif
  
!  call Multispecies_getSum(GAMMA, sim_gamma, massFrac_box)

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

        rho  = (1.-p)*(1.-q)*(1.-r)*sim_densArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_densArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_densArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_densArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_densArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_densArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_densArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_densArr(ixp1, iyp1, izp1)
        pres = (1.-p)*(1.-q)*(1.-r)*sim_presArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_presArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_presArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_presArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_presArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_presArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_presArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_presArr(ixp1, iyp1, izp1)
        gpot = (1.-p)*(1.-q)*(1.-r)*sim_gpotArr(ix  , iy  , iz  ) &
        &    + (1.-p)*(1.-q)*    r *sim_gpotArr(ix  , iy  , izp1) &
        &    + (1.-p)*    q *(1.-r)*sim_gpotArr(ix  , iyp1, iz  ) &
        &    + (1.-p)*    q *    r *sim_gpotArr(ix  , iyp1, izp1) &
        &    +     p *(1.-q)*(1.-r)*sim_gpotArr(ixp1, iy  , iz  ) &
        &    +     p *(1.-q)*    r *sim_gpotArr(ixp1, iy  , izp1) &
        &    +     p *    q *(1.-r)*sim_gpotArr(ixp1, iyp1, iz  ) &
        &    +     p *    q *    r *sim_gpotArr(ixp1, iyp1, izp1)
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

        rho  = max(rho, smlrho)
        pres = max(pres, smallp)

        !if ((rho < 1e-25) .or. (rho > 1e-17)) then
        !    print *, xx, yy, zz, ix, iy, iz, p, q, r, rho, sim_densArr(ix,iy,iz)
        !endif

      

        ! assume gamma-law equation of state
        ek  = 0.5*(vx*vx + vy*vy + vz*vz)
        e   = pres/(sim_gamma-1.)
        e   = e/rho + ek
        e   = max (e, smallP)


        solnData(DENS_VAR,i,j,k)=rho
        solnData(PRES_VAR,i,j,k)=pres
        solnData(ENER_VAR,i,j,k)=e
        solnData(GAME_VAR,i,j,k)=sim_gamma
        solnData(GAMC_VAR,i,j,k)=sim_gamma
        solnData(VELX_VAR,i,j,k)=vx
        solnData(VELY_VAR,i,j,k)=vy
        solnData(VELZ_VAR,i,j,k)=vz
#ifdef TDUS_VAR
        solnData(TDUS_VAR,i,j,k)=sim_tdust
#endif
#ifdef MAGX_VAR
        ! Adding Bfield data to block centers
        solnData(MAGX_VAR,i,j,k)= sim_magx
        solnData(MAGY_VAR,i,j,k)= sim_magy
        solnData(MAGZ_VAR,i,j,k)= sim_magz
#endif

#if NFACE_VARS > 0
  ! Adding Bfield data to block faces - SCL 10/2020
  if (sim_killdivb) then      
     facexData(:,:,:,:)=sim_magx
     faceyData(:,:,:,:)=sim_magy
     if (NDIM == 3) facezData(:,:,:,:)=sim_magz
  endif
#endif

#ifdef IHP_SPEC
! if ionization is excluded from setup, this loop will simply not execute,
! so pre-processor ifdef is not strictly needed.
! but, it is useful to explicitly mark where code hooks into ray-tracing unit
          do n=1,NSPECIES
            solnData(SPECIES_BEGIN+n-1,i,j,k) = massFrac_box(n)
          end do
#endif

      enddo
    enddo
  enddo
  call Grid_releaseBlkPtr(blockID, solnData)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif
     
  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID)
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock
