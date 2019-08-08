!!****f* source/physics/sourceTerms/GridInject/GridInjectMain/GridInject_thermalSN
!!
!! NAME
!!
!!  GridInject_thermalSN
!!
!! SYNOPSIS
!!
!!  call GridInject_thermalSN(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: energy
!!      real(IN)          :: mass
!!      real(IN)          :: r_min
!!      real(IN)          :: r_max
!!      real(IN)          :: nms
!!      real(OUT)         :: r_exp
!!      real(OUT)         :: m_exp
!!      real(OUT)         :: rho_avg
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit SN onto grid as thermal energy
!!
!! ARGUMENTS
!!
!!  xloc    : where to inject
!!  yloc    : where to inject
!!  zloc    : where to inject
!!  energy  : amount of energy to inject
!!  mass    : how much mass should lie in the explosion sphere?
!!  r_min   : minimum SN explosion radius
!!  r_max   : minimum SN explosion radius
!!  nms     : number of mass shells to search between r_min and r_max
!!
!! RETURNS
!!
!!  r_exp : explosion radius
!!  m_exp : mass within explosion sphere
!!  rho_avg : average density within explosion sphere
!!
!!***
subroutine GridInject_thermalSN (xloc, yloc, zloc, energy, mass, &
    r_min, r_max, nms, r_exp, m_exp, rho_avg)

  use Driver_interface, ONLY : Driver_getMype
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr, &
                             Grid_getDeltas, &
                             Grid_getListOfBlocks, &
                             Grid_releaseBlkPtr, &
                             Grid_getCellCoords

  use gi_interface, ONLY : gi_distance

  implicit none

#include "constants.h"
#include "Eos.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real, intent(IN)    :: xloc, yloc, zloc, energy, mass, r_min, r_max
  integer, intent(IN) :: nms
  real, intent(OUT)   :: r_exp, m_exp, rho_avg

  ! block/zone loops
  integer :: injBlockCount
  integer, dimension(MAXBLOCKS) :: injBlockList
  integer         :: blockID, thisBlock
  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord
  real, allocatable, dimension(:)   :: dx, dy, dz
  integer         :: xSizeCoord, ySizeCoord, zSizeCoord
  logical         :: getGuardCells = .true.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, pointLimit
  real, dimension(MDIM)      :: del

  ! communication and mass calculation
  real, dimension(nms)     :: vol_sum_loc, vol_sum_tot
  real, dimension(nms)     :: mass_sum_tot, mass_sum_loc
  integer :: ierr, mype

  ! iterators
  integer :: i, j, k, im, im_exp
  real    :: xx, yy, zz, r

  ! variables for the SN explosion
  real  :: ei, ek, rho, edens

  call Driver_getMype(MESH_COMM, mype)

  call GridInject_getInjBlks(xloc, yloc, zloc, r_max, &
                             injBlockList, injBlockCount)

  !====================================================
  ! find a sphere containing the requested mass
  !====================================================

  mass_sum_loc = 0.e0  ! sum of mass within initial explosion sphere on this proc
  mass_sum_tot = 0.e0  ! sum of local masses from all procs
  vol_sum_loc  = 0.e0  ! sum of zone volumes in initial explosion sphere on this proc
  vol_sum_tot  = 0.e0  ! sum of local zone volumes in initial explosion sphere from all procs

  do thisBlock = 1, injBlockCount
    blockID = injBlockList(thisBlock)

    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,del) !grid spacing dx dz dy

    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    xSizeCoord = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
    ySizeCoord = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
    zSizeCoord = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))
    call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

    allocate(dx(xSizeCoord))
    allocate(dy(ySizeCoord))
    allocate(dz(zSizeCoord))
    dx(:) = del(IAXIS)
    dy(:) = del(JAXIS)
    dz(:) = del(KAXIS)

    ! loop over all zones in block
    ! for 2d k is just 1
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      zz = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

          rho = solnData(DENS_VAR,i,j,k)
          r = gi_distance(xx, yy, zz, xloc, yloc, zloc)
          do im = 1, nms
            if ( r <= r_min*(r_max/r_min)**((im-1)/(nms-1.0)) ) then
              mass_sum_loc(im) = mass_sum_loc(im) + rho*dx(i)*dy(j)*dz(k)
              vol_sum_loc(im) = vol_sum_loc(im) + dx(i)*dy(j)*dz(k)
            endif
          enddo

        enddo !x
      enddo !y
    enddo !z

    call Grid_releaseBlkPtr(blockID,solnData)

    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

    deallocate(dx)
    deallocate(dy)
    deallocate(dz)

  enddo ! block loop

  ! TODO make faster by only talking to neighbouring cpus, or using FLASH4 communication routines
  !  most SN information just has to be send to neighbours

  call MPI_AllReduce (mass_sum_loc, mass_sum_tot, nms, &
      & MPI_Double_Precision, MPI_Sum, MPI_Comm_World, ierr)

  call MPI_AllReduce (vol_sum_loc, vol_sum_tot, nms, &
      & MPI_Double_Precision, MPI_Sum, MPI_Comm_World, ierr)

  ! ==============================
  ! choose explosion radius
  ! ==============================

  im_exp = 0
  do im = 1, nms
    if (mass_sum_tot(im) > mass) then
      im_exp = im
      exit
    endif
  enddo

  if (im_exp == 0) then
    im_exp = nms
  endif

  rho_avg = mass_sum_tot(im_exp)/vol_sum_tot(im_exp)
  r_exp = r_min*(r_max/r_min)**((im_exp-1)/(nms-1.0))
  m_exp = mass_sum_tot(im_exp)
  edens = energy / vol_sum_tot(im_exp)  ! [erg/cm^3]

  !====================================================
  ! insert SN energy
  !====================================================

  do thisBlock = 1, injBlockCount
    blockID = injBlockList(thisBlock)

    call Grid_getBlkPtr(blockID,solnData)

    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    xSizeCoord = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
    ySizeCoord = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
    zSizeCoord = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))
    call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)

    ! loop over all zones in block, no guard cells
    ! for 2d k is just 1
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      zz = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

          r = gi_distance(xx, yy, zz, xloc, yloc, zloc)

          if (r <= r_exp) then

            ei  = solnData(EINT_VAR,i,j,k)
            rho = solnData(DENS_VAR,i,j,k)
            ek  = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
            &            solnData(VELY_VAR,i,j,k)**2 + &
            &            solnData(VELZ_VAR,i,j,k)**2)

            ! apply SN as constant specific heating rate (uniform in mass)
            ! versus constant volumetric heating rate
            ! this will heat high-density zones more, and avoid over-heating
            ! low-density zones, as compared to applying constant heat/volume.
            solnData(DENS_VAR,i,j,k) = rho
            solnData(EINT_VAR,i,j,k) = ei + edens/rho_avg
            solnData(ENER_VAR,i,j,k) = ei + edens/rho_avg + ek
#ifdef PHIO_VAR
            ! ionize all gas if radiation transport is included
            ! this halves mean mass and hence temperature through Eos call
            solnData(IHA_SPEC,i,j,k) = 0d0
            solnData(IHP_SPEC,i,j,k) = 1d0
#endif
            ! call pointwise eos
            pointLimit(:,IAXIS) = i
            pointLimit(:,JAXIS) = j
            pointLimit(:,KAXIS) = k
            call Eos_wrapped(MODE_DENS_EI,pointLimit,blockID)

          endif

        enddo !x
      enddo !y
    enddo !z

    call Grid_releaseBlkPtr(blockID,solnData)

    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

  enddo ! block loop

  call Grid_notifySolnDataUpdate()

  return
end subroutine GridInject_thermalSN
