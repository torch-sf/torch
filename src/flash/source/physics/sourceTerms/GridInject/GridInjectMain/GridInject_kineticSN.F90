!!****f* source/physics/sourceTerms/GridInject/GridInjectMain/GridInject_kineticSN
!!
!! NAME
!!
!!  GridInject_kineticSN
!!
!! SYNOPSIS
!!
!!  call GridInject_kineticSN(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: energy
!!      real(IN)          :: mass
!!      optional,logical(IN) :: snap_to_grid
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit SN onto grid with both kinetic and thermal energy
!!
!! ARGUMENTS
!!
!!  xloc    : where to inject
!!  yloc    : where to inject
!!  zloc    : where to inject
!!  energy  : amount of energy to inject
!!  mass    : amount of mass to inject
!!  snap_to_grid : for testing/debugging
!!
!!***
subroutine GridInject_kineticSN (xloc, yloc, zloc, energy, mass, snap_to_grid)

  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY: Grid_getBlkIDFromPos, &
                            Grid_getBlkIndexLimits, &
                            Grid_getBlkPtr, &
                            Grid_getBlkRefineLevel, &
                            Grid_getCellCoords, &
                            Grid_getDeltas, &
                            Grid_getMinCellSize, &
                            Grid_notifySolnDataUpdate, &
                            Grid_releaseBlkPtr
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get

  use Grid_data, ONLY: gr_meshComm, gr_meshMe
  use tree, ONLY: refine, derefine, stay

  use GridInject_data, only : gi_maxref
  use GridInject_interface, only : GridInject_getInjBlks
  use gi_interface, only : gi_distanceVector, gi_overlap

  implicit none

#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  ! For Simpson+ (2015) kinetic SN prescription
  real, parameter :: pc = 3.086e18
  real, parameter :: mu = 1.3*1.6726e-24

  real, intent(IN)              :: xloc, yloc, zloc, energy, mass
  logical, optional, intent(IN) :: snap_to_grid

  real :: xx, yy, zz, radx, rady, radz, rad, dVol
  real :: loc(3), delta(3), SNdelta(3), cell_top(3), cell_bot(3)

  ! block loops
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, pointLimit
  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: xSizeCoord, ySizeCoord, zSizeCoord
  logical :: getGuardCells = .true.
  integer :: i, j, k, n, ierr
  integer :: blockID, procID

  integer :: injBlkNum, injBlocks(MAXBLOCKS)
  real, allocatable, dimension(:,:,:,:) :: injOverlap
  real, allocatable, dimension(:,:,:,:,:) :: injUnitVel

  logical :: allRefined
  integer :: refineLevel

  real :: overlap_frac, sumOverlap
  real :: injRadius, injKE, injMom
  real :: bgDens, rPDS, tPDS, fracKin

  real :: deltaKE, deltaEint
  real :: oldDens, newDens
  real :: oldVel(3), newVel(3)

  ! ================================================
  ! Which blocks overlap SN injection sphere,
  ! and what is the spatial resolution in those blocks?
  ! ================================================

  ! TODO this is gonna be weird for SN at block refinement boundaries...
  ! we can pre-emptively ask code to refine, before explosion, like the
  ! logic in SN unit of adjusting timestep before explosion.
  ! BUT, code needs to know that SN failed, we should either
  ! force the refine now, or tell the caller "plz refine and call again"

  ! Get and communicate the proper cell size at SN position
  call Grid_getBlkIDFromPos([xloc, yloc, zloc], blockID, procID, gr_meshComm)
  if (gr_meshMe .eq. procID) then
    call Grid_getDeltas(blockID, SNdelta)
  end if
  call MPI_Bcast(SNdelta, MDIM, MPI_DOUBLE_PRECISION, procID, gr_meshComm, ierr)

  ! Place the star in the center of a cell.
  if (present(snap_to_grid) .and. snap_to_grid) then
    loc = (floor([xloc, yloc, zloc]/SNdelta) + 0.5)*SNdelta
  else
    loc = [xloc, yloc, zloc]
  end if

  injRadius = 3.0*minval(SNdelta)

  ! ================================================
  ! count and record blocks partially within injRadius,
  ! check that they are maximally refined
  ! ================================================

  call GridInject_getInjBlks(loc(1), loc(2), loc(3), injRadius, injBlocks, injBlkNum)

  allRefined = .true.

  do n = 1, injBlkNum
    blockID = injBlocks(n)

    call Grid_getBlkRefineLevel(blockID, refineLevel)
    if (refineLevel < gi_maxref) then
      allRefined = .false.
      ! set this block to be refined! - JW
      refine(blockID) = .true.
      derefine(blockID) = .false.
      stay(blockID) = .true.
    end if
  end do

  call MPI_ALLREDUCE(MPI_IN_PLACE, allRefined, 1, MPI_LOGICAL, MPI_LAND, &
                     gr_meshComm, ierr)

  if (.not. allRefined) then
    if (gr_meshMe == 0) then
      print *, "SN blocks not fully refined; skipping SN for now"
    end if
    return
  end if

  ! =============================================================
  ! Compute overlap fractions, velocity unit vector, and bkg dens
  ! =============================================================

  ! "Fortran is cool with 0-sized arrays"
  ! https://stackoverflow.com/a/35251929
#ifdef FIXEDBLOCKSIZE
  allocate(injOverlap(injBlkNum, &
                    & GRID_ILO:GRID_IHI, &
                    & GRID_JLO:GRID_JHI, &
                    & GRID_KLO:GRID_KHI))
  allocate(injUnitVel(injBlkNum, &
                & GRID_ILO:GRID_IHI, &
                & GRID_JLO:GRID_JHI, &
                & GRID_KLO:GRID_KHI, &
                & 1:3))
#else
  call Driver_abortFlash("Non-fixed block size is not supported")
#endif

  injOverlap = 0.0  ! correct default for non-overlap case
  injUnitVel = 0.0

  bgDens = 0.0
  sumOverlap = 0.0

  do n = 1, injBlkNum
    blockID = injBlocks(n)

    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,delta)

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
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      zz = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

          ! compute radx, rady, radz accounting for periodic BCs
          call gi_distanceVector(xx, yy, zz, loc(1), loc(2), loc(3), radx, rady, radz)

          rad = sqrt(radx**2 + rady**2 + radz**2)

          ! Calculate the overlapping areas of the sphere and this cell.
          if (rad .le. injRadius + 0.5*sqrt(sum(delta(:)**2))) then

            cell_bot = [ sign(abs(xx) - 0.5*delta(1), xx), &
                         sign(abs(yy) - 0.5*delta(2), yy), &
                         sign(abs(zz) - 0.5*delta(3), zz) ]

            cell_top = [ sign(abs(xx) + 0.5*delta(1), xx), &
                         sign(abs(yy) + 0.5*delta(2), yy), &
                         sign(abs(zz) + 0.5*delta(3), zz) ]

            ! gi_overlap(...) samples ~nsteps**3 points to do integral; that's
            ! 125,000 points for nsteps=50.  Call on as few cells as possible.
            call gi_overlap(1, injRadius, loc, cell_bot, cell_top, 50, overlap_frac)

            injOverlap(n,i,j,k) = overlap_frac

            ! unit velocity vector
            if (rad .gt. 0.0) then
              injUnitVel(n,i,j,k,1:3) = [radx/rad, rady/rad, radz/rad]
            end if

            bgDens      = bgDens + overlap_frac * solnData(DENS_VAR,i,j,k)
            sumOverlap  = sumOverlap + overlap_frac

          end if

        end do !x
      end do !y
    end do !z

    call Grid_releaseBlkPtr(blockID, solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

  end do ! block loop

  call MPI_AllReduce(MPI_IN_PLACE, bgDens, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, gr_meshComm, ierr)

  call MPI_AllReduce(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, gr_meshComm, ierr)

  bgDens = bgDens / sumOverlap
  injOverlap = injOverlap / sumOverlap
  ! normalize to get fractional overlap wrt entire explosion sphere

  ! ================================================================
  ! Calculate the kinetic energy fraction as in Simpson et al 2015.
  ! ================================================================

  ! Time (units: 10^3 yr) when pressure-dominated snow-plow phase begins.
  ! Z=1 metallicity assumed.
  tPDS = 26.5 * (energy/1e51)**(3.0/14.0) * (bgDens/mu)**(-4.0/7.0)

  ! SNR radius (parsecs) at t=tPDS.
  rPDS = 18.5 * (energy/1e51)**(2.0/7.0) * (bgDens/mu)**(-3.0/7.0)

  ! Calculate the kinetic fraction. Assumes mean molecular weight of 1.0
  fracKin = 3.97e-6 * (bgDens/mu) * rPDS**7 * tPDS**(-2) * (minval(SNdelta)/pc)**(-2) * (energy/1e51)**(-1)
  if (4.5*minval(SNdelta) < rPDS*pc) then
    fracKin = 0.0d0
  end if

  ! =============================
  ! Inject directly onto solnData
  ! =============================

  injKE   = fracKin*energy
  injMom  = sqrt(2.0*injKE*mass)

  do n = 1, injBlkNum
    blockID = injBlocks(n)

    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,delta)
    dVol = product(delta)

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
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

          ! Round off errors in calculations lead to small
          ! changes in thermal energy that have big
          ! consequences in EOS calls. So just skip.
          if (injOverlap(n,i,j,k) .le. 0.0) cycle

          oldDens = solnData(DENS_VAR,i,j,k)
          newDens = solnData(DENS_VAR,i,j,k) + injOverlap(n,i,j,k)*mass/dVol

          oldVel  = solnData(VELX_VAR:VELZ_VAR,i,j,k)
          newVel  = ( injMom*injUnitVel(n,i,j,k,:)*injOverlap(n,i,j,k)/dVol + oldVel*oldDens ) / newDens

          ! Energy is in erg/g, same regardless of refinement.
          deltaKE = 0.5 * (newDens*sum(newVel**2) - oldDens*sum(oldVel**2)) / newDens
          ! Thermal energy is difference b/t injected specific E and KE.
          deltaEint = energy*injOverlap(n,i,j,k)/(newDens*dVol) - deltaKE

          ! must update EINT before ENER
          solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)*oldDens/newDens + deltaEint
          solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + 0.5*sum(newVel**2)
          solnData(VELX_VAR:VELZ_VAR,i,j,k) = newVel
          solnData(DENS_VAR,i,j,k) = newDens

          ! call pointwise eos
          pointLimit(:,IAXIS) = i
          pointLimit(:,JAXIS) = j
          pointLimit(:,KAXIS) = k
          call Eos_wrapped(MODE_DENS_EI,pointLimit,blockID)

        end do !x
      end do !y
    end do !z

    call Grid_releaseBlkPtr(blockID,solndata)

    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

  end do ! block loop

  deallocate(injOverlap)
  deallocate(injUnitVel)

  call Grid_notifySolnDataUpdate()

  return
end subroutine GridInject_kineticSN
