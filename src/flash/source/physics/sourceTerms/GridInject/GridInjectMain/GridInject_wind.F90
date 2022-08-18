!!****f* source/physics/sourceTerms/GridInject/GridInjectMain/GridInject_wind
!!
!! NAME
!!
!!  GridInject_wind
!!
!! SYNOPSIS
!!
!!  call GridInject_wind(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: injectMassIn
!!      real(IN)          :: injectVelocityIn
!!      real(IN)          :: twind
!!      real(IN)          :: dt
!!      real(INOUT)       :: bgDens
!!      optional,logical(IN) :: snap_to_grid
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit a wind onto grid with both kinetic and thermal energy;
!!  winds can be modified (mass loaded) via runtime parameters.
!!
!! ARGUMENTS
!!
!!  xloc              : where to inject
!!  yloc              : where to inject
!!  zloc              : where to inject
!!  injectMassIn      : amount of mass to inject
!!  injectVelocityIn  : velocity of injected mass
!!  twind             : lifetime of injected wind
!!  dt                : timestep for wind injection, to get dm/dt
!!  bgDens            : initial background density of wind.
!!                      if bgDens==0 and runtime param wind_var_radius is set,
!!                      calculate and return to caller, else use input value
!!  snap_to_grid      : for testing/debugging
!!
!!***
subroutine GridInject_wind (xloc, yloc, zloc, injectMassIn, injectVelocityIn, &
                            twind, dt, bgDens, snap_to_grid)

  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, &
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

  use GridInject_data, ONLY: gi_maxref, &
                             wind_cons_quant, wind_add_therm_e, &
                             wind_mass_load, wind_target_temp, &
                             wind_min_radius, wind_var_radius, &
                             wind_perturb, wind_perturb_stdev
  use GridInject_interface, only : GridInject_getInjBlks
  use gi_interface, ONLY: gi_distanceVector, gi_overlap, gi_normal_rand

  implicit none

#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  real, intent(IN)    :: xloc, yloc, zloc
  real, intent(IN)    :: injectMassIn, injectVelocityIn, twind, dt
  real, intent(INOUT) :: bgDens
  logical, optional, intent(IN) :: snap_to_grid

  real :: xx, yy, zz, radx, rady, radz, rad, dVol
  real :: rad2, del2
  real :: loc(3), delta(3), locDelta(3), cell_top(3), cell_bot(3)

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

  real :: overlap_frac, sumOverlap, solidAngle
  real :: injectRadius, injectRadiusMax, injectMass, injectVelocity
  real :: injE, deltaKE, deltaEint
  real :: oldVel(3), injVel(3), newVel(3), newVelSq(3)
  real :: oldDens, dDens, newDens

  real :: mass_load_factor, refVel, randnum
  integer :: refineLevel
  logical :: allRefined, calcBgDens

  ! ================================================
  ! SETUP blagh
  ! ================================================

  ! Get and communicate the proper cell size at wind inject position
  call Grid_getBlkIDFromPos([xloc, yloc, zloc], blockID, procID, gr_meshComm)
  if (gr_meshMe .eq. procID) then
    call Grid_getDeltas(blockID, locDelta)
  end if
  call MPI_Bcast(locDelta, MDIM, MPI_DOUBLE_PRECISION, procID, gr_meshComm, ierr)

  ! Place the star in the center of a cell.
  if (present(snap_to_grid) .and. snap_to_grid) then
    loc = (floor([xloc, yloc, zloc]/locDelta) + 0.5)*locDelta
  else
    loc = [xloc, yloc, zloc]
  end if

  calcBgDens      = .false.
  injectRadius    = 3.5*sqrt(3.0)*minval(locDelta)
  injectRadiusMax = 3.5*sqrt(3.0)*minval(locDelta)

  if (wind_var_radius .and. bgDens /= 0.0) then
    ! Weaver et al 1977 eqn 12
    injectRadius = 0.74296*(injectMassIn/dt/bgDens)**0.3 * injectVelocityIn**0.1 * twind**0.4
    ! floor at sqrt(3)*dx, ceil at 3.5*sqrt(3)*dx
    injectRadius = max(injectRadius, sqrt(3.0)*minval(locDelta))
    injectRadius = min(injectRadius, injectRadiusMax)
  end if

  ! ================================================
  ! count and record blocks partially within injectRadius,
  ! check that they are maximally refined
  ! ================================================

  call GridInject_getInjBlks(loc(1), loc(2), loc(3), injectRadius, injBlocks, injBlkNum)

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

  call MPI_AllReduce(MPI_IN_PLACE, allRefined, 1, MPI_LOGICAL, MPI_LAND, &
                     gr_meshComm, ierr)

  if (.not. allRefined) then
    if (gr_meshMe == 0) then
      print *, "wind injection blocks not fully refined; not injecting"
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
  sumOverlap = 0.0

  if (bgDens == 0.0) then
    calcBgDens = .true.
  end if

  do n = 1, injBlkNum
    blockID = injBlocks(n)

    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,delta)
    del2  = minval(delta)**2

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

          ! compute dx, dy, dz accounting for periodic BCs
          call gi_distanceVector(xx, yy, zz, loc(1), loc(2), loc(3), radx, rady, radz)

          rad2 = radx**2 + rady**2 + radz**2
          rad  = sqrt(rad2)

          if (rad .le. injectRadius + 0.5*sqrt(sum(delta**2))) then

            cell_bot = [ sign(abs(xx) - 0.5*delta(1), xx), &
                         sign(abs(yy) - 0.5*delta(2), yy), &
                         sign(abs(zz) - 0.5*delta(3), zz) ]

            cell_top = [ sign(abs(xx) + 0.5*delta(1), xx), &
                         sign(abs(yy) + 0.5*delta(2), yy), &
                         sign(abs(zz) + 0.5*delta(3), zz) ]

            ! calculate overlap between injection sphere and this cell.
            ! returns value to overlap_frac.
            call gi_overlap(1, injectRadius, loc, cell_bot, cell_top, 20, overlap_frac)

            if (.not. wind_mass_load) then
              ! flash4.2.2-rad commit 88fc47a (2018-08-13), inject winds
              ! such that mass density decreases as r^2.
              ! Calculate the overlapping solid angle of a square at distance rad from the sphere.
              ! Calculation from MPIA: http://www.mpia.de/~mathar/public/mathar20051002.pdf
              solidAngle   = 4.*acos(sqrt((1.+del2/(2.*rad2))/(1.+del2/(2.*rad2)+(del2/rad2)**2.)))
              overlap_frac = overlap_frac*solidAngle
              ! not sure if del2=minval(delta)**2 is correct for non-cube cells
              ! -AT, 2019 Jul 03
            end if
            ! sum all weightings for normalization later

            if (rad .ge. wind_min_radius) then
              sumOverlap = sumOverlap + overlap_frac
              injOverlap(n,i,j,k) = overlap_frac
              if (rad .gt. 0.0) then
                ! normalized components of the star --> cell center vector
                injUnitVel(n,i,j,k,1:3) = [radx/rad, rady/rad, radz/rad]
              end if
              if (calcBgDens) then
                bgDens = bgDens + overlap_frac*solnData(DENS_VAR,i,j,k)
              end if
            end if
          end if  ! rad .le. injectRadius + ...

        end do !x
      end do !y
    end do !z

    call Grid_releaseBlkPtr(blockID, solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

  end do ! injBlkNum

  call MPI_AllReduce(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, gr_meshComm, ierr)

  injOverlap = injOverlap / sumOverlap  ! normalize to get frac overlap wrt entire injection sphere

  if (calcBgDens) then
    call MPI_AllReduce(MPI_IN_PLACE, bgDens, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, gr_meshComm, ierr)
    bgDens = bgDens / sumOverlap
  end if

  ! =============================
  ! Inject directly onto solnData
  ! =============================

  ! mass_load_factor is the ratio between input velocity and the velocity to
  ! get T = wind_target_temp [K] shocked gas.
  ! and use it to set the mass and velocity to get the same
  ! input energy. - JW

  if (wind_mass_load) then
    ! Draine 2011 eqn 36.28. Note this gas is always ionized by
    ! radiation first. - JW
    refVel = sqrt(wind_target_temp/1.38d7)*1e8
    if (wind_cons_quant .eq. "momentum") then
      mass_load_factor = injectVelocityIn/refVel
      injectMass       = injectMassIn*mass_load_factor
      injectVelocity   = injectVelocityIn/mass_load_factor
    else if (wind_cons_quant .eq. "energy") then
      mass_load_factor = (injectVelocityIn/refVel)**2
      injectMass       = injectMassIn*mass_load_factor
      injectVelocity   = injectVelocityIn/sqrt(mass_load_factor)
    endif
  else
    injectMass         = injectMassIn
    injectVelocity     = injectVelocityIn
  end if

  do n = 1, injBlkNum
    blockID = injBlocks(n)

    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,delta)
    dVol = product(delta)

    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

    ! loop over all zones in block, no guard cells
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

          ! Round off errors in calculations lead to small
          ! changes in thermal energy that have big
          ! consequences in EOS calls. So just skip.
          if (injOverlap(n,i,j,k) .le. 0.0) cycle

          oldDens = solnData(DENS_VAR,i,j,k)
          oldVel  = solnData(VELX_VAR:VELZ_VAR,i,j,k)
          injVel = injUnitVel(n,i,j,k,1:3) * injectVelocity
          if (wind_perturb) then
            call gi_normal_rand(1.0, wind_perturb_stdev, randnum)
            injVel = injVel * randnum
          end if

          if (wind_cons_quant .eq. "momentum") then

            ! includes mass-loading adjustment
            dDens   = injOverlap(n,i,j,k)*injectMass/dVol
            newDens = oldDens + dDens ! m + delta_m
            newVel  = (injVel*dDens + oldVel*oldDens) / newDens

            ! To make up for the energy loss from conserving momentum,
            ! use the original velocity from injE to add back enough thermal energy
            ! to get back to the proper total energy at the end, and
            ! NOT the modified velocity from mass loading or velocity
            ! calculated from momentum conservation. - JW
            deltaEint = 0.0
            if (wind_add_therm_e) then
              deltaKE = 0.5*sum(newVel**2)*newDens - 0.5*sum(oldVel**2)*oldDens  ! erg cm^-3
              ! Total mechanical energy injected by the wind. - JW
              injE = 0.5*injectMassIn*injectVelocityIn**2
              ! I think this is wrong, since the "correct" energy being added is the mass
              ! loaded energy, and the mass loading is the stuff falling in from the shell, which
              ! has zero kinetic energy. So we *should lose some kinetic energy here*! - JW
              deltaEint = injE*injOverlap(n,i,j,k)/dVol - deltaKE  ! erg cm^-3
            end if

            solnData(EINT_VAR,i,j,k) = (solnData(EINT_VAR,i,j,k)*oldDens + deltaEint)/newDens
            solnData(ENER_VAR,i,j,k) = 0.5*sum(newVel**2) + solnData(EINT_VAR,i,j,k)
            solnData(VELX_VAR:VELZ_VAR,i,j,k) = newVel
            solnData(DENS_VAR,i,j,k) = newDens

          else if (wind_cons_quant .eq. "energy") then

            dDens   = injOverlap(n,i,j,k)*injectMass/dVol
            newDens = oldDens + dDens

            newVelSq = (oldDens*sign(oldVel**2,oldVel) &
                        + dDens*sign(injVel**2,injVel)) / newDens
            newVel = sign(sqrt(abs(newVelSq)),newVelSq)

            ! must update EINT before ENER
            solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)*oldDens/newDens
            solnData(ENER_VAR,i,j,k) = 0.5*sum(newVel**2) + solnData(EINT_VAR,i,j,k)
            solnData(VELX_VAR:VELZ_VAR,i,j,k) = newVel
            solnData(DENS_VAR,i,j,k) = newDens

          ! method from Simpson et al. 2015. Should be same as normal energy conservation
!          else if (wind_cons_quant .eq. "energy_simpson") then
!
!            dDens = injOverlap(n,i,j,k)*injectMass/dVol
!            newDens = oldDens + dDens
!            newVel = (injVel*dDens + oldVel*oldDens) / newDens
!
!            a_ = 0.5 / newDens
!            b_ = oldDens * oldVel / newDens
!            c_ = (oldDens * oldVel)**2 / ( 2 * newDens ) &
!                 - 0.5 * oldDens * oldVel**2 &
!                 - 0.5 * dDens * injVel**2
!            deltaP = - b_ + sqrt(b_**2 - 4*a_*c_)/(2*a_)
!            ! new velocity is calculated by conserving momentum when dDens moving
!            ! at injectVelocity is added to oldDens moving at original velocity
!            solnData(VELX_VAR:VELZ_VAR,i,j,k) = deltaP / dDens
!            solnData(DENS_VAR,i,j,k) = newDens
!
          end if  ! wind_cons_quant

          ! call pointwise eos
          pointLimit(:,IAXIS) = i
          pointLimit(:,JAXIS) = j
          pointLimit(:,KAXIS) = k
          call Eos_wrapped(MODE_DENS_EI,pointLimit,blockID)

        end do !x
      end do !y
    end do !z

    call Grid_releaseBlkPtr(blockID, solnData)

  end do ! block loop

  deallocate(injOverlap)
  deallocate(injUnitVel)

  call Grid_notifysolnDataUpdate()

  return
end subroutine GridInject_wind
