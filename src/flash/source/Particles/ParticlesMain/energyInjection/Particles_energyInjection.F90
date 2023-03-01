
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! subroutine inject_directSN
!!!
!!! Authors: Joshua Wall and Andrew Pellegrino
!!!          Drexel University
!!!          Summer and Fall 2016
!!!
!!! A routine for injection of SN energy. Adapted from our winds model.
!!!
!!! This routine injects both mass and energy in a 
!!! momentum conserving fashion. Note that we require
!!! both momentum and energy conservation with the
!!! momentum method.
!!!
!!! In the momentum method we chose the injection radius such that
!!! the radius of the free flowing kinetic wind (R_1) is resolved. If it is not resolved
!!! by the radius, we make the injection radius smaller down to dx. If the
!!! injection radius is smaller than R_1, we make the injection radius larger
!!! up to injectionRadiusMax, which is a runtime parameter. The default is
!!! 5 cells at the diagional in 3d (3.5*sqrt(3.0)*dx).
!!!
!!! Note that both this and the further requirement that energy be
!!! conserved by inelasitc collision of the wind mass and the cell mass
!!! enforces that under-resolved winds always end up on the grid as
!!! thermal energy (as they should if the shocked wind region is not
!!! resolved). In this way, one can argue that the shocked hot wind is
!!! a result of inelastic collision. - Joshua Wall
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine Particles_energyInjection(injectEnergyin, fracKin, injectMassIn, xloc, yloc, zloc, dt)

!#define DEBUG
!#define DEBUG2
#define DEBUG_ENERGY

use Grid_data, ONLY: gr_meshComm, gr_meshMe

use Hydro_data, ONLY: hy_cfl

use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getMinCellSize, &
    Grid_notifySolnDataUpdate, Grid_getBlkRefineLevel, Grid_getBlkIDFromPos
     
use Driver_interface, ONLY : Driver_abortFlash
    
use Eos_interface, ONLY : Eos_wrapped

use pt_enerInjInterface, only : overlap, sphere_and_cell_frac

use RuntimeParameters_interface, ONLY: RuntimeParameters_get

use tree, ONLY: nodetype, coord, bsize, lnblocks, refine, derefine, stay

#include "Flash.h"
#include "constants.h"

implicit none

#include "Flash_mpi.h"

integer, parameter :: dp = kind(1.d0)

real(dp), intent(in)    :: injectMassIn, injectEnergyIn, xloc, yloc, zloc
real(dp), intent(out)   :: dt, fracKin

logical, save :: first_call = .true.
logical :: iHaveInjectBlk

logical :: snap_to_grid
logical :: useTimeStep

real(dp) :: mass, energy, injectVelocity, bgDens

real(dp) :: loc(3)

real(dp) :: cell_top(3), cell_bot(3)
real(dp) :: dVol
real(dp) :: x,y,z, dx, dy, dz

real(dp), save :: injectRadius, delta(3), SNdelta(3), min_delta
real(dp) :: overlap_frac,  sumOverlap, oldDens, dDens, newDens
real(dp) :: newVel(3), totP(3), totE, globalDeltaE, globalDeltaP, emech, pmech
real(dp) :: oldE, newE, oldP, newP
real(dp) :: oldVel(3), injVel(3), newVelSq(3)
real(dp) :: idir, jdir, kdir, rad, xvel, yvel, zvel
real(dp) :: xcoll, ycoll, zcoll, d2coll

integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), blockID, procID
integer :: i, j, k, m, n
integer :: ind_array(3)


integer :: ierr

real(dp), pointer, dimension(:,:,:,:) :: solndata


! Indices: blockID, i, j, k, [dimension]
real(dp), allocatable, dimension(:,:,:,:) :: injectDataOverlap
real(dp), allocatable, dimension(:,:,:,:,:) :: injectDataVel
integer, allocatable, dimension(:) :: localInjectBlocks

real(dp) :: cs2(NXB,NYB,NZB), v2(NXB,NYB,NZB), old_dt
real(dp), save :: gamma_


integer :: injBlkNum
real(dp) :: blkCtr(3), blkSize(3)  ! code requires MDIM=3
real(dp), parameter :: yr = (60.0_dp**2.0)*24.0_dp*365.25_dp, pc = 3.086e18
real(dp), parameter :: solarMass = 1.989d33, mu = 1.3*1.6726e-24

real(dp) :: starMass, ThermE, KinE, injE
real(dp) :: injectKinE, injectThermE, injectMom, injectMass
real(dp) :: globalTE, globalKE, largestTE, largestKE
real(dp) :: deltaKinE, deltaThermE, finalThermE, deltaE
real(dp) :: initialKE, initialTE

real(dp) :: sumMass ! Sum of all the mass array.


real(dp) :: rPDS, tPDS, four_point_five_dx


if (first_call) then
    call RuntimeParameters_get("gamma", gamma_)
    !call RuntimeParameters_get("lrefine_max", maxref)
    first_call = .false.
end if


iHaveInjectBlk = .false.
useTimeStep = .true.
bgDens = 0.0_dp
sumOverlap = 0.0_dp
sumMass = 0.0_dp

!!  USE WITH CAUTION! 
!!  This snap_to_grid code does not correctly position star particles!
!!  floor() will move the star to roughly (0,0,0) in parsecs.
!!  See jets code inject_direct.F90 for issues with snap_to_grid and 
!!  improved version of the code.  -SA 20230301

snap_to_grid = .false. !.true. ! For testing / debugging.
globalDeltaE = 0.0
globalDeltaP = 0.0
globalTE     = 0.0
globalKE     = 0.0
largestKE    = 0.0
largestTE    = 0.0
deltaThermE  = 0.0
deltaKinE    = 0.0


mass   = injectMassIn
energy = injectEnergyIn

call Grid_getBlkIDFromPos([xloc, yloc, zloc] ,blockID ,procID, gr_meshComm)

! Note only the proc that owns this block has the proper cell size... other procs
! may have a block with this number but its the wrong block (and maybe a parent, 
! wrong size, etc).

if (gr_meshMe .eq. procID) then
  call Grid_getDeltas(blockID,SNdelta)
end if

! Send the proper cell size to all procs.
call MPI_Bcast(SNdelta, 3, MPI_DOUBLE_PRECISION, procID, gr_meshComm, ierr)

call Grid_getMinCellSize(min_delta)

! Place the star in the center of a cell?
if (snap_to_grid) then
    loc = floor([xloc, yloc, zloc]) + 0.5_dp*SNdelta
else
    loc = [xloc, yloc, zloc]    
end if

injectRadius = 3.0_dp*minval(SNdelta)
four_point_five_dx = 4.5_dp*minval(SNdelta)


#ifdef DEBUG_ENERGY
if (gr_meshMe == 0) &
    write(*, '(A,ES13.3e3)') "SN injection radius / dx = ", injectRadius/SNdelta(1)
#endif


#ifdef DEBUG2
!print *, "inject_direct.F90 is conserving ", conserved_quant
write(*, '(A,ES13.3e3, I4)') " injection radius = ", injectRadius, gr_meshMe
print*, "loc =", loc, gr_meshMe
#endif

! count # of blocks which are at least partially within injectRadius, check that
! they are maximally refined
injBlkNum = 0
do blockID = 1, lnblocks
    if(nodetype(blockID) == LEAF) then
        ! exact collision detection for sphere and rectangular prism
        ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
        call Grid_getBlkCenterCoords(blockID,blkCtr)
        call Grid_getBlkPhysicalSize(blockID,blkSize)
        ! point within block that is closest to SN location
        xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
        ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
        zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))
        if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadius**2) then
            injBlkNum = injBlkNum + 1
            iHaveInjectBlk = .true.
        end if
    end if
end do

! build array of block IDs for all blocks to be injected. Get their center
! distances from the injection star. Then for each cell in each block, calculate
! its overlap with the injection sphere and store the value. Store the
! components of a velocity vector pointing radially outwards from the star
if (iHaveInjectBlk) then

#ifdef DEBUG2
print *, "Found", injBlkNum, "injection blocks on proc ", gr_meshMe
#endif

    allocate(localInjectBlocks(injBlkNum))
    allocate(injectDataOverlap(injBlkNum,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI, &
                                                            GRID_KLO:GRID_KHI))
    allocate(injectDataVel(injBlkNum,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI, &
                                                        GRID_KLO:GRID_KHI,1:3))

    localInjectBlocks = 0
    injectDataOverlap = 0.0d0
    injectDataVel     = 0.0d0


#ifdef DEBUG2
    print *, "Allocations done"
#endif

    n = 1
    do blockID = 1, lnblocks
        if(nodetype(blockID) == LEAF) then
            ! exact collision detection for sphere and rectangular prism
            ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
            call Grid_getBlkCenterCoords(blockID,blkCtr)
            call Grid_getBlkPhysicalSize(blockID,blkSize)
            ! point within block that is closest to SN location
            xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
            ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
            zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))
            if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadius**2) then
                localInjectBlocks(n) = blockID
                n = n + 1
            end if
        end if
    end do

#ifdef DEBUG2
    print *, "Found injection blocks:", localInjectBlocks, "on proc", gr_MeshMe
#endif
 
    do n = 1, injBlkNum
        blockID = localInjectBlocks(n)
        call Grid_getDeltas(blockID,delta)
        call Grid_getBlkPtr(blockID, solndata)
        do k = GRID_KLO, GRID_KHI
            do j = GRID_JLO, GRID_JHI
                do i = GRID_ILO, GRID_IHI
                    ! since we have checked that all cells are refined, use
                    ! mindelta

                    x = (i - NGUARD - NXB/2.0 - 0.5)*delta(1) + coord(1,blockID)
                    y = (j - NGUARD - NYB/2.0 - 0.5)*delta(2) + coord(2,blockID)
                    z = (k - NGUARD - NZB/2.0 - 0.5)*delta(3) + coord(3,blockID)

                    ! exact collision detection for sphere and rectangular prism
                    ! https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection#Sphere_vs._AABB
                    ! point within cell that is closest to SN location
                    xcoll = max(x-0.5*delta(1),min(loc(1),x+0.5*delta(1)))
                    ycoll = max(y-0.5*delta(2),min(loc(2),y+0.5*delta(2)))
                    zcoll = max(z-0.5*delta(3),min(loc(3),z+0.5*delta(3)))
                    d2coll = (xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2

                    ! is cell outside injection sphere?
                    if (d2coll > injectRadius**2) then
                        cycle
                    end if

                    ! get overlapping volume of inject sphere and this cell,
                    ! modified by a tapered center-weighting within overlap(..)
                    cell_bot = [ sign(abs(x) - 0.5*delta(1), x), &
                                 sign(abs(y) - 0.5*delta(2), y), &
                                 sign(abs(z) - 0.5*delta(3), z) ]
                    cell_top = [ sign(abs(x) + 0.5*delta(1), x), &
                                 sign(abs(y) + 0.5*delta(2), y), &
                                 sign(abs(z) + 0.5*delta(3), z) ]
                    call overlap(1, injectRadius, loc, cell_bot, &
                                 cell_top, 10, overlap_frac)

                    dx = x - loc(1)
                    dy = y - loc(2)
                    dz = z - loc(3)

                    rad = sqrt(dx**2 + dy**2 + dz**2)
                    
                    ! normalized components of the star --> cell center vector 
                    if (rad .ne. 0.0_dp) then
                    
                        idir = dx/rad
                        jdir = dy/rad
                        kdir = dz/rad
                        
                    else
                      
                        idir = 0.0_dp
                        jdir = 0.0_dp
                        kdir = 0.0_dp
                      
                    endif

                    ! velocity of the wind points radially outwards
                        
                    xvel = idir !* injectVelocity
                    yvel = jdir !* injectVelocity
                    zvel = kdir !* injectVelocity
                    
!                    if (abs(sqrt(xvel**2.0_dp+yvel**2.0_dp+zvel**2.0_dp) - injectVelocity)/injectVelocity > 1.0d-6) then
!                      write(*,'(A,2ES13.3E3)') "vels don't match: calc vel, inj vel= ",  &
!                           sqrt(xvel**2.0_dp+yvel**2.0_dp+zvel**2.0_dp), injectVelocity
!                      stop
!                    end if
                    
                    ! Get the background density to estimate what the inject
                    ! radius should be physically. - JW

                    bgDens = bgDens + overlap_frac*solndata(DENS_VAR, i, j, k)

                    ! sum all weightings for normalization later

                    if (overlap_frac .gt. 0.0d0) then
                        sumOverlap = sumOverlap + overlap_frac
                        injectDataOverlap(n,i,j,k) = overlap_frac
                        injectDataVel(n,i,j,k,1:3) = [xvel,yvel,zvel]
                    end if
                end do
            end do
        end do
        call Grid_releaseBlkPtr(blockID, solndata)
    end do
#ifdef DEBUG2
    print *, "Calculated overlaps"
#endif
end if

#ifdef DEBUG2
print*, "Proc ", gr_meshMe, " about to call MPI with sumOverlap = ", sumOverlap
#endif
call MPI_ALLREDUCE(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

#ifdef DEBUG2
    print*, "Before MPI background density is", bgDens
#endif

call MPI_ALLREDUCE(MPI_IN_PLACE, bgDens, 1, MPI_DOUBLE_PRECISION, &
                                        MPI_SUM, gr_meshComm, ierr)

bgDens = bgDens / sumOverlap


#ifdef DEBUG2
    if (gr_meshMe == 0) then
        print *, "Total overlap weight is", sumOverlap
        print*, "Background density is", bgDens
    end if
#endif


if (fracKin .lt. 0.0) then

  ! Now calculate the fraction of kinetic energy as in Simpson et al 2015.
  ! First find the pressure dominated snow plow phase time (Z=1 metallicity assumed).
  ! Also note that tPDS is in units of 10^3 years.
  tPDS = 26.5*((energy/1e51)**(3.0/14.0))*((bgDens/mu)**(-4.0/7.0))
#ifdef DEBUG_ENERGY
  if (gr_meshMe == 0) print*, "tPDS =", tPDS
#endif
  ! Then calculate the radius of the SN front at this time.
  ! rPDS is in units of parsecs.
  rPDS = 18.5*((energy/1e51)**(2.0/7.0))*((bgDens/mu)**(-3.0/7.0))
#ifdef DEBUG_ENERGY
  if (gr_meshMe == 0) print*, "rPDS= ", rPDS
#endif
  ! Finally calculate the kinetic fraction. Here we assume mean molecular weight of 1.0
  fracKin = 3.97e-6*(bgDens/mu)*(rPDS**7.0)*(tPDS**(-2.0))*((SNdelta(1)/pc)**(-2.0))*((energy/1e51)**(-1.0))
  if (four_point_five_dx < (rPDS*pc)) fracKin = 0.0d0
#ifdef DEBUG_ENERGY
  if (gr_meshMe == 0) print*, "fracKin", fracKin
#endif
end if

injectMass   = mass 
injectKinE   = (fracKin*energy)
injectThermE = ((1.0-fracKin)*energy)
injectMom    = sqrt(2.0*injectKinE*InjectMass)

#ifdef DEBUG
write(*,'(A,ES13.3E3)') 'InjectMass =', injectMass
write(*,'(A,ES13.3E3)') 'InjectKinE =', injectKinE
write(*,'(A,ES13.3E3)') 'InjectThermE =', injectThermE
write(*,'(A,ES13.3E3)') 'InjectMom =', injectMom
#endif

if (iHaveInjectBlk) then
    !if (conserved_quant .eq. "momentum") then
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getDeltas(blockID,delta)
            dVol = product(delta)
            call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
            call Grid_getBlkPtr(blockID, solndata)

            do k = GRID_KLO, GRID_KHI
                do j = GRID_JLO, GRID_JHI
                    do i = GRID_ILO, GRID_IHI
                      
                      ! Round off errors in calculations lead to small
                      ! changes in thermal energy that have big
                      ! consequences in EOS calls. So just skip.
                      if (injectDataOverlap(n,i,j,k) .le. 0.0_dp) cycle
                        
                      dDens   = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
                      oldDens = solndata(DENS_VAR,i,j,k)
                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
                      oldVel  = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                      injVel  = injectDataVel(n,i,j,k,1:3)
                      totP    = injVel*injectMom*injectDataOverlap(n,i,j,k)/sumOverlap/dVol + oldVel * oldDens
                      newVel  = totP / newDens
                      
                      ! How much mass are we really putting in? - JW
                      sumMass = sumMass + dDens*dVol
                      initialKE = 0.5_dp*sum(oldVel**2.0_dp)
                      initialTE = solndata(EINT_VAR,i,j,k)
!<<<<<<< HEAD
! My code - JW

                      !oldE = oldDens * sum(oldVel**2)
                      oldE = (0.5_dp*sum(oldVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*oldDens*dVol

                      oldP = oldDens * sqrt(sum((oldVel+0.0_dp)**2.0_dp))

		      ! Note in Flash energy is in erg/g so that the energy is the same regardless of refinement. - JW
		      ! Also don't forget the already existing internal energy in this cell. - JW
		      
		      ! Specific KE added to this cell. - JW
                      deltaKinE = (0.5_dp * newDens * sum(newVel**2.0_dp) &
				      - 0.5_dp * oldDens * sum(oldVel**2.0_dp)) &
				      / newDens
                      
                      !deltaKinE = 0.5_dp * sum(newVel**2.0_dp) - 0.5_dp * sum(oldVel**2.0_dp)
                      
		      ! Thermal energy should be the difference b/t the injected specific E
		      ! and the actual specific kinetic energy. - JW

		      ! Final TE in this cell = portion of mech
		      ! energy added to the cell + the existing TE in cell.
              ! Note what we did here. In order to make up for the energy loss from conserving momentum
              ! we used the original velocity from injE to ensure that we added back enough energy
              ! in thermal to get back to the proper total energy at the end, and
              ! NOT the modified velocity from mass loading or what was calculated from momentum
              ! conservation. - JW
		      ! - KE added to cell. - JW !0.5_dp*dDens*sum(injVel**2.0_dp)/newDens
!		              finalThermE = injE*injectDataOverlap(n,i,j,k)/sumOverlap/(newDens*dVol)  &
!		                  + solndata(EINT_VAR, i, j, k)*oldDens/newDens &
!				          - deltaKinE
		      
                      !deltaThermE = finalThermE - solndata(EINT_VAR, i, j, k)
                      deltaThermE = energy*injectDataOverlap(n,i,j,k)/sumOverlap/(newDens*dVol)  &
                                  - deltaKinE
		      
		      !ThermE = (0.5_dp * newDens * sum(newVel**2.0_dp) & 
		      !       + solndata(EINT_VAR, i, j, k) - oldE

		      
		      ! Now add the "missing" kinetic energy as thermal energy. - JW

                      !solndata(EINT_VAR, i, j, k) =  finalThermE
                      solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR,i,j,k) + deltaThermE
                      solndata(ENER_VAR, i, j, k) =  0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k)
!=======
! Andrew code
                      !oldE = oldDens * sum(oldVel**2)
                      !oldP = oldDens * sqrt(sum((oldVel+0.0_dp)**2))

                      !totE = 0.5 * (dDens * sum(injVel**2) + oldDens * sum(oldVel**2))

                      !new thermal energy per unit mass 
                      !solndata(EINT_VAR, i, j, k) = totE/newDens - 0.5 * sum(newVel**2)
!>>>>>>> 8a465d18d19d157e666df203cfdfcf0cdb4d63fe

                      ! new velocity is calculated by conserving momentum when dDens moving
                      ! at injectVelocity is added to oldDens moving at original velocity

                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
                      solndata(DENS_VAR, i, j, k) = newDens

		      
		              if (deltaKinE .gt. largestKE) largestKE = deltaKinE
		              if (deltaThermE .gt. largestTE) largestTE = deltaThermE
		      
		      ! Delta E (not specific!) added to this cell. - JW
                      deltaE = (deltaKinE + deltaThermE)*newDens*dVol !+ oldE

                      !newE = newDens * sum(newVel**2)
                      !newE = sum(newVel**2.0_dp)
                      newE = (0.5_dp * sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*newDens*dVol
		      !newP = newDens * sqrt(sum(newVel**2))
                      newP = sqrt(sum((totP+0.0_dp)**2.0_dp))
                      !globalDeltaE = globalDeltaE + 0.5 * (newE &
		      !             - oldE )*dVol
                           

#ifdef DEBUG
                      if ((injectDataOverlap(n,i,j,k) .le. 0.0_dp) .and. &
                          (dDens .gt. 0.0d0 .or. deltaKinE .gt. 0.0d0 .or. deltaThermE .gt. 0.0 &
                           .or. abs(newE - oldE) .gt. 0.0d0)) then
                           
                           write(*,*) "BS detected! Overlap is zero but someone is up to no good!"
                           write(*,'(A,ES13.3E3)') "Delta KE= ", deltaKinE
                           write(*,'(A,ES13.3E3)') "Delta TE= ", deltaThermE
                           write(*,'(A,ES13.3E3)') "dDens= ", dDens
                           write(*,'(A,2ES13.3E3)') "oldDens, newDens= ", oldDens, newDens
                           write(*,'(A,6ES13.3E3)') "oldVel, newVel= ", oldVel, newVel
                           write(*,'(A,3ES13.3E3)') "injVel =", injVel
                           write(*,'(A,2ES13.3E3)') "old KE, old TE= ", initialKE, initialTE
                           write(*,'(A,6ES13.3E3)') "new KE, new TE= ", 0.5_dp*sum(newVel**2.0_dp), solndata(EINT_VAR, i, j, k)
                           write(*,'(A,ES13.3E3)') "initial Total E= ", oldE
                           write(*,'(A,ES13.3E3)') "final Total E= ", newE
                           write(*,'(A,ES13.3E3)') "(newE-oldE)= ", (newE-oldE)
                           call flush(6)
                      end if

                    if (injectDataOverlap(n,i,j,k) > 0.0) then
                      write(*,'(A,ES13.3E3)') "Delta KE= ", deltaKinE
                      write(*,'(A,ES13.3E3)') "Delta TE= ", deltaThermE
                      write(*,'(A,ES13.3E3)') "dDens= ", dDens
                      write(*,'(A,2ES13.3E3)') "oldDens, newDens= ", oldDens, newDens
                      write(*,'(A,6ES13.3E3)') "oldVel, newVel= ", oldVel, newVel
                      write(*,'(A,3ES13.3E3)') "injVel =", injVel
                      write(*,'(A,2ES13.3E3)') "old KE, old TE= ", initialKE, initialTE
                      write(*,'(A,6ES13.3E3)') "new KE, new TE= ", 0.5_dp*sum(newVel**2.0_dp), solndata(EINT_VAR, i, j, k)
                      write(*,'(A,ES13.3E3)') "final TE= ", solndata(EINT_VAR, i, j, k)
                      write(*,'(A,ES13.3E3)') "initial Total E= ", oldE
                      write(*,'(A,ES13.3E3)') "final Total E= ", 0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR, i, j, k)
                      write(*,'(A,ES13.3E3)') "deltaE - (newE-oldE)= ", deltaE - (newE-oldE)
                      call flush(6)
                    end if
#endif
		              globalKE     = globalKE + deltaKinE*newDens*dVol
                      globalTE     = globalTE + deltaThermE*newDens*dVol
                      globalDeltaE = globalDeltaE +  deltaE !newE - oldE !deltaE
                      globalDeltaP = globalDeltaP + (newP - oldP)*dVol

                    end do
                end do
            end do

            call Grid_releaseBlkPtr(blockID, solndata)
            call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
        end do
!    else if (conserved_quant .eq. "energy") then
!        do n = 1, injBlkNum
!            blockID = localInjectBlocks(n)
!            call Grid_getDeltas(blockID,delta)
!            dVol = product(delta)
!            call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!            call Grid_getBlkPtr(blockID, solndata)

!            do k = GRID_KLO, GRID_KHI
!                do j = GRID_JLO, GRID_JHI
!                    do i = GRID_ILO, GRID_IHI
                        
!                      dDens   = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
!                      oldDens = solndata(DENS_VAR,i,j,k)
!                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
!                      oldVel  = solndata(VELX_VAR:VELZ_VAR, i, j, k)
!                      injVel  = injectDataVel(n,i,j,k,1:3)
!                      oldE    = oldDens * sum(oldVel**2)
!                      oldP    = oldDens * sqrt(sum(oldVel**2))

!                      ! Same as above, but conserving kinetic energy

!                      newVelSq = oldDens/newDens * sign(oldVel**2,oldVel) &
!                               + dDens/newDens * sign(injVel**2,injVel)
!                      newVel = sign(sqrt(abs(newVelSq)),newVelSq)

!                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
!                      solndata(DENS_VAR, i, j, k) = newDens

!                      newE = newDens * sum(newVel**2)
!                      newP = newDens * sqrt(sum(newVel**2))
!                      globalDeltaE = globalDeltaE + 0.5 * (newE - oldE)*dVol
!                      globalDeltaP = globalDeltaP + (newP - oldP)*dVol
!                    end do
!                end do
!            end do

!            call Grid_releaseBlkPtr(blockID, solndata)
!#ifdef DEBUG
!            print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
!                                                   blockID, gr_meshMe, injBlkNum
!#endif
!            call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
!        end do

    ! energy conserving method from Simpson et al. 2015. Should be equivalent to
    ! normal energy conservation
    !else if (conserved_quant .eq. "energy_simpson") then
    !    do n = 1, injBlkNum
    !        blockID = localInjectBlocks(n)
    !        call Grid_getBlkPtr(blockID, solndata)

    !        do k = GRID_KLO, GRID_KHI
    !            do j = GRID_JLO, GRID_JHI
    !                do i = GRID_ILO, GRID_IHI
    !                    
    !                  dDens = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
    !                  oldDens = solndata(DENS_VAR,i,j,k)
    !                  newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
    !                  oldVel = solndata(VELX_VAR:VELZ_VAR, i, j, k)
    !                  injVel = injectDataVel(n,i,j,k,1:3)
    !                  totP = injVel * dDens + oldVel * oldDens
    !                  newVel = totP / newDens

    !                  a_ = 0.5 / newDens
    !                  b_ = oldDens * oldVel / newDens
    !                  c_ = (oldDens * oldVel)**2 / ( 2 * newDens ) &
    !                    - 0.5 * oldDens * oldVel**2 &
    !                    - 0.5 * dDens * injVel**2
    !                  deltaP = - b_ + sqrt(b_**2 - 4*a_*c_)/(2*a_)

    !                  !solndata(EINT_VAR, i, j, k) = totE/newDens - 0.5 * sum(newVel**2)

    !                  ! new velocity is calculated by conserving momentum when dDens moving
    !                  ! at injectVelocity is added to oldDens moving at original velocity

    !                  solndata(VELX_VAR:VELZ_VAR, i, j, k) = deltaP / dDens
    !                  solndata(DENS_VAR, i, j, k) = newDens

    !                end do
    !            end do
    !        end do

    !        call Grid_releaseBlkPtr(blockID, solndata)
    !    end do

    !end if

end if


#ifdef DEBUG_ENERGY
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalTE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestTE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestKE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    
#ifdef DEBUG2
print*, "Proc ", gr_meshMe, " about to call MPI with sumMass = ", sumMass
#endif
call MPI_ALLREDUCE(MPI_IN_PLACE, sumMass, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

    if (gr_meshMe == 0) then
        emech = InjectEnergyIn
        pmech = sqrt(2.0 * injectMassIn * fracKin* InjectEnergyIn)
        write(*,'(A,ES10.3,A)') "We got ", globalDeltaE, "ergs"
        write(*,'(A,F10.3,A)') "Error in injected SN total E is ", abs(globalDeltaE - emech)/emech*100, "%"
        write(*,'(A,F10.3,A)') "Percentage of energy that is thermal is ", globalTE/globalDeltaE*100.0, "%"
        write(*,'(A,F10.3,A)') "Percentage of energy that is kinetic is ", globalKE/globalDeltaE*100.0, "%"
        if (pmech .ne. 0.0d0) write(*,'(A,F10.3,A)') "Error in injected SN pressure is ", abs(globalDeltaP - pmech)/pmech*100, "%"
        write(*,'(A,ES10.3e2,A)') "Largest TE is ", largestTE
        write(*,'(A,ES10.3e2,A)') "Largest KE is ", largestKE
        write(*,'(A,ES10.3e2,X,A,ES10.3e2)') "Total mass from SN is ", injectMass, "Injected mass is ", sumMass
    end if
#endif


call Grid_notifySolnDataUpdate()

dt = 1d99

if (iHaveInjectBlk) then

    ! Update the EOS variables if we have an affected block for only those blocks.

    do n=1, injBlkNum
        blockID = localInjectBlocks(n)
!#ifdef DEBUG
!        print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
!                                                   blockID, gr_meshMe, injBlkNum
!#endif

    ! Moved into the actual loops where values were modified to be
    ! absolutely certain the correct blocks get this call. - JW
        !call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        !call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
	

    	! calculate min. crossing time of wind velocity on all injection blocks
    	if (useTimeStep) then

            call Grid_getBlkPtr(blockID, solndata)

	! square of the max sound speed in the block
	!c_sound^2 = gamma * k * T / m_H

            cs2 = gamma_ * 1.3807d-16 * &
                solndata(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) &
                / mu 

	! square of magnitude of velocities in all cells of block
            v2 = sum(solndata(VELX_VAR:VELZ_VAR,GRID_ILO:GRID_IHI, &
                GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0,1)
                
    ! NOTE NOTE NOTE ! We need to also calculate the Alfven velocity here
    ! or otherwise take it into account in Hydro_computeDt (<--- this should be done now. JW) - JW
    
	! minimum dt to capture wind movement in this block
    ! Lets use the hydro CFL here. - JW
            old_dt = hy_cfl * min_delta/sqrt( maxval(cs2 + v2) )
            if (old_dt .lt. dt) then
                dt = old_dt
            end if
#ifdef DEBUG2
            print*, "Max temp on proc, blk", gr_meshMe, blockID, "is", maxval(solndata(TEMP_VAR,:,:,:))
#endif
            call Grid_releaseBlkPtr(blockID, solndata)

        end if


    end do

    deallocate(localInjectBlocks)
    deallocate(injectDataOverlap)
    deallocate(injectDataVel)
#ifdef DEBUG2
    print *, "Deallocating done for proc", gr_meshMe
#endif
    !call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
!#ifdef DEBUG
!    print *, "Guard cells filled on proc", gr_meshMe
!#endif
end if

call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_MIN, gr_meshComm, ierr)
                                            
if (useTimeStep .and. gr_meshMe == 0) then
#ifdef DEBUG_ENERGY
    write(*,'(A,ES10.3)') "Timestep set by Particles_energyInjection = ", dt
#endif
end if

call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)

#ifdef DEBUG2
print *, "Exiting Particles_energyInjection for proc", gr_meshMe
#endif
end subroutine Particles_energyInjection
