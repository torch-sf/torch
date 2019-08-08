
subroutine inject_direct(loc_in, injectMassIn, injectVelocityIn, star_mass, twind, dt)

!#define DEBUG

use Grid_data, ONLY: gr_meshComm, gr_meshMe

use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getMinCellSize, &
    Grid_notifySolnDataUpdate, Grid_getBlkRefineLevel
    
use Driver_interface, ONLY : Driver_abortFlash
    
use Eos_interface, ONLY : Eos_wrapped

use pt_windInterface, only : overlap, sphere_and_cell_frac

use Particles_sinkData, ONLY: min_wind_dt

use RuntimeParameters_interface, ONLY: RuntimeParameters_get

use PhysicalConstants_interface, ONLY: PhysicalConstants_get

use tree, ONLY: nodetype, coord, bsize, lnblocks

#include "Flash.h"
#include "constants.h"

implicit none

#include "Flash_mpi.h"

integer, parameter :: dp = kind(1.d0)

real(dp), intent(in)    :: loc_in(3)
real(dp), intent(in)    :: injectMassIn, injectVelocityIn, twind, star_mass
real(dp), intent(inout) :: dt

logical, save :: first_call = .true.
logical :: iHaveInjectBlk
logical :: iHaveUnrefined

logical :: snap_to_grid
logical :: rampVelocity
logical, save :: isRestart
logical :: useTimeStep

real(dp) :: injectMass, injectVelocity

real(dp) :: star_x, star_y, star_z, loc(3)
real(dp) :: delta(3)
real(dp) :: cell_top(3), cell_bot(3)
real(dp) :: dVol
real(dp) :: x,y,z, dx, dy, dz

real(dp) :: injectRadius
real(dp) :: overlap_frac,  sumOverlap, oldDens, dDens, newDens
real(dp) :: oldVel(3), injVel(3), newVelSq(3)
real(dp) :: idir, jdir, kdir, rad, xvel, yvel, zvel

integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), blockID
integer :: i, j, k, m, n
integer :: ind_array(3)

!integer, parameter :: trans_tag = 103180
!integer, parameter :: mess_tag  = 801031
integer :: messages, trans_status(MPI_STATUS_SIZE), ierr

real(dp), pointer, dimension(:,:,:,:) :: solndata


! Indices: blockID, i, j, k, [dimension]
real(dp), allocatable, dimension(:,:,:,:) :: injectDataOverlap
real(dp), allocatable, dimension(:,:,:,:,:) :: injectDataVel
integer, allocatable, dimension(:) :: localInjectBlocks

real(dp) :: maxdir
real(dp) :: rem, dist, mult
integer :: d
real(dp) :: cs2(NXB,NYB,NZB), v2(NXB,NYB,NZB), old_dt
real(dp), save :: gamma_
real(dp) :: mass_load_factor
integer :: refineLevel
integer, save :: maxref
real(dp), save :: rampVel = 0d0
real(dp), save :: rampVelStep = 3.0d7 !3.0d6 ! cm/s
real(dp) :: injectVelocityMax
character(len=10), save :: conserved_quant
integer :: injBlkNum
real(dp) :: blkCenterDist(3), blkCornerDist(3)
real(dp), parameter :: yr = 60.0_dp**2.0*24.0_dp*365.25_dp
real(dp), parameter :: solarMass = 1.989d33
real(dp) :: v_esc
real(dp), save :: newton


if (first_call) then
    call RuntimeParameters_get("gamma", gamma_)
    call RuntimeParameters_get("lrefine_max", maxref)
    call RuntimeParameters_get("restart", isRestart)
    call RuntimeParameters_get("cons_quant", conserved_quant)
    call PhysicalConstants_get("Newton", newton)
    first_call = .false.
end if

#ifdef DEBUG
print *, "inject_direct.F90 is conserving ", conserved_quant
#endif

iHaveUnrefined = .false.
iHaveInjectBlk = .false.
useTimeStep = .true.
rampVelocity = .true.
sumOverlap = 0.0_dp
mass_load_factor = 1.0_dp
snap_to_grid = .false.

call Grid_getMinCellSize(delta(1))
delta=delta(1)
dVol = product(delta)
injectRadius = 3.5_dp*sqrt(3.0_dp)*delta(1)

! Place the star in the center of a cell. Gives best results
if (snap_to_grid) then
    do i=1,3
        loc(i) = floor(loc_in(i)) + 0.5_dp*delta(i)
    end do
end if

! count # of blocks which are at least partially within injectRadius, check that
! they are maximally refined
injBlkNum = 0
do blockID = 1, lnblocks
    if(nodetype(blockID) == LEAF) then
        blkCenterDist(:) = coord(:,blockID) - loc(:)
        blkCornerDist(:) = abs(blkCenterDist(:)) - 0.5_dp*bsize(:,blockID)

        if (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + &
            blkCornerDist(3)**2.0 < injectRadius**2.0) then

            injBlkNum = injBlkNum + 1

            iHaveInjectBlk = .true.
            call Grid_getBlkRefineLevel(blockID, refineLevel)
            
            iHaveUnrefined = (refineLevel < maxref)
            if ( iHaveUnrefined ) then
#ifdef DEBUG
                print *, "Block ", blockID, " on proc ", gr_meshMe, " is &
                    refined to level ", refineLevel, ", should be ", maxref
#endif DEBUG
                exit
            end if
        end if
    end if
end do

! if there are any not maximally refined, skip winds for now. Check across all
! procs
call MPI_ALLREDUCE(MPI_IN_PLACE, iHaveUnrefined, 1, MPI_LOGICAL, MPI_LOR, &
gr_meshComm, ierr)

if (iHaveUnrefined) then
    if (gr_meshMe == 0) then
        print *, "Not all wind injection blocks refined. Holding off until &
                    they are"
    end if
    return
end if

! at this point we know we are ready to inject wind

if (conserved_quant .eq. "momentum") then
    injectVelocityMax = injectVelocityIn / mass_load_factor
else if (conserved_quant .eq. "energy") then
    injectVelocityMax = injectVelocityIn / sqrt(mass_load_factor)
end if

if (rampVel .lt. injectVelocityMax) then
    if(injectVelocityMax - rampVel .lt. rampVelStep .or. isRestart) then
        rampVel = injectVelocityMax
    else
        rampVel = rampVel + rampVelStep
    end if
end if
injectVelocity = rampVel

if (conserved_quant .eq. "momentum") then
    injectMass = injectMassIn * injectVelocityIn / injectVelocity
else if (conserved_quant .eq. "energy") then
    injectMass = injectMassIn * (injectVelocityIn / injectVelocity)**2.0
end if


if (gr_meshMe == 0) then
    write(*,'(A,ES10.3,A,ES10.3,A)') "Injecting", injectMass/dt/solarMass*yr, &
                                 " solar masses/yr at ", injectVelocity, "cm/s"
end if

#ifdef DEBUG
print *, "Found", injBlkNum, "injection blocks on proc ", gr_meshMe
#endif

! build array of block IDs for all blocks to be injected. Get their center
! distances from the injection star. Then for each cell in each block, calculate
! its overlap with the injection sphere and store the value. Store the
! components of a velocity vector pointing radially outwards from the star
if (iHaveInjectBlk) then
    allocate(localInjectBlocks(injBlkNum))
    allocate(injectDataOverlap(injBlkNum,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI, &
                                                            GRID_KLO:GRID_KHI))
    allocate(injectDataVel(injBlkNum,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI, &
                                                        GRID_KLO:GRID_KHI,1:3))

#ifdef DEBUG
    print *, "Allocations done"
#endif

    n = 1
    do blockID = 1, lnblocks
        if(nodetype(blockID) == LEAF) then
            blkCenterDist(:) = coord(:,blockID) - loc(:)
            blkCornerDist(:) = abs(blkCenterDist(:)) - 0.5_dp*bsize(:,blockID)

            if (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + &
                blkCornerDist(3)**2.0 < injectRadius**2.0) then

                localInjectBlocks(n) = blockID
                n = n + 1

            end if
        end if
    end do

#ifdef DEBUG
    print *, "Found injection blocks:", localInjectBlocks, "on proc", gr_MeshMe
#endif
 
    do n = 1, injBlkNum
        blockID = localInjectBlocks(n)

        do i = GRID_ILO, GRID_IHI
            do j = GRID_JLO, GRID_JHI
                do k = GRID_KLO, GRID_KHI
                    ! since we have checked that all cells are refined, use
                    ! mindelta

                    x = (i - NGUARD - NXB/2.0 - 0.5)*delta(1) + coord(1,blockID)
                    y = (j - NGUARD - NYB/2.0 - 0.5)*delta(2) + coord(2,blockID)
                    z = (k - NGUARD - NZB/2.0 - 0.5)*delta(3) + coord(3,blockID)

                    dx = x - loc(1)
                    dy = y - loc(2)
                    dz = z - loc(3)

                    cell_bot = [ sign(abs(x) - 0.5*delta(1), x), &
                                 sign(abs(y) - 0.5*delta(2), y), &
                                 sign(abs(z) - 0.5*delta(3), z) ]

                    cell_top = [ sign(abs(x) + 0.5*delta(1), x), &
                                 sign(abs(y) + 0.5*delta(2), y), &
                                 sign(abs(z) + 0.5*delta(3), z) ]

                    rad = sqrt(dx**2.0 + dy**2.0 + dz**2.0)
                    
                    ! calculate the escape velocity at this point - JW
                    
                    v_esc = sqrt(2.0_dp * newton * star_mass / rad)
                    
                    ! add the escape velocity to the injection velocity - JW
                    
                    injectVelocity = sqrt(injectVelocity**2.0_dp + v_esc**2.0_dp)

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
                        
                    xvel = idir * injectVelocity
                    yvel = jdir * injectVelocity
                    zvel = kdir * injectVelocity
        
                    overlap_frac = 0.0

                    ! Now calculate the overlapping areas of the sphere and this cell.
                    call overlap(1, injectRadius, loc, cell_bot, &
                                    cell_top, 20, overlap_frac)

                    ! distribute mass differently based on whether star --> cell vector is
                    ! orthogonal or diagonal, to achieve spherical bubble

                    ! get the length of the longest component of the star --> cell vector,
                    ! normalized (i.e. values go from 1/sqrt(3) to 1)
                    maxdir = max(abs(idir),abs(jdir),abs(kdir))
                    ! if pointing diagonally, distribute extra mass here
                    !if (maxdir < .85_dp) then
                        ! linear relationship where maxdir = .65, multiplier = 1.8
                        !                           maxdir = 1.0, multiplier = 1.0
                    !  overlap_frac = overlap_frac * ((.85_dp - maxdir) * 3.4_dp + 1.0_dp)
                        ! if pointing orthogonally, distribute mass only if star is close to
                        ! a cell corner. Values decided empirically
                    !else if (maxdir > .96_dp) then
                        ! Find distance from star to the center of its own cell
                    !    dist = 0.0_dp
                    !    mult = 1.0
                    !    do d=1,3
                          ! distance in each direction
                    !      rem = abs(abs(mod(loc(d)/delta(d),1.0_dp)) - 0.5_dp)
                          !mult = mult + (.05 * rem / 0.5)
                    !      dist = dist + min( rem, 1.0_dp - rem )**2.0_dp
                    !      dist = dist + rem**2.0_dp
                    !    end do
                          ! total distance
                    !    dist = dist**0.5_dp
                          ! Give multiplier to overlap_frac from 1.0 if at center to 1.15 if
                          ! at corner
                    !    overlap_frac = overlap_frac * ( 1.0_dp + 0.15_dp * dist / &
                    !                                                (0.5_dp * sqrt(3.0_dp)))
                        !overlap_frac = overlap_frac * mult
                    !end if

                    !overlap_frac = overlap_frac * exp(18*(maxdir-1))+1

                    ! more simplistic method
                    !if (max(abs(idir),abs(jdir),abs(kdir)) < .65) then
                    !  overlap_frac = overlap_frac * 1.8
                    !else if (max(abs(idir),abs(jdir),abs(kdir)) < .8) then
                    !  overlap_frac = overlap_frac * 1.5
                    !end if

                    ! sum all weightings for normalization later

                    sumOverlap = sumOverlap + overlap_frac
                    injectDataOverlap(n,i,j,k) = overlap_frac
                    injectDataVel(n,i,j,k,1:3) = [xvel,yvel,zvel]
                end do
            end do
        end do
    end do
#ifdef DEBUG
    print *, "Calculated overlaps"
#endif
end if

#ifdef DEBUG
print*, "Proc ", gr_meshMe, " about to call MPI with sumOverlap = ", sumOverlap
#endif
call MPI_ALLREDUCE(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

#ifdef DEBUG
if (gr_meshMe == 0) then
    print *, "Total overlap weight is", sumOverlap
end if
#endif

if (iHaveInjectBlk) then
    if (conserved_quant .eq. "momentum") then
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkPtr(blockID, solndata)

            do i = GRID_ILO,GRID_IHI
                do j = GRID_JLO, GRID_JHI
                    do k = GRID_KLO, GRID_KHI
                        
                      dDens = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
                      oldDens = solndata(DENS_VAR,i,j,k)
                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
                      oldVel = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                      injVel = injectDataVel(n,i,j,k,1:3)

                      ! new velocity is calculated by conserving momentum when dDens moving
                      ! at injectVelocity is added to oldDens moving at original velocity

                      !solndata(VELX_VAR:VELZ_VAR, i, j, k) = &
                      !solndata(VELX_VAR:VELZ_VAR, i, j, k) * (oldDens/(oldDens+dDens)) &
                      !+ injectDataVel(n,i,j,k,1:3) * (dDens/(oldDens+dDens))
                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = &
                        injVel * dDens/newDens + oldVel * oldDens/newDens

                      solndata(DENS_VAR, i, j, k) = newDens
		      
                      if (with_thermal) then
		      
	                    KinE   = 0.5_dp*injVel**2.0_dp
			            ThermE = 0.5_dp*injectVelocityMax**2.0_dp - KinE
			            solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR, i, j, k) + ThermE
			
                      end if

                    end do
                end do
            end do

            call Grid_releaseBlkPtr(blockID, solndata)
        end do
    else if (conserved_quant .eq. "energy") then
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkPtr(blockID, solndata)

            do i = GRID_ILO,GRID_IHI
                do j = GRID_JLO, GRID_JHI
                    do k = GRID_KLO, GRID_KHI
                        
                      dDens = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
                      oldDens = solndata(DENS_VAR,i,j,k)
                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
                      oldVel = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                      injVel = injectDataVel(n,i,j,k,1:3)

                      ! Same as above, but conserving kinetic energy

                      newVelSq = oldDens/newDens * oldVel * abs(oldVel) &
                                + dDens/newDens * injVel * abs(injVel)

                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVelSq/sqrt(abs(newVelSq))

                      solndata(DENS_VAR, i, j, k) = newDens

                    end do
                end do
            end do

            call Grid_releaseBlkPtr(blockID, solndata)
        end do

    end if

    ! calculate min. crossing time of wind velocity on all injection blocks
    if (useTimeStep) then
        !min_wind_dt = 1d99
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkPtr(blockID, solndata)

            ! square of the max sound speed in the block
            !c_sound^2 = gamma * k * T / m_H

            cs2 = gamma_ * 1.3807d-16 * &
            solndata(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) &
            / 1.6726d-24 
            !v2 = solndata(VELX_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp + &
            !     solndata(VELY_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp + &
            !     solndata(VELZ_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp

            ! square of magnitude of velocities in all cells of block
            v2 = sum(solndata(VELX_VAR:VELZ_VAR,GRID_ILO:GRID_IHI, &
                                GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0,1)
            ! minimum dt to capture wind movement in this block
            !old_dt = minval(delta)/sqrt(cs2 + injectVelocity**2.0_dp)
            old_dt = minval(delta)/sqrt( maxval(cs2 + v2) )
            if (old_dt .lt. min_wind_dt) then
                min_wind_dt = old_dt
            end if

            call Grid_releaseBlkPtr(blockID, solndata)
        end do
    end if
end if

call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_MIN, gr_meshComm, ierr)

if (useTimeStep .and. gr_meshMe == 0) then
    write(*,'(A,ES10.3)') "Timestep set by inject_direct = ", min_wind_dt
end if

if (iHaveInjectBlk) then
    ! Update the EOS variables if we have an affected block for only those blocks.
    call Grid_notifySolnDataUpdate()

    do n=1, injBlkNum
        blockID = localInjectBlocks(n)
#ifdef DEBUG
        print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
                                                   blockID, gr_meshMe, injBlkNum
#endif

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
    end do

    deallocate(localInjectBlocks)
    deallocate(injectDataOverlap)
    deallocate(injectDataVel)
#ifdef DEBUG
    print *, "Deallocating done for proc", gr_meshMe
#endif
    !call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
#ifdef DEBUG
    print *, "Guard cells filled on proc", gr_meshMe
#endif
end if
call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
#ifdef DEBUG
print *, "Exiting inject_direct for proc", gr_meshMe
#endif
end subroutine inject_direct
