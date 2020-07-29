
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! subroutine inject_direct
!!!
!!! Authors: Joshua Wall and Andrew Pellegrino
!!!          Drexel University
!!!          Summer and Fall 2016
!!!
!!! A routine for injection of stellar winds (radiative type for massive
!!! stars).
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




subroutine inject_direct(loc_in, injectMassIn, injectVelocityIn, starMass, twind, dt, bgDens)

!#define DEBUG
#define DEBUG_ENERGY

use Grid_data, ONLY: gr_meshComm, gr_meshMe

use Hydro_data, ONLY: hy_cfl, hy_eswitch

use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getMinCellSize, &
    Grid_notifySolnDataUpdate, Grid_getBlkRefineLevel
    
use Driver_interface, ONLY : Driver_abortFlash
    
use Eos_interface, ONLY : Eos_wrapped

use pt_windInterface, only : overlap, sphere_and_cell_frac
use normal_rand

use Particles_windData, ONLY: min_wind_dt, ref_radius, mass_load, add_therm_e, &
                              wind_target_temp, var_radius, min_radius, &
                              perturb_velocity, perturb_std_dev, use_wind_compute_dt

use RuntimeParameters_interface, ONLY: RuntimeParameters_get

use tree, ONLY: nodetype, coord, bsize, lnblocks, refine, derefine, stay

#include "Flash.h"
#include "constants.h"

implicit none

#include "Flash_mpi.h"

integer, parameter :: dp = kind(1.d0)

real(dp), intent(in)    :: loc_in(3)
real(dp), intent(in)    :: injectMassIn, injectVelocityIn, twind, dt
real(dp), intent(inout) :: bgDens

logical, save :: first_call = .true.
logical :: iHaveInjectBlk
logical :: iHaveUnrefined

logical :: snap_to_grid
logical :: rampVelocity
logical, save :: isRestart
logical :: useTimeStep

real(dp) :: injectMass, injectVelocity, pertVelocity

real(dp) :: star_x, star_y, star_z, loc(3)

real(dp) :: cell_top(3), cell_bot(3)
real(dp) :: dVol
real(dp) :: x,y,z, dx, dy, dz

real(dp), save :: injectRadius, delta(3)
real(dp) :: overlap_frac,  sumOverlap, oldDens, dDens, newDens, solidAngle
real(dp) :: newVel(3), totP(3), totE, globalDeltaE, globalDeltaP, emech, pmech
real(dp) :: oldE, newE, oldP, newP
real(dp) :: oldVel(3), injVel(3), newVelSq(3)
real(dp) :: idir, jdir, kdir, rad, rad2, del2, xvel, yvel, zvel

integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), blockID
integer :: i, j, k, m, n
integer :: ind_array(3)


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
real(dp), save :: rampVelStep = 3.0d6 ! cm/s
real(dp) :: injectVelocityMax
character(len=10), save :: conserved_quant
integer :: injBlkNum
real(dp) :: blkCenterDist(3), blkRadDist , blkCornerDist(3)
real(dp), parameter :: yr = (60.0_dp**2.0)*24.0_dp*365.25_dp
real(dp), parameter :: solarMass = 1.989d33

real(dp) :: starMass, ThermE, KinE, injE, oldThermE
real(dp) :: globalTE, globalKE, largestTE, largestKE
real(dp) :: deltaKinE, deltaThermE, finalThermE, deltaE, addThermE, addKinE
real(dp) :: initialKE, initialTE
real(dp) :: delr ! distance between bottom of cell and either 
                 ! top of cell or max injectRadius, whichever is closer
real(dp) :: r1, r2 ! distance from star to top and bottom of cell
real(dp) :: sumMass ! Sum of all the mass array.

! Reference velocity to use for mass loading.
! Need a switch here for ionized or neutral.
!real(dp), parameter :: refVel = 1.86d7 ! For neutral gas to be = 10^6 K, but ours should be ionizied...
!real(dp), parameter :: refVel = 2.6919d7 ! For ionized gas to be = 10^6 K.
!real(dp), parameter :: refVel = 6.0d7 ! For ionized gas to be ~ 5x10^6 K.
real(dp), save :: refVel

! Here, save R_1 and the background denisty. On the *next* loop,
! compare R_1 to the current injection radius. If its larger, enlarge
! the injection radius up to the max we allow. This lets us start
! injecting in a region one dx in radius and expand as we go. - JW
real(dp) :: R_1

real, save :: injectRadiusMax
logical  :: calcBgDens, star_in_block


!integer  :: blkStar, iStar, jStar, kStar
!logical  :: hostCell

if (first_call) then
    call RuntimeParameters_get("gamma", gamma_)
    call RuntimeParameters_get("lrefine_max", maxref)
    call RuntimeParameters_get("restart", isRestart)
    call RuntimeParameters_get("cons_quant", conserved_quant)
    call RuntimeParameters_get("add_therm_e", add_therm_e)
    call RuntimeParameters_get("ref_radius", ref_radius)
    call RuntimeParameters_get("min_radius", min_radius)
    call RuntimeParameters_get("wind_target_temp", wind_target_temp)
    call RuntimeParameters_get("mass_load", mass_load)
    call RuntimeParameters_get("var_radius", var_radius)
    call RuntimeParameters_get("perturb_velocity", perturb_velocity)
    call RuntimeParameters_get("perturb_std_dev", perturb_std_dev)
    call RuntimeParameters_get("use_wind_compute_dt", use_wind_compute_dt)
    call Grid_getMinCellSize(delta(1))
    delta=delta(1)
    
    ! If we are setting ref_radius to -1, this means we want to
    ! vary the injection radius up to 
    if (ref_radius == -1.0) then
      ref_radius = 3.5_dp*sqrt(3.0_dp)*delta(1)
    end if
    
    injectRadius = ref_radius
    injectRadiusMax = injectRadius
      
    if (gr_meshMe == 0) then  
        print*, " Max injection radius / dx = ", injectRadiusMax/delta(1)
        print *, "inject_direct.F90 is conserving ", conserved_quant
    end if
    
    if (mass_load) then
        ! Draine 2011 eqn 36.28. Note this gas is always ionized by
        ! radiation first. - JW
        refVel = sqrt(wind_target_temp/1.38d7)*1e8
        if (gr_meshMe == 0) print*, "[inject_direct]: Reference velocity for wind is", &
                                     refVel, "for reference temp", wind_target_temp
     end if
    first_call = .false.
    
end if


! Mechanical Energy injected by the wind. - JW
injE = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp

!write(*,'(A,X,ES13.3e3,X,A)') "Injected specific E =", injE, "ergs."

iHaveUnrefined = .false.
iHaveInjectBlk = .false.
useTimeStep = .true.
rampVelocity = .false. !.true.
sumOverlap = 0.0_dp
sumMass = 0.0_dp

calcBgDens   = .false.
mass_load_factor = 0.0d0

snap_to_grid = .false. !.true.
globalDeltaE = 0.0
globalDeltaP = 0.0
globalTE     = 0.0
globalKE     = 0.0
largestKE    = 0.0
largestTE    = 0.0
deltaThermE  = 0.0
deltaKinE    = 0.0

injectMass = injectMassIn
injectVelocity = injectVelocityIn

call Grid_getMinCellSize(delta(1))
delta=delta(1)
del2 =delta(1)**2.
dVol = product(delta)

! If we've never calculated the background density for this region,
! just use the largest injection radius to be safe. Note that if its
! too large, the gas is denser and the injected amount will be more
! spread out, so the effect should be minor.

if (bgDens == 0.0 .and. var_radius) then 
    calcBgDens = .true.
#ifdef DEBUG_ENERGY
    if (gr_meshMe == 0) print*, "Calculating background density for this star for the first time."
#endif
end if

injectRadius = injectRadiusMax

if (var_radius) then

    if (calcBgDens) then
        ! If this is the first loop for this wind ever, calculate the
        ! average background density over the max injection range.
        injectRadius = injectRadiusMax !minval(delta)
        
        R_1 = injectRadius
    
        ! If we're going to recalculate density, zero it out here.
        bgDens = 0.0
    else

    ! Here we check if R_1 is inside of injectRadius. If it is,
    ! we recalculate overlap with this new radius. - JW

    ! Weaver et al 77 eqn 12
        R_1 = .74296_dp*(1.0_dp / bgDens * injectMassIn / dt)**(3.0_dp/10.0_dp) &
            * injectVelocityIn**(1.0_dp/10.0_dp) * twind**(2.0_dp/5.0_dp)
            
        !print*, R_1
        
!        if (R_1 < (injectRadius-2.0*minval(delta))) then !
!            ! Step down by two dx
!            injectRadius = max(R_1, minval(delta), (injectRadius-2.0*minval(delta)))
            
!        else if (R_1 > (injectRadius+minval(delta))) then !
!            ! Step up by one dx
!            injectRadius = min(injectRadiusMax, R_1, (injectRadius+minval(delta)))
            
!        end if

! If we are just calculating R_1 based on the original amb density, then just
! set the injection radius = R_1 each step (should only increase here...).
        if (R_1 < injectRadius) then !
            ! Step down to R_1
            injectRadius = max(R_1, sqrt(3.0d0)*minval(delta)) !, (injectRadius-2.0*minval(delta)))
            
        else if (R_1 > injectRadius) then !
            ! Step up to R_1
            injectRadius = min(injectRadiusMax, R_1) !, (injectRadius+minval(delta)))
            
        end if
        
    end if


    if (gr_meshMe == 0) then
#ifdef DEBUG_ENERGY
        write(*, '(A,ES13.3e3)') "R_1 from Weaver 77 = ", R_1
        write(*, '(A,ES13.3e3)') "# cells from free expansion edge R_1 from Weaver 77 = ", R_1/delta(1)
#endif
    end if

else

    injectRadius = injectRadiusMax

end if ! var_radius

! Always at least one cell in size.
injectRadius = max(injectRadius, sqrt(3.0d0)*minval(delta))

#ifdef DEBUG_ENERGY
if (gr_meshMe == 0) &
    write(*, '(A,ES13.3e3)') " injection radius / dx = ", injectRadius/delta(1)
#endif

! Place the star in the center of a cell. Gives best results
if (snap_to_grid) then
    do i=1,3
        loc(i) = floor(loc_in(i)) + 0.5_dp*delta(i)
    end do

else

    loc(:) = loc_in(:)    

end if

#ifdef DEBUG
if (gr_meshMe == 0) then
    print *, "inject_direct.F90 is conserving ", conserved_quant
    write(*, '(A,ES13.3e3, I4)') " injection radius = ", injectRadius, gr_meshMe
    print*, "loc =", loc, gr_meshMe
end if
#endif

! count # of blocks which are at least partially within injectRadius, check that
! they are maximally refined
injBlkNum = 0
blkCenterDist = 0.0d0
do blockID = 1, lnblocks
    if(nodetype(blockID) == LEAF) then
        star_in_block = .false.
        blkCenterDist(:) = abs(coord(:,blockID) - loc(:))
!        blkCornerDist    = blkCenterDist(:)-abs(0.5*bsize(:,blockID))
!        print*, "Hi ", (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + &
!            blkCornerDist(3)**2.0),  injectRadiusMax**2.0, gr_meshMe
            
!        if (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + & 
!            blkCornerDist(3)**2.0 <= (injectRadiusMax)**2.0) then ! + sqrt(3.0_dp)*delta(1)
! If all three directions of the block have been pierced by the injection sphere
! then this is an injection block. Note the previous method missed some cases
! where the sphere was totally contained in the block, such as:
! blk_low_left_corner = (0.0, 0.0, 0.0)
! blk_cent            = (1.0, 1.0, 1.0)
! blk_size            = (2.0, 2.0, 2.0)
! star_loc            = (0.5, 0.5, 0.0)
! inj_rad             = 0.1
!
! sum(blkCornerDist**2.0) = 0.5
! inj_rad**2.0            = 0.01
! sum(blkCornerDist**2.0) < inj_rad**2.0 = False in this case! - JW
        if (blkCenterDist(1) .le. 0.5*bsize(1,blockID) .and. &
            blkCenterDist(2) .le. 0.5*bsize(2,blockID) .and. &
            blkCenterDist(3) .le. 0.5*bsize(3,blockID)) then
            ! Then the star is already in this block.
            star_in_block = .true.
        else
            blkCornerDist    = blkCenterDist(:)-abs(0.5*bsize(:,blockID))        

            if ((blkCornerDist(1) .le. injectRadiusMax) .and. &
                 (blkCornerDist(2) .le. injectRadiusMax) .and. &
                 (blkCornerDist(3) .le. injectRadiusMax)) then
                ! Then the star is in a different block but the sphere hits this block.
                star_in_block = .true.
            end if
        end if
                
        if (star_in_block) then
            injBlkNum = injBlkNum + 1
            
#ifdef DEBUG_MPI
            print*, "injBlkNum =", injBlkNum
#endif

            iHaveInjectBlk = .true.
            call Grid_getBlkRefineLevel(blockID, refineLevel)
            
            if (refineLevel < maxref) then
                iHaveUnrefined = .true.
                ! Also set this block to be refined! - JW
                refine(blockID) = .true.
                derefine(blockID) = .false.
                stay(blockID) = .true.
#ifdef DEBUG_MPI
                print *, "Block ", blockID, " on proc ", gr_meshMe, " is &
                    refined to level ", refineLevel, ", should be ", maxref
#endif
            end if
        end if
    end if
end do

! if there are any not maximally refined, skip winds for now. Check across all
! procs
call MPI_ALLREDUCE(MPI_IN_PLACE, iHaveUnrefined, 1, MPI_LOGICAL, MPI_LOR, &
gr_meshComm, ierr)

if (iHaveUnrefined) then !  .and. ((.not. calcBgDens) .and. var_radius)
    if (gr_meshMe == 0) then
        print *, "Not all wind injection blocks refined. Holding off until &
                    they are"
    end if
    return
end if

! at this point we know we are ready to inject wind

!if (conserved_quant .eq. "momentum") then
!    injectVelocityMax = injectVelocityIn / mass_load_factor
!else if (conserved_quant .eq. "energy") then
!    injectVelocityMax = injectVelocityIn / sqrt(mass_load_factor)
!end if

! Lets make this simpler. Here we take the mass_load_factor, which
! is the ratio between the input velocity and the velocity to get 10^6 K
! shocked gas, and use it to set the mass and velocity to get the same
! input energy. - JW

if (mass_load) then  
    if (conserved_quant .eq. "momentum") then
        mass_load_factor = injectVelocityIn/refVel - 1.0d0
        injectVelocity   = injectVelocityIn / (1.0d0+mass_load_factor)
        injectVelocityMax= injectVelocity
        injectMass       = injectMassIn*(1.0d0+mass_load_factor)
    else if (conserved_quant .eq. "energy") then
        mass_load_factor = injectVelocityIn**2.0d0/refVel**2.0d0 - 1.0d0
        injectVelocity   = injectVelocityIn / dsqrt(1.0d0+mass_load_factor)
        injectVelocityMax= injectVelocity
        injectMass       = injectMassIn*(1.0d0+mass_load_factor)
    endif
end if

if (perturb_velocity) then
    ! all procs compute to keep random stream synchronized
    pertVelocity = norm_rand(injectVelocity, injectVelocity*perturb_std_dev)
end if

#ifdef DEBUG_ENERGY
if (gr_meshMe == 0) then
    write(*,'(A,ES10.3)') "Mass load factor = ", mass_load_factor
    write(*,'(A,ES10.3,A,ES10.3,A)') "Injecting", injectMassIn/dt/solarMass*yr, &
                                 " solar masses/yr at ", injectVelocityIn, " cm/s"
    !write(*,'(A,ES10.3,A,ES10.3,A)') "(", injectMass, " g over ", dt, " seconds)"
    write(*,'(A,ES10.3)') "Star mass =", starMass
    write(*,'(A,3ES10.3)') "Star loc  =", loc
    write(*,'(A,ES10.3,A)') "For a total of ", &
            0.5 * injectMassIn * injectVelocityIn**2/dt, " ergs/s"
end if
#endif


! build array of block IDs for all blocks to be injected. Get their center
! distances from the injection star. Then for each cell in each block, calculate
! its overlap with the injection sphere and store the value. Store the
! components of a velocity vector pointing radially outwards from the star
if (iHaveInjectBlk) then

#ifdef DEBUG_MPI
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

#ifdef DEBUG_MPI
    print *, "Allocations done"
#endif

    n = 1
    do blockID = 1, lnblocks
        if(nodetype(blockID) == LEAF) then
            star_in_block = .false.
            blkCenterDist(:) = abs(coord(:,blockID) - loc(:))
    !        blkCornerDist    = blkCenterDist(:)-abs(0.5*bsize(:,blockID))
    !        print*, "Hi ", (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + &
    !            blkCornerDist(3)**2.0),  injectRadiusMax**2.0, gr_meshMe
                
    !        if (blkCornerDist(1)**2.0 + blkCornerDist(2)**2.0 + & 
    !            blkCornerDist(3)**2.0 <= (injectRadiusMax)**2.0) then ! + sqrt(3.0_dp)*delta(1)
    ! If all three directions of the block have been pierced by the injection sphere
    ! then this is an injection block. Note the previous method missed some cases
    ! where the sphere was totally contained in the block, such as:
    ! blk_low_left_corner = (0.0, 0.0, 0.0)
    ! blk_cent            = (1.0, 1.0, 1.0)
    ! blk_size            = (2.0, 2.0, 2.0)
    ! star_loc            = (0.5, 0.5, 0.0)
    ! inj_rad             = 0.1
    !
    ! sum(blkCornerDist**2.0) = 0.5
    ! inj_rad**2.0            = 0.01
    ! sum(blkCornerDist**2.0) < inj_rad**2.0 = False in this case! - JW
            if (blkCenterDist(1) .le. 0.5*bsize(1,blockID) .and. &
                blkCenterDist(2) .le. 0.5*bsize(2,blockID) .and. &
                blkCenterDist(3) .le. 0.5*bsize(3,blockID)) then
                ! Then the star is already in this block.
                star_in_block = .true.
            else
                blkCornerDist    = blkCenterDist(:)-abs(0.5*bsize(:,blockID))        

                if ((blkCornerDist(1) .le. injectRadiusMax) .and. &
                     (blkCornerDist(2) .le. injectRadiusMax) .and. &
                     (blkCornerDist(3) .le. injectRadiusMax)) then
                    ! Then the star is in a different block but the sphere hits this block.
                    star_in_block = .true.
                end if
            end if
                    
            if (star_in_block) then
             
                localInjectBlocks(n) = blockID
                n = n + 1

            end if
        end if
    end do

#ifdef DEBUG_MPI
    print *, "Found injection blocks:", localInjectBlocks, "on proc", gr_MeshMe
#endif
 
    do n = 1, injBlkNum
        blockID = localInjectBlocks(n)
        call Grid_getBlkPtr(blockID, solndata)
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

                    rad2 = dx**2.0 + dy**2.0 + dz**2.0
                    rad  = sqrt(rad2)

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
                    if (perturb_velocity) then ! Perturb radial velocity?
                      xvel = idir * pertVelocity
                      yvel = jdir * pertVelocity
                      zvel = kdir * pertVelocity
                    else                        
                      xvel = idir * injectVelocity
                      yvel = jdir * injectVelocity
                      zvel = kdir * injectVelocity
                    end if
                    
                    if (abs(sqrt(xvel**2.0_dp+yvel**2.0_dp+zvel**2.0_dp) &
                        - injectVelocity*(1.0+5.0*perturb_std_dev)) &
                        /(injectVelocity*(1.0+5.0*perturb_std_dev))> 5.0d-1) then
                      write(*,'(A,2ES13.3E3)') "vels don't match: calc vel, inj vel= ",  &
                           sqrt(xvel**2.0_dp+yvel**2.0_dp+zvel**2.0_dp), injectVelocity
                      write(*,'(A,3ES13.3E3)') "xvel, yvel zvel =", xvel, yvel, zvel
                      stop
                    end if
        
                    overlap_frac = 0.0

                      ! Now calculate the overlapping areas of the sphere and this cell.
                      call overlap(1, injectRadius, loc, cell_bot, &
                                    cell_top, 20, overlap_frac)

#ifdef DEBUG
                    !write(*,*) "Calling sphere_and_cell_frac."
#endif
                    !call sphere_and_cell_frac(overlap_frac,injectRadius,dx,dy,dz,delta(1))
                    
                    ! Get the background density to estimate what the inject
                    ! radius should be physically. - JW
                    if (calcBgDens) &
                    bgDens = bgDens + overlap_frac*solndata(DENS_VAR, i, j, k)

                    if (.not. mass_load) then
                        ! Calculate the overlapping solid angle of a square at distance rad from the sphere.
                        ! Calculation from MPIA: http://www.mpia.de/~mathar/public/mathar20051002.pdf
                        solidAngle   = 4.*acos(sqrt((1.+del2/(2.*rad2))/(1.+del2/(2.*rad2)+(del2/rad2)**2.)))
                        overlap_frac = overlap_frac*solidAngle
                    end if
                    ! sum all weightings for normalization later
                    
                    if (rad .ge. min_radius) then
                        if (overlap_frac .gt. 0.0d0) then
                            sumOverlap = sumOverlap + overlap_frac
                            injectDataOverlap(n,i,j,k) = overlap_frac
                            injectDataVel(n,i,j,k,1:3) = [xvel,yvel,zvel]
                        end if
                    end if
                end do
            end do
        end do
        call Grid_releaseBlkPtr(blockID, solndata)
    end do
#ifdef DEBUG_MPI
    print *, "Calculated overlaps", gr_meshMe
#endif
end if

#ifdef DEBUG_MPI
print*, "Proc ", gr_meshMe, " about to call MPI with sumOverlap = ", sumOverlap
#endif
call MPI_ALLREDUCE(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

if (var_radius .and. calcBgDens) then

#ifdef DEBUG_MPI
    print*, "Before MPI background density is", bgDens
#endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, bgDens, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

    bgDens = bgDens / sumOverlap

end if

#ifdef DEBUG
    if (gr_meshMe == 0) then
        print *, "Total overlap weight is", sumOverlap
        print*, "Background density is", bgDens
    end if
#endif

!#ifdef DEBUG
!print*, "Proc ", gr_meshMe, " about to call MPI with sumMass = ", sumMass
!#endif
!call MPI_ALLREDUCE(MPI_IN_PLACE, sumMass, 1, MPI_DOUBLE_PRECISION, &
!                                            MPI_SUM, gr_meshComm, ierr)


if (iHaveInjectBlk) then
    if (conserved_quant .eq. "momentum") then
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
            call Grid_getBlkPtr(blockID, solndata)

            do i = GRID_ILO,GRID_IHI
                do j = GRID_JLO, GRID_JHI
                    do k = GRID_KLO, GRID_KHI
                      
                      ! Round off errors in calculations lead to small
                      ! changes in thermal energy that have big
                      ! consequences in EOS calls. So just skip.
                      if (injectDataOverlap(n,i,j,k) .le. 0.0_dp) cycle
                      
                      ! Note: This is includes the mass loaded stuff.
                      dDens   = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
                      oldDens = solndata(DENS_VAR,i,j,k)
                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
                      oldVel  = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                      injVel  = injectDataVel(n,i,j,k,1:3)
                      totP    = injVel * dDens + oldVel * oldDens
                      newVel  = totP / newDens
                      
                      ! How much mass are we really putting in? - JW
                      sumMass = sumMass + dDens*dVol
                      initialKE = 0.5_dp*sum(oldVel**2.0_dp) ! Per particle per unit mass, KE/M = KE / (N * mu * m_H).
                      initialTE = solndata(EINT_VAR,i,j,k)   ! Per particle per unit mass, TE/M = TE / (N * mu * m_H).

                      !oldE = oldDens * sum(oldVel**2)
                      oldE = (0.5_dp*sum(oldVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*oldDens
                      oldP = oldDens * sqrt(sum((oldVel+0.0_dp)**2.0_dp))

		      ! Note in Flash energy is in erg/g so that the energy is the same regardless of refinement. - JW
		      ! Also don't forget the already existing internal energy in this cell. - JW
		      
		      ! Specific KE added to this cell. - JW
!                      deltaKinE = (0.5_dp * newDens * sum(newVel**2.0_dp) &
!				      - 0.5_dp * oldDens * sum(oldVel**2.0_dp)) &
!				      / newDens
                      ! Didn't we already account for mass ratios etc when we inverted the
                      ! momentum eqn (i.e. isnt that the point of doing it this way)? - JW
                      !addKinE   = 0.5_dp * sum(injVel**2.0_dp)*dDens
                      deltaKinE = 0.5_dp * sum(newVel**2.0_dp)*newDens - 0.5_dp * sum(oldVel**2.0_dp)*oldDens
                      addKinE = deltaKinE
		      ! Thermal energy should be the difference b/t the injected specific E
		      ! and the actual specific kinetic energy. - JW

		      ! Final TE in this cell = portion of mech
		      ! energy added to the cell + the existing TE in cell.
              ! Note what we did here. In order to make up for the energy loss from conserving momentum
              ! we used the original velocity from injE to ensure that we added back enough energy
              ! in thermal to get back to the proper total energy at the end, and
              ! NOT the modified velocity from mass loading or what was calculated from momentum
              ! conservation. - JW
              if (add_therm_e) then
		      ! - KE added to cell. - JW !0.5_dp*dDens*sum(injVel**2.0_dp)/newDens
!		              finalThermE = injE*injectDataOverlap(n,i,j,k)/sumOverlap/(newDens*dVol)  &
!		                  + solndata(EINT_VAR, i, j, k)*oldDens/newDens &
!				          - deltaKinE
		      ! Note I think this is wrong, since the "correct" energy being added is the mass
              ! loaded energy, and the mass loading is the stuff falling in from the shell, which
              ! has zero kinetic energy. So we *should lose some kinetic energy here*!
                      !addThermE = finalThermE - solndata(EINT_VAR, i, j, k)
                      addThermE = (injE*injectDataOverlap(n,i,j,k)/sumOverlap/dVol  &
                                  - deltaKinE)
                      
              ! Lets add in the thermal energy of the mass loaded material to be the temperature
              ! the ionized shell material should be on the inner edge, which is 1e4 K.
              ! e = 1/(gamma-1)*k*T/mu - JW
                      !addThermE = (1.d0/(gamma_-1.0d0)*1.3807d-16*1e4/(1.67e-24*0.61))*dDens
              else
                      addThermE = 0.0d0
              end if
		      !ThermE = (0.5_dp * newDens * sum(newVel**2.0_dp) & 
		      !       + solndata(EINT_VAR, i, j, k) - oldE
                      ! Record the change in TE for debugging.
                      deltaThermE = addThermE
		      ! Now add the "missing" kinetic energy as thermal energy. - JW

                      !solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR,i,j,k) + deltaThermE
                      ! Shouldn't this have a factor that includes how much more mass we spread the old energy
                      ! around to by adding mass to the cell? In other words, this is the internal energy
                      ! per particle. So if we had 2 particles at the old TE and 1 particle at the TE of the
                      ! mass added to the cell, the new TE = (2*oldTE + 1*new_massTE) / 3 - JW
                      solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR,i,j,k)*(oldDens/newDens) + addThermE/newDens
                      !if (solndata(EINT_VAR,i,j,k) .lt. oldThermE) then 
                      !  print*, "We lost thermal energy in this cell. oldTE=, newTE=", &
                      !  oldThermE, solndata(EINT_VAR,i,j,k)
                      !end if

                      solndata(ENER_VAR, i, j, k) =  0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k)
                      if (solndata(EINT_VAR,i,j,k)/(0.5_dp*sum(newVel**2.0_dp)) < hy_eswitch) then
                        print*, "[inject_direct]: Warning, eintswitch should be kicking on here."
                        print*, "eint/ke=",solndata(EINT_VAR,i,j,k)/(0.5_dp*sum(newVel**2.0_dp))
                        print*, "eintSwitch=", hy_eswitch
                      end if

                      ! new velocity is calculated by conserving momentum when dDens moving
                      ! at injectVelocity is added to oldDens moving at original velocity

                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
                      solndata(DENS_VAR, i, j, k) = newDens

		      
		              if (addKinE .gt. largestKE) largestKE = addKinE*dVol
		              if (addThermE .gt. largestTE) largestTE = addThermE*dVol
		      
		      ! Delta E (not specific!) added to this cell. - JW
                      deltaE = (deltaKinE + deltaThermE)
                      
                      !newE = newDens * sum(newVel**2)
                      !newE = sum(newVel**2.0_dp)
                      newE = (0.5_dp * sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*newDens
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
                      write(*,'(A,ES13.3E3)') "Delta KE= ", addKinE
                      write(*,'(A,ES13.3E3)') "Delta TE= ", addThermE
                      write(*,'(A,ES13.3E3)') "dDens= ", dDens
                      write(*,'(A,2ES13.3E3)') "oldDens, newDens= ", oldDens, newDens
                      write(*,'(A,6ES13.3E3)') "oldVel, newVel= ", oldVel, newVel
                      write(*,'(A,3ES13.3E3)') "injVel =", injVel
                      write(*,'(A,2ES13.3E3)') "old KE, old TE= ", initialKE, initialTE
                      write(*,'(A,6ES13.3E3)') "new KE, new TE= ", 0.5_dp*sum(newVel**2.0_dp), solndata(EINT_VAR, i, j, k)
                      write(*,'(A,ES13.3E3)') "final TE= ", solndata(EINT_VAR, i, j, k)
                      write(*,'(A,ES13.3E3)') "initial Total E= ", oldE
                      write(*,'(A,ES13.3E3)') "final Total E= ", newE
                      write(*,'(A,ES13.3E3)') "deltaE= ", deltaE
                      write(*,'(A,ES13.3E3)') "newE-oldE= ", (newE-oldE)
                      write(*,'(A,ES13.3E3)') "deltaE - (newE-oldE)= ", deltaE - (newE-oldE)
                      write(*,'(A,ES13.3E3)') "globalTE= ", globalTE
                      write(*,'(A,ES13.3E3)') "globalKE= ", globalKE
                      write(*,'(A,ES13.3E3)') "globalDeltaE= ", globalDeltaE
                      call flush(6)
                    end if
#endif
		              globalKE     = globalKE + addKinE*dVol
                      globalTE     = globalTE + addThermE*dVol
                      globalDeltaE = globalDeltaE +  (addKinE + addThermE)*dVol !newE - oldE !deltaE
                      globalDeltaP = globalDeltaP + dDens * sqrt(sum((injVel)**2.0_dp))*dVol !(newP - oldP)

                    end do
                end do
            end do

            call Grid_releaseBlkPtr(blockID, solndata)
            call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
        end do
    else if (conserved_quant .eq. "energy") then
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
            call Grid_getBlkPtr(blockID, solndata)

            do i = GRID_ILO,GRID_IHI
                do j = GRID_JLO, GRID_JHI
                    do k = GRID_KLO, GRID_KHI
                        
                      dDens   = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
                      oldDens = solndata(DENS_VAR,i,j,k)
                      newDens = solndata(DENS_VAR,i,j,k) + dDens ! m + delta_m
                      oldVel  = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                      injVel  = injectDataVel(n,i,j,k,1:3)
                      oldE    = oldDens * sum(oldVel**2)
                      oldP    = oldDens * sqrt(sum(oldVel**2))

                      ! Same as above, but conserving kinetic energy

                      newVelSq = oldDens/newDens * sign(oldVel**2,oldVel) &
                               + dDens/newDens * sign(injVel**2,injVel)
                      newVel = sign(sqrt(abs(newVelSq)),newVelSq)

                      solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
                      solndata(DENS_VAR, i, j, k) = newDens
                      solndata(ENER_VAR, i, j, k) =  0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k)
                      
                      newE = newDens * sum(newVel**2)
                      newP = newDens * sqrt(sum(newVel**2))
                      !globalDeltaE = globalDeltaE + 0.5 * (newE - oldE)*dVol
                      !globalDeltaP = globalDeltaP + (newP - oldP)*dVol
                      
                      globalDeltaE = globalDeltaE + 0.5 * dDens *sum(injVel**2) *dVol
                      globalDeltaP = globalDeltaP + dDens * sqrt(sum((injVel)**2.0_dp))*dVol 
                      
                      ! How much mass are we really putting in? - JW
                      sumMass = sumMass + dDens*dVol
                    end do
                end do
            end do

            call Grid_releaseBlkPtr(blockID, solndata)
#ifdef DEBUG_MPI
            print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
                                                   blockID, gr_meshMe, injBlkNum
#endif
            call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
        end do

    ! energy conserving method from Simpson et al. 2015. Should be equivalent to
    ! normal energy conservation
    !else if (conserved_quant .eq. "energy_simpson") then
    !    do n = 1, injBlkNum
    !        blockID = localInjectBlocks(n)
    !        call Grid_getBlkPtr(blockID, solndata)

    !        do i = GRID_ILO,GRID_IHI
    !            do j = GRID_JLO, GRID_JHI
    !                do k = GRID_KLO, GRID_KHI
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

    end if

end if


#ifdef DEBUG_ENERGY
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalTE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestTE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestKE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    
#ifdef DEBUG_MPI
print*, "Proc ", gr_meshMe, " about to call MPI with sumMass = ", sumMass
#endif
call MPI_ALLREDUCE(MPI_IN_PLACE, sumMass, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_SUM, gr_meshComm, ierr)

    if (gr_meshMe == 0) then
        emech = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp
        pmech = injectMassIn * injectVelocityIn
        write(*,'(A,ES10.3,A)') "We got ", globalDeltaE/dt, "ergs/s"
        write(*,'(A,ES10.3,A)') "We got global TE=", globalTE, "ergs"
        write(*,'(A,ES10.3,A)') "We got global KE=", globalKE, "ergs"        
        write(*,'(A,F10.3,A)') "Error in injected wind total E is ", abs(globalDeltaE - emech)/emech*100, "%"
        write(*,'(A,F10.3,A)') "Percentage of energy that is thermal is ", globalTE/abs(globalDeltaE)*100.0, "%"
        write(*,'(A,F10.3,A)') "Percentage of energy that is kinetic is ", globalKE/abs(globalDeltaE)*100.0, "%"
        write(*,'(A,F10.3,A)') "Error in injected wind P is ", abs(globalDeltaP - pmech)/pmech*100, "%"
        write(*,'(A,ES10.3e2,A)') "Largest TE is ", largestTE
        write(*,'(A,ES10.3e2,A)') "Largest KE is ", largestKE
        write(*,'(A,ES10.3e2,X,A,ES10.3e2)') "Total mass from wind is ", injectMass, "Injected mass is ", sumMass
    end if
#endif

if (use_wind_compute_dt) then
  call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_MIN, gr_meshComm, ierr)
else
  min_wind_dt = 1d99
endif

call Grid_notifySolnDataUpdate()

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
    	if (use_wind_compute_dt) then

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
                
    ! NOTE NOTE NOTE ! We need to also calculate the Alfven velocity here
    ! or otherwise take it into account in Hydro_computeDt (<--- this should be done now. JW) - JW
    
	! minimum dt to capture wind movement in this block
    ! Lets use the hydro CFL here. - JW
	!old_dt = minval(delta)/sqrt(cs2 + injectVelocity**2.0_dp)
            old_dt = hy_cfl * minval(delta)/sqrt( maxval(cs2 + v2) )
            if (old_dt .lt. min_wind_dt) then
                min_wind_dt = old_dt
            end if
#ifdef DEBUG
            print*, "Max temp on proc, blk", gr_meshMe, blockID, "is", maxval(solndata(TEMP_VAR,:,:,:))
#endif
            call Grid_releaseBlkPtr(blockID, solndata)

        end if


    end do

    deallocate(localInjectBlocks)
    deallocate(injectDataOverlap)
    deallocate(injectDataVel)
#ifdef DEBUG_MPI
    print *, "Deallocating done for proc", gr_meshMe
#endif
    !call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
!#ifdef DEBUG
!    print *, "Guard cells filled on proc", gr_meshMe
!#endif
end if

                                            
if (use_wind_compute_dt) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, &
                                            MPI_MIN, gr_meshComm, ierr)
#ifdef DEBUG_ENERGY
    if (gr_meshMe == 0) write(*,'(A,ES10.3)') "Timestep set by inject_direct = ", min_wind_dt
#endif
end if

! Moving to outside the inject_direct call and into Particles_wind.
!call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)

#ifdef DEBUG_MPI
print *, "Exiting inject_direct for proc", gr_meshMe
#endif

end subroutine inject_direct
