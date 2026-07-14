
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

! Available compiler flags:
!   WIND_VERBOSE : Human-readable per-call/per-star diagnostics.
!   WIND_DEBUG   : Detailed per-cell/per-block diagnostics. Implies WIND_VERBOSE.
#ifdef WIND_DEBUG
#define WIND_VERBOSE
#endif

subroutine inject_direct(loc_in, injectMassIn, injectVelocityIn, injectYieldIn, twind, dt, bgDens)

!#define DEBUG
!#define DEBUG_ENERGY

#include "Flash.h"
#include "constants.h"

#ifdef WIND_DEBUG
    use iso_fortran_env, only : error_unit
#endif
    
    use Grid_data, ONLY: gr_globalNumProcs, gr_meshComm, gr_meshMe

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

#ifdef TRACER_FIELDS
    use Particles_windData, ONLY: mass_load_yields, ism_loading
#endif

    implicit none

#include "Flash_mpi.h"

    ! =========================================================================
    ! Precision Parameters & Constants
    ! =========================================================================
    integer, parameter  :: dp = kind(1.d0)          ! Double precision 
    real(dp), parameter :: yr = 3.1557600d7         ! Year
    real(dp), parameter :: solarMass = 1.989d33     ! Solar mass
    real(dp), parameter :: kB = 1.3807d-16          ! Boltzmann constant
    real(dp), parameter :: mH = 1.6726d-24          ! Proton mass
    real(dp), parameter :: mu = 1.3                 ! Mean molecular weight
    ! =========================================================================
    ! Persitent parameters
    ! =========================================================================
    logical, save           :: first_call = .true.
    real(dp), save          :: injectRadiusMax
    ! TODO: These runtime paramters should be imported from elsewhere.
    character(len=10), save :: conserved_quant
    integer, save           :: maxref
    real(dp), save          :: gamma_
    real(dp), save          :: eint_floor

    ! =========================================================================
    ! Subroutine Arguments
    ! =========================================================================
    real(dp), intent(in)    :: loc_in(3)        ! Injection origin
    real(dp), intent(in)    :: injectMassIn     ! Injected mass
    real(dp), intent(in)    :: injectVelocityIn ! Injected velocity
    real(dp), intent(in)    :: twind            ! Time since onset of wind
    real(dp), intent(in)    :: dt               ! Current timestep
    real(dp), intent(inout) :: bgDens           ! Background density

    ! =========================================================================
    ! Configuration & Control Flags
    ! =========================================================================
    logical :: iHaveInjectBlk                   ! Local MPI rank has injection block
    logical :: iHaveUnrefined                   ! Local MPI rank has unrefined injection block
    logical :: snap_to_grid                     ! Inject at center of grid
    logical :: calcBgDens                       ! Background has not been calculated for source

    ! =========================================================================
    ! Injected quantities
    ! =========================================================================
    real(dp) :: injectMass                      ! Mass
    real(dp) :: injectVelocity                  ! Velocity (magnitude)
    real(dp) :: injectEkin                      ! Kinetic energy

    ! =========================================================================
    ! Injection region and wind modification
    ! =========================================================================
    integer  :: ip                              ! Index used for velocity pertubation
    integer  :: nOverlap                        ! Local number of overlapping cells
    integer  :: nOverlapTot                     ! Total number of overlapping cells
    integer  :: nOverlapArr(gr_globalNumProcs)  ! Number of overlapping cells on each MPI rank
    real(dp) :: injectRadius                    ! Injection radius
    real(dp) :: refVel                          ! Target velocity if mass-loading
    real(dp) :: R_1                             ! Injection region from Weaver et al. (1977)
    real(dp) :: mass_load_factor                ! Mass loading factor
    real(dp) :: overlap_frac                    ! Fraction of cell which overlap with injection
    real(dp) :: sumOverlap                      ! Total overlap
    real(dp) :: weight                          ! Injection weight from overlap fraction
    real(dp) :: solidAngle                      ! Solid angle of cell face as seen by star
    real(dp) :: perturbEkin                     ! Kinetic energy of velocity pertubation
    real(dp) :: perturbNorm                     ! Normalization to limit energy added by pertubations.
    real(dp), allocatable, dimension(:) :: perturbScale ! Wind velocity perturbation

    ! =========================================================================
    ! Cell update
    ! =========================================================================
    real(dp) :: oldDens                         ! Density before injection
    real(dp) :: oldVel(3)                       ! Velocity before injection 
    real(dp) :: oldEkinDens                     ! Kinetic energy density before injection 
    real(dp) :: oldEtheDens                     ! Thermal energy density before injection 

    real(dp) :: injDens                         ! Density of injected material
    real(dp) :: injVel(3)                       ! Velocity of injected material
    real(dp) :: injEkinDens                     ! Kinetic energy density of injected material
    real(dp) :: injEtheDens                     ! Thermal energy density of injected material
    
    real(dp) :: newDens                         ! Density after injection
    real(dp) :: newVel(3)                       ! Velocity after injection
    real(dp) :: newVelSq(3)                     ! Square of velocity after injection 
    real(dp) :: newMom(3)                       ! Momentum density after injection
    real(dp) :: newEkinDens                     ! Kinetic energy density after injection
    real(dp) :: newEkin                         ! Kinetic energy after injection
    real(dp) :: newEthe                         ! Thermal energy after injection

    ! =========================================================================
    ! Timestep Calculation Parameters
    ! =========================================================================
    real(dp) :: cs2(NXB,NYB,NZB)                ! Sound speed squared 
    real(dp) :: v2(NXB,NYB,NZB)                 ! Gas velocity squared

    ! =========================================================================
    ! Geometry & Cell Coordinates
    ! =========================================================================
    real(dp) :: loc(3)
    real(dp) :: cell_top(3)
    real(dp) :: cell_bot(3)
    real(dp) :: delta(3)
    real(dp) :: dVol
    real(dp) :: x, y, z
    real(dp) :: dx, dy, dz
    real(dp) :: xcoll, ycoll, zcoll
    real(dp) :: d2coll
    real(dp) :: rad, rad2, del2
    real(dp) :: xvel, yvel, zvel

    ! =========================================================================
    ! Iteration Loop Indices
    ! =========================================================================
    integer :: i, j, k, n

    ! =========================================================================
    ! Grid Block Indexing & Structure
    ! =========================================================================
    integer :: blockID
    integer :: injBlkNum
    integer :: refineLevel
    integer :: blkLimits(2,MDIM)
    integer :: blkLimitsGC(2,MDIM)
    real(dp) :: blkCtr(3)
    real(dp) :: blkSize(3)

    ! =========================================================================
    ! Solver Buffers
    ! =========================================================================
    real(dp), pointer, dimension(:,:,:,:)       :: solndata
    real(dp), allocatable, dimension(:,:,:,:)   :: injectDataOverlap
    real(dp), allocatable, dimension(:,:,:,:,:) :: injectDataVel
    integer, allocatable, dimension(:)          :: localInjectBlocks

    ! =========================================================================
    ! Error handling & MPI
    ! =========================================================================
    logical :: perturbUsedFullChunk
    integer :: ierr

    ! =========================================================================
    ! Tracer Fields (if enabled)
    ! =========================================================================
#ifdef TRACER_FIELDS
    real(dp), intent(in)    :: injectYieldIn(NMASS_SCALARS)
    real(dp)                :: injectYield(NMASS_SCALARS)
    real(dp)                :: oldTracerField(NMASS_SCALARS)
    real(dp)                :: injTracerField(NMASS_SCALARS)
    real(dp)                :: newTracerField(NMASS_SCALARS)
    real(dp)                :: ism_mass
    integer                 :: itracer
#else
    ! If TRACER_FIELDS is disabled, injectYield is declared for signature consistency
    real(dp), intent(in)    :: injectYieldIn
    real(dp)                :: injectYield
#endif
    
#ifdef WIND_VERBOSE
    ! =========================================================================
    ! Diagnostics & Debug Statistics
    ! =========================================================================
    real(dp) :: emech                           ! Total mechanical energy
    real(dp) :: pmech                           ! Total momentum
    real(dp) :: globalTE                        ! Total thermal energy 
    real(dp) :: globalKE                        ! Total kinetic energy
    real(dp) :: globalDeltaE                    ! Total change in energy
    real(dp) :: globalDeltaP                    ! Total change in momentum
    real(dp) :: largestTE                       ! Largest thermal energy injected in cell
    real(dp) :: largestKE                       ! Largest kinetic energy injected in cell
    real(dp) :: sumMass                         ! Total injected mass
    real(dp) :: largestDelta                    ! Largest pertubation
    real(dp) :: largestDeltaVel                 ! Largest velocity pertubation
    real(dp) :: sumNegativeTherm = 0.0_dp       ! Check how much energy was added.
    integer :: nNegativeTherm = 0               ! Check if thermal energy was limited
#endif

!    ! =========================================================================
!    ! Unused / Obsolete Variables (Preserved for reference & investigation)
!    ! =========================================================================
!    logical :: rampVelocity = .false.
!    logical :: useTimeStep = .true.
!    logical, save :: isRestart
!    real(dp) :: injectVelocityMax
!    real(dp), save :: rampVel = 0d0          ! Leftover velocity ramping state
!    real(dp), save :: rampVelStep = 3.0d6    ! Leftover velocity ramping step parameter
!    real(dp) :: star_x, star_y, star_z       ! Intended for star coordinates; replaced by loc(3)
!    real(dp) :: totE                         ! Obsolete total energy; now using separate oldE/newE
!    real(dp) :: ThermE, KinE                 ! Obsolete cell energy components
!    real(dp) :: oldThermE                    ! Obsolete baseline thermal energy
!    real(dp) :: finalThermE                  ! Commented out energy calculation buffer
!    real(dp) :: delr                         ! Unused cell coordinate delta (distance calculation)
!    real(dp) :: r1, r2                       ! Unused spherical/cell radial coordinates
!    real(dp) :: maxdir                       ! Unused vector coordinate helper
!    real(dp) :: rem, dist, mult              ! Unused loop/interpolation variables
!    integer :: d                             ! Unused dimension index
!    integer :: m                             ! Unused loop index
!    integer :: ind_array(3)                  ! Unused index vector buffer
!    integer :: messages                      ! Unused MPI message counter
!    integer :: trans_status(MPI_STATUS_SIZE) ! Unused MPI status tracker
!    real(dp) :: old_dt
!    real(dp) :: oldE                         ! From old injection 
!    real(dp) :: idir, jdir, kdir <--- No longer in use.
!    real(dp) :: newE
!    real(dp) :: oldP
!    real(dp) :: newP
!    real(dp) :: initialKE
!    real(dp) :: initialTE
!    real(dp) :: addKinE
!    real(dp) :: deltaKinE = 0.0_dp
!    real(dp) :: deltaThermE = 0.0_dp
!    real(dp) :: deltaE

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "Entered subroutine."
#endif

    ! =========================================================================
    ! Initialize variables, and perform initial checks
    ! =========================================================================
    
    ! First, check that the input mass and velocity are greater than 0.
    ! Do we want to allow 0 velocity winds?
    if ((injectMassIn .le. tiny(0.0_dp)) .or. (injectVelocityIn .le. tiny(0.0_dp))) then
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "Wind mass or velocity is zero, returning."
#endif
        return
    end if

    ! Store input in separate variables. 
    loc(:) = loc_in(:)
    injectMass = injectMassIn
    injectVelocity = injectVelocityIn
    injectEkin = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp
#ifdef TRACER_FIELDS
    injectYield = injectYieldIn
#endif
    
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) then
            write(*, 900) "Injection properties"
            write(*, 921) "injectMass = ", injectMass
            write(*, 921) "injectVelocity = ", injectVelocity
            write(*, 921) "twind = ", twind
            write(*, 921) "dt = ", dt
            write(*, 921) "bgDens = ", bgDens
#ifdef TRACER_FIELDS
            do itracer = 1, NMASS_SCALARS
                write(*, 921) "injectYield(i) = ", injectYield(itracer)
            enddo
#endif
        endif
#endif

    ! Get the minimum cell size and set identical in all directions.
    ! Also calculate face area and volume of cells, which are used repeatedly.
    call Grid_getMinCellSize(delta(1))
    delta=delta(1)
    del2 =delta(1)**2.
    dVol = product(delta)
    
    ! Read in and set all parameters if this is the first time inject_direct is being called.
    if (first_call) then

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "First call, loading runtime parameters"
#endif

        call RuntimeParameters_get("gamma", gamma_)
        call RuntimeParameters_get("lrefine_max", maxref)
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
#ifdef TRACER_FIELDS
        call RuntimeParameters_get("mass_load_yields", mass_load_yields)
        call RuntimeParameters_get("ism_loading", ism_loading)
#endif

        if (ref_radius == -1.0) then
            ! If ref_radius is -1, use a default injection radius equal to the distance 
            ! from the center of the central cell to a corner of a 7x7x7 cube of cells.
            injectRadiusMax = 3.5_dp*sqrt(3.0_dp)*delta(1)
        else
            if (gr_meshMe == 0 .and. ref_radius < (sqrt(3.0_dp)*minval(delta))) then
                write(*, 900) "VALUE ERROR: Injection radius is being set to less than the cell size."
                write(*, 922) "ref_radius, sqrt(3)*delta = ", ref_radius, sqrt(3.0_dp)*minval(delta)
                stop
            end if
            injectRadiusMax = max(ref_radius, (sqrt(3.0_dp)*minval(delta)))
        end if
        
        ! Internal energy floor. 10 K for now, add this as a parameter.
        eint_floor = kB * 10.0_dp / ((gamma_ - 1.0_dp) * mu * mH)
    
!        if (gr_meshMe == 0) then
!            print*, " Max injection radius / dx = ", injectRadiusMax/delta(1)
!            print *, "inject_direct.F90 is conserving ", conserved_quant
!        end if

        first_call = .false.
    end if

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) then
            write(*, 900) "Parameters being used:"
            write(*, 921) "gamma = ", gamma_
            write(*, 911) "lrefine_max = ", maxref
            write(*, 941) "cons_quant = ", conserved_quant
            write(*, 931) "add_therm_e = ", add_therm_e
            write(*, 921) "ref_radius = ", ref_radius
            write(*, 921) "min_radius = ", min_radius
            write(*, 921) "wind_target_temp = ", wind_target_temp
            write(*, 931) "mass_load = ", mass_load
            write(*, 931) "var_radius = ", var_radius
            write(*, 931) "perturb_velocity = ", perturb_velocity
            write(*, 921) "perturb_std_dev = ", perturb_std_dev
            write(*, 931) "use_wind_compute_dt = ", use_wind_compute_dt
#ifdef TRACER_FIELDS
            write(*, 931) "mass_load_yields = ", mass_load_yields
            write(*, 931) "ism_loading = ", ism_loading
#endif
        endif
#endif

    ! =========================================================================
    ! Check injection region.
    ! =========================================================================
#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Setting up injection region"
#endif

    ! Place the star in the center of a cell for cleaner injection.
    ! This is hardcoded to .false. - we should consider adding a parameter.
    snap_to_grid = .false.
    if (snap_to_grid) then
        do i=1,3
            loc(i) = floor(loc_in(i)) + 0.5_dp*delta(i)
        end do
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) then
            write(*, 900) "Injection moved to center of grid."
            write(*, 923) "(loc - loc_in)", loc(1)-loc_in(1), loc(2)-loc_in(2), loc(3)-loc_in(3)
        endif
#endif
    end if

    ! First pass over local leaf blocks: count blocks whose bounding boxes overlap 
    ! injectRadiusMax and flag any that still require refinement.
    ! Notes 
    !    - iHaveInjectBlk is local to this MPI rank; it is not globally reduced.
    !    - The refinement check uses injectRadiusMax, not the current injectRadius.
    !      This conservatively refines all blocks that could ever be used by the
    !      wind injection region. If var_radius makes injectRadius smaller, some
    !      refined blocks may lie outside the current injection sphere, but they are
    !      kept refined so the radius can expand later without injecting into
    !      under-refined blocks.

#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Counting injection blocks and checking refinement"
#endif

    injBlkNum = 0
    iHaveInjectBlk = .false.
    iHaveUnrefined = .false.
    do blockID = 1, lnblocks
        if(nodetype(blockID) == LEAF) then
            ! Determine if injection region overlap with leaf block bounding box owned 
            ! by this MPI rank.
            call Grid_getBlkCenterCoords(blockID,blkCtr)
            call Grid_getBlkPhysicalSize(blockID,blkSize)
            
            ! The algorithm checks for overlap by testing if the spherical region intersects 
            ! the block AABB (axis-aligned bounding box). The closest point on the AABB to 
            ! the sphere center is found by clamping each coordinate to the block boundaries. 
            ! If the distance from this point to the sphere center is less than the sphere
            ! radius, then the sphere and AABB overlap.

            ! Coordinates for box which is closest to the injection sphere center.
            xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
            ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
            zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))

            if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadiusMax**2) then
                ! This block is overlapping with the injection region. Add to count and flag
                ! this MPI rank as one that has injection in it.
                iHaveInjectBlk = .true.
                injBlkNum = injBlkNum + 1

!#ifdef DEBUG_MPI
!                print*, "injBlkNum =", injBlkNum
!#endif
                ! Check if this block is maximally refined. If not, it will be flagged for 
                ! refinement and injection will be aborted.
                call Grid_getBlkRefineLevel(blockID, refineLevel)
                
                if (refineLevel < maxref) then
                    iHaveUnrefined = .true.
                    refine(blockID) = .true.
                    derefine(blockID) = .false.
                    stay(blockID) = .true.

#ifdef WIND_DEBUG
                    write(*, 913) "Unrefined block for rank, blockID, refineLevel", gr_meshMe, blockID, refineLevel
#endif
!#ifdef DEBUG_MPI
!                    print *, "Block ", blockID, " on proc ", gr_meshMe, " is &
!                        refined to level ", refineLevel, ", should be ", maxref
!#endif
                end if
            end if
        end if
    end do
    
!#ifdef WIND_VERBOSE
!    if (iHaveUnrefined) then
!        write(*, 912) "WARNING: unrefined injection blocks, gr_meshMe, refineLevel = ", gr_meshMe, refineLevel
!    end if
!#endif

    ! Check if any MPI rank has unrefined blocks in injection region. If so, return 
    ! for this timestep. These blocks have been flagged for refinement and should 
    ! eventually reach maximum refinement unless this is blocked elsewhere.
    call MPI_ALLREDUCE(MPI_IN_PLACE, iHaveUnrefined, 1, MPI_LOGICAL, MPI_LOR, &
        gr_meshComm, ierr)

#ifdef WIND_DEBUG
    if (iHaveInjectBlk .and. injBlkNum <= 0) then
        write(*, 921) "Rank has iHaveInjectBlk but injBlkNum <= 0, injBlkNum =", injBlkNum
    end if
#endif

    if (iHaveUnrefined) then
        if (gr_meshMe == 0) write(*, 900) "Unrefined blocks, will exit injection step."
        return
    end if

    ! =========================================================================
    ! Injection region is refined, derive wind properties.
    ! =========================================================================
#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Deriving wind properties"
#endif

    ! Start with max radius as default for each time step.
    injectRadius = injectRadiusMax

!    ! Calculate mechanical energy injected by the wind.
!    ! I'm moving this to later for clarity
!    injEkinDens = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp
!    write(*,'(A,X,ES13.3e3,X,A)') "Injected specific E =", injEkinDens, "ergs."
    
    calcBgDens = .false.
    if (var_radius) then
        ! If using a variable injection radius, inject only within the freely
        ! expanding wind region interior to the first shock (R_1; Weaver et al. 1977).
        ! Computing R_1 requires the ambient density, which is measured only during
        ! the first injection step and assumed to remain constant thereafter.
        !
        ! Since the background density is not yet known during the first step, the
        ! full injection radius (injectRadiusMax) is used initially.
        !
        ! Notes:
        !   - bgDens is only used when var_radius == .true.
        !
        ! Legacy comments (to be removed):
        !   If we've never calculated the background density for this region,
        !   just use the largest injection radius to be safe. Note that if its
        !   too large, the gas is denser and the injected amount will be more
        !   spread out, so the effect should be minor.

        if (bgDens == 0.0) then
!#ifdef DEBUG_ENERGY
!            if (gr_meshMe == 0) print*, "Calculating background density for this star for the first time."
!#endif

            ! During the first injection step, the background density must be 
            ! calculated. Note that for this first step the injection radius 
            ! will always be injectRadiusMax.
            calcBgDens = .true.    

#ifdef WIND_VERBOSE
            if (gr_meshMe == 0) then
                write(*, 900) "First injection, will calculate background density"
                write(*, 921) "Instead using max injectRadius = ", injectRadius
            endif
#endif
        else
            ! Adapt the injection radius to the radius of the first shock, R_1,
            ! separating the freely expanding supersonic wind region from the
            ! shocked wind region (Weaver et al. 1977, eqn 12). Note that we
            ! assume alpha = 0.88.
            R_1 = 0.74296_dp * (injectMassIn / (bgDens * dt))**0.3_dp &
                * injectVelocityIn**0.1_dp &
                * twind**0.4_dp
            
            ! We limit the injection radius between the diagonal of the cell and 
            ! maximum allowed injection radius.            
            injectRadius = min(injectRadiusMax, max(R_1, sqrt(3.0_dp) * minval(delta)))
#ifdef WIND_VERBOSE
            if (gr_meshMe == 0) then
                write(*, 921) "Updating injectRadius, calculated W77 R1 = ", R_1
                write(*, 921) "Injection limited to injectRadius = ", injectRadius
            endif
#endif
        end if

!        if (gr_meshMe == 0) then
!#ifdef DEBUG_ENERGY
!            write(*, '(A,ES13.3e3)') "R_1 from Weaver 77 = ", R_1
!            write(*, '(A,ES13.3e3)') "# cells from free expansion edge R_1 from Weaver 77 = ", R_1/delta(1)
!#endif
!        end if

    end if
    
!#ifdef DEBUG_ENERGY
!    if (gr_meshMe == 0) &
!        write(*, '(A,ES13.3e3)') " injection radius / dx = ", injectRadius/delta(1)
!#endif


!#ifdef DEBUG
!    if (gr_meshMe == 0) then
!        print *, "inject_direct.F90 is conserving ", conserved_quant
!        write(*, '(A,ES13.3e3, I4)') " injection radius = ", injectRadius, gr_meshMe
!        print*, "loc =", loc, gr_meshMe
!    end if
!#endif

    if (mass_load) then
        ! Mass loading artificially increases the injected mass to lower the wind 
        ! velocity to the reference velocity refVel. The mass loading factor is 
        ! calculated from either input momentum or kinetic energy, depending on 
        ! which is conserved (set by cons_quant parameter).
        ! Notes:
        !   - Assumes injectVelocityIn >= refVel so that mass_load_factor >= 0.
        !     For very slow winds (~300 km/s at 1e6 K), this can results in
        !     mass unloading. This becomes increasingly problematic for higher 
        !     temperatures (e.g. ~850km/s at 1e7 K). 
            
        ! Calculate reference velocity from target temperature (Draine, 2011, eqn 36.28).
        ! Note that this assumes that the gas is ionized.
        refVel = sqrt(wind_target_temp/1.38d7)*1e8
!        if (gr_meshMe == 0) print*, "[inject_direct]: Reference velocity for wind is", &
!            refVel, "for reference temp", wind_target_temp

        mass_load_factor = 0.0_dp
        if (conserved_quant .eq. "momentum") then
            mass_load_factor = injectVelocityIn / refVel - 1.0_dp
        else if (conserved_quant .eq. "energy") then
            mass_load_factor = (injectVelocityIn / refVel)**2.0_dp - 1.0_dp
        end if

        injectVelocity    = refVel
        ! injectVelocityMax = injectVelocity    <--- This is never used
        injectMass        = injectMassIn * (1.0_dp + mass_load_factor)

#ifdef TRACER_FIELDS
        ! If mass_load_yields is enabled, decide how to assign tracer mass to the
        ! added mass. There are two interpretations:
        !
        !   ism_loading = .false.
        !       The added mass is treated as wind material. The wind yield is scaled
        !       with the total loaded wind mass.
        !
        !   ism_loading = .true.
        !       The added mass is treated as swept-up ambient gas. Its tracer content
        !       is added later using the local ISM abundance.
        
        ism_mass = 0.0_dp
        if (mass_load_yields) then
            if (ism_loading) then
                ! Swept-up ISM mass to be multiplied by the local ISM abundance later.
                ism_mass = injectMassIn*mass_load_factor
            else
                ! Treat loaded mass as wind material with the same composition as the
                ! input wind.
                injectYield = injectYieldIn*(1.0d0+mass_load_factor)
            endif
        endif
#endif

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) then
            write(*, 921) "Limiting velocity, mass_loading_factor = ", mass_load_factor
            write(*, 922) "Mass and velocity after loading, injectMass, injectVelocity = ", injectMass, injectVelocity
        endif
#endif

!
!#ifdef DEBUG_ENERGY
!        if (gr_meshMe == 0) then
!            write(*,'(A,ES10.3)') "Mass load factor = ", mass_load_factor
!            write(*,'(A,ES10.3,A,ES10.3,A)') "Injecting", injectMassIn/dt/solarMass*yr, &
!                " solar masses/yr at ", injectVelocityIn, " cm/s"
!            !write(*,'(A,ES10.3,A,ES10.3,A)') "(", injectMass, " g over ", dt, " seconds)"
!            write(*,'(A,3ES10.3)') "Star loc  =", loc
!            write(*,'(A,ES10.3,A)') "For a total of ", &
!                0.5 * injectMassIn * injectVelocityIn**2/dt, " ergs/s"
!        end if
!#endif

    end if

    ! =========================================================================
    ! Prepare for injection.
    ! =========================================================================
#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Preparing injection region"
#endif

    ! Build array of block IDs for all blocks to be injected. Get their center
    ! distances from the injection star. Then for each cell in each block, calculate
    ! its overlap with the injection sphere and store the value. Store the
    ! components of a velocity vector pointing radially outwards from the star

    if (iHaveInjectBlk) then
        ! This MPI rank has blocks overlapping with injection region. We will:
        !   1. Allocate memory and initialize arrays for injection region.
        !   2. Save list of block IDs where injection will occur.
        !   3. Calculate the fraction of each cell in each block that is overlapping 
        !      with injection region and a vector pointing radially outwards from the 
        !      injection star.

!#ifdef DEBUG_MPI
!        print *, "Found", injBlkNum, "injection blocks on proc ", gr_meshMe
!#endif

        ! Initialize arrays local to this MPI rank.
        allocate(localInjectBlocks(injBlkNum))
        allocate(injectDataOverlap(injBlkNum, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI))
        allocate(injectDataVel(injBlkNum, GRID_ILO:GRID_IHI, GRID_JLO:GRID_JHI, GRID_KLO:GRID_KHI, 1:3))
        localInjectBlocks = 0
        injectDataOverlap = 0.0d0
        injectDataVel     = 0.0d0

!#ifdef DEBUG_MPI
!        print *, "Allocations done"
!#endif
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "Finding indices of all inject blocks on each MPI rank."
#endif

        ! Loop over all blocks and save IDs of those overlapping with injection region.
        n = 1
        do blockID = 1, lnblocks
            if(nodetype(blockID) == LEAF) then
                ! Determine if injection region overlap with leaf block bounding box owned 
                ! by this MPI rank. The algorithm is identical to earlier when counting 
                ! blocks with overlap. See earlier code for more details.
                call Grid_getBlkCenterCoords(blockID,blkCtr)
                call Grid_getBlkPhysicalSize(blockID,blkSize)

                ! Coordinates for box which is closest to the injection sphere center.
                xcoll = max(blkCtr(1)-0.5*blkSize(1),min(loc(1),blkCtr(1)+0.5*blkSize(1)))
                ycoll = max(blkCtr(2)-0.5*blkSize(2),min(loc(2),blkCtr(2)+0.5*blkSize(2)))
                zcoll = max(blkCtr(3)-0.5*blkSize(3),min(loc(3),blkCtr(3)+0.5*blkSize(3)))

                if ((xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2<injectRadiusMax**2) then
                    ! This block is overlapping with the injection region. Save block ID.
                    localInjectBlocks(n) = blockID
                    n = n + 1
                end if
            end if
        end do

!#ifdef DEBUG_MPI
!        print *, "Found injection blocks:", localInjectBlocks, "on proc", gr_MeshMe
!#endif
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "Calculating overlap fraction in each cell receiving material."
#endif

        nOverlap = 0
        sumOverlap = 0.0_dp
        ! Loop over all blocks that were identified as within the injection region.
        do n = 1, injBlkNum
            blockID = localInjectBlocks(n)
            call Grid_getBlkPtr(blockID, solndata)
            ! Loop over all cells in the block.
            do k = GRID_KLO, GRID_KHI
                do j = GRID_JLO, GRID_JHI
                    do i = GRID_ILO, GRID_IHI
                        ! Compute physical cell-center coordinates from the block center.
                        ! For NXB=8, the interior-cell offsets relative to the block center are
                        !
                        !   cell #:      1    2    3    4    5    6    7    8
                        !   offset/dx: -3.5 -2.5 -1.5 -0.5 +0.5 +1.5 +2.5 +3.5
                        !
                        ! The array index includes NGUARD guard cells, so the interior-cell number is
                        ! i - NGUARD. Therefore offset/dx = i - NGUARD - 0.5*NXB - 0.5
                        x = coord(1,blockID) + delta(1) * (i - NGUARD - 0.5_dp * (NXB + 1))
                        y = coord(2,blockID) + delta(2) * (j - NGUARD - 0.5_dp * (NYB + 1))
                        z = coord(3,blockID) + delta(3) * (k - NGUARD - 0.5_dp * (NZB + 1))

                        ! Coordinates for box which is closest to the injection sphere center.
                        xcoll = max(x-0.5*delta(1),min(loc(1),x+0.5*delta(1)))
                        ycoll = max(y-0.5*delta(2),min(loc(2),y+0.5*delta(2)))
                        zcoll = max(z-0.5*delta(3),min(loc(3),z+0.5*delta(3)))

                        ! As before, calculate distance from this point to the injection region center.
                        d2coll = (xcoll-loc(1))**2+(ycoll-loc(2))**2+(zcoll-loc(3))**2

                        ! Cycle to next cell if outside injection region.
                        if (d2coll > injectRadius**2) cycle
                        
                        ! Calculate the overlapping volume of injection region and this cell.
                        ! The subroutine overlap(ishp, rad, center, cell_bot, cell_top, nsteps, overlap_vol)
                        ! does a Monte Carlo-like integration between the cell and a sphere (ishp = 1)
                        ! modified with tapered center-weighting using nsteps=10 points in each dimension.
                        ! to sample the overlap.
                        ! Determine the lower and upper cell bounds for the overlap calculation.
                        ! TO DISCUSS: The original implementation computes the cell bounds using
                        ! sign(abs(x) ± dx/2, x) instead of the more obvious x ± dx/2.
                        ! For negative coordinates this reverses the ordering of cell_bot and
                        ! cell_top; e.g. x=-10 gives cell_bot=-9.5 and cell_top=-10.5.
                        !
                        ! This is probably okay for the current overlap() routine because it
                        ! samples uniformly between cell_bot and cell_top, so reversing the order
                        ! only changes the sign of the sampling step and visits the same points in
                        ! reverse order. However, the naming is misleading and this convention
                        ! should be revisited if overlap() is modified or replaced.
                        cell_bot = [ sign(abs(x) - 0.5*delta(1), x), &
                                     sign(abs(y) - 0.5*delta(2), y), &
                                     sign(abs(z) - 0.5*delta(3), z) ]
                        cell_top = [ sign(abs(x) + 0.5*delta(1), x), &
                                     sign(abs(y) + 0.5*delta(2), y), &
                                     sign(abs(z) + 0.5*delta(3), z) ]
                        call overlap(1, injectRadius, loc, cell_bot, cell_top, 10, overlap_frac)

                        ! Determine the vector from the star to the cell center. Then use it to
                        ! set the direction of the wind velocity field.
                        dx = x - loc(1)
                        dy = y - loc(2)
                        dz = z - loc(3)
                        rad2 = dx**2 + dy**2 + dz**2
                        rad  = sqrt(rad2)
                        
                        ! Skip cells inside the minimum radius. This avoids division by zero for a
                        ! cell centered exactly on the star and preserves the original min_radius cut.
                        if (rad .lt. min_radius) cycle
                        
                        ! Set velocity of the wind to point radially outwards.
                        xvel = dx/rad * injectVelocity
                        yvel = dy/rad * injectVelocity
                        zvel = dz/rad * injectVelocity

                        ! If working with variable injection radius we need to calculate the background
                        ! density to determine R_1 (see earlier code). Note that this will only happen
                        ! the first time a star is injecting wind.
                        if (calcBgDens) bgDens = bgDens + overlap_frac*solndata(DENS_VAR, i, j, k)

                        if (.not. mass_load) then
                            ! Unless mass-loading winds, add an additional weight to the volume overlap
                            ! using solid angle calculated from the apparent angular size of the cell face
                            ! as seen by the star. In this way, the injection weighting is more like a 
                            ! radial wind flux
                            ! -> cells covering a larger angular patch around the star receive more weight.
                            !
                            ! The expression for the solid angle calculation is a special case of eqn 27 
                            ! in Mathar (2022, vixra.org/abs/2001.0603, see also Khadjavi, 1968) where the 
                            ! plate is a square, i.e., alpha = beta = delta/(2*rad).
                            ! 
                            ! Notes: 
                            !   - This should be considered an approximation because it assumes the cell face 
                            !     is perpendicular to the radial direction which is typically not the case.
                            !   - It is not clear why this correction is only applied when mass_load is false. 
                            !     Arguably, mass loading can be viewed as sweeping up ambient ISM in which
                            !     case it should only be volume weighted (which is the case in a uniform medium).
                            !   - FIXED! The original expression appears to be missing a factor 4**2 in the final term
                            !     of the denominator. It is kept unchanged here to preserve previous results.
                            solidAngle   = 4.0_dp*acos(sqrt((1.0_dp+del2/(2.0_dp*rad2)) &
                                       & / (1.0_dp+ del2/(2.0_dp*rad2) + (del2/rad2/4.0_dp)**2.0_dp)))
                            overlap_frac = overlap_frac*solidAngle
                        end if
                        
                        ! Save overlap fraction and injection velocity of each cell if they exceed a minimum radius.
                        ! Also store number of overlapping cell. The total number across all MPI ranks will be
                        ! gathered later.
                        if (overlap_frac .gt. 0.0d0) then
                            sumOverlap = sumOverlap + overlap_frac
                            nOverlap = nOverlap + 1
                            injectDataOverlap(n,i,j,k) = overlap_frac
                            injectDataVel(n,i,j,k,1:3) = [xvel,yvel,zvel]
                        end if
                    end do
                end do
            end do
            call Grid_releaseBlkPtr(blockID, solndata)
        end do

!#ifdef DEBUG_MPI
!        print *, "Calculated overlaps", gr_meshMe
!#endif
!
    end if

!#ifdef DEBUG_MPI
!    print*, "Proc ", gr_meshMe, " about to call MPI with sumOverlap = ", sumOverlap
!#endif

    ! Sum up overlap fraction from all MPI ranks. Later, we will use the total overlap
    ! to rescale all quantities we inject.
    call MPI_ALLREDUCE(MPI_IN_PLACE, sumOverlap, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0)  write(*, 921) "Summed overlap fraction across all MPI ranks, sumOverlap = ", sumOverlap
#endif

    if (sumOverlap <= 0.0_dp) then
        call Driver_abortFlash("[inject_direct] sumOverlap reduces to negative/zero after overlap calculation.")
    endif

    if (var_radius .and. calcBgDens) then
        ! If this is the first variable-radius injection step, convert the
        ! overlap-weighted density sum into a global overlap-weighted mean density.
!#ifdef DEBUG_MPI
!        print*, "Before MPI background density is", bgDens
!#endif

        call MPI_ALLREDUCE(MPI_IN_PLACE, bgDens, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
        bgDens = bgDens / sumOverlap
#ifdef WIND_VERBOSE
        if (gr_meshMe == 0)  write(*, 921) "Global overlap-weighted mean density, bgDens = ", bgDens 
#endif
    end if
 

!#ifdef DEBUG
!    if (gr_meshMe == 0) then
!        print *, "Total overlap weight is", sumOverlap
!        print*, "Background density is", bgDens
!    end if
!#endif

    if (perturb_velocity) then
        ! Perturb the magnitude of the wind velocity with a normal distribution of
        ! width perturb_std_dev. This introduces small cell-by-cell differences that
        ! reduce numerical artifacts from perfectly spherical injection.

#ifdef WIND_VERBOSE
        if (gr_meshMe == 0) write(*, 900) "Adding velocity perturbations"
        largestDelta = 0.0_dp
        largestDeltaVel = 0.0_dp
#endif

        ! Calculate velocity perturbations on every MPI rank, including ranks without
        ! injection cells, to maintain a synchronized random number stream.
        call MPI_ALLGATHER(nOverlap, 1, MPI_INTEGER, nOverlapArr, 1, MPI_INTEGER, gr_meshcomm, ierr)
        nOverlapTot = sum(nOverlapArr)

        allocate(perturbScale(nOverlapTot))
        do ip = 1, nOverlapTot
            perturbScale(ip) = norm_rand(1.0_dp, perturb_std_dev)
            ! Limit pertubation to 5 sigma.
            perturbScale(ip) = min(perturbScale(ip), 1.0_dp + 5.0_dp * perturb_std_dev)
            perturbScale(ip) = max(perturbScale(ip), 1.0_dp - 5.0_dp * perturb_std_dev)
            ! Prevent sign flip in case of large perturb_std_dev.
            perturbScale(ip) = max(perturbScale(ip), 1.0e-6_dp)
#ifdef WIND_VERBOSE
            largestDelta = max(largestDelta, abs(perturbScale(ip) - 1.0_dp))
#endif
        end do

        ! Calculate starting index ip in the global perturbScale array for this rank.
        ! nOverlapArr is indexed by MPI rank + 1, since MPI ranks are zero-based
        ! while Fortran arrays are one-based. Thus this rank starts after all cells
        ! owned by lower-rank MPI tasks.
        !
        ! For example, if nOverlapArr = [1, 10, 0, 1] for ranks 0, 1, 2, and 3, then
        ! perturbScale has 12 entries assigned as:
        !   rank 0 -> perturbScale(1)
        !   rank 1 -> perturbScale(2:11)
        !   rank 2 -> no entries
        !   rank 3 -> perturbScale(12)
        ip = 1
        if (gr_meshMe > 0) then

            ! TO DISCUSS: In standard Fortran, sum(nOverlapArr(1:0)) is a zero-sized sum and
            ! evaluates to 0. The explicit gr_meshMe > 0 conditional is redundant. 
            ! Consider removing this check and Arons comment.

            ! gr_meshMe starts at 0, sum(nOverlapArr(1:0)) should give 0, so we may not
            ! need this case logic.  But I'm not sure if the behavior of (1:0) slice is
            ! specified by Fortran standard. -ATr,2020aug27
            ip = 1 + sum( nOverlapArr(1:gr_meshMe) )
        end if

        ! Loop over all cells in all blocks and apply the velocity perturbation.
        perturbEkin = 0.0_dp
        do n = 1, injBlkNum
            do k = GRID_KLO, GRID_KHI
                do j = GRID_JLO, GRID_JHI
                    do i = GRID_ILO, GRID_IHI
                        if (injectDataOverlap(n,i,j,k) > 0.0_dp) then

                            injectDataVel(n,i,j,k,1:3) = injectDataVel(n,i,j,k,1:3) * perturbScale(ip)
                            weight = injectDataOverlap(n,i,j,k) / sumOverlap
                            perturbEkin = perturbEkin + &
                                        weight * 0.5_dp * injectMass * sum(injectDataVel(n,i,j,k,1:3)**2)

#ifdef WIND_VERBOSE
                            largestDeltaVel = max(largestDeltaVel, sqrt(sum(injectDataVel(n,i,j,k,:)**2)))
#endif
                            ip = ip + 1
                            
!                            ! This condition checks whether the final velocity is within 50% of the
!                            ! 5-sigma velocity scale. This is a somewhat unusual stopping criterion.
!                            ! Consider replacing it with a bounded normal distribution or an explicit
!                            ! check for pathological cases, such as a negative perturbation factor that
!                            ! reverses the wind direction.
!                            if (abs(sqrt(sum(injectDataVel(n,i,j,k,:)**2)) &
!                                    - injectVelocity*(1.0+5.0*perturb_std_dev)) &
!                                    /(injectVelocity*(1.0+5.0*perturb_std_dev))> 5.0d-1) then
!                                write(*,'(A,2ES13.3E3)') "vels don't match: calc vel, inj vel= ",  &
!                                    sqrt(sum(injectDataVel(n,i,j,k,:)**2)), injectVelocity
!                                write(*,'(A,3ES13.3E3)') "xvel, yvel zvel =", &
!                                    injectDataVel(n,i,j,k,1), injectDataVel(n,i,j,k,2), injectDataVel(n,i,j,k,3)
!                                stop
!                            end if
                        end if
                    end do
                end do
            end do
        end do

        ! Prevent the global kinetic energy from velocity perturbations from being
        ! too large by renormalizing to the injected kinetic energy. Large perturbations
        ! can otherwise cause negative thermal energies from developing.
        ! Note that this normalization conserves the total injected kinetic energy, but does
        ! not guarantee that every individual cell has non-negative thermal
        ! compensation.

        call MPI_ALLREDUCE(MPI_IN_PLACE, perturbEkin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
        if (perturbEkin > 0.0_dp) then
            perturbNorm = sqrt(injectEkin / perturbEkin)
        else
            call Driver_abortFlash("[inject_direct] perturbed wind kinetic energy <= 0")
        end if

#ifdef WIND_VERBOSE
        call MPI_ALLREDUCE(MPI_IN_PLACE, largestDelta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, largestDeltaVel, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
        if (gr_meshMe == 0) then
            write(*, 922) "Perturbed kinetic energy before renorm: perturbEkin, injectEkin =", &
                perturbEkin, injectEkin
            write(*, 921) "Velocity perturbation renormalization: perturbNorm =", perturbNorm
            write(*, 922) "Largest perturbation: |scale-1|, max speed before renorm =", &
                largestDelta, largestDeltaVel
        end if
        largestDeltaVel = 0.0_dp
#endif

        do n = 1, injBlkNum
            do k = GRID_KLO, GRID_KHI
                do j = GRID_JLO, GRID_JHI
                    do i = GRID_ILO, GRID_IHI
                        if (injectDataOverlap(n,i,j,k) > 0.0_dp) then
                            injectDataVel(n,i,j,k,1:3) = injectDataVel(n,i,j,k,1:3) * perturbNorm

#ifdef WIND_VERBOSE
                            largestDeltaVel = max(largestDeltaVel, sqrt(sum(injectDataVel(n,i,j,k,:)**2)))
#endif
                        end if
                    end do
                end do
            end do
        end do

#ifdef WIND_VERBOSE
        call MPI_ALLREDUCE(MPI_IN_PLACE, largestDeltaVel, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
        if (gr_meshMe == 0) write(*, 921) &
            "Maximum perturbed wind speed after renorm =", largestDeltaVel
#endif

        ! Sanity check to ensure that each MPI rank used exactly its assigned
        ! perturbation values. Every value in perturbScale should be used exactly
        ! once, i.e. ranks must take disjoint chunks from perturbScale that together
        ! cover the entire array.
        perturbUsedFullChunk = .true.
        if (ip /= 1 + sum(nOverlapArr(1:gr_meshMe+1))) then
            perturbUsedFullChunk = .false.
        end if
        call MPI_ALLREDUCE(MPI_IN_PLACE, perturbUsedFullChunk, 1, MPI_LOGICAL, &
            MPI_LAND, gr_meshComm, ierr)

        if (.not. perturbUsedFullChunk) then
!            if (gr_meshMe == 0) then
!                print*, "Error in wind velocity perturbation, RNG stream sampled wrongly"
!            end if
            call Driver_abortFlash("[inject_direct] Error in wind velocity perturbation, RNG stream sampled wrongly")
        end if
        
        deallocate(perturbScale)
    end if

    ! =========================================================================
    ! Inject material.
    ! =========================================================================
    ! Inject mass, momentum/energy, and tracers into each overlapping cell.
    ! Two modes exists and differ in how the post-injection velocity is chosen:
    !
    !   conserved_quant == "momentum":
    !       Momentum-conserving injection with optional thermal compensation.
    !       We first mix the injected material with the existing cell gas by conserving
    !       vector momentum. This fixes the post-injection velocity and therefore the
    !       kinetic-energy change in the cell. If add_therm_e is enabled, any difference
    !       between the target injected mechanical energy and the kinetic-energy change
    !       is added as thermal energy. This represents unresolved thermalization of the
    !       wind kinetic energy in an inelastic collision/shock.
    !
    !   conserved_quant == "energy":
    !       New velocity is computed by mixing signed velocity-squared
    !       components. 
    !       TO DISCUSS: This looks like a legacy energy-conserving prescription 
    !       which has fairly severe algorithmic limitations. Notably, the
    !       equivalent implementation for SNe has this version commented out.
    !       We should consider deprecating this verion and remove it from the
    !       code. The more physical way of conserving energy is using momentum
    !       conservation together with add_therm_e.
    
#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Injection step"
    
    ! Reset diagnostics
    globalTE = 0.0_dp    
    globalKE = 0.0_dp    
    globalDeltaE = 0.0_dp
    globalDeltaP = 0.0_dp
    largestTE = 0.0_dp   
    largestKE = 0.0_dp   
    sumMass = 0.0_dp     
    nNegativeTherm = 0
    sumNegativeTherm = 0.0_dp
#endif

    if (iHaveInjectBlk) then
        if (conserved_quant .eq. "momentum") then
            do n = 1, injBlkNum
                blockID = localInjectBlocks(n)
                call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
                call Grid_getBlkPtr(blockID, solndata)

                do k = GRID_KLO, GRID_KHI
                    do j = GRID_JLO, GRID_JHI
                        do i = GRID_ILO, GRID_IHI

                            ! Skip cells that receive no injected material. Even when no mass is added,
                            ! recomputing the conserved variables can introduce small roundoff errors
                            ! that propagate through the EOS update.
                            if (injectDataOverlap(n,i,j,k) .le. 0.0_dp) cycle

                            ! State before injection.
                            oldDens     = solndata(DENS_VAR,i,j,k)
                            oldVel      = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                            oldEkinDens = 0.5_dp * oldDens * sum(oldVel**2)
                            oldEtheDens = oldDens * solndata(EINT_VAR,i,j,k)

                            ! Injected material.
                            ! Note that injectEkin, which is the target mechanical energy density 
                            ! assigned to this cell, is based on the original input wind energy before 
                            ! any mass-loading modification.
                            weight      = injectDataOverlap(n,i,j,k) / sumOverlap
                            injDens     = weight * injectMass/dVol
                            injVel      = injectDataVel(n,i,j,k,1:3)
                            injEkinDens = weight * injectEkin / dVol
                            injEtheDens = 0.0_dp
                            
                            ! State after injection.
                            newDens     = oldDens + injDens
                            newMom      = injVel * injDens + oldVel * oldDens
                            newVel      = newMom / newDens
                            newEkin     = 0.5_dp * sum(newVel**2)
                            newEkinDens = newEkin * newDens
                            if (add_therm_e) then
                                ! Update thermal energy to conserve total energy.
                                injEtheDens = injEkinDens  - (newEkinDens - oldEkinDens)
                            
                                if (injEtheDens < 0.0_dp) then
#ifdef WIND_VERBOSE
                                    nNegativeTherm = nNegativeTherm + 1
                                    sumNegativeTherm = sumNegativeTherm + injEtheDens
#endif
#ifdef WIND_DEBUG
                                    write(error_unit,921) "Negative thermal compensation in wind injection: twind = ", twind
                                    write(error_unit,912) "rank, blockID", gr_meshMe, blockID
                                    write(error_unit,923) "oldEkinDens, injEkinDens, newEkinDens =", &
                                        oldEkinDens, injEkinDens, newEkinDens
                                    write(error_unit,922) "oldEtheDens, injEtheDens =", &
                                        oldEtheDens, injEtheDens
                                    write(error_unit,923) "oldVel =", oldVel(1), oldVel(2), oldVel(3)
                                    write(error_unit,923) "injVel =", injVel(1), injVel(2), injVel(3)
#endif
                                    injEtheDens = 0.0_dp
                                end if
                            endif
                            newEthe     = (oldEtheDens + injEtheDens) / newDens

#ifdef WIND_DEBUG
                            if (newEkin > 0.0_dp .and. newEthe / newEkin < hy_eswitch) then
                                write(*, 922) "Cell is kinetic-energy dominated: eint/ekin, hy_eswitch =" &
                                    , newEthe / newEkin, hy_eswitch
                            end if
#endif

                            ! Limit the internal energy to the temperature floor.
                            ! This is primarily a problem for runs where momentum
                            ! is conserved but no thermal energy is added to 
                            ! compensate. This cause energy to be diluted and the
                            ! internal energy starts approaching zero.
                            newEthe = max(newEthe, eint_floor)
                            
                            ! Update the hydrodynamic solution
                            solndata(DENS_VAR, i, j, k) = newDens
                            solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
                            solndata(EINT_VAR, i, j, k) = newEthe
                            solndata(ENER_VAR, i, j, k) =  newEkin + newEthe

#ifdef WIND_VERBOSE
                            ! Summary statistics for reporting and sanity checks
                            globalKE     = globalKE + (newEkinDens - oldEkinDens) * dVol
                            globalTE     = globalTE + injEtheDens*dVol
                            globalDeltaE = globalDeltaE + ((newEkinDens - oldEkinDens) + injEtheDens)*dVol
                            globalDeltaP = globalDeltaP + injDens * sqrt(sum((injVel)**2.0_dp))*dVol
                            largestTE    = max(largestTE, injEtheDens*dVol)
                            largestKE    = max(largestKE, injEkinDens*dVol)
                            sumMass = sumMass + injDens * dVol
#endif

#ifdef TRACER_FIELDS
                            ! Finally, update the tracer fields. Note that we work in conserved quantities 
                            ! during the injection, and convert to fraction after the update.
                            do itracer = 1, NMASS_SCALARS
                                ! Tracer density before the update.
                                if(mass_load_yields .and. ism_loading) then
                                    ! Here, ism_mass is the additional mass due to mass loading if applied.
                                    ! It assumes additional mass is from swept-up material (old metallicity)
                                    oldTracerField(itracer) = solndata(MASS_SCALARS_BEGIN + (itracer-1), i, j, k) &
                                                            * (oldDens + weight * ism_mass / dVol)
                                else
                                    oldTracerField(itracer) = solndata(MASS_SCALARS_BEGIN+(itracer-1),i,j,k) * oldDens
                                endif

                                ! Tracers density to be injected
                                if (injectYield(itracer) < 0.0_dp) then
                                    injTracerField(itracer) = solndata(MASS_SCALARS_BEGIN + (itracer-1), i, j, k) * injDens
                                else
                                    injTracerField(itracer) = weight * injectYield(itracer) / dVol
                                endif
                                ! Tracers density after the injection
                                newTracerField(itracer) = oldTracerField(itracer) + injTracerField(itracer)
                            enddo

                            ! Update scalar field to new metallicity after wind mass injection
                            do itracer = 1, NMASS_SCALARS
                                solndata(MASS_SCALARS_BEGIN+(itracer-1), i, j, k) = newTracerField(itracer) / newDens
                            enddo
#endif

#ifdef WIND_DEBUG
                            if (oldDens <= 0.0_dp .or. injDens < 0.0_dp .or. newDens <= 0.0_dp .or. &
                                oldEtheDens <= 0.0_dp .or. newEthe <= 0.0_dp .or. &
                                oldEkinDens > 1.0e100_dp .or. newEkinDens > 1.0e100_dp .or. newEthe > 1.0e100_dp .or. &
                                injEtheDens > 1.0e100_dp .or. injEtheDens < -1.0e100_dp) then
                            
                                write(error_unit,900) "Bad state detected during wind injection."
                                write(error_unit,912) "rank, blockID", gr_meshMe, blockID
                                write(error_unit,923) "twind, dt, weight =", twind, dt, weight
                            
                                write(error_unit,923) "oldDens, injDens, newDens =", oldDens, injDens, newDens
                                write(error_unit,923) "oldEkinDens, injEkinDens, newEkinDens =", &
                                    oldEkinDens, injEkinDens, newEkinDens
                                write(error_unit,923) "oldEtheDens, injEtheDens, newEthe =", &
                                    oldEtheDens, injEtheDens, newEthe
                            
                                write(error_unit,923) "oldVel =", oldVel(1), oldVel(2), oldVel(3)
                                write(error_unit,923) "injVel =", injVel(1), injVel(2), injVel(3)
                                write(error_unit,923) "newVel =", newVel(1), newVel(2), newVel(3)
                                call flush(error_unit) 
                            end if
#endif

!                            Summary statistics for reporting and sanity checks
!                            sumMass = sumMass + injDens*dVol         ! Total injected mass 
!                            initialKE = 0.5_dp * sum(oldVel**2.0_dp) ! Per particle per unit mass, KE/M = KE / (N * mu * m_H).
!                            initialTE = solndata(EINT_VAR,i,j,k)     ! Per particle per unit mass, TE/M = TE / (N * mu * m_H).
!
!                            !oldE = oldDens * sum(oldVel**2)
!                            oldE = (0.5_dp*sum(oldVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*oldDens
!                            oldP = oldDens * sqrt(sum((oldVel+0.0_dp)**2.0_dp))
!
!                            ! Note in Flash energy is in erg/g so that the energy is the same regardless of refinement. - JW
!                            ! Also don't forget the already existing internal energy in this cell. - JW
!
!                            ! Specific KE added to this cell. - JW
!                            !                      deltaKinE = (0.5_dp * newDens * sum(newVel**2.0_dp) &
!                            !				      - 0.5_dp * oldDens * sum(oldVel**2.0_dp)) &
!                            !				      / newDens
!                            ! Didn't we already account for mass ratios etc when we inverted the
!                            ! momentum eqn (i.e. isnt that the point of doing it this way)? - JW
!                            !addKinE   = 0.5_dp * sum(injVel**2.0_dp)*injDens
!                            deltaKinE = 0.5_dp * sum(newVel**2.0_dp)*newDens - 0.5_dp * sum(oldVel**2.0_dp)*oldDens
!                            addKinE = deltaKinE
!                            ! Thermal energy should be the difference b/t the injected specific E
!                            ! and the actual specific kinetic energy. - JW
!
!                            ! Final TE in this cell = portion of mech
!                            ! energy added to the cell + the existing TE in cell.
!                            ! Note what we did here. In order to make up for the energy loss from conserving momentum
!                            ! we used the original velocity from injEkinDens to ensure that we added back enough energy
!                            ! in thermal to get back to the proper total energy at the end, and
!                            ! NOT the modified velocity from mass loading or what was calculated from momentum
!                            ! conservation. - JW
!                            if (add_therm_e) then
!                                ! - KE added to cell. - JW !0.5_dp*injDens*sum(injVel**2.0_dp)/newDens
!                                !		              finalThermE = injEkinDens*injectDataOverlap(n,i,j,k)/sumOverlap/(newDens*dVol)  &
!                                !		                  + solndata(EINT_VAR, i, j, k)*oldDens/newDens &
!                                !				          - deltaKinE
!                                ! Note I think this is wrong, since the "correct" energy being added is the mass
!                                ! loaded energy, and the mass loading is the stuff falling in from the shell, which
!                                ! has zero kinetic energy. So we *should lose some kinetic energy here*!
!                                !injEtheDens = finalThermE - solndata(EINT_VAR, i, j, k)
!                                injEtheDens = (injEkinDens*injectDataOverlap(n,i,j,k)/sumOverlap/dVol  &
!                                    - deltaKinE)
!
!                                ! Lets add in the thermal energy of the mass loaded material to be the temperature
!                                ! the ionized shell material should be on the inner edge, which is 1e4 K.
!                                ! e = 1/(gamma-1)*k*T/mu - JW
!                                !injEtheDens = (1.d0/(gamma_-1.0d0)*1.3807d-16*1e4/(1.67e-24*0.61))*injDens
!                            else
!                                injEtheDens = 0.0d0
!                            end if
!                            !ThermE = (0.5_dp * newDens * sum(newVel**2.0_dp) &
!                            !       + solndata(EINT_VAR, i, j, k) - oldE
!                            ! Record the change in TE for debugging.
!                            deltaThermE = injEtheDens
!                            ! Now add the "missing" kinetic energy as thermal energy. - JW
!
!                            !solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR,i,j,k) + deltaThermE
!                            ! Shouldn't this have a factor that includes how much more mass we spread the old energy
!                            ! around to by adding mass to the cell? In other words, this is the internal energy
!                            ! per particle. So if we had 2 particles at the old TE and 1 particle at the TE of the
!                            ! mass added to the cell, the new TE = (2*oldTE + 1*new_massTE) / 3 - JW
!                            solndata(EINT_VAR, i, j, k) = solndata(EINT_VAR,i,j,k)*(oldDens/newDens) + injEtheDens/newDens
!                            !if (solndata(EINT_VAR,i,j,k) .lt. oldThermE) then
!                            !  print*, "We lost thermal energy in this cell. oldTE=, newTE=", &
!                            !  oldThermE, solndata(EINT_VAR,i,j,k)
!                            !end if
!
!                            solndata(ENER_VAR, i, j, k) =  0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k)
!                            if (solndata(EINT_VAR,i,j,k)/(0.5_dp*sum(newVel**2.0_dp)) < hy_eswitch) then
!                                print*, "[inject_direct]: Warning, eintswitch should be kicking on here."
!                                print*, "eint/ke=",solndata(EINT_VAR,i,j,k)/(0.5_dp*sum(newVel**2.0_dp))
!                                print*, "eintSwitch=", hy_eswitch
!                            end if
!
!                            ! new velocity is calculated by conserving momentum when injDens moving
!                            ! at injectVelocity is added to oldDens moving at original velocity
!
!                            solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
!                            solndata(DENS_VAR, i, j, k) = newDens
!
!                            if (addKinE .gt. largestKE) largestKE = addKinE*dVol
!                            if (injEtheDens .gt. largestTE) largestTE = injEtheDens*dVol
!
!                            ! Delta E (not specific!) added to this cell. - JW
!                            deltaE = (deltaKinE + deltaThermE)
!
!                            !newE = newDens * sum(newVel**2)
!                            !newE = sum(newVel**2.0_dp)
!                            newE = (0.5_dp * sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k))*newDens
!                            !newP = newDens * sqrt(sum(newVel**2))
!                            newP = sqrt(sum((newMom+0.0_dp)**2.0_dp))
!                            !globalDeltaE = globalDeltaE + 0.5 * (newE &
!                            !             - oldE )*dVol
!
!
!#ifdef DEBUG
!                            if ((injectDataOverlap(n,i,j,k) .le. 0.0_dp) .and. &
!                                (injDens .gt. 0.0d0 .or. deltaKinE .gt. 0.0d0 .or. deltaThermE .gt. 0.0 &
!                                .or. abs(newE - oldE) .gt. 0.0d0)) then
!
!                                write(*,*) "BS detected! Overlap is zero but someone is up to no good!"
!                                write(*,'(A,ES13.3E3)') "Delta KE= ", deltaKinE
!                                write(*,'(A,ES13.3E3)') "Delta TE= ", deltaThermE
!                                write(*,'(A,ES13.3E3)') "injDens= ", injDens
!                                write(*,'(A,2ES13.3E3)') "oldDens, newDens= ", oldDens, newDens
!                                write(*,'(A,6ES13.3E3)') "oldVel, newVel= ", oldVel, newVel
!                                write(*,'(A,3ES13.3E3)') "injVel =", injVel
!                                write(*,'(A,2ES13.3E3)') "old KE, old TE= ", initialKE, initialTE
!                                write(*,'(A,6ES13.3E3)') "new KE, new TE= ", 0.5_dp*sum(newVel**2.0_dp), solndata(EINT_VAR, i, j, k)
!                                write(*,'(A,ES13.3E3)') "initial Total E= ", oldE
!                                write(*,'(A,ES13.3E3)') "final Total E= ", newE
!                                write(*,'(A,ES13.3E3)') "(newE-oldE)= ", (newE-oldE)
!                                call flush(6)
!                            end if
!
!                            if (injectDataOverlap(n,i,j,k) > 0.0) then
!                                write(*,'(A,ES13.3E3)') "Delta KE= ", deltaKinE
!                                write(*,'(A,ES13.3E3)') "Delta TE= ", deltaThermE
!                                write(*,'(A,ES13.3E3)') "Delta KE= ", addKinE
!                                write(*,'(A,ES13.3E3)') "Delta TE= ", injEtheDens
!                                write(*,'(A,ES13.3E3)') "injDens= ", injDens
!                                write(*,'(A,2ES13.3E3)') "oldDens, newDens= ", oldDens, newDens
!                                write(*,'(A,6ES13.3E3)') "oldVel, newVel= ", oldVel, newVel
!                                write(*,'(A,3ES13.3E3)') "injVel =", injVel
!                                write(*,'(A,2ES13.3E3)') "old KE, old TE= ", initialKE, initialTE
!                                write(*,'(A,6ES13.3E3)') "new KE, new TE= ", 0.5_dp*sum(newVel**2.0_dp), solndata(EINT_VAR, i, j, k)
!                                write(*,'(A,ES13.3E3)') "final TE= ", solndata(EINT_VAR, i, j, k)
!                                write(*,'(A,ES13.3E3)') "initial Total E= ", oldE
!                                write(*,'(A,ES13.3E3)') "final Total E= ", newE
!                                write(*,'(A,ES13.3E3)') "deltaE= ", deltaE
!                                write(*,'(A,ES13.3E3)') "newE-oldE= ", (newE-oldE)
!                                write(*,'(A,ES13.3E3)') "deltaE - (newE-oldE)= ", deltaE - (newE-oldE)
!                                write(*,'(A,ES13.3E3)') "globalTE= ", globalTE
!                                write(*,'(A,ES13.3E3)') "globalKE= ", globalKE
!                                write(*,'(A,ES13.3E3)') "globalDeltaE= ", globalDeltaE
!                                call flush(6)
!                            end if
!#endif
!                            globalKE     = globalKE + addKinE*dVol
!                            globalTE     = globalTE + injEtheDens*dVol
!                            globalDeltaE = globalDeltaE +  (addKinE + injEtheDens)*dVol !newE - oldE !deltaE
!                            globalDeltaP = globalDeltaP + injDens * sqrt(sum((injVel)**2.0_dp))*dVol !(newP - oldP)

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

                do k = GRID_KLO, GRID_KHI
                    do j = GRID_JLO, GRID_JHI
                        do i = GRID_ILO, GRID_IHI
                            
                            ! As for momentum conserving, we should probably add 
                            ! here to not spend time updating variables that do
                            ! not change.
                            if (injectDataOverlap(n,i,j,k) .le. 0.0_dp) cycle
                            
                            ! State before injection.
                            oldDens = solndata(DENS_VAR, i, j, k)
                            oldVel  = solndata(VELX_VAR:VELZ_VAR, i, j, k)
                            ! oldE    = oldDens * sum(oldVel**2)         <--- not used, also missing factor 0.5
                            ! oldP    = oldDens * sqrt(sum(oldVel**2))   <--- not used
                            
                            ! Injected material.
                            weight  = injectDataOverlap(n,i,j,k) / sumOverlap
                            injDens = weight * injectMass / dVol
                            injVel  = injectDataVel(n, i, j, k, 1:3)
                            
                            ! State after injection.
                            ! Notes:
                            !   - The resolution is uniform (maximally refined), therefore
                            !     dens + delta dens is proportional to mass + delta mass
                            !   - Velocity is done component-by-component. This does not
                            !     follow any conservation law, e.g. if oldVel=(0,100,0) 
                            !     and injVel=(100,0,0), the scheme preserves |v| but does
                            !     not correspond to vector momentum conservation.
                            !   - Revisit this implementation as it is clearly legacy at 
                            !     this point. A pure thermal dump may be a cleaner alter-
                            !     native if energy injection without momentum conservation 
                            !     is desired.

                            newDens = oldDens + injDens
                            newVelSq = oldDens / newDens * sign(oldVel**2, oldVel) &
                                     + injDens / newDens * sign(injVel**2, injVel)
                            newVel = sign(sqrt(abs(newVelSq)), newVelSq)
!                            ! newE = newDens * sum(newVel**2)       <--- not used, missing factor 0.5
!                            ! newP = newDens * sqrt(sum(newVel**2)) <--- not used

                            ! Update the hydrodynamic solution
                            solndata(VELX_VAR:VELZ_VAR, i, j, k) = newVel
                            solndata(DENS_VAR, i, j, k) = newDens
                            solndata(ENER_VAR, i, j, k) =  0.5_dp*sum(newVel**2.0_dp) + solndata(EINT_VAR,i,j,k)

#ifdef WIND_VERBOSE
                            ! Summary statistics for reporting and sanity checks
                            globalDeltaE = globalDeltaE + 0.5_dp * injDens *sum(injVel**2.0_dp) * dVol
                            globalDeltaP = globalDeltaP + injDens * sqrt(sum((injVel)**2.0_dp)) * dVol
                            sumMass = sumMass + injDens * dVol
#endif

#ifdef TRACER_FIELDS
                            ! Finally, update the tracer fields. Note that we work in conserved quantities 
                            ! during the injection, and convert to fraction after the update.
                            do itracer = 1, NMASS_SCALARS
                                ! Tracer density before the update.
                                if(mass_load_yields .and. ism_loading) then
                                    ! Here, ism_mass is the additional mass due to mass loading if applied.
                                    ! It assumes additional mass is from swept-up material (old metallicity)
                                    oldTracerField(itracer) = solndata(MASS_SCALARS_BEGIN + (itracer-1), i, j, k) &
                                                            * (oldDens + weight * ism_mass / dVol)
                                else
                                    oldTracerField(itracer) = solndata(MASS_SCALARS_BEGIN+(itracer-1),i,j,k) * oldDens
                                endif

                                ! Tracers density to be injected
                                if(injectYield(itracer) < 0.0_dp) then
                                    injTracerField(itracer) = solndata(MASS_SCALARS_BEGIN + (itracer-1), i, j, k) * injDens
                                else
                                    injTracerField(itracer) = weight * injectYield(itracer) / dVol
                                endif
                                ! Tracers density after the injection
                                newTracerField(itracer) = oldTracerField(itracer) + injTracerField(itracer)
                            enddo

                            ! Update scalar field to new metallicity after wind mass injection
                            do itracer = 1, NMASS_SCALARS
                                solndata(MASS_SCALARS_BEGIN+(itracer-1), i, j, k) = newTracerField(itracer) / newDens
                            enddo
#endif
                        end do
                    end do
                end do

                call Grid_releaseBlkPtr(blockID, solndata)
!#ifdef DEBUG_MPI
!                print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
!                    blockID, gr_meshMe, injBlkNum
!#endif
                call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
            end do

!            
!            !else if (conserved_quant .eq. "energy_simpson") then
!            
!                ! TO DISCUSS: This appears to implement an energy-conserving update inspired 
!                ! by Simpson et al. (2015). The method attempts increase the kinetic and
!                ! thermal energy to the specified amount by adding momentum isotropically
!                ! based on the cloud-in-cell approach. As written, the solve is component-
!                ! wise and assigns deltaP/injDens as the final velocity, so it is not clear 
!                ! that it consistently conserves either vector momentum or total energy. 
!                ! Unless we want to pursue a proper and tested implementation, we should
!                ! remove this.
!            
!            !    do n = 1, injBlkNum
!            !        blockID = localInjectBlocks(n)
!            !        call Grid_getBlkPtr(blockID, solndata)
!
!            !        do i = GRID_ILO,GRID_IHI
!            !            do j = GRID_JLO, GRID_JHI
!            !                do k = GRID_KLO, GRID_KHI
!            !
!            !                  injDens = injectDataOverlap(n,i,j,k)/sumOverlap*injectMass/dVol
!            !                  oldDens = solndata(DENS_VAR,i,j,k)
!            !                  newDens = solndata(DENS_VAR,i,j,k) + injDens ! m + delta_m
!            !                  oldVel = solndata(VELX_VAR:VELZ_VAR, i, j, k)
!            !                  injVel = injectDataVel(n,i,j,k,1:3)
!            !                  newMom = injVel * injDens + oldVel * oldDens
!            !                  newVel = newMom / newDens
!
!            !                  a_ = 0.5 / newDens
!            !                  b_ = oldDens * oldVel / newDens
!            !                  c_ = (oldDens * oldVel)**2 / ( 2 * newDens ) &
!            !                    - 0.5 * oldDens * oldVel**2 &
!            !                    - 0.5 * injDens * injVel**2
!            !                  deltaP = - b_ + sqrt(b_**2 - 4*a_*c_)/(2*a_)
!
!            !                  !solndata(EINT_VAR, i, j, k) = totE/newDens - 0.5 * sum(newVel**2)
!
!            !                  ! new velocity is calculated by conserving momentum when injDens moving
!            !                  ! at injectVelocity is added to oldDens moving at original velocity
!
!            !                  solndata(VELX_VAR:VELZ_VAR, i, j, k) = deltaP / injDens
!            !                  solndata(DENS_VAR, i, j, k) = newDens
!
!            !                end do
!            !            end do
!            !        end do
!
!            !        call Grid_releaseBlkPtr(blockID, solndata)
!            !    end do

        end if

    end if

!#ifdef DEBUG_ENERGY
!    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
!    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
!    call MPI_ALLREDUCE(MPI_IN_PLACE, globalTE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
!    call MPI_ALLREDUCE(MPI_IN_PLACE, globalKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
!    call MPI_ALLREDUCE(MPI_IN_PLACE, largestTE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
!    call MPI_ALLREDUCE(MPI_IN_PLACE, largestKE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
!
!#ifdef DEBUG_MPI
!    print*, "Proc ", gr_meshMe, " about to call MPI with sumMass = ", sumMass
!#endif
!
!    call MPI_ALLREDUCE(MPI_IN_PLACE, sumMass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
!
!    if (gr_meshMe == 0) then
!        emech = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp
!        pmech = injectMassIn * injectVelocityIn
!        write(*,'(A,ES10.3,A)') "We got ", globalDeltaE/dt, "ergs/s"
!        write(*,'(A,ES10.3,A)') "We got global TE=", globalTE, "ergs"
!        write(*,'(A,ES10.3,A)') "We got global KE=", globalKE, "ergs"
!        write(*,'(A,F10.3,A)') "Error in injected wind total E is ", abs(globalDeltaE - emech)/emech*100, "%"
!        write(*,'(A,F10.3,A)') "Percentage of energy that is thermal is ", globalTE/abs(globalDeltaE)*100.0, "%"
!        write(*,'(A,F10.3,A)') "Percentage of energy that is kinetic is ", globalKE/abs(globalDeltaE)*100.0, "%"
!        write(*,'(A,F10.3,A)') "Error in injected wind P is ", abs(globalDeltaP - pmech)/pmech*100, "%"
!        write(*,'(A,ES10.3e2,A)') "Largest TE is ", largestTE
!        write(*,'(A,ES10.3e2,A)') "Largest KE is ", largestKE
!        write(*,'(A,ES10.3e2,X,A,ES10.3e2)') "Total mass from wind is ", injectMass, "Injected mass is ", sumMass
!    end if
!#endif

#ifdef WIND_VERBOSE
    call MPI_ALLREDUCE(MPI_IN_PLACE, nNegativeTherm, 1, MPI_INTEGER, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sumNegativeTherm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalDeltaP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalTE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, globalKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestTE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, largestKE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sumMass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
    
    if (gr_meshMe == 0) then
        emech = 0.5_dp * injectMassIn * injectVelocityIn**2.0_dp
        pmech = injectMassIn * injectVelocityIn

        write(*, 900) "Injection completed, hydrodynamical solution has been updated."
        write(*, 922) "Mechanical energy and momentum of wind, emech, pmech = ", emech, pmech
        write(*, 921) "Energy rate of injected material, globalDeltaE/dt = ", globalDeltaE/dt
        write(*, 922) "Kinetic and thermal energy, globalKE, globalTE = ", globalKE, globalTE
        write(*, 921) "Percentage in kinetic energy, globalKE/abs(globalDeltaE)", globalKE/abs(globalDeltaE)
        write(*, 921) "Percentage in thermal energy, globalTE/abs(globalDeltaE)", globalTE/abs(globalDeltaE)
        write(*, 921) "Error in total energy, abs(globalDeltaE - emech)/emech", abs(globalDeltaE - emech)/emech
        write(*, 921) "Error in momentum, abs(globalDeltaP - pmech)/pmech", abs(globalDeltaP - pmech)/pmech
        write(*, 921) "Largest kinetic energy in any cell is ", largestKE
        write(*, 921) "Largest thermal energy in any cell is ", largestTE
        write(*, 922) "Total and injected mass, injectMass, sumMass = ", injectMass, sumMass
        if (nNegativeTherm > 0) then
            write(error_unit,921) "WARNING: Clamped negative thermal compensation: sumNegativeTherm =", sumNegativeTherm
        end if
    endif
#endif

    ! =========================================================================
    ! Injection completed, now finalize.
    ! =========================================================================

#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Finalizing injection"
#endif

    ! Notify FLASH that the hydrodynamic solution has been modified.
    call Grid_notifySolnDataUpdate()

    ! Update timestep.
    min_wind_dt = 1d99
    if (use_wind_compute_dt) then

        if (iHaveInjectBlk) then
            do n=1, injBlkNum
                blockID = localInjectBlocks(n)
                call Grid_getBlkPtr(blockID, solndata)

                ! Note that EOS must have been called after updating hydro-variables
                ! otherwise sound speed calculation might use stale temperatures.
                cs2 = gamma_ * kB / mH * &
                    solndata(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)
                v2 = sum(solndata(VELX_VAR:VELZ_VAR,GRID_ILO:GRID_IHI, &
                    GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp, 1)

                ! While ignored here, the alfven speed is taken into account in Hydro_ComputeDt
                min_wind_dt = min(min_wind_dt, hy_cfl * minval(delta) / sqrt( maxval(cs2 + v2) ))
                
                call Grid_releaseBlkPtr(blockID, solndata)
            end do
        end if

        call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, gr_meshComm, ierr)

!#ifdef DEBUG_ENERGY
!        if (gr_meshMe == 0) write(*, 901) "Timestep limited by wind injection: dt = ", min_wind_dt
!#endif
    endif

#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 921) "Timestep limit from wind injection, min_wind_dt = ", min_wind_dt
#endif

    ! Deallocate local memory
    if (iHaveInjectBlk) then
        deallocate(localInjectBlocks)
        deallocate(injectDataOverlap)
        deallocate(injectDataVel)

!#ifdef DEBUG_MPI
!        print *, "Deallocating done for proc", gr_meshMe
!#endif   
    end if

#ifdef WIND_VERBOSE
    if (gr_meshMe == 0) write(*, 900) "Memory deallocated, exiting inject_direct"
#endif

!    if (use_wind_compute_dt) then
!        call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, &
!            MPI_MIN, gr_meshComm, ierr)
!    else
!        min_wind_dt = 1d99
!    endif
!
!
!if (iHaveInjectBlk) then
!
!    ! Update the EOS variables if we have an affected block for only those blocks.
!
!    do n=1, injBlkNum
!        blockID = localInjectBlocks(n)
!        !#ifdef DEBUG
!        !        print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ",&
!        !                                                   blockID, gr_meshMe, injBlkNum
!        !#endif
!
!        ! Moved into the actual loops where values were modified to be
!        ! absolutely certain the correct blocks get this call. - JW
!        !call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!        !call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
!
!
!        ! calculate min. crossing time of wind velocity on all injection blocks
!        if (use_wind_compute_dt) then
!
!            call Grid_getBlkPtr(blockID, solndata)
!
!            ! square of the max sound speed in the block
!            !c_sound^2 = gamma * k * T / m_H
!
!            cs2 = gamma_ * 1.3807d-16 * &
!                solndata(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI) &
!                / 1.6726d-24
!            !v2 = solndata(VELX_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp + &
!            !     solndata(VELY_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp + &
!            !     solndata(VELZ_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0_dp
!
!            ! square of magnitude of velocities in all cells of block
!            v2 = sum(solndata(VELX_VAR:VELZ_VAR,GRID_ILO:GRID_IHI, &
!                GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)**2.0,1)
!
!            ! NOTE NOTE NOTE ! We need to also calculate the Alfven velocity here
!            ! or otherwise take it into account in Hydro_computeDt (<--- this should be done now. JW) - JW
!
!            ! minimum dt to capture wind movement in this block
!            ! Lets use the hydro CFL here. - JW
!            !old_dt = minval(delta)/sqrt(cs2 + injectVelocity**2.0_dp)
!            old_dt = hy_cfl * minval(delta)/sqrt( maxval(cs2 + v2) )
!            if (old_dt .lt. min_wind_dt) then
!                min_wind_dt = old_dt
!            end if
!#ifdef DEBUG
!            print*, "Max temp on proc, blk", gr_meshMe, blockID, "is", maxval(solndata(TEMP_VAR,:,:,:))
!#endif
!            call Grid_releaseBlkPtr(blockID, solndata)
!
!        end if
!
!
!    end do
!
!    deallocate(localInjectBlocks)
!    deallocate(injectDataOverlap)
!    deallocate(injectDataVel)
!#ifdef DEBUG_MPI
!    print *, "Deallocating done for proc", gr_meshMe
!#endif
!    !call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
!    !#ifdef DEBUG
!    !    print *, "Guard cells filled on proc", gr_meshMe
!    !#endif
!end if
!
!
!if (use_wind_compute_dt) then
!    call MPI_ALLREDUCE(MPI_IN_PLACE, min_wind_dt, 1, MPI_DOUBLE_PRECISION, &
!        MPI_MIN, gr_meshComm, ierr)
!#ifdef DEBUG_ENERGY
!    if (gr_meshMe == 0) write(*,'(A,ES10.3)') "Timestep set by inject_direct = ", min_wind_dt
!#endif
!end if
!
!! Moving to outside the inject_direct call and into Particles_wind.
!!call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)
!
!#ifdef DEBUG_MPI
!    print *, "Exiting inject_direct for proc", gr_meshMe
!#endif

! The last digit indicates the number of values printed.
900 format("[inject_direct] ",A)
! Integer
911 format("[inject_direct] ",A,1X,I0)
912 format("[inject_direct] ",A,2(1X,I0))
913 format("[inject_direct] ",A,3(1X,I0))
! Real
921 format("[inject_direct] ",A,1X,ES13.5E3)
922 format("[inject_direct] ",A,2(1X,ES13.5E3))
923 format("[inject_direct] ",A,3(1X,ES13.5E3))
! Logical
931 format("[inject_direct] ",A,1X,L1)
! Strings
941 format("[inject_direct] ",A,1X,A)

end subroutine inject_direct
