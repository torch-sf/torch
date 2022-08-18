! Star particle formation routine.

! Adapted from the sink creation routine from Federrath et al 2010.

subroutine Particles_starForm()

use Driver_data, ONLY : dr_simTime
use Grid_data, ONLY: gr_meshMe, gr_meshComm, gr_meshNumProcs
use Particles_sinkData, ONLY: particles_local

use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use PhysicalConstants_interface, ONLY : PhysicalConstants_get
use Cosmology_interface, ONLY : Cosmology_getRedshift
use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getSingleCellCoords, &
                           Grid_getCellCoords, Grid_getBlkIndexLimits, &
                           Grid_getBlkPtr, Grid_releaseBlkPtr, &
                           Grid_putPointData, Grid_getBlkRefineLevel, &
                           Grid_getCellIndexFromPosition, &
                           Grid_notifySolnDataUpdate, Grid_getMinCellSize

use Eos_interface, ONLY : Eos_wrapped

use pt_sinkInterface, ONLY : pt_sinkCreateParticle, pt_sinkGatherGlobal

use sort
use sampleIMF
use gaus

implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

real, pointer, dimension(:,:,:,:)   :: solnData

! Properties for both the dense cell and star arrays. 
integer, save :: MASS=1, XLOC=2, YLOC=3, ZLOC=4, &
                 PROC=5, BLCK=6

! For just the dense cell array.
integer, save :: DENS=7, LAMJ=8

! For just the star array.
integer, save :: VELX=7, VELY=8, VELZ=9

! Here we assume some upper limit
! on the number of dense cells that meet the threshold for the star
! particle creation routine, say 20.

integer, parameter :: maxNumDenseCells = 20, NumDenseCellProps = 8, &
                      maxNumStars = 1000, NumStarProps = 7

real, parameter    :: eff = 1.0

real   :: denseCells(NumDenseCellProps,maxNumDenseCells), &
          recDenseCells(NumDenseCellProps,maxNumDenseCells), &
          stars(NumStarProps, maxNumStars), &
          recStars(NumStarProps, maxNumStars), numJeansCells

integer :: localNumDenseCells, totalNumDenseCells, n_stars, totalNumStars

integer :: ii, n, p, ierr, ind(3), lrefine_max, n_loops, pno
integer, dimension(gr_meshNumProcs) :: disp, recNumDenseCellsArray, &
                                       recNumStars

logical :: dense_cell, accept, place_by_dens, accept_v

real    :: cs, high_dens_thres, low_dens_thres, min_mass, r1, r_fac, &
           removed_dens, xr, yr, zr, blk_size(3)

real, dimension(:), allocatable    :: mass_sorted
integer, dimension(:), allocatable :: QSindex

integer :: returnProc, returnBlk

! Stuff for saving cell velocities and getting average and dispersions
real          :: min_cell_dx, jeans_rad, avg_vx, avg_vy, avg_vz, &
                 sig_vx, sig_vy, sig_vz, &
                 gaus_r, gaus_vx, gaus_vy, gaus_vz, vx_r, vy_r, vz_r
integer, save :: inside_count
real, dimension(:,:), allocatable  :: cell_vx, cell_vy, cell_vz, &
                                      recCellVx, recCellVy, recCellVz, &
                                      recNumVelCells
integer, dimension(:), allocatable :: numVelCells, totalNumVelCells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Stuff from Particles_sinkCreateAccrete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical, save       :: first_call=.true.
real, dimension(:), allocatable     :: xc, yc, zc

integer             :: blkRefine

real                :: redshift, onePlusRedshift, onePlusRedshift2, onePlusRedshift3
real                :: dVol, dx_block, dy_block, dz_block, rad


integer, save       :: MyPE, MasterPE, comm, idens, ipres, ivelx, ively, ivelz, igpot, ieint
integer, save       :: Nblockx, Nblocky, Nblockz, izn
integer             :: lb, llb, blk, blockID, i, j, ip, jp, kp, blockCount
integer             :: size_x, size_y, size_z
integer             :: i1, j1, k1, ii1, jj1, kk1, ncells_in_vol
integer             :: i_val, j_val, loc_i, iXcoord, iYcoord, iZcoord 

real                :: egrav_gas, egrav_par, vrad
logical             :: within_inner_r_accr, affected_block, create_sink


real, save          :: accretion_radius
real                :: accretion_radius_comoving, gpot, absgpot, comovingCellDens
real                :: density, density_thresh_comoving, gas_dens
real                :: distx, disty, distz

character(len=80), save :: grav_boundary_type

logical, save       :: GasAccretionChecks, convergingFlowCheck, potentialMinCheck, &
                       jeansCheck, negativeEtotCheck
integer, dimension(MAXBLOCKS) :: blockList
integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
integer, parameter  :: ngc_sink_creation_check = 2
real, parameter     :: ngc_sink_creation_check_radius_sqr = (ngc_sink_creation_check+1.0)**2

real, save          :: Newton
real                :: vxcm_in_vol, vycm_in_vol, vzcm_in_vol, ekindisp_in_vol, etherm_in_vol, emag_in_vol
real                :: r_search, mass_in_vol, maxgpot_in_vol, egravdeltapot_in_vol

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  real, save          :: mu_zero
  character(4), save  :: units
  integer, save       :: imagx, imagy, imagz
#endif



  if (first_call) then

    call RuntimeParameters_get("low_dens_thres", low_dens_thres)
    call RuntimeParameters_get("high_dens_thres", high_dens_thres)
    call RuntimeParameters_get("min_mass", min_mass)
    call RuntimeParameters_get("sink_accretion_radius", accretion_radius)
    call RuntimeParameters_get("sink_GasAccretionChecks", GasAccretionChecks)
    call RuntimeParameters_get("sink_convergingFlowCheck", convergingFlowCheck)
    call RuntimeParameters_get("sink_potentialMinCheck", potentialMinCheck)
    call RuntimeParameters_get("sink_jeansCheck", jeansCheck)
    call RuntimeParameters_get("sink_negativeEtotCheck", negativeEtotCheck)
    call RuntimeParameters_get("lrefine_max", lrefine_max)
    
    call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
    
    MyPE = gr_meshMe
    MasterPE = MASTER_PE

    comm = gr_meshComm

    call PhysicalConstants_get("Newton", Newton)

    idens = DENS_VAR
    ipres = PRES_VAR
    ivelx = VELX_VAR
    ively = VELY_VAR
    ivelz = VELZ_VAR
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
    imagx = MAGX_VAR
    imagy = MAGY_VAR
    imagz = MAGZ_VAR
#endif
    igpot = GPOT_VAR
    ieint = EINT_VAR

    iXcoord  = IAXIS
    iYcoord  = JAXIS
    iZcoord  = KAXIS
    izn = CENTER
    
    call RuntimeParameters_get("Nblockx", Nblockx)
    call RuntimeParameters_get("Nblocky", Nblocky)
    call RuntimeParameters_get("Nblockz", Nblockz)
    

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
    call RuntimeParameters_get("UnitSystem", units)
    if ( units == "SI" .or. units == "si" ) then
      mu_zero = 4.0*PI*1.e-7
    else if ( units == "CGS" .or. units == "cgs" ) then
      mu_zero = 4.0*PI
    else
      mu_zero = 1.0
    end if
#endif

    ! Calculate the maximum number of cells inside the Jeans' sphere.

    ! How many cells across is the Jeans length.
    call RuntimeParameters_get("jeans_ncells_ref", numJeansCells)
    call Grid_getMinCellSize(min_cell_dx)

    ! Minimum munber of cells to refine the Jeans length plus a little
    ! wiggle room.
    jeans_rad = numJeansCells*min_cell_dx*1.2

    inside_count = 0 
    call cells_in_sphere(jeans_rad, min_cell_dx, inside_count)

  end if

! Initialize to zero.
denseCells = 0.0
recDenseCells = 0.0
recNumDenseCellsArray = 0
localNumDenseCells = 0
totalNumDenseCells = 0

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  call Cosmology_getRedshift(redshift)
  onePlusRedshift = 1.0 + redshift
  onePlusRedshift2 = onePlusRedshift * onePlusRedshift
  onePlusRedshift3 = onePlusRedshift2 * onePlusRedshift

  ! Convert accretion_radius and density_thresh to comoving coordinates
  accretion_radius_comoving = accretion_radius * onePlusRedshift
  density_thresh_comoving = high_dens_thres / onePlusRedshift3
  
  call Grid_fillGuardCells(CENTER, ALLDIR, unitReadsMeshDataOnly=.true.)

! Jeans radius is set by the low density threshold from the parameter
! file and should only have to be calculated once if we make assumptions
! about the temperature of the gas.

! Could consider variable Jeans radii in the future for more flexible
! star formation.

  ! loop over leaf blocks (note that passed blockList only contains leafs)
  do lb = 1, blockCount

        blockID = blockList(lb)
        
        call Grid_getBlkRefineLevel(blockID, blkRefine)
        
        if (blkRefine == lrefine_max) then

        call Grid_getBlkPtr(blockID,solnData)

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

        allocate(xc(size_x))
        allocate(yc(size_y))
        allocate(zc(size_z))

        affected_block = .false.
        call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xc, size_x)
        call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yc, size_y)
        call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zc, size_z)
        
        call Grid_getBlkPhysicalSize(blockID,blk_size)
        dx_block = blk_size(1)/real(NXB)
        dy_block = blk_size(2)/real(NYB)
        dz_block = blk_size(3)/real(NZB)
        dVol = dx_block*dy_block*dz_block

        ! loop over cells (not including guard cells)
        do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef ALFL_MSCALAR
                 ! in case we are using the Alfven limiter
                 comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(idens,ip,jp,kp)
#else
                 comovingCellDens = solnData(DENS_VAR,ip,jp,kp)
#endif
                 if (comovingCellDens .gt. high_dens_thres) then

                    dense_cell = .true.   ! for now...
#ifdef debug_form                    
                    print*, "Proc", gr_meshMe, "found an over dense cell."
#endif
                    ! Is there an existing particle in range?
                    ! Looping over global particles
!                    do pno = 1, localnpf
!                       distx = xc(ip) - particles_global(ipx,pno)
!                       disty = yc(jp) - particles_global(ipy,pno)
!                       distz = zc(kp) - particles_global(ipz,pno)
!                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
!                       rad = sqrt(distx**2 + disty**2 + distz**2 )
!                       ! Does this position fall within 2x accretion radius of existing sink?
!                       if (rad .le. 2.0*accretion_radius_comoving) then
!                          dense_cell = .false.
!                       end if
!                    end do

                    ! Is there an over dense cell in range?
                    ! Looping over dense cells.
                    do n = 1, localNumDenseCells
                       distx = xc(ip) - denseCells(XLOC,n)
                       disty = yc(jp) - denseCells(YLOC,n)
                       distz = zc(kp) - denseCells(ZLOC,n)
                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
                       rad = sqrt(distx**2 + disty**2 + distz**2 )
                       ! Does this position fall within 2x accretion radius of existing sink?
                       if (rad .le. 2.0*accretion_radius_comoving) then
                          dense_cell = .false.
                       end if
                    end do

                    if (dense_cell) then

                       ! just use the isothermal sound speed as an estimate for the v_rad-check
                       cs = sqrt(solnData(PRES_VAR,ip,jp,kp)/comovingCellDens)

                       ! check for converging flow in all surrounding cells
                       if (convergingFlowCheck) then
                          do k1 = -1, 1
                             do j1 = -1, 1
                                do i1 = -1, 1

                                   rad = sqrt(real(i1*i1+j1*j1+k1*k1))
                                   if (rad .GT. 0.) then
                                      vrad = ( i1*(solnData(VELX_VAR, ip+i1, jp+j1, kp+k1)-solnData(VELX_VAR, ip, jp, kp)) + &
                                               j1*(solnData(VELY_VAR, ip+i1, jp+j1, kp+k1)-solnData(VELY_VAR, ip, jp, kp)) + &
                                               k1*(solnData(VELZ_VAR, ip+i1, jp+j1, kp+k1)-solnData(VELZ_VAR, ip, jp, kp)) ) / rad
                                      if (vrad .GT. 1.e-5*cs) then ! a surrounding cell diverges, so do not create sink
                                         dense_cell = .false.
#ifdef debug_form
                                         print*, "Proc", gr_meshMe, "rejected cell because no converging flow."
#endif
                                      endif
                                   endif ! rad > 0

                                enddo
                             enddo
                          enddo
                       end if

                       ! check for potential minimum
                       if (potentialMinCheck) then
                          if (dense_cell) then

                             ncells_in_vol = 1
                             gpot = solnData(GPOT_VAR,ip,jp,kp)
                             absgpot = abs(gpot)
                             do k1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = -ngc_sink_creation_check, ngc_sink_creation_check

                                      if (i1**2 + j1**2 + k1**2 .le. ngc_sink_creation_check_radius_sqr) then

                                         if ( ((solnData(GPOT_VAR,ip+i1,jp+j1,kp+k1)-gpot)/absgpot) .lt. -1.e-5) then
                                            dense_cell = .false.
#ifdef debug_form
                                            print*, "Proc", gr_meshMe, "rejected cell because no potential min."
#endif
                                         end if

                                         ncells_in_vol = ncells_in_vol + 1

                                      end if

                                   end do
                                end do
                             end do

                          end if
                       end if

                       ! check for Jeans condition and total energies
                       if (jeansCheck .or. negativeEtotCheck) then
                          if (dense_cell) then

                             etherm_in_vol = 0.0
                             vxcm_in_vol = 0.0
                             vycm_in_vol = 0.0
                             vzcm_in_vol = 0.0
                             mass_in_vol = 0.0
                             maxgpot_in_vol = solnData(GPOT_VAR,ip,jp,kp)
                             do k1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = -ngc_sink_creation_check, ngc_sink_creation_check

                                      if (i1**2 + j1**2 + k1**2 .le. ngc_sink_creation_check_radius_sqr) then

                                         ii1 = ip+i1
                                         jj1 = jp+j1
                                         kk1 = kp+k1
#ifdef ALFL_MSCALAR
                                         ! in case we are using the Alfven limiter
                                         gas_dens = (1.0-solnData(ALFL_MSCALAR,ii1,jj1,kk1)) * &
                                                    solnData(DENS_VAR,ii1,jj1,kk1)
#else
                                         gas_dens = solnData(DENS_VAR,ii1,jj1,kk1)
#endif
                                         etherm_in_vol = etherm_in_vol + solnData(EINT_VAR,ii1,jj1,kk1)*gas_dens 

                                         vxcm_in_vol = vxcm_in_vol + solnData(VELX_VAR,ii1,jj1,kk1)*gas_dens
                                         vycm_in_vol = vycm_in_vol + solnData(VELY_VAR,ii1,jj1,kk1)*gas_dens
                                         vzcm_in_vol = vzcm_in_vol + solnData(VELZ_VAR,ii1,jj1,kk1)*gas_dens

                                         mass_in_vol = mass_in_vol + gas_dens

                                         if (solnData(igpot,ii1,jj1,kk1) .gt. maxgpot_in_vol) & 
                                              maxgpot_in_vol = solnData(GPOT_VAR,ii1,jj1,kk1)

                                      end if

                                   end do
                                end do
                             end do

                             etherm_in_vol = etherm_in_vol*dVol
                             vxcm_in_vol = vxcm_in_vol/mass_in_vol
                             vycm_in_vol = vycm_in_vol/mass_in_vol
                             vzcm_in_vol = vzcm_in_vol/mass_in_vol
                             mass_in_vol = mass_in_vol*dVol

                             ekindisp_in_vol = 0.0
                             egravdeltapot_in_vol = 0.0
                             emag_in_vol = 0.0

                             do k1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = -ngc_sink_creation_check, ngc_sink_creation_check

                                      if (i1**2 + j1**2 + k1**2 .LE. ngc_sink_creation_check_radius_sqr) then

                                         ii1 = ip+i1
                                         jj1 = jp+j1
                                         kk1 = kp+k1
#ifdef ALFL_MSCALAR
                                         ! in case we are using the Alfven limiter
                                         gas_dens = (1.0-solnData(ALFL_MSCALAR,ii1,jj1,kk1)) * &
                                                    solnData(DENS_VAR,ii1,jj1,kk1)
#else
                                         gas_dens = solnData(DENS_VAR,ii1,jj1,kk1)
#endif
                                         ekindisp_in_vol = ekindisp_in_vol + gas_dens * &
                                              ( (solnData(VELX_VAR, ii1, jj1, kk1) - vxcm_in_vol)**2 + &
                                                (solnData(VELY_VAR, ii1, jj1, kk1) - vycm_in_vol)**2 + &
                                                (solnData(VELZ_VAR, ii1, jj1, kk1) - vzcm_in_vol)**2  )

                                         egravdeltapot_in_vol = egravdeltapot_in_vol + & 
                                              (solnData(GPOT_VAR, ii1, jj1, kk1) - maxgpot_in_vol)*gas_dens

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
                                         emag_in_vol = emag_in_vol + solnData(MAGX_VAR, ii1, jj1, kk1)**2 + &
                                                                     solnData(MAGY_VAR, ii1, jj1, kk1)**2 + &
                                                                     solnData(MAGZ_VAR, ii1, jj1, kk1)**2
#endif
                                      end if

                                   end do
                                end do
                             end do

                             ekindisp_in_vol = 0.5*ekindisp_in_vol*dVol
                             egravdeltapot_in_vol = -egravdeltapot_in_vol*dVol
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
                             emag_in_vol = 0.5/mu_zero*emag_in_vol*dVol
#endif
                             ! Jeans mass virial argument (see e.g., Bate Bonnell Price 1995)
                             if (jeansCheck) then
                                if (2.0*etherm_in_vol + emag_in_vol .GT. egravdeltapot_in_vol) then
                                   dense_cell = .false.
#ifdef debug_form
                                   print*, "Proc", gr_meshMe, "rejected cell because failed Jeans mass virial check."
#endif
                                endif
                             end if

                             ! total energy should be negative (see e.g., Bate Bonnell Price 1995)
                             if (negativeEtotCheck) then
                                if (dense_cell) then
                                   if (etherm_in_vol + ekindisp_in_vol + emag_in_vol .GT. egravdeltapot_in_vol) then
                                      dense_cell = .false.
#ifdef debug_form
                                      print*, "Proc", gr_meshMe, "rejected cell because failed total E negative check."
#endif
                                   endif
                                endif
                             endif

                          end if ! energy check

                       end if

                       ! finally create the sink in the cell centre
                       if (dense_cell) then
                          print*, "Proc", gr_meshMe, "kept a dense cell. Storing." 
                                                
                          ! Make sure we update the variables with the new density.
                          
                          ! Now make the list of high_dens cells locs and densities.
                          
                          localNumDenseCells = localNumDenseCells + 1
                          
                          denseCells(DENS, localNumDenseCells) = comovingCellDens
                          denseCells(XLOC, localNumDenseCells) = xc(ip)
                          denseCells(YLOC, localNumDenseCells) = yc(jp)
                          denseCells(ZLOC, localNumDenseCells) = zc(kp)
                          
                          ! Calculate Jeans diameter.
                          denseCells(LAMJ, localNumDenseCells) = sqrt((PI * cs**2.0) / (Newton * low_dens_thres))
                          denseCells(PROC, localNumDenseCells) = gr_meshMe
                          denseCells(BLCK, localNumDenseCells) = blockID
                          
                       end if

                    end if

                    

                 end if    ! cell density > max gas density

              end do
           end do
        end do

!        if (affected_block) then
!           llb = llb+1
!           block_list(llb) = blockID
!        end if

        call Grid_releaseBlkPtr(blockID, solnData)

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)
        
        end if ! If block at max refinement.

  end do ! block loop
  
  ! ======== END find high dens cells loop ========
  
  ! Now gather the denseCells array on all processors using the vector
  ! gather commands of MPI.
  
  ! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array.
  
  ! Make an array of the # of dense cells from each processor.
  call MPI_AllGather(localNumDenseCells, 1, MPI_INTEGER, &
                     recNumDenseCellsArray, 1, MPI_INTEGER, &
                     gr_meshComm, ierr)
                  
  !ret_cnt(:dr_globalNumProcs) = rec_count      
            
  ! Set the displacement for the incoming data based on how many
  ! cells are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.
  
  !write(*,*) "received count on proc ", gr_meshMe, " = ", recNumDenseCellsArray

  totalNumDenseCells = sum(recNumDenseCellsArray)
  
  ! Are there any cells that met the star formation criterion? If not
  ! return.
  if (totalNumDenseCells == 0) then
    if (gr_meshMe == 0) write(*,*) "Total number of dense cells = 0, no &
                                    star formation."
    return
  end if

  disp = 0
  do i=1, gr_meshNumProcs-1
  
    disp(i+1) = disp(i) + recNumDenseCellsArray(i)
    
  end do

  
  ! Now actually gather the cells on all procs using the variable length array
  ! AllGatherV command in MPI. This command can gather different length
  ! arrays together into one array on all processes safely, and it can
  ! do it in place. Nice! -JW
  
  do n=1, NumDenseCellProps
    
    call MPI_AllGatherV(denseCells(n,:), localNumDenseCells, MPI_DOUBLE_PRECISION, &
                        recDenseCells(n,:), recNumDenseCellsArray, disp, &
                        MPI_DOUBLE_PRECISION, gr_meshComm, ierr)
                        
  end do

  denseCells = recDenseCells
  
  ! Now lets gather up cell mass to make stars with.
  
  ! We have to keep track of cell velocities to use to set star
  ! velocities later.
  
  allocate(cell_vx(totalNumDenseCells,inside_count))
  allocate(cell_vy(totalNumDenseCells,inside_count))
  allocate(cell_vz(totalNumDenseCells,inside_count))

  allocate(recCellVx(totalNumDenseCells,inside_count))
  allocate(recCellVy(totalNumDenseCells,inside_count))
  allocate(recCellVz(totalNumDenseCells,inside_count))

  allocate(numVelCells(totalNumDenseCells))
  allocate(recNumVelCells(totalNumDenseCells, gr_meshNumProcs))
  allocate(totalNumVelCells(totalNumDenseCells))
  
  cell_vx=0.0; cell_vy=0.0; cell_vz=0.0
  recCellVx=0.0; recCellVy=0.0; recCellVz=0.0
  numVelCells=0; recNumVelCells=0; totalNumVelCells=0 

  ! loop over leaf blocks (note that passed blockList only contains leafs)
  do lb = 1, blockCount

    blockID = blockList(lb)
        
    call Grid_getBlkRefineLevel(blockID, blkRefine)
        
    if (blkRefine == lrefine_max) then
        
        call Grid_getBlkPhysicalSize(blockID,blk_size)
        
        dx_block = blk_size(1)/real(NXB)
        dy_block = blk_size(2)/real(NYB)
        dz_block = blk_size(3)/real(NZB)
        dVol = dx_block*dy_block*dz_block

        call Grid_getBlkPtr(blockID,solnData)

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

        allocate(xc(size_x))
        allocate(yc(size_y))
        allocate(zc(size_z))

        affected_block = .false.
        call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xc, size_x)
        call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yc, size_y)
        call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zc, size_z)



        ! loop over cells (not including guard cells)
        do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
            do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                
              ! Set the probability scratch variable to zero for all cells.
              ind = [ip, jp, kp]
              !call Grid_putPointData(blockID, SCRATCH_CTR, PROB_SCRATCH_CENTER_VAR, &
              !                       INTERIOR, ind, 0.0) 
              
              do n=1, totalNumDenseCells


#ifdef ALFL_MSCALAR
                 ! in case we are using the Alfven limiter
                 comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(DENS_VAR,ip,jp,kp)
#else
                 comovingCellDens = solnData(DENS_VAR,ip,jp,kp)
#endif
                  
!                  if (comovingCellDens .gt. low_dens_thres) then
                  
                    !!! NOTE: This should be corrected for PBCs at some point.
                    r1 = (xc(ip)-denseCells(XLOC,n))**2.0 &
                       + (yc(jp)-denseCells(YLOC,n))**2.0 &
                       + (zc(kp)-denseCells(ZLOC,n))**2.0
                    r1 = sqrt(r1)
                    
                    if (r1 .le. denseCells(LAMJ, n)) then
                    
                    !print*, "Inside lambda_j check."
                         
                    ! Note this is basically volume weighted harvesting of the gas. So we might
                    ! consider doing volume weighted placement of the stars, i.e. the most
                    ! massive 15% get placed in the inner 0.4*jeans_rad, then the next
                    ! 0.15*max_mass < mass < 0.55*max_mass goes between 0.4*jeans_rad < rad < 0.8*jeans_rad,
                    ! finally the rest is outside of 0.8*jeans_rad. - JW
                    
#ifdef ALFL_MSCALAR
                        ! in case we are using the Alfven limiter
                        alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(idens,ip,jp,kp)
                        
                        if (r1 .lt. 0.4*denseCells(LAMJ, n)) then
                        
                            removed_dens = eff*0.85*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.85)*(comovingCellDens - low_dens_thres) + low_dens_thres + alfl_dens)
                            
                        else if ((r1 .ge. 0.4*denseCells(LAMJ, n)) &
                                 .and. (r1 .lt. 0.8*denseCells(LAMJ, n))) then
                                 
                            removed_dens = eff*0.55*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.45)*(comovingCellDens - low_dens_thres) + low_dens_thres + alfl_dens)
                            
                        else
                        
                            removed_dens = eff*0.15*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.85)*(comovingCellDens - low_dens_thres) + low_dens_thres + alfl_dens)
                            
                        end if
                                                
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving + alfl_dens
                        solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(idens,ip,jp,kp)
#else
                        if (r1 .lt. 0.4*denseCells(LAMJ, n)) then
                        
                            removed_dens = eff*0.85*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.15)*(comovingCellDens - low_dens_thres) + low_dens_thres
                            
                        else if ((r1 .ge. 0.4*denseCells(LAMJ, n)) &
                                 .and. (r1 .lt. 0.8*denseCells(LAMJ, n))) then
                                 
                            removed_dens = eff*0.55*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.45)*(comovingCellDens - low_dens_thres) + low_dens_thres
                        
                        else
                        
                            removed_dens = eff*0.15*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.85)*(comovingCellDens - low_dens_thres) + low_dens_thres
                            
                        end if
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving
#endif

                        denseCells(MASS,n) = denseCells(MASS,n) + removed_dens*dVol
                        
                        !print*, "After mass was removed."

                      affected_block = .true. ! If we changed dens, we need to call EOS.
                      
                      ! Get the cell velocitites so we can calculate the
                      ! average velocity and the dispersion.
                      
                      numVelCells(n) = numVelCells(n) + 1
                      cell_vx(n, numVelCells(n)) = solnData(VELX_VAR,ip,jp,kp)
                      cell_vy(n, numVelCells(n)) = solnData(VELY_VAR,ip,jp,kp)
                      cell_vz(n, numVelCells(n)) = solnData(VELZ_VAR,ip,jp,kp)
                      
                        
                    end if ! low_jeans radius criterion
                      
!                  end if !  low thres criterion
                  
                end do  ! high threshold dens loop  

            end do  ! i loop
          end do  ! j look
        end do  ! k loop
        
        !print*, "Outside cell loops."
        
        call Grid_releaseBlkPtr(lb, solnData)
      
        deallocate(xc)
        deallocate(yc)
        deallocate(zc)
      
        if(affected_block) then
          call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)
        endif
      
        !print*, "After EOS_wrapped."
    end if ! if block at max refinement

  end do  ! block loop
  
  ! Reduce the gathered mass from each processor to its proper array slot
  ! to make stars on the processor that the original high dense cell was on.
  
  !print*, "Before AllReduce denseCells(MASS,:)."
  
  call MPI_ALLREDUCE(MPI_IN_PLACE, denseCells(MASS,:), maxNumDenseCells, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)

  
  ! Make an array of the # of cell velocities from each processor.
  do n=1, totalNumDenseCells

    call MPI_AllGather(numVelCells(n), 1, MPI_INTEGER, &
                       recNumVelCells(n,:), 1, MPI_INTEGER, &
                       gr_meshComm, ierr)

    disp = 0
    do i=1, gr_meshNumProcs-1
  
      disp(i+1) = disp(i) + recNumVelCells(n,i)
    
    end do
    
    totalNumVelCells(n) = sum(recNumVelCells(n,:))

  ! Now actually gather the cell velocities on all procs using the variable length array
  ! AllGatherV command in MPI.
  
    call MPI_AllGatherV(cell_vx(n,:), numVelCells(n), MPI_DOUBLE_PRECISION, &
                        recCellVx(n,:), recNumVelCells(n,:), disp, &
                        MPI_DOUBLE_PRECISION, gr_meshComm, ierr)
    
    call MPI_AllGatherV(cell_vy(n,:), numVelCells(n), MPI_DOUBLE_PRECISION, &
                        recCellVy(n,:), recNumVelCells(n,:), disp, &
                        MPI_DOUBLE_PRECISION, gr_meshComm, ierr)
    
    call MPI_AllGatherV(cell_vz(n,:), numVelCells(n), MPI_DOUBLE_PRECISION, &
                        recCellVz(n,:), recNumVelCells(n,:), disp, &
                        MPI_DOUBLE_PRECISION, gr_meshComm, ierr)
                        
  end do


  print*, "Going to try and make stars on proc ", gr_meshMe


  do n=1, totalNumDenseCells
  
    stars         = 0.0
    n_stars       = 0
    totalNumStars = 0
    recStars      = 0.0
    recNumStars   = 0
                          
      ! convert to solar masses to sample the IMF.
      
      denseCells(MASS,n) = denseCells(MASS,n) / 1.989e33
                
    if (denseCells(MASS,n) .ge. min_mass) then
      
        if (denseCells(PROC,n) == gr_meshMe) then
      
          print*, "Sampling the IMF with ", denseCells(MASS,n), &
                  " solar masses on proc ", gr_meshMe
        
          call sample_IMF(denseCells(MASS,n), min_mass, &
                          stars(MASS, :), n_stars)
                          
          totalNumStars = totalNumStars + n_stars
          
          print*, "After star formation, total mass left is =", denseCells(MASS,n)
          
          ! Convert back to cgs.
          ! Do I need this anymore?! - JW
          denseCells(MASS,n) = denseCells(MASS,n) * 1.989e33
          
          ! Sort the star masses from smallest to largest.
          stars(MASS,:) = stars(MASS,:) * 1.989e33 ! Put back in grams.
                    
          print*, "We got ", n_stars, " with total mass ", & 
                   sum(stars(MASS,:))/1.989e33, "on proc ", gr_meshMe
                   
          ! Calculate the average velocities and velocity dispersions
          ! of the gas.
          
          print*, "recCellVx=", recCellVx(n,:)
          print*, "totalNumVelCells=", totalNumVelCells(n)
          
          avg_vx = sum(recCellVx(n,:)) / totalNumVelCells(n)
          avg_vy = sum(recCellVy(n,:)) / totalNumVelCells(n)
          avg_vz = sum(recCellVz(n,:)) / totalNumVelCells(n)
          
          sig_vx = sqrt( sum( (recCellVx(n,:) - avg_vx)**2.0 ) / totalNumVelCells(n) )
          sig_vy = sqrt( sum( (recCellVy(n,:) - avg_vy)**2.0 ) / totalNumVelCells(n) )
          sig_vz = sqrt( sum( (recCellVz(n,:) - avg_vz)**2.0 ) / totalNumVelCells(n) )
          
          print*, "avg_vx=", avg_vx
          print*, "sig_vx", sig_vx
          
          ! Get a dictionary to the stars sorted by mass.
          
          allocate(mass_sorted(n_stars))
          allocate(QSindex(n_stars))
          
          do p=1, n_stars
            mass_sorted(p) = stars(MASS,p)
          end do
          
          call sortInd(mass_sorted, QSindex)
        
          do i=1, n_stars
      
            accept = .false.
            n_loops = 0
          
            !dist_from_mean = stars(i) / mean_stars(i_val)
      
            do while ((.not. accept) .and. (n_loops .lt. 1e4))
          
              if (place_by_dens) then
            
              else ! Place randomly.
            
                ! For random placement, we "loosely" follow the highest
                ! density with the most massive stars by expanding the 
                ! the radius of placement with decreasing mass. - JW
                
                ! Note stars are sorted by ascending mass.
                
                if (real(i)/real(n_stars) .lt. 1.0/3.0) then
                    r_fac = 1.0
                else if (real(i)/real(n_stars) .gt. 2.0/3.0) then
                    r_fac = 0.4
                else
                    r_fac = 0.8
                end if
          
                ! Get a random location within the Jeans radius.
                call random_number(xr)
                !print*, "Random x =", xr
                xr = (xr - 0.5)*2.0*r_fac*denseCells(LAMJ,n)
                call random_number(yr)
                yr = (yr - 0.5)*2.0*r_fac*denseCells(LAMJ,n)
                call random_number(zr)
                zr = (zr - 0.5)*2.0*r_fac*denseCells(LAMJ,n)
                
                if ( sqrt(xr**2.0 + yr**2.0 + zr**2.0) .gt. denseCells(LAMJ,n)) cycle
                
                stars(XLOC, QSindex(i)) = xr + denseCells(XLOC,n)
                stars(YLOC, QSindex(i)) = yr + denseCells(YLOC,n)
                stars(ZLOC, QSindex(i)) = zr + denseCells(ZLOC,n)
                
                ! Choose and compare random velocities within one sigma.
                
                accept_v = .false.
                do while (.not. accept_v)
                    call random_number(vx_r)
                    vx_r = vx_r*sig_vx
                    call gaussian(vx_r, gaus_vx, avg_vx, sig_vx)
                    call random_number(gaus_r)
                    if (gaus_vx/gaus_r .gt. 1.0) then
                        stars(VELX, QSindex(i)) = vx_r
                        accept_v = .true.
                    end if
                    
                    print*, "vx_r=", vx_r
                    print*, "avg_vx=", avg_vx
                    print*, "gaus_r=", gaus_r
                    print*, "gaus_vx=", gaus_vx
                    print*, "gaus_vx/gaus_r=", gaus_vx/gaus_r
                    call flush(6)
                end do
                
                accept_v = .false.
                do while (.not. accept_v)
                    call random_number(vy_r)
                    vy_r = vy_r*sig_vy
                    call gaussian(vy_r, gaus_vy, avg_vy, sig_vy)
                    call random_number(gaus_r)
                    if (gaus_vy/gaus_r .gt. 1.0) then
                        stars(VELY, QSindex(i)) = vy_r
                        accept_v = .true.
                    end if
                end do
                
                accept_v = .false.
                do while (.not. accept_v)
                    call random_number(vz_r)
                    vz_r = vz_r*sig_vz
                    call gaussian(vz_r, gaus_vz, avg_vz, sig_vz)
                    call random_number(gaus_r)
                    if (gaus_vz/gaus_r .gt. 1.0) then
                        stars(VELZ, QSindex(i)) = vx_r
                        accept_v = .true.
                    end if
                end do
                
                accept = .true.
                n_loops = n_loops + 1

              end if ! If place by density.
            
            end do !end do while .not. accept
          
          end do ! loop over n_stars
        
          deallocate(mass_sorted)
          deallocate(QSindex)

        end if ! If denseCells(PROC == gr_meshMe) ! Because everyone has to call
        
        ! Now lets place and exchange stars.
        
        ! Now gather the stars array on all processors using the vector
        ! gather commands of MPI.
        
        ! Here's how many total stars we made for this dense cell across
        ! all processors. We'll use this to place stars.
        call MPI_AllReduce(MPI_IN_PLACE, totalNumStars, 1, MPI_INTEGER, MPI_SUM, &
                           gr_meshComm, ierr)
          
        ! Gather the array on the root process. Note that we require the
        ! user to pass the proper length of the final array.
          
        ! Make an array of the # of stars from each processor.
        call MPI_AllGather(n_stars, 1, MPI_INTEGER, &
                           recNumStars, 1, MPI_INTEGER, &
                           gr_meshComm, ierr)
                          
        ! Set the displacement for the incoming data based on how many
        ! cells are coming in from each processor. Note the displacement
        ! for the root process is zero, for rank 1 disp = num on root,
        ! for rank 2 disp = num on root + num on 1, etc etc.
          
        write(*,*) "received count on proc ", gr_meshMe, " = ", recNumStars         
         
        disp = 0
        
        do i=1, gr_meshNumProcs-1
          
          disp(i+1) = disp(i) + recNumstars(i)
            
        end do

          
          ! Now actually gather the stars on all procs using the variable length array
          ! AllGatherV command in MPI. This command can gather different length
          ! arrays together into one array on all processes safely, and it can
          ! do it in place. Nice! -JW
          
        do i=1, NumDenseCellProps
            
          call MPI_AllGatherV(stars(i,:), n_stars, MPI_DOUBLE_PRECISION, &
                              recStars(i,:), recNumStars, disp, &
                              MPI_DOUBLE_PRECISION, gr_meshComm, ierr)
                                
        end do
          
          
        totalNumStars = sum(recNumStars)
        stars = recStars
        
        print*, "Total number of stars =", totalNumStars, "on proc", gr_meshMe
        
        !if (gr_meshMe == 0) print*, "Star masses = ", stars(MASS, 1:totalNumStars)/1.989e33
            
        ! Now that all procs have all the stars, including there locations,
        ! we can safely place them on each processor correctly.
        
        ! Now lets count how many each proc places. Should sum to the
        ! original total.
        n_stars = 0
        
        do i=1, totalNumStars
            
          call Grid_getCellIndexFromPosition(returnProc, returnBlk, &
                                             ip, jp, kp, &
                                             stars(XLOC,i), &
                                             stars(YLOC,i), &
                                             stars(ZLOC,i))
          stars(PROC,i) = real(returnProc)
          stars(BLCK,i) = real(returnBlk)
          if (int(stars(PROC,i)) == gr_meshMe) then

              n_stars = n_stars + 1
              call Grid_getBlkPtr(int(stars(BLCK,i)), solnData)
              pno = pt_sinkCreateParticle(stars(XLOC,i), stars(YLOC,i), stars(ZLOC,i), &
                                        dr_simTime, int(stars(BLCK,i)), int(stars(PROC,i)))
              particles_local(MASS_PART_PROP, pno) = stars(MASS,i)
              particles_local(VELX_PART_PROP, pno) = stars(VELX, i)
              particles_local(VELY_PART_PROP, pno) = stars(VELY, i)
              particles_local(VELZ_PART_PROP, pno) = stars(VELZ, i)
            
              call Grid_releaseBlkPtr(blockID, solnData)
            end if ! If proc == gr_meshMe
          
        end do ! Loop over all stars.
        
        print*, "Proc ", gr_meshMe, " placed ", n_stars, " stars."

    end if ! if mass > min_mass
    
  end do  ! Loop over all n_dens_cells
  
  deallocate(cell_vx)
  deallocate(cell_vy)
  deallocate(cell_vz)

  deallocate(recCellVx)
  deallocate(recCellVy)
  deallocate(recCellVz)

  deallocate(numVelCells)
  deallocate(totalNumVelCells)
  
  call pt_sinkGatherGlobal()


contains

subroutine cells_in_sphere(rad, dx, inside_count)

! Get the number of cells inside a sphere assuming
! dx=dy=dz and the center of a cell and the sphere
! coincide at x=y=z=0.

  implicit none

  integer, parameter          :: dp = kind(1.d0)

  real(dp), intent(in)        :: rad, dx
  integer, intent(inout)      :: inside_count

  real(dp)                    :: r, factor
  real(dp), allocatable       :: xsq(:), ysq(:), zsq(:)
  integer                     :: i, j, k, nsteps

! Get the deltas and positions to sample across the sphere.
  nsteps = int(2.0*rad/dx)
  
  allocate(xsq(nsteps))
  allocate(ysq(nsteps))
  allocate(zsq(nsteps))
  
  do i=1, nsteps

    xsq(i) = (-real(nsteps/2)*dx + (0.5_dp+real(i-1))*dx)**2.0_dp
    
  end do
  
  ysq = xsq
  zsq = xsq
  
! Now sample these to see if they are inside the object.
  inside_count = 0

  do i=1, nsteps
    do j=1, nsteps
      do k=1, nsteps
      
        r = sqrt(xsq(i) + ysq(j) + zsq(k))
        if (r .le. rad) inside_count = inside_count + 1

      end do
    end do
  end do
  
  deallocate(xsq)
  deallocate(ysq)
  deallocate(zsq)

end subroutine cells_in_sphere


end subroutine Particles_starForm
