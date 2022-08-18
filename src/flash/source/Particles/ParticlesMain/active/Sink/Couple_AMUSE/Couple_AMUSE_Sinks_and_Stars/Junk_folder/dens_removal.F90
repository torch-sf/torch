

!#define debug_form



subroutine dens_removal(made_stars)

use Driver_data, ONLY : dr_simTime
use Grid_data, ONLY   : gr_meshNumProcs, gr_meshMe, gr_meshComm
use Driver_interface, ONLY : Driver_abortFlash
use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use PhysicalConstants_interface, ONLY : PhysicalConstants_get
use Cosmology_interface, ONLY : Cosmology_getRedshift
use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getSingleCellCoords, &
                           Grid_getCellCoords, Grid_getBlkIndexLimits, &
                           Grid_getBlkPtr, Grid_releaseBlkPtr, &
                           Grid_putPointData, Grid_getBlkRefineLevel, &
                           Grid_getCellIndexFromPosition, Grid_notifySolnDataUpdate
use Eos_interface, ONLY : Eos_wrapped

use pt_sinkInterface, ONLY : pt_sinkCreateParticle, pt_sinkGatherGlobal

use ut_qsortInterface, ONLY : ut_qsort

use pt_sinkSort, ONLY : NewQsort_IN

use Particles_sinkData

use sampleIMF
use gaus
use qsort2D

implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

logical, intent(out) :: made_stars

real, dimension(:,:), allocatable   :: dens_matrix
real, dimension(:), allocatable     :: xc, yc, zc, proc_array
real, dimension(:), allocatable     :: x_high, y_high, z_high, dens_high ! For the high dens thres cells
real, pointer, dimension(:,:,:,:)   :: solnData

real, dimension(:), allocatable     :: stars, stars_offproc, starsproc_offproc, &
                                       starsx_offproc, starsy_offproc, starsz_offproc !, dens_gaus_pdf, star_gaus_pdf
real, dimension(:), allocatable     :: mass, jeans_rad, tot_star_mass
integer, dimension(:), allocatable  :: n_save, QSindex
integer                             :: starsnum_offproc, idest, isrc, comm_stat
real, dimension(:), allocatable     :: mean_dens, var_dens, mean_stars, &
                                       var_stars 

logical, parameter :: place_by_dens = .false.

! For storing the information on all cells where mass is gathered.
integer, dimension(:), allocatable :: cell_count
real, dimension(:,:), allocatable  :: dens_cells, dens_blk, dens_x, dens_y, dens_z

real                :: time, loc !, low_jeans_diameter !, mass
real                :: rnum, r1, prob, mu, var, cs
real                :: xr, yr, zr, dist_from_mean
real                :: coords(3)
integer             :: ind(3)
real, save          :: low_dens_thres, high_dens_thres, min_mass, eff
logical, save       :: first_call=.true.
integer, save       :: array_length, array_start, start_ind, end_ind, nproc
integer, save       :: lrefine_max
integer             :: tot_n_cells, n_stars, n_loops, n !, n_save
integer             :: n_cells_high(gr_meshNumProcs), ierr, dens_loc, pno
integer             :: rind, blkRefine, proc

logical             :: overlap, accept, offproc

character(len=12)   :: fout

real                :: blk_size(3)
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
    
    eff = 0.5

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
    !call Driver_getComm(GLOBAL_COMM, comm)
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

    
!    allocate(proc_array(nproc))
    
!    do i=1, nproc
!      start_ind = (i-1)*10+1
!      end_ind   = i*10
!      proc_array(start_ind:end_ind) = i-1
!    end do
    
    if (MyPE == MASTER_PE) then
      print*, "Star formation turned on!"
      print*, "High density thresh = ", high_dens_thres
      print*, "Low density thresh = ", low_dens_thres
    end if

!    call init_random_seed()

    first_call = .false.
    
  end if

  made_stars = .false.

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  call Cosmology_getRedshift(redshift)
  onePlusRedshift = 1.0 + redshift
  onePlusRedshift2 = onePlusRedshift * onePlusRedshift
  onePlusRedshift3 = onePlusRedshift2 * onePlusRedshift

  ! Convert accretion_radius and density_thresh to comoving coordinates
  accretion_radius_comoving = accretion_radius * onePlusRedshift
  density_thresh_comoving = high_dens_thres / onePlusRedshift3
  
  call Grid_fillGuardCells(CENTER, ALLDIR, unitReadsMeshDataOnly=.true.)
  
  ! ======== START find high dens cells loop ========

  nproc = gr_meshNumProcs
  array_length = (nproc)*10 ! These arrays allow for 10 dense cells on each proc
  array_start  = MyPE*10 + 1 !(nproc-1)*10 + 1 

  ! Allocate for storing the information on all cells that mass is
  ! gathered from.
  
  if (place_by_dens) then
  allocate(dens_cells(nproc, 10000))
  !allocate(dens_blk(nproc, 10000))
  allocate(dens_x(nproc, 10000))
  allocate(dens_y(nproc, 10000))
  allocate(dens_z(nproc, 10000))
  end if

  ! Allocate for storing the information on the cells that exceed the
  ! high density threshold.

  allocate(x_high(array_length))
  allocate(y_high(array_length))
  allocate(z_high(array_length))
  allocate(dens_high(array_length))
  allocate(jeans_rad(array_length))
  
  n_cells_high = 0
  x_high       = 0.0
  y_high       = 0.0
  z_high       = 0.0
  dens_high    = 0.0
  jeans_rad    = 0.0
  
  !mass  = 0.0
  time  = dr_simTime

  llb   = 0

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

                    create_sink = .true.   ! for now...
#ifdef debug_form                    
                    print*, "Proc", MyPE, "found an over dense cell."
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
!                          create_sink = .false.
!                       end if
!                    end do

                    ! Is there an over dense cell in range?
                    ! Looping over dense cells.
                    do n = 1, n_cells_high(MyPe+1)
                       distx = xc(ip) - x_high(n)
                       disty = yc(jp) - y_high(n)
                       distz = zc(kp) - z_high(n)
                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
                       rad = sqrt(distx**2 + disty**2 + distz**2 )
                       ! Does this position fall within 2x accretion radius of existing sink?
                       if (rad .le. 2.0*accretion_radius_comoving) then
                          create_sink = .false.
                       end if
                    end do

                    if (create_sink) then

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
                                         create_sink = .false.
#ifdef debug_form
                                         print*, "Proc", MyPE, "rejected cell because no converging flow."
#endif
                                      endif
                                   endif ! rad > 0

                                enddo
                             enddo
                          enddo
                       end if

                       ! check for potential minimum
                       if (potentialMinCheck) then
                          if (create_sink) then

                             ncells_in_vol = 1
                             gpot = solnData(GPOT_VAR,ip,jp,kp)
                             absgpot = abs(gpot)
                             do k1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = -ngc_sink_creation_check, ngc_sink_creation_check

                                      if (i1**2 + j1**2 + k1**2 .le. ngc_sink_creation_check_radius_sqr) then

                                         if ( ((solnData(GPOT_VAR,ip+i1,jp+j1,kp+k1)-gpot)/absgpot) .lt. -1.e-5) then
                                            create_sink = .false.
#ifdef debug_form
                                            print*, "Proc", MyPE, "rejected cell because no potential min."
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
                          if (create_sink) then

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
                                   create_sink = .false.
#ifdef debug_form
                                   print*, "Proc", MyPE, "rejected cell because failed Jeans criterion."
#endif
                                endif
                             end if

                             ! total energy should be negative (see e.g., Bate Bonnell Price 1995)
                             if (negativeEtotCheck) then
                                if (create_sink) then
                                   if (etherm_in_vol + ekindisp_in_vol + emag_in_vol .GT. egravdeltapot_in_vol) then
                                      create_sink = .false.
#ifdef debug_form
                                      print*, "Proc", MyPE, "rejected cell because failed total E negative check."
#endif
                                   endif
                                endif
                             endif

                          end if ! energy check

                       end if

                       ! finally create the sink in the cell centre
                       if (create_sink) then
                          print*, "Proc", MyPE, "kept a dense cell. Storing." 
                        
                          ! Calculate Jeans diameter. Then loop over blocks to
                          ! see if any of the block is contained within Jeans
                          ! diameter. Then loop over those cells in those blocks.
                        
                          ! Make sure we update the variables with the new density.
                          
                          ! Now make the list of high_dens cells locs and densities.
                          
                          n_cells_high(MyPE+1) = n_cells_high(MyPE+1) + 1
                          loc_i = array_start + n_cells_high(MyPE+1) - 1
                          
                          if (n_cells_high(MyPE+1) > 10) then
                            call Driver_abortFlash("Number of high density cells on one proc exceeds 10!")
                          end if
                          
                          x_high(loc_i) = xc(ip)
                          y_high(loc_i) = yc(jp)
                          z_high(loc_i) = zc(kp)
                          dens_high(loc_i) = comovingCellDens
                          !jeans_rad(loc_i) = sqrt((PI * cs**2.0) / (Newton * low_dens_thres)) / 2.0
                          jeans_rad(loc_i) = sqrt((PI * cs**2.0) / (Newton * low_dens_thres)) / 2.0
!                          x = xc(ip)
!                          y = yc(jp)
!                          z = zc(kp)
!                          pt = time

!                          pno = pt_sinkCreateParticle(x, y, z, pt, blockID, MyPE)

!                          write(*,'(A,4(1X,ES16.9),3I8)') "sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
!                             & x, y, z, pt, blockID, MyPE, int(particles_local(iptag,pno))

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
  
  ! Now send the information about where the cells that are over the
  ! high dens threshold are so that we can find cells near them on
  ! all processors.
  
  !print*, "On proc ", MyPE, " number of high dens cells = ", n_cells_high(MyPE+1)
  
!  if (MyPE .eq. MasterPE) then
!  print*, "size of n_cells_high =", size(n_cells_high)
!  print*, "size of x_high =", size(x_high)
!  print*, "size of dens_high =", size(dens_high)
!  end if

  call MPI_ALLREDUCE(MPI_IN_PLACE, x_high, array_length, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, y_high, array_length, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, z_high, array_length, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, dens_high, array_length, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, n_cells_high, nproc, MPI_INTEGER, &
                     MPI_SUM, comm, ierr)

  if (MyPE .eq. MasterPE) print*, "MPI_ALLREDUCE success."
  
  tot_n_cells = sum(n_cells_high)
!  print*, "I'm ", MyPE, "and total n cells =", tot_n_cells
!  print*, "I'm ", MyPE, "and n_cells_high(MyPE+1) =", n_cells_high(MyPE+1)
  ! If there are no cells over high dens threshold, return here.
  if (tot_n_cells .eq. 0) then
  
    if (MyPE .eq. MasterPE) print*, "No cells hit the high density threshold this step."
    return

  end if

  if (MyPE .eq. MasterPE) then 
  
  print*, "Total number of over dense regions = ", tot_n_cells
  end if
  !print*, "Jeans radius in parsecs = "
  !write(*,'(ES10.3e3)') jeans_rad / (2.0 * 3.086e18)
  ! Do any two dense regions overlap?
!  overlap = .false.
  ! This loop needs to skip values in x_high, y_high, etc that = 0.
!  do n=0, nproc-1
!    do i=1, n_cells_high(n+1)-1
!      i_val = n*10 + i
!      do j=i+1, n_cells_high(n+1)
!        j_val = n*10 + j
        !!! NOTE: This should be corrected for PBCs at some point.
!        r1 = (x_high(i_val)-x_high(j_val))**2.0 + (y_high(i_val)-y_high(j_val))**2.0 & 
!           + (z_high(i_val)-z_high(j_val))**2.0
!        r1 = sqrt(r1)
      
!        if (r1 .le. jeans_rad(i_val)) then
!          overlap = .true.
          !print*, "Zoiks! Two regions overlap!"
!        end if
!      end do
!    end do
!  end do
  
!  end if
  
  ! Now we know how many mean and var calculations we'll have to do.
  ! Note its faster to loop over all cells on proc once and then loop
  ! over the few dense cells, than to loop over the dense cells then over
  ! all the cells on the proc tot_n_cells times. So we make some arrays to
  ! do this just once.
  
  allocate(mean_dens(array_length)) ! Originally done with tot_n_cells
  allocate(var_dens(array_length))  ! for these arrays, and I'd like to
  allocate(mass(array_length))      ! go back to that at some point.
  allocate(n_save(array_length))
  allocate(tot_star_mass(array_length))
  
  mean_dens = 0.0
  var_dens  = 0.0
  mass      = 0.0
  n_save    = 0
  
  !print*, size(x_high), size(y_high), size(z_high)
  
  ! Check to see if any gas is within the region initially defined on
  ! a different processor. Note if this is NOT the case (on any proc),
  !  we can just sample the IMF and distribute stars locally without 
  ! any MPI communication.
  offproc = .false.

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
              
              do n=0, nproc-1
                  
                do i=1, n_cells_high(n+1)
                    
                  i_val = n*10 + i
                 
                  !print*, "n=", n
                  !print*, "i_val=", i_val
                  !print*, "n_cells_high(n+1)=", n_cells_high(n+1) 
                  !!! NOTE: This should be corrected for PBCs at some point.
                  
                  r1 = (xc(ip)-x_high(i_val))**2.0 &
                     + (yc(jp)-y_high(i_val))**2.0 &
                     + (zc(kp)-z_high(i_val))**2.0
                  r1 = sqrt(r1)
                  
#ifdef ALFL_MSCALAR
                 ! in case we are using the Alfven limiter
                 comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(DENS_VAR,ip,jp,kp)
#else
                 comovingCellDens = solnData(DENS_VAR,ip,jp,kp)
#endif
                  !cs = sqrt(solnData(PRES_VAR,ip,jp,kp)/comovingCellDens)
                  !low_jeans_diameter = sqrt(PI*cs**2.0 / (Newton * low_dens_thres))
                  !print*, "The sound speed is ", cs
                  !print*, "Low density Jeans radius = ", low_jeans_diameter/(2.0*3.086e18), " parsecs."
                    
                  !if (r1 .le. (low_jeans_diameter/2.0)) then
                  if (r1 .le. jeans_rad(i_val)) then
                    
                    if (comovingCellDens .gt. low_dens_thres) then
                        if (n .ne. MyPE) offproc = .true.
                         
                    ! Note this is basically volume weighted harvesting of the gas. So we might
                    ! consider doing volume weighted placement of the stars, i.e. the most
                    ! massive 15% get placed in the inner 0.4*jeans_rad, then the next
                    ! 0.15*max_mass < mass < 0.55*max_mass goes between 0.4*jeans_rad < rad < 0.8*jeans_rad,
                    ! finally the rest is outside of 0.8*jeans_rad. - JW
                    
#ifdef ALFL_MSCALAR
                        ! in case we are using the Alfven limiter
                        !alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(idens,ip,jp,kp)
                        
                        if (r1 .lt. 0.4*jeans_rad(i_val)) then
                        
                            mass(i_val) = mass(i_val) +  eff*0.85*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = (density_thresh_comoving + 0.15*(comovingCellDens - low_dens_thres) + alfl_dens)
                            
                        else if ((r1 .ge. 0.4*jeans_rad(i_val)) &
                                 .and. (r1 .lt. 0.8*jeans_rad(i_val))) then
                                 
                            mass(i_val) = mass(i_val) + eff*0.55*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = (density_thresh_comoving + 0.45*(comovingCellDens - low_dens_thres) + alfl_dens)
                            
                        else
                        
                            mass(i_val) = mass(i_val) + eff*0.15*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = (density_thresh_comoving + 0.85*(comovingCellDens - low_dens_thres) + alfl_dens)
                            
                        end if
                        
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving + alfl_dens
                        !solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(idens,ip,jp,kp)
#else
                        if (r1 .lt. 0.4*jeans_rad(i_val)) then
                        
                            mass(i_val) = mass(i_val) + eff*0.85*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = 0.15*(comovingCellDens - low_dens_thres) + low_dens_thres
                            
                        else if ((r1 .ge. 0.4*jeans_rad(i_val)) &
                                 .and. (r1 .lt. 0.8*jeans_rad(i_val))) then
                                 
                            mass(i_val) = mass(i_val) + eff*0.55*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = 0.45*(comovingCellDens - low_dens_thres) + low_dens_thres
                        
                        else
                        
                            mass(i_val) = mass(i_val) + eff*0.15*(comovingCellDens - low_dens_thres)*dVol
                            !solnData(idens,ip,jp,kp) = 0.85*(comovingCellDens - low_dens_thres) + low_dens_thres
                            
                        end if
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving
#endif
                        
!                      if (n .ne. MyPE) then  ! Done higher up.
!                        offproc = .true.
!                      end if
                        
                      !mass(i_val) = mass(i_val) + (comovingCellDens - low_dens_thres) * dVol ! Total mass gathered on this proc.
                      
                      !!!!! This wont work, n_save is the total number of all saved cells,
                      !!!!! it should be the number for this over dens cell.

					  ! Note here we intend to use this number later to remove any "extra" mass from
					  ! oversampling the IMF to get the proper # of massive stars. This allows us to
                      ! borrow more mass than initially granted and then just take the extra from all
					  ! the Jeans region.

                      n_save(i_val) = n_save(i_val) + 1  ! Number of cells we have stored on this proc.
                        
                      !if (n_save .gt. n_dens_max) then  ! Have we stored more than the array limit?
                      !   print*, "n_save = ", n_save, " n_dens_max =", n_dens_max
                      !   call Driver_abortFlash("Overflow for dense cell array, n_save > n_dens_max!")
                      !end if
                      
                      
                      ! Store the (original) dens, proc, blk and position.
                      ! We'll use this to place the stars we make.
                      
                      if (place_by_dens) then
                      
                      cell_count(n) = cell_count(n) + 1
                      
                      dens_cells(n,cell_count(n)) = solnData(DENS_VAR,ip,jp,kp)  
                      dens_blk(n,cell_count(n)) = blockID
                      dens_x(n,cell_count(n)) = xc(ip)
                      dens_y(n,cell_count(n)) = yc(jp)
                      dens_z(n,cell_count(n)) = zc(kp)
                      
                      ! Don't store anything. Lets just calculate the mean and var.
                      
!                      mean_dens(i_val) = mean_dens(i_val) + comovingCellDens
!                      var_dens(i_val)  = var_dens(i_val) + comovingCellDens**2.0
                      
                      end if
                      
                      ! Now update the density of the cell.  
!#ifdef ALFL_MSCALAR
!                      ! in case we are using the Alfven limiter
!                      alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(DENS_VAR,ip,jp,kp)
!                      solnData(DENS_VAR,ip,jp,kp) = low_dens_thres + alfl_dens
!                      solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(DENS_VAR,ip,jp,kp)
!#else
!                      solnData(DENS_VAR,ip,jp,kp) = low_dens_thres
!#endif
                      !affected_block = .true. ! If we changed dens, we need to call EOS.
                        
                    end if !  low thres criterion
                      
                  end if ! low_jeans radius criterion
                end do  ! high threshold dens loop  
              end do  ! processor loop
            end do  ! i loop
          end do  ! j look
        end do  ! k loop
        
      call Grid_releaseBlkPtr(lb, solnData)
      
      deallocate(xc)
      deallocate(yc)
      deallocate(zc)
      
!      if(affected_block) then
!         call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)
!      endif
      
      end if ! if block at max refinement

  end do  ! block loop
  
  ! Here we need to sum up the mass collected on each processor...
  ! Well, maybe. Or do we want to just try to make stars on each
  ! processor that gathered mass, regardless of whether its also
  ! the processor with the central high dense cell?
  
  mass = mass/1.989e33
  
  ! First check to see if we need to do global or local star formation.
  print*, "Do we have off proc data?", offproc, " on ", MyPE
  
  call MPI_ALLREDUCE(MPI_IN_PLACE, offproc, 1, MPI_LOGICAL, MPI_LOR, &
                     comm, ierr)
                     
  print*, "Does anyone have off proc data?", offproc
  
  ! Here only processors that have stored a cell participate in the MPI
  ! gather calls. So only processors that have cells inside the jeans
  ! radius are in the calls here.
  
  if (offproc) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, mass, array_length, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, n_save, array_length, MPI_INTEGER, &
                     MPI_SUM, comm, ierr)
    
    if (place_by_dens) then
                     
    do n=0, nproc
    
      if (n_save(n) .ne. 0) then
    
        call MPI_ALLGATHER(MPI_IN_PLACE, cell_count(n), MPI_DOUBLE_PRECISION, &
                         dens_cells(n,:), n_save(n), MPI_DOUBLE_PRECISION, &
                         comm, ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE, cell_count(n), MPI_DOUBLE_PRECISION, &
                         dens_x(n,:), n_save(n), MPI_DOUBLE_PRECISION, &
                         comm, ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE, cell_count(n), MPI_DOUBLE_PRECISION, &
                         dens_y(n,:), n_save(n), MPI_DOUBLE_PRECISION, &
                         comm, ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE, cell_count(n), MPI_DOUBLE_PRECISION, &
                         dens_z(n,:), n_save(n), MPI_DOUBLE_PRECISION, &
                         comm, ierr)
                         
      end if
    end do
    
    end if ! place by density.
    
  end if
  
  if (maxval(mass) .gt. 0.0) then
  
    !print*, "Current mass collected is ", mass, " on ", MyPE
    print*, "Current mass collected is ", mass(10*MyPE+1:10*(MyPE+1)), " on ", MyPE
    print*, "n_save = ", n_save
  
  end if 
  
  do n=0, nproc-1
      
    !if (n .ne. MyPE) offproc = .true.
    do i=1, n_cells_high(n+1)
        
      i_val = n*10 + i
        
      if (mass(i_val) .gt. min_mass) then
		made_stars=.true.
	  !else
        !n_save(i_val) = 0 ! Not making stars here, so reset the count to 0 for this one.
	  end if
        
    end do
  
  end do
  
  if (made_stars) then
  
  print*, "Going to try and make stars now."
  

!  mean_dens = mean_dens/n_save
!  var_dens  = var_dens/n_save + mean_dens**2.0
  
!  deallocate(x_high)
!  deallocate(y_high)
!  deallocate(z_high)
!  deallocate(dens_high)

  offproc = .false.
  
  !print*, "n_save = ", n_save

    ! Make stars from all the mass that was collected for a local high
    ! density cell. We will then place them locally followed by a call
    ! to Particles_moveAndSort to put them on the correct processor.
    
    !mass = mass / 1.989e33 ! Put in solar masses.
    !min_mass = 0.5 ! Smallest star we can make.
    allocate(stars(1000)) ! Should make an estimate here.
    allocate(mean_stars(1000))
    allocate(var_stars(1000))
    
    allocate(stars_offproc(1000))
    allocate(starsproc_offproc(1000))
    allocate(starsx_offproc(1000))
    allocate(starsy_offproc(1000))
    allocate(starsz_offproc(1000))
    
    print*, "stars arrays allocated."
    
  do n=1, n_cells_high(MyPE+1)
                    
      i_val = MyPE*10 + n
      !i_val = n*10 + i
      stars   = 0.0
      n_stars = 0
      
      print*, "i_val =", i_val
      print*, "mass(i_val) =", mass(i_val)
      print*, "min_mass =", min_mass
      
      if (mass(i_val) .le. min_mass) cycle
      
      print*, "Sampling the IMF with ", mass(i_val), " solar masses."
    
      call sample_IMF(mass(i_val), min_mass, stars, n_stars)
      
      print*, "After star formation, total mass left is =", mass(i_val)
      
      ! Convert back to cgs.
      mass(i_val) = mass(i_val) * 1.989e33
      
    ! Sort the star masses from smallest to largest.
      stars = stars * 1.989e33 ! Put back in grams.
      call ut_qsort(stars, n_stars, .true.)

	  tot_star_mass(i_val) = sum(stars)
      
      print*, "We got ", n_stars, " with total mass ", tot_star_mass(i_val)/1.989e33
      
      ! Get the mean and variance of the star masses.
      mean_stars(i_val) = 0.0
      var_stars(i_val)  = 0.0
      do i=1, n_stars
        mean_stars(i_val) = mean_stars(i_val) + stars(i)
        var_stars(i_val)  = var_stars(i_val) + stars(i)**2.0
      end do
      
      mean_stars(i_val) = mean_stars(i_val)/n_stars
      var_stars(i_val)  = var_stars(i_val)/n_stars + mean_stars(i_val)**2.0
    
    ! Calculate the probabilty distribution of the dense gas.
    ! Note this step might not be necessary in the final calculation,
    ! but I want to see the data output for plotting for now.

      !call ut_qsort(dens_cells(MyPe,:), n_save(MyPE), .true.)
      
      if (place_by_dens) then
      
        allocate(QSindex(n_save(MyPE)))
        !allocate(mean_dens(n_save(MyPE)))
        !allocate(var_dens(n_save(MyPE)))
        
      
        ! Get the mean and variance of the star masses.
        mean_dens(i_val) = 0.0
        var_dens(i_val)  = 0.0
        do i=1, n_save(MyPE)
        
        ! This I think is set up wrong. I think dens_cells will have
        ! to be dens_cells(num_cells_high*nproc, num_dens_cells) and not
        ! dens_cells(nproc, num_dens_cells) to work.
            mean_dens(i_val) = mean_dens(i_val) + dens_cells(MyPE,i)
            var_dens(i_val)  = var_dens(i_val) + dens_cells(MyPE,i)**2.0
        end do
      
        mean_stars(i_val) = mean_stars(i_val)/n_stars
        var_stars(i_val)  = var_stars(i_val)/n_stars + mean_stars(i_val)**2.0
      
        !call NewQsort_IN(dens_cells(MyPE,:), QSindex)

      end if

      !allocate(dens_gaus_pdf(n_save))
    
      !call gaussian(dens_matrix(:,1), dens_gaus_pdf, mu, var)

    !fout = 'gaus_out.txt'
    !open(unit=10, file=fout)
    
    !write(10, '(2A15)') "density", "prob"
    
    !do i=1, n_save
    
    !  write(10, '(2ES15.3)') dens_matrix(i,1), dens_gaus_pdf(i)
      
    !end do
    
    !close(10)
    

    
    !loc = real(i) / real(n_stars) ! Location (index) in the stars array
    !dens_loc = ceiling(loc*real(n_save)) ! Turned into loc (index) in the dens array
  
    ! Center the gaussian distribution on the location in the dense array.
    !mu = dens_matrix(dens_loc,1)
    ! Get the probability distribution centered on this point and 
    ! with HALF the standard deviation of the original density distribution
    ! (note this is a free parameter which we can change to match observations).
    !allocate(star_gaus_pdf(n_stars))

    !call gaussian(dens_matrix(:,1), star_gaus_pdf, mu, var)
    
    ! How far from the mean mass is this star?
  
!        if (i==2) then
  
!          do j=1, n_save
!           blk=3, i=4, j=5, k=6
!            blk = int(dens_matrix(j, 3))
!            ind = [int(dens_matrix(j, 4)),int(dens_matrix(j, 5)),int(dens_matrix(j, 6))]
!            call Grid_putPointData(blk, SCRATCH_CTR, &
!                                   PROB_SCRATCH_CENTER_VAR, INTERIOR, &
!                                   ind, star_gaus_pdf(j))
!          end do
    
!        end if
    
    stars_offproc = 0.0
    starsproc_offproc = -1
    starsx_offproc = 0.0
    starsy_offproc = 0.0
    starsz_offproc = 0.0
    starsnum_offproc = 0
    
    do i=1, n_stars
  
      accept = .false.
      n_loops = 0
      
      dist_from_mean = stars(i) / mean_stars(i_val)
  
      do while ((.not. accept) .and. (n_loops .lt. 1e4))
      
        if (place_by_dens) then
        
        else ! Place randomly.
      
            ! Get a random location within the Jeans radius.
            call random_number(xr)
            !print*, "Random x =", xr
            xr = (xr - 0.5)*2.0*jeans_rad(n) !+ x_high(n)
            call random_number(yr)
            yr = (yr - 0.5)*2.0*jeans_rad(n) !+ y_high(n)
            call random_number(zr)
            zr = (zr - 0.5)*2.0*jeans_rad(n) !+ z_high(n)
            
            if ( sqrt(xr**2.0 + yr**2.0 + zr**2.0) .gt. jeans_rad(n)) cycle
            
            xr = xr + x_high(n)
            yr = yr + y_high(n)
            zr = zr + z_high(n)
            
        
        end if
        
    
      ! Get the cell density at this location
    
        call Grid_getCellIndexFromPosition(proc, blockID, ip, jp, kp, xr, yr, zr)
        !print*, "Random position ", xr, yr, zr
        !print*, "Indices are ", ip, jp, kp
    
        if (proc /= MyPE) then
          !offproc = .true.
          stars_offproc(i)  = stars(i)
          starsproc_offproc = proc
          starsx_offproc(i) = xr
          starsy_offproc(i) = yr
          starsz_offproc(i) = zr

          starsnum_offproc = starsnum_offproc + 1
          
          !print*, "Cell is off proc for star placement, need to write that code!"
          !call Driver_abortFlash("Cell off proc in star placement")
        else
          !call Grid_getPointData(blockID, CENTER, DENS_VAR, INTERIOR, [ip, jp, kp], ComovingCellDens)
          
          n_loops = n_loops + 1
      
        ! Does this cell meet the over low dens threshold? If not, pick
        ! a new one.
        ! Can't do this ^^^, because we already took the over dense gas out.
        ! Just check and see if it EQUALS the low dens thres, cause
        ! that means we took gas out from here.
        !if (solnData(DENS_VAR, ip, jp, kp) .lt. low_dens_thres) then 
        !  cycle
        !else
        
        ! Here's where I think we try something like volume weighted
        ! placement described above where we harvested the gas. - JW

        ! Otherwise lets continue...
      
        ! If just using random placement (vol equally weighted) then just
        ! place the star.
      
        !call Grid_getSingleCellCoords([ip, jp, kp], blockID, CENTER, &
        !                              INTERIOR, coords)
          call Grid_getBlkPtr(blockID, solnData)
          pno = pt_sinkCreateParticle(xr, yr, zr, &
                                    time, blockID, MyPE)
          particles_local(MASS_PART_PROP, pno) = stars(i)
          particles_local(VELX_PART_PROP, pno) = solnData(VELX_VAR, ip, jp, kp)
          particles_local(VELY_PART_PROP, pno) = solnData(VELY_VAR, ip, jp, kp)
          particles_local(VELZ_PART_PROP, pno) = solnData(VELZ_VAR, ip, jp, kp)
        
          call Grid_releaseBlkPtr(blockID, solnData)
        end if ! If offproc
        accept = .true.
        !print*, "Star placed!"
        
        !end if
        ! Now the Monte Carlo magic. Draw a random number...
      
        !!! NOTE NOTE NOTE !!!
        ! With the new (single return) gaussian prob function, we should
        ! make sure to normalize our probabilities so that the center of 
        ! the probability function has a perfect (1.0) chance of accepting
        ! if it got pulled (since no internal normalization occurs anymore).
      
        !call random_number(rnum)
      
        !prob = gaussian(solnData(DENS_VAR, ip, jp, kp), 
!          ! And then use that to get a density from the array...
!          rind = ceiling(rnum*n_save)
!          ! And use it to calculate a probability for that location to be
!          ! accepted to place the star in, normalized by the central probability.
!          prob = star_gaus_pdf(rind) / star_gaus_pdf(dens_loc)
!          ! Get your second random number...
!          call random_number(rnum)
!          ! And compare versus your previous probability. If greater than 1.0, accept!
!          if ( (prob/rnum) .ge. 1.0) then
      
!            accept = .true.
!            ! Now place your star at this cell location.
!            ! Sink particle creation input looks like this:
!            ! call pt_sinkCreateParticle(x, y, z, create_sim_time, blk, proc) 
!            ! which translates into the following for us:
!            blk = int(dens_matrix(rind,3))
!            ind(1) = int(dens_matrix(rind, 4))
!            ind(2) = int(dens_matrix(rind, 5))
!            ind(3) = int(dens_matrix(rind, 6))
        
!            call Grid_getSingleCellCoords(ind, blk, CENTER, &
!                                          INTERIOR, coords)
!            pno = pt_sinkCreateParticle(coords(1), coords(2), coords(3), &
!                                        time, blk, MyPE)
!            particles_local(MASS_PART_PROP, pno) = stars(i)
!          end if  ! Did we accept?
      
        end do ! Did we accept the placement location?
      
      end do ! Loop to place stars.
      
      ! Now lets place the off processor stars.
      
      do i=1, nproc-1
      
        idest = MyPe + 1
        if (idest > nproc-1) idest = 0
        isrc  = MyPe - 1
        if (isrc < 0) isrc = nproc -1
        
        call MPI_SENDRECV_REPLACE(stars_offproc, 1000, MPI_DOUBLE_PRECISION, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
        call MPI_SENDRECV_REPLACE(starsnum_offproc, 1, MPI_INTEGER, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
        call MPI_SENDRECV_REPLACE(starsproc_offproc, 1000, MPI_INTEGER, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
        call MPI_SENDRECV_REPLACE(starsx_offproc, 1000, MPI_DOUBLE_PRECISION, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
        call MPI_SENDRECV_REPLACE(starsy_offproc, 1000, MPI_DOUBLE_PRECISION, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
        call MPI_SENDRECV_REPLACE(starsz_offproc, 1000, MPI_DOUBLE_PRECISION, &
                 idest, MyPE, isrc, isrc, comm, comm_stat, ierr)
                 
        do j=1, starsnum_offproc
        
          if (starsproc_offproc(j) .eq. MyPe) then
          
            call Grid_getCellIndexFromPosition(proc, blockID, ip, jp, kp, &
                               starsx_offproc(j),starsy_offproc(j),starsz_offproc(j))
            call Grid_getBlkPtr(blockID, solnData)
            pno = pt_sinkCreateParticle(starsx_offproc(j),starsy_offproc(j),starsz_offproc(j), &
                                    time, blockID, MyPE)
            particles_local(MASS_PART_PROP, pno) = stars_offproc(j)
            particles_local(VELX_PART_PROP, pno) = solnData(VELX_VAR, ip, jp, kp)
            particles_local(VELY_PART_PROP, pno) = solnData(VELY_VAR, ip, jp, kp)
            particles_local(VELZ_PART_PROP, pno) = solnData(VELZ_VAR, ip, jp, kp)
            call Grid_releaseBlkPtr(blockID, solnData)
            
          end if
          
        end do ! Loop over stars_offproc array.
        
      end do ! Loop over # of procs, rotating the stars_offproc arrays
             ! around the processors in turn.

    end do  ! Loop over all n_dens_cells
    
!  else
  
!    print*, "Well damn we got some off proc stuff, better write that off proc code!"
!    call Driver_abortFlash("And we're aborting due to off proc information needed.")
  
!  end if ! if off_proc
    deallocate(stars)
    deallocate(stars_offproc)
    deallocate(starsproc_offproc)
    deallocate(starsx_offproc)
    deallocate(starsy_offproc)
    deallocate(starsz_offproc)
    deallocate(mean_stars)
    deallocate(var_stars)

! Now that we actually know how much mass we used to make stars (which might
! be more than we caculated from the cells) lets actually remove what we
! used from the gas. Note that any extra taken or left over density should be
! equal to (tot_star_mass(i_val) - mass(i_val)) / n_save(i_val) for each cell.


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
              
              do n=0, nproc-1
                  
                if (n .ne. MyPE) offproc = .true.
                do i=1, n_cells_high(n+1)
                    
                  i_val = n*10 + i
                  
                  !!! NOTE: This should be corrected for PBCs at some point.
                  
                  r1 = (xc(ip)-x_high(i_val))**2.0 &
                     + (yc(jp)-y_high(i_val))**2.0 &
                     + (zc(kp)-z_high(i_val))**2.0
                  r1 = sqrt(r1)
                  
#ifdef ALFL_MSCALAR
                 ! in case we are using the Alfven limiter
                 comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(DENS_VAR,ip,jp,kp)
#else
                 comovingCellDens = solnData(DENS_VAR,ip,jp,kp)
#endif
                  !cs = sqrt(solnData(PRES_VAR,ip,jp,kp)/comovingCellDens)
                  !low_jeans_diameter = sqrt(PI*cs**2.0 / (Newton * low_dens_thres))
                  !print*, "The sound speed is ", cs
                  !print*, "Low density Jeans radius = ", low_jeans_diameter/(2.0*3.086e18), " parsecs."
                    
                  !if (r1 .le. (low_jeans_diameter/2.0)) then
                  if (r1 .le. jeans_rad(i_val)) then
                    
                    if (comovingCellDens .gt. low_dens_thres) then
                    ! Note here we take out the "extra" mass from oversampling the IMF.
                    
#ifdef ALFL_MSCALAR
                        ! in case we are using the Alfven limiter
                        alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(idens,ip,jp,kp)
                        
                        if (r1 .lt. 0.4*jeans_rad(i_val)) then
                        
                            !mass(i) = mass(i) +  0.85*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.85)*(comovingCellDens - low_dens_thres) &
                                                     + alfl_dens + mass(i_val)/n_save(i_val)/dVol
                            
                        else if ((r1 .ge. 0.4*jeans_rad(i_val)) &
                                 .and. (r1 .lt. 0.8*jeans_rad(i_val))) then
                                 
                            !mass(i) = mass(i) + 0.55*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.55)*(comovingCellDens - low_dens_thres) &
                                                     + alfl_dens + mass(i_val)/n_save(i_val)/dVol
                            
                        else
                        
                            !mass(i) = mass(i) + 0.15*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.15)*(comovingCellDens - low_dens_thres) &
                                                     + alfl_dens) + mass(i_val)/n_save(i_val)/dVol
                            
                        end if
                        
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving + alfl_dens
                        solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(idens,ip,jp,kp)
#else
                        if (r1 .lt. 0.4*jeans_rad(i_val)) then
                        
                            !mass(i) = mass(i) + 0.85*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.85)*(comovingCellDens - low_dens_thres) & 
                                                     + low_dens_thres + mass(i_val)/n_save(i_val)/dVol
                            
                        else if ((r1 .ge. 0.4*jeans_rad(i_val)) &
                                 .and. (r1 .lt. 0.8*jeans_rad(i_val))) then
                                 
                            !mass(i) = mass(i) + 0.55*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.55)*(comovingCellDens - low_dens_thres) &
                                                     + low_dens_thres + mass(i_val)/n_save(i_val)/dVol
                        
                        else
                        
                            !mass(i) = mass(i) + 0.15*(comovingCellDens - low_dens_thres)
                            solnData(idens,ip,jp,kp) = (1.0-eff*0.15)*(comovingCellDens - low_dens_thres) &
                                                     + low_dens_thres + mass(i_val)/n_save(i_val)/dVol
                            
                        end if
                        !solnData(idens,ip,jp,kp) = density_thresh_comoving
#endif
                        
!                      if (n .ne. MyPE) then  ! Done higher up.
!                        offproc = .true.
!                      end if
                        
                      !mass(i_val) = mass(i_val) + (comovingCellDens - low_dens_thres) * dVol ! Total mass gathered on this proc.
                      
                      !!!!! This wont work, n_save is the total number of all saved cells,
                      !!!!! it should be the number for this over dens cell.
                      !n_save(i_val) = n_save(i_val) + 1  ! Number of cells we have stored on this proc.
                        
                      !if (n_save .gt. n_dens_max) then  ! Have we stored more than the array limit?
                      !   print*, "n_save = ", n_save, " n_dens_max =", n_dens_max
                      !   call Driver_abortFlash("Overflow for dense cell array, n_save > n_dens_max!")
                      !end if
                      
                      ! Don't store anything. Lets just calculate the mean and var.
                      
!                      mean_dens(i_val) = mean_dens(i_val) + comovingCellDens
!                      var_dens(i_val)  = var_dens(i_val) + comovingCellDens**2.0
                      
                      ! Store the (original) dens, proc, blk and position.
                      ! We'll use this to place the stars we make.
!                      dens_cells(n_save) = solnData(DENS_VAR,ip,jp,kp)  
!                      dens_blk(n_save) = blockID
!                      dens_x(n_save) = ip !xc(ip)
!                      dens_y(n_save) = jp !yc(jp)
!                      dens_z(n_save) = kp !zc(kp)
                      
                      ! Now update the density of the cell.  
!#ifdef ALFL_MSCALAR
!                      ! in case we are using the Alfven limiter
!                      alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(DENS_VAR,ip,jp,kp)
!                      solnData(DENS_VAR,ip,jp,kp) = low_dens_thres + alfl_dens
!                      solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(DENS_VAR,ip,jp,kp)
!#else
!                      solnData(DENS_VAR,ip,jp,kp) = low_dens_thres
!#endif
                      affected_block = .true. ! If we changed dens, we need to call EOS.
                      
                    end if !  low thres criterion
                      
                  end if ! low_jeans radius criterion
                end do  ! high threshold dens loop  
              end do  ! processor loop
            end do  ! i loop
          end do  ! j look
        end do  ! k loop
        
      call Grid_releaseBlkPtr(lb, solnData)
      
      deallocate(xc)
      deallocate(yc)
      deallocate(zc)
      
      if(affected_block) then
         ! Changed from MODE_DENS_EI because if we took mass out
         ! how can energy be the same? - JW
         ! Because eint is specific mass in Flash. - JW
         call Eos_wrapped(MODE_DENS_TEMP, blkLimits, lb)
      endif
      
      end if ! if block at max refinement

  end do  ! block loop
  
  call Grid_notifySolnDataUpdate()
  

  end if ! if make_stars true
  
  deallocate(x_high)
  deallocate(y_high)
  deallocate(z_high)
  if (place_by_dens) then
    deallocate(dens_high)
    deallocate(dens_x)
    deallocate(dens_y)
    deallocate(dens_z)
    deallocate(QSindex)
  end if
  deallocate(n_save)
  deallocate(mass)
  deallocate(mean_dens)
  deallocate(var_dens)



  deallocate(jeans_rad)
  
  call pt_sinkGatherGlobal()
  
end subroutine dens_removal
