

!#define debug_form



subroutine stochastic(dt)

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

real, pointer, dimension(:,:,:,:)   :: solnData

real, dimension(:), allocatable     :: xc, yc, zc
real                :: coords(3)
integer             :: ind(3)

logical, save       :: first_call=.true.
integer, save       :: array_length, array_start, start_ind, end_ind, nproc
integer, save       :: lrefine_max
integer             :: tot_n_cells, n_stars, n_loops, n !, n_save
integer             :: pno
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

! My variables

real, save           :: low_dens_thres, high_dens_thres, min_mass, f_eff
real                 :: mean_dens, var_dens, mean_stars, &
                                       var_stars 
real                 :: time, loc !, low_jeans_diameter !, mass
real                 :: rnum, r1, prob, mu, var, cs
real                 :: xr, yr, zr, dist_from_mean, max_prob, local_dt
real                 :: pt, rn, t_ff, rnd_mass, cell_mass, timeSinceLastSF
real                 :: max_mass, max_dens, max_star, dm, norm, shift
real, dimension(2)   :: sendbuf, recvbuf
integer              :: max_prob_proc, lp

real, intent(in)     :: dt
integer              :: ierr, num_stars
logical              :: made_stars

! Note these 3 defines let us get the proper number of tags
! previously handed out if this is a restart, so that tags
! never get duplicated. Its important! - JW

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  if (first_call) then
    
    f_eff = 1.0 ! Stellar formation effeciency... less useful now that
                ! we are normalizing things. - JW

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

    call pt_sinkGatherGlobal()

    local_tag_number = 0
    do lp = 1, localnpf
        if (get_ppe(int(particles_global(iptag,lp))) .EQ. MyPE) then
           local_tag_number = max(local_tag_number, get_pno(int(particles_global(iptag,lp))))
        endif
    enddo

    if (MyPE == MASTER_PE) then
      print*, "Star formation turned on!"
      print*, "High density thresh = ", high_dens_thres
      print*, "Low density thresh = ", low_dens_thres
    end if
    
    call Grid_getListOfBlocks(LEAF,blockList,blockCount)
    
    do lb=1, blockCount
    
      call Grid_getBlkPtr(lb,solnData)
      solnData(LSFT_VAR, :, :, :) = 0.0
      call Grid_releaseBlkPtr(lb,solnData)
      
    end do
    

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

  time  = dr_simTime
  
  num_stars = 0
  
  max_prob = 0.0
  max_mass = 0.0
  max_star = 0.0
  max_dens = 0.0

  
  ! Here we define the shift for zero probability to be probability to
  ! make a star greater than 100 Msun, so that the probability terminates
  ! at the correct location. Basically we evaulate it at 100 Msun and
  ! subtract off that value. - JW
  shift = 0.0 !0.109139138961
  
  ! Here we define the normalization for the probability density, which
  ! we integrated with quad from 1 Msun to 100 Msun in python. - JW
  
  ! This number gives a probability of 1 for 0.2 Msun at rho=1e-16 g/cm**3
  ! (basically the IMF turns over at 0.2 at high densities) and for 1-100
  ! Msun we get a perfect power law.
  norm = 39.4556561298 !22.6388363319 <- newer #, but makes higher prob... - JW
  
  ! Here we define the mass interval we sample over (and that we multiply
  ! the probability density by). - JW
  dm = 0.1

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
!                    do n = 1, n_cells_high(MyPe+1)
!                       distx = xc(ip) - x_high(n)
!                       disty = yc(jp) - y_high(n)
!                       distz = zc(kp) - z_high(n)
!                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
!                       rad = sqrt(distx**2 + disty**2 + distz**2 )
!                       ! Does this position fall within 2x accretion radius of existing sink?
!                       if (rad .le. 2.0*accretion_radius_comoving) then
!                          create_sink = .false.
!                       end if
!                    end do

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
                             do k1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check

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
                                         
                                         ! Here we'd probably want to calculate the velocity variance using <u^2> - <u>^2
                                         ! For now we'll just skip it. - JW

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

                             do k1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = 0, 0 !-ngc_sink_creation_check, ngc_sink_creation_check

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
                       
                       
                      ! Now check the probability of making a star here
                      ! using the mass in the sink radius from the calculation
                      ! above. Note this means the sink radius is tied to this routine. - JW
                      
                      ! t_ff is from the cell in question, weighting star creation
                      ! towards the most massive cell instead of spreading it around
                      ! the sink radius region. - JW
                      rnd_mass  = 0.0
                      cell_mass = solnData(DENS_VAR, ip, jp, kp) * dVol
                      t_ff = sqrt((3.0*PI) / (32.0 * Newton * ComovingCellDens))
                      
                      ! Local timestep from isothermal sound speed. - JW
                      local_dt = min(dx_block, dy_block, dz_block) / cs
                                                
                      ! Pick a random number from 1 to 100 solar masses
                      call random_number(rnd_mass)
                      rnd_mass = (rnd_mass)*100 + 1.0
                      rnd_mass = NINT(rnd_mass*10.0)*0.1 ! Round to the nearest tenth.
                      rnd_mass = rnd_mass * 1.989e33 ! Convert to grams.
                      
                      timeSinceLastSF = dt + solnData(LSFT_VAR, ip, jp, kp)

                      
                      ! Note we have to be a bit careful here. We can't
                      ! take all the mass out of a cell obiviously,
                      ! but what is a good fraction? Let's start
                      ! with 1/2 the mass I guess (as in Vazquez-Semademi 2016)? - JW
                                            
                      ! From Goldbaum et. al. 2015. - JW
                      !prob = f_eff * (cell_mass / (rnd_mass)) * (timeSinceLastSF / t_ff)
                      prob = (f_eff*(cell_mass/(rnd_mass)) * (local_dt / t_ff) - shift) &
                             / norm * dm
                      
                      ! Compare probability to a random number
                      
                      if (prob .gt. max_prob) then
                          max_prob = prob
                          max_mass = cell_mass/1.989e33
                          max_dens = ComovingCellDens
                          max_star = rnd_mass/1.989e33
                      end if
                      
                      call random_number(rn)
#ifdef debug_form                      
                      print*, "Prob is ", prob 
                      print*, "Prob ratio is ", prob/rn
#endif
                      if ( ((prob/rn) .lt. 1.0) &
                           .or. (rnd_mass .gt. 0.99*cell_mass) ) then
                          create_sink = .false.
                          solnData(LSFT_VAR, ip, jp, kp) = timeSinceLastSF
                      end if
                      
                      ! If probability is greater than 1, make
                      ! enough stars to equal up the probability.
                      if (prob .gt. 1.0) then
                          num_stars = floor(prob)
                      else
                          num_stars = 1
                      end if
#ifdef debug_form                      
                      print*, "Create sink =", create_sink
                      call flush(6)
#endif
                       ! finally create the sink in the cell
                      if (create_sink) then
                        
                        print*, "Proc", MyPE, "making", num_stars, "stars."
                        
                        do i=1, num_stars

                          print*, "Proc", MyPE, "making a star of mass", rnd_mass/1.989e33
                                          
                          call random_number(xr)
                          call random_number(yr)
                          call random_number(zr)

                          xr = xc(ip) + (xr-0.5)*dx_block
                          
                          yr = yc(jp) + (yr-0.5)*dy_block
                          
                          zr = zc(kp) + (zr-0.5)*dz_block
                          pt = time

                          pno = pt_sinkCreateParticle(xr, yr, zr, &
                                                    time, blockID, MyPE)
                                                    
                          ! Randomly perturb the velocity by a max of 10%
                          ! just for testing. To be replaced by a real sigma_v
                          ! soon. - JW.
                          
                          call random_number(xr)
                          call random_number(yr)
                          call random_number(zr)
                          
                          xr = ((xr-0.5)/5.0)*solnData(VELX_VAR, ip, jp, kp)
                          yr = ((yr-0.5)/5.0)*solnData(VELY_VAR, ip, jp, kp)
                          zr = ((zr-0.5)/5.0)*solnData(VELZ_VAR, ip, jp, kp)
                          
                          particles_local(MASS_PART_PROP, pno) = rnd_mass
                          particles_local(VELX_PART_PROP, pno) = solnData(VELX_VAR, ip, jp, kp) + xr
                          particles_local(VELY_PART_PROP, pno) = solnData(VELY_VAR, ip, jp, kp) + yr
                          particles_local(VELZ_PART_PROP, pno) = solnData(VELZ_VAR, ip, jp, kp) + zr
                        
                          write(*,'(A,4(1X,ES16.9),3I8)') "sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
                             & xr, yr, zr, pt, blockID, MyPE, int(particles_local(iptag,pno))

                        end do
                        
                             
                        ! Take that mass out of the cell here.
                      
                        solnData(DENS_VAR, ip, jp, kp) = (cell_mass - num_stars*rnd_mass) / dVol
                        solnData(LSFT_VAR, ip, jp, kp) = 0.0 
                        affected_block = .true. ! If we changed dens, we need to call EOS.
                        
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

      
        if(affected_block) then
           ! Changed from MODE_DENS_EI because if we took mass out
           ! how can energy be the same? - JW
           ! Because eint is specific mass in Flash. - JW
           call Eos_wrapped(MODE_DENS_TEMP, blkLimits, lb)
        endif

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)
        
        end if ! If block at max refinement.

  end do ! block loop
  
  ! ======== END find high dens cells loop ========
  
  sendbuf  = 0.0
  recvbuf  = 0.0
  max_prob_proc = -1
  sendbuf(1) = max_prob
  sendbuf(2) = gr_meshMe
  
  ! Reduce both the maximum probability and its location.
  call MPI_allReduce(sendbuf, recvbuf, 1, MPI_DOUBLE_INT, &
                     MPI_MAXLOC, gr_meshComm, ierr)
  
  max_prob = recvbuf(1)
  max_prob_proc = int(recvbuf(2))
                     
  call MPI_Bcast(max_mass, 1, MPI_DOUBLE_PRECISION, &
                     max_prob_proc, gr_meshComm, ierr)
  call MPI_Bcast(max_dens, 1, MPI_DOUBLE_PRECISION, &
                     max_prob_proc, gr_meshComm, ierr)
  call MPI_Bcast(max_star, 1, MPI_DOUBLE_PRECISION, &
                     max_prob_proc, gr_meshComm, ierr)

  if (gr_meshMe == 0) print*, "Recvbuf = ", recvbuf
  if (gr_meshMe == 0) print*, "Max probability = ", max_prob, "on proc", max_prob_proc
  if (gr_meshMe == 0) print*, "Max mass = ", max_mass
  if (gr_meshMe == 0) print*, "Max star = ", max_star
  if (gr_meshMe == 0) print*, "Max density = ", max_dens

  call Grid_notifySolnDataUpdate()

  call pt_sinkGatherGlobal()
  
end subroutine stochastic
