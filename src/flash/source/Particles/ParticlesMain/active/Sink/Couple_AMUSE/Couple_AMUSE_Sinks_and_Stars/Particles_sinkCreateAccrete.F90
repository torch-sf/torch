!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkCreateAccrete
!!
!! NAME
!!
!!  Particles_sinkCreateAccrete
!!
!! SYNOPSIS
!!
!!  call Particles_sinkCreateAccrete(real, intent(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Handles creation of sink particles and accretion of gas onto sink particles.
!!  It also calls routines for sink particle merging and for dumping sink
!!  particle properties every time step (pt_sinkDumpParticles).
!!
!! ARGUMENTS
!!
!!   dt - the current simulation time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2014
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   renamed and cleaned by Christoph Federrath, 2013
!!
!!***

!#define debug1

subroutine Particles_sinkCreateAccrete(dt)

  use Particles_sinkData
  use Particles_data, ONLY: pt_posInitialized, pt_velInitialized, pt_reduceGcellFills, &
                            pt_numLocal, allproc_particles
  use Particles_interface, ONLY: pt_gatherGlobal
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal, pt_sinkMergingAfterCreation, &
      pt_sinkFindList, pt_sinkParticleMerging, pt_sinkDumpParticles, pt_sinkCreateParticle, &
      pt_sinkCorrectForPeriodicBCs
  use Driver_interface, ONLY : Driver_abortFlash
  use Driver_data, ONLY : dr_globalMe, dr_simTime
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_data, ONLY : gr_maxRefine
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getCellCoords, & 
       Grid_getBlkPhysicalSize, Grid_getBlkPtr, Grid_releaseBlkPtr,  & 
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkRefineLevel
  use Eos_interface, ONLY : Eos_wrapped
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use Logfile_interface, ONLY : Logfile_stamp

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "Eos.h"
#include "GridParticles.h"

  implicit none

  real, intent(IN) :: dt
  
  integer :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList
  
  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(:), allocatable :: xc, yc, zc
  real, dimension(maxsinks) :: tot_mass, cm_x, cm_y, cm_z, vel_cm_x, vel_cm_y, vel_cm_z, &
                             & ang_x, ang_y, ang_z, etot, vr, radius
  integer, dimension(maxsinks) :: pindex_found
  integer             :: np_found, npf, pno_to_accrete
  real                :: egrav_gas, egrav_part, egrav, ekin, etot_min_inner_r_accr, etot_min
  logical             :: within_inner_r_accr
  real                :: size(3)
  real                :: dx_block, dy_block, dz_block, dVol, inner_r_accr
  real                :: mass, gpot, absgpot
  real                :: rad, pt, time, dvx, dvy, dvz
  real                :: x, y, z, cs, vrad, distx, disty, distz, px_old, py_old, pz_old
  real, save          :: delta_at_lrefmax
  integer             :: ip, jp, kp, lp, nlp, npart
  real, save          :: density_thresh, accretion_radius, Newton, xmin, xmax, ymin, ymax, zmin, zmax

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  real, save          :: mu_zero
  character(4), save  :: units
  integer, save       :: imagx, imagy, imagz
#endif

  logical, parameter  :: write_accretion_checks_info = .false.
  logical, parameter  :: print_creation_info = .false.
  logical, parameter  :: debug = .false.
  integer, parameter  :: funit_accretion_checks = 43
  integer, parameter  :: funit_accretion        = 44
  integer, parameter  :: ngc_sink_creation_check = 2
  real, parameter     :: ngc_sink_creation_check_radius_sqr = (ngc_sink_creation_check+1.0)**2
  integer             :: i1, j1, k1, ii1, jj1, kk1, ncells_in_vol
  real                :: vxcm_in_vol, vycm_in_vol, vzcm_in_vol, ekindisp_in_vol, etherm_in_vol, emag_in_vol
  real                :: r_search, mass_in_vol, maxgpot_in_vol, egravdeltapot_in_vol
  logical             :: create_sink, affected_block

  integer             :: lb, llb, llnblocks, blockID, pno, old_localnp

  integer, dimension(MAXBLOCKS) :: block_list

  character(len=80), save :: grav_boundary_type

  logical, save       :: convergingFlowCheck
  logical, save       :: negativeEtotCheck
  logical, save       :: jeansCheck
  logical, save       :: potentialMinCheck
  logical, save       :: GasAccretionChecks
  logical, save       :: sink_merging
  logical, save       :: first_call = .true.

  integer, save       :: iXcoord, iYcoord, iZcoord, izn, lrefine_max
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, save       :: MyPE, MasterPE, idens, ipres, ivelx, ively, ivelz, igpot, ieint
  integer, save       :: Nblockx, Nblocky, Nblockz

  integer             :: size_x, size_y, size_z
  real                :: redshift, comovingCellDens, gas_dens
  real                :: onePlusRedshift, onePlusRedshift2, onePlusRedshift3
  real                :: accretion_radius_comoving, density_thresh_comoving

#ifdef ALFL_MSCALAR
  real                :: alfl_dens
#endif

  integer, parameter :: gather_nprops_pos = 3
  integer, dimension(gather_nprops_pos), save :: gather_propinds_pos = &
    (/ integer :: POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP /)

  integer, parameter :: gather_nprops = 11
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: TAG_PART_PROP, MASS_PART_PROP, &
                  POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, &
                  VELX_PART_PROP, VELY_PART_PROP, VELZ_PART_PROP, &
                  X_ANG_PART_PROP, Y_ANG_PART_PROP, Z_ANG_PART_PROP /)
                  
  integer :: blk_refine_level
  
  real :: gas_mean_vel_x, gas_mean_vel_y, gas_mean_vel_z, &
          gas_var_vel_x, gas_var_vel_y, gas_var_vel_z, gas_mean_cs
  integer :: counter

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  if (MDIM .ne. 3) call Driver_abortFlash('Sink particles only work in three spatial dimensions')

  if (first_call) then

    MyPE = dr_globalMe
    MasterPE = MASTER_PE

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

    call RuntimeParameters_get("sink_density_thresh", density_thresh)
    call RuntimeParameters_get("sink_accretion_radius", accretion_radius)
    call RuntimeParameters_get("sink_GasAccretionChecks", GasAccretionChecks)
    call RuntimeParameters_get("sink_convergingFlowCheck", convergingFlowCheck)
    call RuntimeParameters_get("sink_potentialMinCheck", potentialMinCheck)
    call RuntimeParameters_get("sink_jeansCheck", jeansCheck)
    call RuntimeParameters_get("sink_negativeEtotCheck", negativeEtotCheck)

    call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

    lrefine_max = gr_maxRefine

    call RuntimeParameters_get("Nblockx", Nblockx)
    call RuntimeParameters_get("Nblocky", Nblocky)
    call RuntimeParameters_get("Nblockz", Nblockz)

    call RuntimeParameters_get("xmin", xmin)
    call RuntimeParameters_get("xmax", xmax)
    call RuntimeParameters_get("ymin", ymin)
    call RuntimeParameters_get("ymax", ymax)
    call RuntimeParameters_get("zmin", zmin)
    call RuntimeParameters_get("zmax", zmax)

    delta_at_lrefmax = max ( &
            (xmax-xmin)/real(Nblockx)/(2**(real(lrefine_max)+(log(real(NXB))/log(2.)-1.))), &
            (ymax-ymin)/real(Nblocky)/(2**(real(lrefine_max)+(log(real(NYB))/log(2.)-1.))), &
            (zmax-zmin)/real(Nblockz)/(2**(real(lrefine_max)+(log(real(NZB))/log(2.)-1.))) )

    if (MyPE .eq. MasterPE) write(*,'(A,F6.2,A)') 'SinkParticles: You have set the sink particle accretion radius to ', &
                            & accretion_radius/delta_at_lrefmax, ' * (1+redshift) cells at the highest level of refinement.'

    if (accretion_radius/delta_at_lrefmax .LT. 2 .OR. accretion_radius/delta_at_lrefmax .GT. 3) then
       if (MyPE .eq. MasterPE) write(*,'(A)') &
             & 'CAUTION: Sink particle accretion radius is not within the recommended range (2-3 cells)!'
       if (MyPE .eq. MasterPE) write(*,'(A)') '         Sink particle creation checks might fail!'
    endif

    call RuntimeParameters_get("sink_merging", sink_merging)
    if (sink_merging .and. (MyPE .eq. MasterPE)) &
          write(*,'(A)') 'SinkParticles: Sink particles are allowed to merge.'

    call pt_sinkGatherGlobal()
!    call pt_gatherGlobal()

!    In our version where all particles share this tag method, we can't
!    call this only on sinks, since for restarts or starts with particles it messes
!    up the tag tracking for all particles. Making this called on all particles array. - JW
!    Here this is called on particles_global... why? Is the fix to just make a global array
!    for regular particles array?
!    local_tag_number = 0
!    do lp = 1, pt_numGlobal
!        if (get_ppe(int(allproc_particles(iptag,lp))) .EQ. MyPE) then
!           local_tag_number = max(local_tag_number, get_pno(int(allproc_particles(iptag,lp))))
!        endif
!    enddo

    if (write_accretion_checks_info) then
      open(funit_accretion_checks, file='sinks_accretion_checks_info.dat', position='APPEND')
      if (MyPE == MasterPE) write(funit_accretion_checks,'(8(1X,A14))') 'part_tag', 'time', 'dmass', &
        & 'distance', 'v_rad_of_dmass', 'etot_of_dmass', 'egrav_of_dmass', 'ekin_of_dmass'
      close(funit_accretion_checks)
      open(funit_accretion, file='sinks_accretion_info.dat', position='APPEND')
      if (MyPE == MasterPE) write(funit_accretion,'(6(1X,A14))') 'part_tag', 'time', 'dmass', &
        & 'distance', 'v_rad_of_dmass', 'etot_of_dmass'
      close(funit_accretion)
    endif

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

    first_call = .false.

  end if

   if (debug .and. (dr_globalMe .eq. MASTER_PE)) then
      print*, "Particles_sinkCreateAccrete: entering."
   end if

  ! blockList,blockCount used to be passed in as args but not anymore
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  
  call Cosmology_getRedshift(redshift)
  onePlusRedshift = 1.0 + redshift
  onePlusRedshift2 = onePlusRedshift * onePlusRedshift
  onePlusRedshift3 = onePlusRedshift2 * onePlusRedshift

  ! Convert accretion_radius and density_thresh to comoving coordinates
  accretion_radius_comoving = accretion_radius * onePlusRedshift
  density_thresh_comoving = density_thresh / onePlusRedshift3

  call Logfile_stamp(localnpf, "[SinkParticles]: localnpf now")

  if (pt_reduceGcellFills) then
     call Grid_fillGuardCells(CENTER, ALLDIR,unitReadsMeshDataOnly=.true.)
     ! We shall call Grid_notifySolnDataUpdate below if unitReadsMeshDataOnly turns
     ! out not to be appropriate, after all. - KW
  else
     call Grid_fillGuardCells(CENTER, ALLDIR)
  end if

  ! update particle's cpu info
  particles_local(ipcpu, 1:localnp) = MyPE

  call pt_sinkGatherGlobal(gather_propinds_pos, gather_nprops_pos)

  ! ======== START creation loop ========

  mass = 0.0

  time = dr_simTime

  llb = 0

  ! loop over leaf blocks (note that passed blockList only contains leafs)
  do lb = 1, blockCount

        blockID = blockList(lb)
        
        call Grid_getBlkRefineLevel(blockID, blk_refine_level)
        ! if this block isn't max refined, cycle to the next one. - JW
        if (blk_refine_level .ne. lrefine_max) cycle

        call Grid_getBlkPtr(blockID,solnData)

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

        allocate(xc(size_x))
        allocate(yc(size_y))
        allocate(zc(size_z))

        affected_block = .false.
        call Grid_getCellCoords(iXcoord, blockID, izn, .true., xc, size_x)
        call Grid_getCellCoords(iYcoord, blockID, izn, .true., yc, size_y)
        call Grid_getCellCoords(iZcoord, blockID, izn, .true., zc, size_z)

        call Grid_getBlkPhysicalSize(blockID,size)
        dx_block = size(1)/real(NXB)
        dy_block = size(2)/real(NYB)
        dz_block = size(3)/real(NZB)
        dVol = dx_block*dy_block*dz_block

        ! loop over cells (not including guard cells)
        do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef ALFL_MSCALAR
                 ! in case we are using the Alfven limiter
                 comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(idens,ip,jp,kp)
#else
                 comovingCellDens = solnData(idens,ip,jp,kp)
#endif
                 if (comovingCellDens .gt. density_thresh_comoving) then

                    create_sink = .true.   ! for now...

                    ! Is there an existing particle in range?
                    ! Looping over global particles
                    do pno = 1, localnpf
                       distx = xc(ip) - particles_global(ipx,pno)
                       disty = yc(jp) - particles_global(ipy,pno)
                       distz = zc(kp) - particles_global(ipz,pno)
                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
                       rad = sqrt(distx**2 + disty**2 + distz**2 )
                       ! Does this position fall within 2x accretion radius of existing sink?
                       if (rad .le. 2.0*accretion_radius_comoving) then
                          create_sink = .false.
#ifdef debug1
                          if(abs(xc(ip)) .lt. 2.0*dx_block .and. &
                             abs(yc(jp)) .lt. 2.0*dy_block .and. &
                             abs(zc(kp)) .lt. 2.0*dz_block) &
                            print*, "Failed neighbor sink check."
#endif
                       end if
                    end do

                    if (create_sink) then

                       ! just use the isothermal sound speed as an estimate for the v_rad-check
                       cs = sqrt(solnData(ipres,ip,jp,kp)/comovingCellDens)

                       ! check for converging flow in all surrounding cells
                       if (convergingFlowCheck) then
                          do k1 = -1, 1
                             do j1 = -1, 1
                                do i1 = -1, 1

                                   rad = sqrt(real(i1*i1+j1*j1+k1*k1))
                                   if (rad .GT. 0.) then
                                      vrad = ( i1*(solnData(ivelx, ip+i1, jp+j1, kp+k1)-solnData(ivelx, ip, jp, kp)) + &
                                               j1*(solnData(ively, ip+i1, jp+j1, kp+k1)-solnData(ively, ip, jp, kp)) + &
                                               k1*(solnData(ivelz, ip+i1, jp+j1, kp+k1)-solnData(ivelz, ip, jp, kp)) ) / rad
                                      if (vrad .GT. 1.e-5*cs) then ! a surrounding cell diverges, so do not create sink
                                         create_sink = .false.
#ifdef debug1
                                         if(abs(xc(ip)) .lt. 2.0*dx_block .and. &
                                            abs(yc(jp)) .lt. 2.0*dy_block .and. &
                                            abs(zc(kp)) .lt. 2.0*dz_block) &
                                         print*, "Failed converging flow check."
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
                             gpot = solnData(igpot,ip,jp,kp)
                             absgpot = abs(gpot)
                             do k1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                do j1 = -ngc_sink_creation_check, ngc_sink_creation_check
                                   do i1 = -ngc_sink_creation_check, ngc_sink_creation_check

                                      if (i1**2 + j1**2 + k1**2 .le. ngc_sink_creation_check_radius_sqr) then

                                         if ( ((solnData(igpot,ip+i1,jp+j1,kp+k1)-gpot)/absgpot) .lt. -1.e-5) then
                                            create_sink = .false.
#ifdef debug1
                                         if(abs(xc(ip)) .lt. 2.0*dx_block .and. &
                                            abs(yc(jp)) .lt. 2.0*dy_block .and. &
                                            abs(zc(kp)) .lt. 2.0*dz_block) &
                                         print*, "Failed potential minimun check."
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
                             maxgpot_in_vol = solnData(igpot,ip,jp,kp)
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
                                                    solnData(idens,ii1,jj1,kk1)
#else
                                         gas_dens = solnData(idens,ii1,jj1,kk1)
#endif
                                         etherm_in_vol = etherm_in_vol + solnData(ieint,ii1,jj1,kk1)*gas_dens 

                                         vxcm_in_vol = vxcm_in_vol + solnData(ivelx,ii1,jj1,kk1)*gas_dens
                                         vycm_in_vol = vycm_in_vol + solnData(ively,ii1,jj1,kk1)*gas_dens
                                         vzcm_in_vol = vzcm_in_vol + solnData(ivelz,ii1,jj1,kk1)*gas_dens

                                         mass_in_vol = mass_in_vol + gas_dens

                                         if (solnData(igpot,ii1,jj1,kk1) .gt. maxgpot_in_vol) & 
                                              maxgpot_in_vol = solnData(igpot,ii1,jj1,kk1)

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
                                                    solnData(idens,ii1,jj1,kk1)
#else
                                         gas_dens = solnData(idens,ii1,jj1,kk1)
#endif
                                         ekindisp_in_vol = ekindisp_in_vol + gas_dens * &
                                              ( (solnData(ivelx, ii1, jj1, kk1) - vxcm_in_vol)**2 + &
                                                (solnData(ively, ii1, jj1, kk1) - vycm_in_vol)**2 + &
                                                (solnData(ivelz, ii1, jj1, kk1) - vzcm_in_vol)**2  )

                                         egravdeltapot_in_vol = egravdeltapot_in_vol + & 
                                              (solnData(igpot, ii1, jj1, kk1) - maxgpot_in_vol)*gas_dens

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
                                         emag_in_vol = emag_in_vol + solnData(imagx, ii1, jj1, kk1)**2 + &
                                                                     solnData(imagy, ii1, jj1, kk1)**2 + &
                                                                     solnData(imagz, ii1, jj1, kk1)**2
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
#ifdef debug1
                                if(abs(xc(ip)) .lt. 2.0*dx_block .and. &
                                   abs(yc(jp)) .lt. 2.0*dy_block .and. &
                                   abs(zc(kp)) .lt. 2.0*dz_block) &
                                print*, "Failed Jeans check."
#endif
                                endif
                             end if

                             ! total energy should be negative (see e.g., Bate Bonnell Price 1995)
                             if (negativeEtotCheck) then
                                if (create_sink) then
                                   if (etherm_in_vol + ekindisp_in_vol + emag_in_vol .GT. egravdeltapot_in_vol) then
                                      create_sink = .false.
#ifdef debug1
                                   if(abs(xc(ip)) .lt. 2.0*dx_block .and. &
                                      abs(yc(jp)) .lt. 2.0*dy_block .and. &
                                      abs(zc(kp)) .lt. 2.0*dz_block) &
                                   print*, "Failed total energy check."
#endif
                                   endif
                                endif
                             endif

                          end if ! energy check

                       end if

                       ! finally create the sink in the cell centre
                       if (create_sink) then
                       
#ifdef debug1
                          print*, "Success!"
#endif

                          x = xc(ip)
                          y = yc(jp)
                          z = zc(kp)
                          pt = time

                          pno = pt_sinkCreateParticle(x, y, z, pt, blockID, MyPE)

                          write(*,'(A,4(1X,ES16.9),3I8)') "sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
                             & x, y, z, pt, blockID, MyPE, int(particles_local(iptag,pno))

                       end if

                    end if

                    affected_block = .true.
#ifdef debug2
                 else
                   print*, "Failed density check."
#endif

                 end if    ! cell density > max gas density

              end do
           end do
        end do

        if (affected_block) then
           llb = llb+1
           block_list(llb) = blockID
        end if

        call Grid_releaseBlkPtr(blockID, solnData)

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)

  end do ! block loop

  ! ======== END creation loop ========

  ! call this here instead of in pt_sinkMergingAfterCreation(), in preparation for pt_sinkFindList()
  ! and accumulation below, involving particles_global;
  ! NOTE THAT ALL particle properties must be communicated here.
  call pt_sinkGatherGlobal()

  ! Merges sink particles that were created close to one another
  call pt_sinkMergingAfterCreation(delta_at_lrefmax)

  ! ======== START accretion loop ========

  llnblocks = llb

  old_localnp = localnp

  ! clear mass & velocity
  tot_mass(:)   = 0.
  cm_x(:)       = 0.
  cm_y(:)       = 0.
  cm_z(:)       = 0.
  vel_cm_x(:)   = 0.
  vel_cm_y(:)   = 0.
  vel_cm_z(:)   = 0.
  ang_x(:)      = 0.
  ang_y(:)      = 0.
  ang_z(:)      = 0.

  ! do it again, but only loop over affected blocks and
  ! add mass to particles

  do llb = 1, llnblocks

     lb = block_list(llb)

     call Grid_getBlkPtr(lb, solnData)

     call Grid_getBlkPhysicalSize(lb, size)
     dx_block = size(1)/real(NXB)
     dy_block = size(2)/real(NYB)
     dz_block = size(3)/real(NZB)
     dVol = dx_block*dy_block*dz_block

     call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

     size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

     allocate(xc(size_x))
     allocate(yc(size_y))
     allocate(zc(size_z))

     call Grid_getCellCoords(iXcoord, lb, izn, .true., xc, size_x)
     call Grid_getCellCoords(iYcoord, lb, izn, .true., yc, size_y)
     call Grid_getCellCoords(iZcoord, lb, izn, .true., zc, size_z)

     affected_block = .false.

     do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! Let's start to accrete mass

#ifdef ALFL_MSCALAR
              ! in case we are using the Alfven limiter
              comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(idens,ip,jp,kp)
#else
              comovingCellDens = solnData(idens,ip,jp,kp)
#endif

              if (comovingCellDens .gt. density_thresh_comoving) then

                 mass = (comovingCellDens - density_thresh_comoving) * dVol
                 cs = sqrt(solnData(ipres,ip,jp,kp) / comovingCellDens)

                 ! return a list 'pindex_found' containing all particles found within r_accr
                 ! all affected particles are in local list 'particles'
                 ! if necessary extend the list with a dummy particle (TRUE switch)
                 ! if TRUE, then particles_local also contains the global particles which fall
                 ! within accretion_radius

                 r_search = accretion_radius_comoving

                 call pt_sinkFindList(xc(ip), yc(jp), zc(kp), r_search, .true., pindex_found, np_found)

                 if (np_found .gt. 0) then
                    ! there is a particle within accretion_radius of this cell

                    do npf = 1, np_found
                       ! loop over all particles within accretion_radius of this cell

                       pno = pindex_found(npf)

                       dvx = solnData(ivelx,ip,jp,kp) - particles_local(ipvx, pno)
                       dvy = solnData(ively,ip,jp,kp) - particles_local(ipvy, pno)
                       dvz = solnData(ivelz,ip,jp,kp) - particles_local(ipvz, pno)

                       distx = xc(ip) - particles_local(ipx, pno)
                       disty = yc(jp) - particles_local(ipy, pno)
                       distz = zc(kp) - particles_local(ipz, pno)

                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)

                       radius(npf) = sqrt(distx**2 + disty**2 + distz**2)

                       egrav_Gas = -Newton*2.0 * PI / 3.0 * density_thresh_comoving * &
                          & (accretion_radius_comoving**2 - radius(npf)**2)*mass

                       if (radius(npf) .gt. 0) then
                          egrav_part = -Newton*particles_local(ipm, pno)*mass * &
                                        (1.0/radius(npf)-1./accretion_radius_comoving)
                       else
                          egrav_part = -1.0e99
                       end if

                       egrav = egrav_gas + egrav_part
                       ! CTSS - factor of (1+z)^3 needed for correct comoving potential
                       egrav = egrav * onePlusRedshift3
                       ekin = 0.5 * mass * ( dvx**2 + dvy**2 + dvz**2 )

                       etot(npf) = egrav + ekin

                       ! calculate the radial velocity wrt each particle found
                       if (radius(npf) .gt. 0.) then
                          vr(npf) = ( distx*dvx + disty*dvy + distz*dvz ) / radius(npf)
                       else
                          vr(npf) = 0.
                       endif

                    end do

                    pno_to_accrete = 0

                    inner_r_accr = max( 0.2*accretion_radius_comoving, (dVol**(1.0/3.0)) )
                    within_inner_r_accr = .false.
                    etot_min_inner_r_accr = 1.0e99
                    etot_min = 1.0e99

                    do npf = 1, np_found

                       if (radius(npf) .lt. inner_r_accr) then
                          if (etot(npf) .lt. etot_min_inner_r_accr) then
                             pno_to_accrete = pindex_found(npf)
                             pno = npf
                             etot_min_inner_r_accr = etot(npf)
                             within_inner_r_accr = .true.
                          end if
                       else

!                         bug fix by Andrea Gatto 150515
                          if (GasAccretionChecks)then

                            if (vr(npf) .lt. 1.0e-5*cs .and. etot(npf) .lt. 0. .and. etot(npf) &
                                .lt. etot_min_inner_r_accr) then

                              pno_to_accrete = pindex_found(npf)
                              pno = npf
                              etot_min_inner_r_accr = etot(npf)
                              within_inner_r_accr = .true.

                            end if

                          else

                            if (etot(npf) .lt. etot_min_inner_r_accr) then

                              pno_to_accrete = pindex_found(npf)
                              pno = npf
                              etot_min_inner_r_accr = etot(npf)
                              within_inner_r_accr = .true.
                              end if

                          endif
!                          if (GasAccretionChecks) then
!
!                             if (.not. within_inner_r_accr .and. vr(npf) .lt. 1.0e-5*cs &
!                                  .and. etot(npf) .lt. 0. .and. etot(npf) .lt. etot_min) then
!                                pno_to_accrete = pindex_found(npf)
!                                pno = npf
!                                etot_min = etot(npf)
!                             end if
!
!                          else
!
!                             if (.not. within_inner_r_accr .and. (etot(npf) .lt. etot_min)) then
!                                pno_to_accrete = pindex_found(npf)
!                                pno = npf
!                                etot_min = etot(npf)
!                             end if
!
!                          end if    ! perform gas accretion checks?

                       end if   ! inner accretion?

                    end do   ! potential sinks

                    if (pno_to_accrete .gt. 0) then

                        dvx = solnData(ivelx,ip,jp,kp) - particles_local(ipvx, pno_to_accrete)
                        dvy = solnData(ively,ip,jp,kp) - particles_local(ipvy, pno_to_accrete)
                        dvz = solnData(ivelz,ip,jp,kp) - particles_local(ipvz, pno_to_accrete)

                        distx = xc(ip) - particles_local(ipx, pno_to_accrete)
                        disty = yc(jp) - particles_local(ipy, pno_to_accrete)
                        distz = zc(kp) - particles_local(ipz, pno_to_accrete)

                        if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
#ifdef ALFL_MSCALAR
                        ! in case we are using the Alfven limiter
                        alfl_dens = solnData(ALFL_MSCALAR,ip,jp,kp) * solnData(idens,ip,jp,kp)
                        solnData(idens,ip,jp,kp) = density_thresh_comoving + alfl_dens
                        solnData(ALFL_MSCALAR,ip,jp,kp) = alfl_dens / solnData(idens,ip,jp,kp)
#else
                        solnData(idens,ip,jp,kp) = density_thresh_comoving
#endif
                        affected_block = .true.

                        tot_mass(pno_to_accrete) = tot_mass(pno_to_accrete) + mass
                        cm_x(pno_to_accrete)     = cm_x(pno_to_accrete) + distx*mass
                        cm_y(pno_to_accrete)     = cm_y(pno_to_accrete) + disty*mass
                        cm_z(pno_to_accrete)     = cm_z(pno_to_accrete) + distz*mass
                        vel_cm_x(pno_to_accrete) = vel_cm_x(pno_to_accrete) + dvx*mass
                        vel_cm_y(pno_to_accrete) = vel_cm_y(pno_to_accrete) + dvy*mass
                        vel_cm_z(pno_to_accrete) = vel_cm_z(pno_to_accrete) + dvz*mass
                        ang_x(pno_to_accrete)    = ang_x(pno_to_accrete) + &
                             & (disty*solnData(ivelz,ip,jp,kp)-distz*solnData(ively,ip,jp,kp))*mass
                        ang_y(pno_to_accrete)    = ang_y(pno_to_accrete) + &
                             & (distz*solnData(ivelx,ip,jp,kp)-distx*solnData(ivelz,ip,jp,kp))*mass
                        ang_z(pno_to_accrete)    = ang_z(pno_to_accrete) + &
                             & (distx*solnData(ively,ip,jp,kp)-disty*solnData(ivelx,ip,jp,kp))*mass

                    end if

                  end if

               end if

               ! End of mass accretion

            end do      ! i
         end do      ! j
      end do      ! k

      if(affected_block) then
         call Eos_wrapped(MODE_DENS_EI, blkLimits, lb)
      end if

      call Grid_releaseBlkPtr(lb, solnData)

      deallocate(xc)
      deallocate(yc)
      deallocate(zc)

   end do

   ! ======== END accretion loop ========

   ! Copy grid info to dummy particle for data exchange
   do lp = old_localnp+1, localnp
      particles_local(ipm,lp) = tot_mass(lp)
      particles_local(ipx,lp) = cm_x(lp)
      particles_local(ipy,lp) = cm_y(lp)
      particles_local(ipz,lp) = cm_z(lp)
      particles_local(ipvx,lp) = vel_cm_x(lp)
      particles_local(ipvy,lp) = vel_cm_y(lp)
      particles_local(ipvz,lp) = vel_cm_z(lp)
      particles_local(iplx,lp) = ang_x(lp)
      particles_local(iply,lp) = ang_y(lp)
      particles_local(iplz,lp) = ang_z(lp)
   end do

   ! Exchange information across CPUs
   call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

   npart = localnp

   ! delete dummy (non-local) particles from list
   ! do this because pt_sinkFindList with .true. call raises localnp
   do lp = 1, localnp
      if (int(particles_local(ipcpu,lp)) .ne. MyPE) then
         npart = npart-1
      end if
   end do
   localnp = npart

   do lp = 1, localnp

      ! check if local particle is affected by regions on other CPUs
      do nlp = localnp+1, localnpf

         if(int(particles_local(iptag,lp)) .eq. int(particles_global(iptag,nlp))) then
            tot_mass(lp) = tot_mass(lp) + particles_global(ipm,nlp)
            cm_x(lp) = cm_x(lp) + particles_global(ipx,nlp)
            cm_y(lp) = cm_y(lp) + particles_global(ipy,nlp)
            cm_z(lp) = cm_z(lp) + particles_global(ipz,nlp)
            vel_cm_x(lp) = vel_cm_x(lp) + particles_global(ipvx,nlp)
            vel_cm_y(lp) = vel_cm_y(lp) + particles_global(ipvy,nlp)
            vel_cm_z(lp) = vel_cm_z(lp) + particles_global(ipvz,nlp)
            ang_x(lp) = ang_x(lp) + particles_global(iplx,nlp)
            ang_y(lp) = ang_y(lp) + particles_global(iply,nlp)
            ang_z(lp) = ang_z(lp) + particles_global(iplz,nlp)
         end if
      end do

      ! update particle properties (conservation laws)

      particles_local(iold_pmass,lp) = particles_local(ipm,lp)

      if (tot_mass(lp) .ne. 0.0) then

         ! mass update
         particles_local(ipm,lp) = particles_local(ipm,lp) + tot_mass(lp)

         ! position update
         px_old = particles_local(ipx,lp)
         py_old = particles_local(ipy,lp)
         pz_old = particles_local(ipz,lp)
         particles_local(ipx,lp) = particles_local(ipx,lp) + cm_x(lp)/particles_local(ipm,lp)
         particles_local(ipy,lp) = particles_local(ipy,lp) + cm_y(lp)/particles_local(ipm,lp)
         particles_local(ipz,lp) = particles_local(ipz,lp) + cm_z(lp)/particles_local(ipm,lp)

         ! velocity update
         particles_local(ipvx,lp) = particles_local(ipvx,lp) + vel_cm_x(lp)/particles_local(ipm,lp)
         particles_local(ipvy,lp) = particles_local(ipvy,lp) + vel_cm_y(lp)/particles_local(ipm,lp)
         particles_local(ipvz,lp) = particles_local(ipvz,lp) + vel_cm_z(lp)/particles_local(ipm,lp)

         ! spin update
         particles_local(iplx_old,lp) = particles_local(iplx,lp)
         particles_local(iply_old,lp) = particles_local(iply,lp)
         particles_local(iplz_old,lp) = particles_local(iplz,lp)
         particles_local(iplx,lp) = particles_local(iplx,lp) + ang_x(lp) - particles_local(ipm,lp) * &
                                    ( (particles_local(ipy,lp)-py_old)*particles_local(ipvz,lp) - &
                                      (particles_local(ipz,lp)-pz_old)*particles_local(ipvy,lp)  )
         particles_local(iply,lp) = particles_local(iply,lp) + ang_y(lp) - particles_local(ipm,lp) * &
                                    ( (particles_local(ipz,lp)-pz_old)*particles_local(ipvx,lp) - &
                                      (particles_local(ipx,lp)-px_old)*particles_local(ipvz,lp)  )
         particles_local(iplz,lp) = particles_local(iplz,lp) + ang_z(lp) - particles_local(ipm,lp) * &
                                    ( (particles_local(ipx,lp)-px_old)*particles_local(ipvy,lp) - &
                                      (particles_local(ipy,lp)-py_old)*particles_local(ipvx,lp)  )
      end if ! mass was accreted

      particles_local(ipmdot,lp) = (particles_local(ipm,lp) - particles_local(iold_pmass,lp)) / dt

   end do

   lp = 1
   do while (lp .le. localnp)
      if (particles_local(ipm,lp) .le. 0.0) then
         print*, "SinkParticles: deleted particle due to zero mass"
         particles_local(:,lp) = particles_local(:,localnp)
         particles_local(ipblk,localnp) = NONEXISTENT
         n_empty = n_empty + 1
         localnp = localnp - 1
         lp = lp - 1
      end if
      lp = lp + 1
   end do

   if (sink_merging) call pt_sinkParticleMerging(dt)

   ! write sink particle data to sinks_evol.dat
   call pt_sinkDumpParticles(time, "sinks_evol.dat")

   ! This is needed to signal the particle unit that particle
   ! positions have been initialized, which is important for
   ! Particles_updateRefinement() to actually move particle to
   ! to the right blocks and processors.
   if (localnpf .gt. 0) then
      pt_posInitialized = .true.
      pt_velInitialized = .true.
   endif
  
! Calculate the sample mean velocity and sample velocity variance of the
! surrounding gas cells that are under 100 K. We will use this to set 
! stellar velocities when we make stars from the sink mass. - JW

  do lp = 1, localnp
  
    gas_mean_cs    = 0.0
    
    gas_mean_vel_x = 0.0
    gas_mean_vel_y = 0.0
    gas_mean_vel_z = 0.0
    
    gas_var_vel_x = 0.0
    gas_var_vel_x = 0.0
    gas_var_vel_z = 0.0

    blockID = particles_local(BLK_PART_PROP, lp)

    call Grid_getBlkPtr(blockID, solnData)
    call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
    
    size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
    size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
    size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

    allocate(xc(size_x))
    allocate(yc(size_y))
    allocate(zc(size_z))

    call Grid_getCellCoords(iXcoord, blockID, izn, .true., xc, size_x)
    call Grid_getCellCoords(iYcoord, blockID, izn, .true., yc, size_y)
    call Grid_getCellCoords(iZcoord, blockID, izn, .true., zc, size_z)
    
    counter = 0
    ! Loop over cells (including guard cells, since here we start with
    ! the particle location and look for cells, unlike above).
     do kp = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do jp = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do ip = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

             distx = xc(ip) - particles_local(ipx, lp)
             disty = yc(jp) - particles_local(ipy, lp)
             distz = zc(kp) - particles_local(ipz, lp)

             if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)

             rad = sqrt(distx**2 + disty**2 + distz**2 )
             
             ! Look within double the accr radius to calculate gas velocity variance.
             if (rad .le. 2.0*accretion_radius_comoving) then

#ifdef ALFL_MSCALAR
              ! in case we are using the Alfven limiter
              comovingCellDens = (1.0-solnData(ALFL_MSCALAR,ip,jp,kp))*solnData(idens,ip,jp,kp)
#else
              comovingCellDens = solnData(idens,ip,jp,kp)
#endif

!                 do k1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
!                    do j1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
!                       do i1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
                   
                         if (solnData(TEMP_VAR,ip,jp,kp) .le. 100.0) then 
                   
                         counter = counter + 1
                   
                         gas_mean_vel_x = gas_mean_vel_x + solnData(ivelx,ip,jp,kp)
                         gas_mean_vel_y = gas_mean_vel_y + solnData(ively,ip,jp,kp)
                         gas_mean_vel_z = gas_mean_vel_z + solnData(ivelz,ip,jp,kp)
                         gas_mean_cs    = gas_mean_cs + sqrt((5.0 / 3.0) * solnData(ipres,ip,jp,kp) / comovingCellDens)
                         end if
                     
!                       end do
!                    end do
!                 end do
             end if
           end do ! cell i in block
        end do ! cell j in block
     end do ! cell k in block
     
     gas_mean_vel_x = gas_mean_vel_x / counter
     gas_mean_vel_y = gas_mean_vel_y / counter
     gas_mean_vel_z = gas_mean_vel_z / counter
     gas_mean_cs    = gas_mean_cs / counter
     
     counter = 0
                 
     do kp = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do jp = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do ip = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

             distx = xc(ip) - particles_local(ipx, lp)
             disty = yc(jp) - particles_local(ipy, lp)
             distz = zc(kp) - particles_local(ipz, lp)

             if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)

             rad = sqrt(distx**2 + disty**2 + distz**2)

             if (rad .le. 2.0*accretion_radius_comoving) then
                                  
!                 do k1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
!                    do j1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
!                       do i1 = -4,4 !-ngc_sink_creation_check, ngc_sink_creation_check
                   
                         if (solnData(TEMP_VAR,ip,jp,kp) .le. 100.0) then 
                          
                          counter = counter + 1
                          
                          gas_var_vel_x = gas_var_vel_x + (solnData(ivelx,ip,jp,kp) - gas_mean_vel_x)**2.0
                          gas_var_vel_y = gas_var_vel_y + (solnData(ively,ip,jp,kp) - gas_mean_vel_y)**2.0
                          gas_var_vel_z = gas_var_vel_z + (solnData(ivelz,ip,jp,kp) - gas_mean_vel_z)**2.0
                     
                         end if
                     
!                       end do
!                    end do
!                 end do
             end if
           end do ! cell i in block
        end do ! cell j in block
     end do ! cell k in block
     
     gas_var_vel_x = gas_var_vel_x / (counter - 1)
     gas_var_vel_y = gas_var_vel_y / (counter - 1)
     gas_var_vel_z = gas_var_vel_z / (counter - 1)

     call Grid_releaseBlkPtr(lb, solnData)

     deallocate(xc)
     deallocate(yc)
     deallocate(zc)
     
     ! Set particle mean gas sound speed.
     particles_local(CSGM_PART_PROP, lp) = gas_mean_cs
     ! Set particle velocity gas mean.
     particles_local(VXGM_PART_PROP, lp) = gas_mean_vel_x
     particles_local(VYGM_PART_PROP, lp) = gas_mean_vel_y
     particles_local(VZGM_PART_PROP, lp) = gas_mean_vel_z
     ! Set particle velocity gas variance.
     particles_local(VXGV_PART_PROP, lp) = gas_var_vel_x
     particles_local(VYGV_PART_PROP, lp) = gas_var_vel_y
     particles_local(VZGV_PART_PROP, lp) = gas_var_vel_z

  end do ! local particle
  
   ! We don't have to leave with guardcells updated.
   ! call Grid_fillGuardCells(CENTER, ALLDIR)

   if (debug .and. (dr_globalMe .eq. MASTER_PE)) then
      print*, "Particles_sinkCreateAccrete: exiting, localnpf = ", localnpf
   end if

   return

end subroutine Particles_sinkCreateAccrete

