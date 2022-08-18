!!****if* source/physics/RadTrans/RadTransMain/RayRad/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!
!!      Modified to properly cycle over the ray tracing and ionization
!!      calculation here, instead of running the entirety of Flash on
!!      the very short fractional ionization change timesteps. Also,
!!      the ionization solver was made implicit so as to take larger
!!      timesteps overall. Finally, coupling to sinks was added. - Josh Wall 5-15-16
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!   pass   : reverses solve direction
!!
!!***

! this corresponds to ionization from FLASH2.5 without the raytracing 
!#ifdef ONE_CELL_TESTING
!#define DEBUG_RADTRANS
!#endif
#define timing

subroutine RadTrans3(nblk, blklst, dt, pass)

#include "Flash.h"

! extra unit specific data
  use rt_data, only : rt_protonMass, abu_c, rt_abar, rt_idealgas, rt_dt, rt_dt_pos, rt_maxHchange, &
                      rt_rayTrace, rt_vary_atomic_frac, rt_ion_threshold

! general radiation transport data
  use RadTrans_data,  ONLY : rt_meshMe, rt_globalComm
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getDeltas

  use Driver_interface, ONLY : Driver_abortFlash
  use Diffuse_interface, ONLY: Diffuse_solveScalar, Diffuse_fluxLimiter

  use RadTrans_data, ONLY: rt_useRadTrans
  use Eos_interface, ONLY: Eos_wrapped
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_rayAdvance, Particles_getGlobalNum
  use Driver_data, ONLY: dr_simTime, dr_globalcomm, dr_globalMe, dr_nStep

  use Logfile_interface, ONLY: Logfile_stamp ! Lets start recording some info.

! for subcycling internal energy
  use Heat_interface, ONLY : RadHeat
! for background CR ionization
  use Heat_data !, ONLY : he_crIonRate, he_crIonNH, he_crIonExp, he_use_cr_heating, he_subfactor, he_stratifyHeat, he_h_UV
  use heatCool !, ONLY : approx_column_dens, dei_dt
  use cool_vars
  use calc_ion
  use get_ion

  implicit none

#include "Flash_mpi.h"
#include "Eos.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  integer :: j, k, i, l
  real    :: xx, yy, zz
!  solndata
  real, pointer, dimension(:,:,:,:) :: solnData, solnDataCtr
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord
  integer                           :: xSizeCoord, ySizeCoord, zSizeCoord
  integer, dimension(2,MDIM)        :: blkLimits, blkLimitsGC
  logical                           :: getGuardCells = .true.

  integer :: tmpID

! for solving 
  real  :: numdens, xH1, xH0 !, tmp_dt <-Changing this name, its a bit confusing
  real  :: phih, dens, hvphih ! since its not an actual timestep. -JW
  real  :: store, eldens, temp, fac, check, sub_dt, ion_dt, ei, ek
  real  :: hc_dt, frac_dt
  real  :: hFracNew(0:1), xh(0:1)

! for source sink check
  real  :: del(3), x_dis, y_dis, z_dis, frac_change, fac_old ! frac_change is new tmp_dt 
  real  :: send_buff, rec_buff, pos_save(3) !send_buff(2), rec_buff(2) <- repurposed.
  real  :: dt_old, frac_old_save(2), frac_new_save(2)
  
! For new solver. - JW  
  real  :: total_timestep, xion_save, global_dt_save, largest_xion, &
           NumH0, phoio_old, phoio_new, del_phoio, phoio_save, &
           mu_old, xHp_old
  real, save :: dt_save=1d99
  real  :: t_start, t_stop, ion_time, heat_time, raytrace_time
  logical :: first_loop, converged, early_exit, all_done, global_all_done, fully_ionized, fileexists
  integer :: nsteps, global_nsteps, step, heatcool_nsteps
  integer :: p, ierr, stat(MPI_STATUS_SIZE), sendproc, dt_request, done_request
  character(len=MAX_STRING_LENGTH*2) :: strbuff
  integer :: globalNumParticles
  
!  real, parameter :: ion_threshold=1d-10, ion_frac_change_accept=1d-1

  !=========================================================================

! Don't ever trust a function written by someone else to ensure guardcells
! are filled. Always fill on entry to a routine.

  call Grid_fillGuardCells(CENTER, ALLDIR)

!  reset radiation timestep

!  rt_dt = min(0.1*rt_maxHchange*dt,dt_save) 
 rt_dt = dt

  total_timestep = 0.0
  pos_save = 0.0
!  fac_old = 1.0
  dt_old = 0.0
  dt_save = 1e99
  xion_save = 1.0d99
  del_phoio = 1.0d99
  phoio_old = 0.0d0
  phoio_new = 0.0d0
  phoio_save = 1.0d99
  nsteps = 0
!  global_nsteps = 0
  all_done  = .false.


#ifdef timing
ion_time = 0.0
heat_time = 0.0
raytrace_time = 0.0
#endif

  if(.not. rt_useRadTrans) return

call Logfile_stamp("Entering RadRay.", "[RadTrans]")

if (rt_rayTrace) call Particles_getGlobalNum(globalNumParticles)

#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  print*,'entering rad. solver'
#endif

  call Timers_start("RadTrans")

! Subcycle the timestep over the ray tracing method
! and the ionization calculation until we hit the
! hydro timestep or we hit convergence everywhere on the
! grid.

! This loop needs to cycle over the heating and cooling solver on timesteps
! based on the change in internal energy.
do while ((dt-total_timestep) .gt. (1d-6*dt))
  xion_save  = 1.0d99
  converged  = .false.
  first_loop = .true.

! Lets just run til convergence, or til the cows come home, to solve ionization fraction.
! Note all_done lets us make one more loop where we store the actual values.
  do while (.not. all_done) ! .or. phoio_save > 1d-3 .or. nsteps < 2)
  
! Clear previous fluxes and heating from ray tracing.
    do l = 1, nblk
      tmpID = blklst(l)
      call Grid_getBlkPtr(tmpID,solnData)
      solnData(UVFL_VAR,:,:,:) = 0d0
      solnData(FUFL_VAR,:,:,:) = 0d0
      ! Store the old photoionization heating rates for convergence check.
      !solnData(OPHH_VAR,:,:,:) = solnData(PHHE_VAR,:,:,:)
      solnData(PHHE_VAR,:,:,:) = 0d0
#ifdef PE_HEAT
      solnData(PEFL_VAR,:,:,:) = 0d0
#endif
      solnData(PHIO_VAR,:,:,:) = 0d0
      call Grid_releaseBlkPtr(tmpID,solnData)
    end do
  
! First step we try to go the whole dt.
!    converged = .false. ! If you converge, you get to go home!
!    early_exit = .false. ! If you're bad, you have to start over!
    xion_save = 0.0d0
    phoio_save = 0.0d0
    largest_xion = 0.0
    nsteps = nsteps + 1

! If we reset here, it allows the timestep to get bigger as the
! solution converges.
    sub_dt = 1d99
#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  print*,'entering raytracing'
#endif

! raytracing, should go in Particles_advance.F90 but then
! order in Driver is screwed up, as source term is calculated after call to hydro solver
! as long as the raytracing manipulates the data structure orderly then nothing should happen

! Actually we'd like to do this here, lets keep this together with
! solving dx_ion. - JW

! Actual actual radiation transport. -JW
  if (rt_rayTrace .and. (globalNumParticles > 0)) then
#ifdef DEBUG_RADTRANS    
    if (rt_Meshme == MASTER_PE) print*, "Calling ray tracing with dt=.", rt_dt
#endif

#ifdef timing
  t_start = MPI_Wtime()
#endif

    call Timers_start("raytracing")
    call Particles_rayAdvance(rt_dt) ! Just make this dt_hydro, not rt_dt
    call Timers_stop("raytracing")
    
#ifdef timing
  t_stop = MPI_Wtime()
  raytrace_time = raytrace_time + t_stop - t_start
#endif
    
  endif

#ifdef DEBUG_RADTRANS
  print*,'leaving raytracing'
#endif

!  call Timers_start("solving_ionization")

#ifdef timing
  t_start = MPI_Wtime()
#endif

!===========================================
! actual ionization calc... -JW
!===========================================
  block: do l = 1, nblk
    tmpID = blklst(l)

! allocate space for dimensions
    call Grid_getBlkPtr(tmpID,solnData)
!		call Grid_getBlkPtr(tmpID,solnDataCtr,SCRATCH_CTR)
    call Grid_getBlkIndexLimits(tmpID,blkLimits,blkLimitsGC)

    xSizeCoord = blkLimitsGC(HIGH,IAXIS)
    ySizeCoord = blkLimitsGC(HIGH,JAXIS)
    zSizeCoord = blkLimitsGC(HIGH,KAXIS)

    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))

    call Grid_getCellCoords(IAXIS,tmpID,CENTER,getGuardCells,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,tmpID,CENTER,getGuardCells,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,tmpID,CENTER,getGuardCells,zCoord,zSizeCoord)    
    
    call Grid_getDeltas(tmpID,del)

  ! loop over all zones in block
    ! for 2d k is just 1
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      zz = zCoord(k)
      if (he_stratifyHeat) then
        strat_factor = exp(-abs(zz)/(he_h_UV))
      else
        strat_factor = 1.0d0
      end if
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

! phot ionization rate. 
          phih   = solnData(PHIO_VAR,i,j,k)

! Set vars to approx rate of change in energy.
          rho    = solnData(DENS_VAR,i,j,k)
#ifdef VARY_ATM_FRAC
          mu_mol = solnData(ATMU_VAR,i,j,k) 
#endif
          ndens  = rho/(mu_mol*rt_protonMass)
          temp   = solnData(TEMP_VAR,i,j,k)
          ei     = solnData(EINT_VAR,i,j,k)
          tdust  = solnData(TDUS_VAR,i,j,k)
          xHp    = solnData(IHP_SPEC,i,j,k)
          ephen  = solnData(PHHE_VAR,i,j,k)
          TtoEI   = ei / temp
          conf    = ndens*ndens/rho
          
          ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
                    & solnData(VELY_VAR,i,j,k)**2 + &
                    & solnData(VELZ_VAR,i,j,k)**2)
                    
! convert to erg /(g s)
          ephen  = ephen/rho
#ifdef PE_HEAT
! Get the photoelectric flux from radiation sources. -JW
          Gflux = solnData(PEFL_VAR,i,j,k)
#else
          Gflux = 0.0
#endif
          xHp_old = xHp
          mu_old = mu_mol

!!! Now do the actual physics things.
          ! Ionization.
          xHp = get_ionization(xHp, ndens, rho, temp, phih, rt_dt)
          ! Atomic weight.
          ! Define mean molecular weight from temperature and ionization fraction
          ! since we don't actually follow the elements of the gas. - JW
#ifdef VARY_ATM_FRAC
          if (rt_vary_atomic_frac) then
              mu_mol = get_mu(xHp, ndens)
              mu_mol = (mu_mol+mu_old)/2.0 ! Beat down oscillations by averaging.
              solnData(ATMU_VAR,i,j,k) = mu_mol
          end if
#endif
!!! Done with physics things.
          if(xHp .gt. 1.01) then
            print*,phih,dens,temp,xHp
            print*,'ion wrong'
            call flush(6)
            stop
          endif
          
          frac_change = abs(xHp - xHp_old) !/min(xHp,xHp_old)

! Look I only care about the dangerous instability due to rapidly varying radiation fields
! at very high ionization. Therefore, I'm going to try and ignore timestepping criterion for
! everyone else and I think this is safe (and will be much faster).
        if (ephen .gt. 0d0) then
! Now get the estimated heating and cooling timestep from ei /  d(ei)/dt.
          hc_dt   = he_subfactor * ei / max(abs(dei_dt(0.0,ei)), 1d-50)          
          frac_dt = rt_dt * (rt_maxHchange / frac_change)
        
#ifdef ONE_CELL_TESTING
          write(*,'(A,ES18.6E3,I4)') 'ei / de_dt   =', hc_dt, dr_globalMe
          write(*,'(A,ES18.6E3,I4)') 'frac_dt      =', frac_dt, dr_globalMe
          write(*,'(A,ES18.6E3,I4)') 'first sub_dt =', sub_dt, dr_globalMe
#endif
          sub_dt  = min(sub_dt,frac_dt)
          !sub_dt  = min(sub_dt,ion_dt,frac_dt)
          !sub_dt  = min(sub_dt,hc_dt,frac_dt)
          !sub_dt  = min(sub_dt,ion_dt,hc_dt,frac_dt)
        else
          ion_dt  = 1d99
          hc_dt   = 1d99
          frac_dt = 1d99
        end if

          largest_xion = max(largest_xion, xHp)
          xion_save = max(xion_save, frac_change)
          
! update the ionisation state
          solnData(IHA_SPEC,i,j,k) = 1.0d0 - xHp
          solnData(IHP_SPEC,i,j,k) = xHp
          
! Only store the new state solutions if we are actually converged.

          if (converged .and. .not. first_loop) then
            ! Heating and cooling.
            call heating_and_cooling(ei, rt_dt, temp, heatcool_nsteps)
            solnData(ENER_VAR,i,j,k) = ei + ek
            solnData(EINT_VAR,i,j,k) = ei
            solnData(TEMP_VAR,i,j,k) = temp
            solnData(TDUS_VAR,i,j,k) = tdust
            all_done = .true.
          end if
          
        enddo ! coord loops
      enddo
    enddo

!  clean up memory 
    call Grid_releaseBlkPtr(tmpID,solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)

!===========================================
! call to EOS to set zone
!===========================================
! no call to eos hydrogen does not partake in hydrodynamics it is too good for that.
    call Eos_wrapped(MODE_DENS_EI, blkLimits, tmpID)

  enddo block ! block

#ifdef timing
  t_stop = MPI_Wtime()
  ion_time = ion_time + t_stop - t_start
#endif

!  call Timers_stop("solving_ionization")
  
  ! Unfortunately we can't get around a blocking allreduce here. Its
  ! important that if the ionization changes enough to allow radiation
  ! to travel from one processor to another that it is caught, and therefore
  ! all processors are going to have to take the same sized steps. Honestly
  ! though, this is still much cheaper than calling all of Driver_evolveFlash
  ! on small steps. - JW

#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  write(*,*) "[RadTrans]: Reducing sub_dt.", dr_globalMe
#endif

  call MPI_allReduce(MPI_IN_PLACE, xion_save, 1, FLASH_REAL, &
                   MPI_MAX, dr_globalcomm, ierr)
!  call MPI_allReduce(MPI_IN_PLACE, phoio_save, 1, FLASH_REAL, &
!                   MPI_MAX, dr_globalcomm, ierr)

  if (xion_save .lt. rt_ion_threshold) converged = .true.

#ifdef DEBUG_RADTRANS
  call MPI_allReduce(MPI_IN_PLACE, largest_xion, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)
    write(*,'(A,ES18.6E3,I4)') 'rt_dt =', rt_dt, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'sub_dt =', sub_dt, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'total_timestep =', total_timestep, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'hydro timestep =', dt, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'highest xion =', largest_xion, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'final xion =', hFracNew(1), dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'final neutral =', hFracNew(0), dr_globalMe
    write(*,'(A,ES18.6E3,I4)') 'highest xion change =', xion_save, dr_globalMe
    !write(*,'(A,ES18.6E3,I4)') 'highest phio change =', phoio_save, dr_globalMe
    write(*,'(A,ES18.6E3,I4)') "nsteps =", real(nsteps), dr_globalMe
    call flush(6)
#endif

  if (first_loop) first_loop = .false.

  end do ! End ionization convergence loop. 

  call MPI_allReduce(MPI_IN_PLACE, sub_dt, 1, FLASH_REAL, &
                   MPI_MIN, dr_globalcomm, ierr)
                   
!! Now subcycle heating and cooling. - JW
!#ifdef DEBUG_RADTRANS
!!if (rt_meshMe .eq. MASTER_PE) &
!  write(*,'(A,ES12.3E3,I4)') "[RadTrans]: Calling RadHeat with dt =", rt_dt, dr_globalMe
!#endif

!#ifdef timing
!  t_start = MPI_Wtime()
!#endif
!  call RadHeat(nblk, blklst, rt_dt, total_timestep)
!#ifdef timing
!  t_stop = MPI_Wtime()
!  heat_time = heat_time + t_stop - t_start
!#endif

!#ifdef DEBUG_RADTRANS
!!if (rt_meshMe .eq. MASTER_PE) &
!  write(*,*) "[RadTrans]: Leaving RadHeat.", dr_globalMe
!#endif

  total_timestep = total_timestep + rt_dt
  rt_dt = min(dt-total_timestep, sub_dt) !, rt_dt*1.5)

end do ! subcycle over heating and cooling ! do while (dt-rt_dt .lt. 1e-6)

dt_save = sub_dt

#ifdef DEBUG_RADTRANS
!if (rt_meshMe .eq. MASTER_PE) then
  write(*,*) "[RadTrans]: After subcycle loop.", dr_globalMe
  write(*,'(A,ES18.6E3,I4)') 'rt_dt =', rt_dt, dr_globalMe
  write(*,'(A,ES18.6E3,I4)') 'dt_save =', dt_save, dr_globalMe
!end if
#endif

#ifdef timing
call MPI_allReduce(MPI_IN_PLACE, ion_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

call MPI_allReduce(MPI_IN_PLACE, heat_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)
                     
call MPI_allReduce(MPI_IN_PLACE, raytrace_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)
                     
call MPI_allReduce(nsteps, global_nsteps, 1, MPI_INT, &
                     MPI_MAX, dr_globalcomm, ierr)
                     
if (rt_meshMe .eq. MASTER_PE) then
    write(*, '(A,ES12.3E3, X, A)') "Ionization solver took", ion_time, 'seconds.'
    write(*, '(A,ES12.3E3, X, A)') "Heating and cooling solver took", heat_time, 'seconds.'
    write(*, '(A,ES12.3E3, X, A)') "Ray tracing took", raytrace_time, 'seconds.'
    write(*, '(A,I4, X, A)') "Ray tracing took a max", global_nsteps, 'steps.'
    write(*,'(A,ES12.3E3)') 'Final highest xion =', largest_xion
    write(*,'(A,ES12.3E3)') 'Final highest xion change =', xion_save
end if
#endif

rt_dt = 1e99


#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  print*,'leaving rad. solver'
#endif

  call Grid_notifySolnDataUpdate( (/ EINT_VAR, ENER_VAR, TEMP_VAR, TDUS_VAR, IHP_SPEC, IHA_SPEC /) )

  call Grid_fillGuardCells(CENTER, ALLDIR)
  !call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.true., eosMode=MODE_DENS_EI, selectBlockType=ACTIVE_BLKS)

  call Timers_stop("RadTrans")

  call Logfile_stamp("Leaving RadRay.", "[RadTrans]")

  return
end subroutine RadTrans3
