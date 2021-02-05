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
#define DEBUG_RADTRANS
#define timing

subroutine RadTrans(nblk, blklst, dt, pass)

#include "Flash.h"

! extra unit specific data
  use rt_data, only : rt_protonMass, abu_c, rt_abar, rt_idealgas, rt_dt, rt_dt_pos, rt_maxHchange, &
                      rt_rayTrace

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
  use Driver_data, ONLY: dr_simTime, dr_globalcomm

  use Logfile_interface, ONLY: Logfile_stamp ! Lets start recording some info.

! for subcycling internal energy
  use Heat_interface, ONLY : RadHeat
! for background CR ionization
  use Heat_data, ONLY : he_crIonRate, he_crIonNH, he_crIonExp, he_use_cr_heating
  use heatCool, ONLY : approx_column_dens

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
  real  :: ndens, xH1, xH0 !, tmp_dt <-Changing this name, its a bit confusing
  real  :: phih, dens, hvphih ! since its not an actual timestep. -JW
  real  :: store, eldens, temp, fac, check, sub_dt
  real  :: hFracNew(0:1), xh(0:1)

! for source sink check
  real  :: del(3), x_dis, y_dis, z_dis, frac_change, fac_old ! frac_change is new tmp_dt 
  real  :: send_buff, rec_buff, pos_save(3) !send_buff(2), rec_buff(2) <- repurposed.
  real  :: dt_old, frac_old_save(2), frac_new_save(2)
  
! For new solver. - JW  
  real  :: total_timestep, xion_save, global_dt_save, largest_xion
  real, save :: dt_save=1d99
  real  :: t_start, t_stop, ion_time, heat_time, raytrace_time
  logical :: first_loop, converged, early_exit, all_done, global_all_done
  integer :: nsteps, global_nsteps, step
  integer :: p, ierr, stat(MPI_STATUS_SIZE), sendproc, dt_request, done_request
  character(len=MAX_STRING_LENGTH*2) :: strbuff
  integer :: globalNumParticles

  !=========================================================================

! Don't ever trust a function written by someone else to ensure guardcells
! are filled. Always fill on entry to a routine.

!  call Grid_fillGuardCells(CENTER, ALLDIR)

!  reset radiation timestep

  rt_dt = min(dt,dt_save) 


  total_timestep = 0.0
  pos_save = 0.0
!  fac_old = 1.0
  dt_old = 0.0
!  dt_save = 1e99
  xion_save = 0.0
  first_loop = .true.
  nsteps = 0
!  global_nsteps = 0
!  all_done = .false.
!  global_all_done = .false.

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
  
  do while ((dt-total_timestep) .gt. (1d-6*dt))
  
! Clear previous fluxes.
    do l = 1, nblk
      tmpID = blklst(l)
      call Grid_getBlkPtr(tmpID,solnData)
      solnData(UVFL_VAR,:,:,:) = 0.0
      solnData(FUFL_VAR,:,:,:) = 0.0
      solnData(AUVF_VAR,:,:,:) = 0.0
      solnData(AFUF_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(tmpID,solnData)
    end do
  
! First step we try to go the whole dt.
!    converged = .false. ! If you converge, you get to go home!
!    early_exit = .false. ! If you're bad, you have to start over!
    xion_save = 0.0
    largest_xion = 0.0
    dt_save = 1d99
    nsteps = nsteps + 1

! If we reset here, it allows the timestep to get bigger as the
! solution converges.
    
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
    call Particles_rayAdvance(rt_dt)
    call Timers_stop("raytracing")
    
#ifdef timing
  t_stop = MPI_Wtime()
  raytrace_time = raytrace_time + t_stop - t_start
#endif
    
  endif

#ifdef DEBUG_RADTRANS
  print*,'leaving raytracing'
#endif

  call Timers_start("solving_ionization")

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
        strat_factor = 1d0
      end if
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

! not sure if I need the blk pointers or just a solnData array
! temp and density calls for multiple sources here
! additional state variables
          rho    = solnData(DENS_VAR,i,j,k)
          temp = solnData(TEMP_VAR,i,j,k)
          ei     = solnData(EINT_VAR,i,j,k)
          tdust  = solnData(TDUS_VAR,i,j,k)
          
! neutral and ionised hydrogen fraction
          xH0  = solnData(IHA_SPEC,i,j,k)
          xHp  = solnData(IHP_SPEC,i,j,k)

! Convert mass density to number density
          ndens = dens /(rt_abar* rt_protonMass)
          xh = (/xH0,xHp/)

! phot ionization/ heating rates
          phih   = solnData(PHIO_VAR,i,j,k)
! get heating from feedback processes, i.e. radiation, SN
! this is in [erg/(s cm^3)], as radiation heating only sees atomic Hydrogen
          ephen  = solnData(PHHE_VAR,i,j,k)
          !if (phih .gt. 0.0) print*, "phih=", phih
#ifdef PE_HEAT
! Get the photoelectric flux from radiation sources. -JW
          Gflux = solnData(PEFL_VAR,i,j,k)
#else
          Gflux = 0.0
#endif

! Define mean molecular weight from temperature and ionization fraction
! since we don't actually follow the elements of the gas. - JW

          !if (temp .gt. 100.0) then ! Probably not molecular
          !    if (xHp .gt. 0.5) then ! More than half ionized.
          !       mu_mol = 0.61  ! Ionized
          !    else
                 mu_mol = 24.0d0/11.0d0   ! Atomic
          !    end if
          !else
          !       mu_mol = 2.4   ! Molecular
          !end if


          eiold  = ei
          tmpold = temp

! not sure ek is needed
          ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
                  & solnData(VELY_VAR,i,j,k)**2 + &
                  & solnData(VELZ_VAR,i,j,k)**2)

          
! Add background cosmic ray ionization rate to the photoionization rate.
          if (he_use_cr_heating) then
            if (he_crIonNH == 1.0 .and. he_crIonExp == 1.0) then
              phih = phih + he_crIonRate ! Uniform background ionization from CRs
            else ! Use eqn 27 from Padovani et. al. 2009 
              phih = phih + he_crIonRate*(min(he_crIonNH,approx_column_dens(rho, ndens, temp)) / he_crIonNH)**(-he_crIonExp)
            end if
          end if


#ifdef one_cell_testing
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin energy =", ei
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin temp =", temp
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin density =", rho
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin dust temp =", tdust
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin PEFL =", solnData(PEFL_VAR,i,j,k)
            write(*,'(A,ES13.3E3)') "[RadHeat]: Begin PHHE =", solnData(PHHE_VAR,i,j,k)
            call flush(6)
        !stop
#endif

! Integrate over the change in energy.

#ifdef timing
  nstep   = 0
  t_start = MPI_Wtime()
#endif

               call heating_and_cooling(ei, dt, temp, nstep)

#ifdef timing
                t_stop = MPI_Wtime()
                if (t_stop-t_start > heat_time) then
                !if (nstep > tnstep) then
                    heat_time  = t_stop-t_start
                    ttempstart = tmpold
                    ttempend   = temp
                    tnstep     = nstep
                end if
#endif

! Initialize new ionization fractions and new temperature to initial ones.
! Needed for convergence and flip-flop check below.
! Hopefully no more flip-flopping with implicit solver! - JW

          hFracNew = xh

! Here we call the new implicit solver for ionization. It will return
! the subcycling timestep if it is smaller than the hydro timestep,
! otherwise it returns the hydro timestep. - JW
!          call calc_ionization(dt, temp, eldens, ndens, hFracNew, xh, phih )
          
          call calc_ionization(rt_dt, sub_dt, temp, ndens, hFracNew, xh, phih)
          ! Now set rt_dt for comparison in the while loop. If  current 
          ! rt_dt+sub_dt > dt-rt_dt, we'll just finish the last loop on
          ! dt-rt_dt. Note that if the first subcycle sub_dt > dt, then
          ! rt_dt = dt which ends the while loop. 
         
          if(hFracNew(0) + hFracNew(1) .gt. 1.01) then
            print*,phih,dens,temp,hFracNew,xh
            print*,'ion wrong'
            call flush(6)
            stop
          endif
          
          frac_change = abs(hFracNew(1) - xh(1))
          
! Take steps that are at least 0.1 % in size. Much better than the 10% before,
! which was the cause of the flip-flopping. - JW

! NOTE NOTE NOTE: If this number is zero (for instance if you are in a cell
!                 with no ionization) you're about to divide by zero.
!                 So we add a min number to not do that.
          !if (frac_change .gt. 0.1) then
          !  sub_dt = min(sub_dt,(rt_dt * 0.1 / frac_change))
          !else if (xh(1) .gt. 0.99d0) then
          !else if (frac_change .lt. 0.001) then 
          !  sub_dt = min(rt_dt*2.0,sub_dt * 0.001 / max(frac_change, 1d-50))
            if (frac_change .lt. 0.001) sub_dt = sub_dt * 0.001 / max(frac_change, 1d-50)
          !end if
          ! Set for global comparision on this processor. - JW
          dt_save = min(dt_save, sub_dt, 1.0d1*rt_dt) !, dt)
!#ifdef DEBUG_RADTRANS
!          write(*,'(A,ES12.3E3)') 'dt_save after ionization =', dt_save
!#endif

! Don't let neutral fraction ever be less than some given value.
! 1. This should be an input runtime parameter and
! 2. With proper CR heating and ionization this would not be needed. - JW

          if (hFracNew(1) .gt. 1.0d0 .or. hFracNew(0) .lt. 0.0d0) then
            hFracNew(0)=1.0d-50
            hFracNew(1)=1.0d0 - hFracNew(0)
          end if 


          largest_xion = max(largest_xion, hFracNew(1))
          xion_save = max(xion_save, frac_change)
          

! update the ionisation state
! Unless your the first loop, in which just calculate the
! timestep that's proper for the next loop. - JW

          if (.not. first_loop) then
            solnData(IHA_SPEC,i,j,k) = hFracNew(0)
            solnData(IHP_SPEC,i,j,k) = hFracNew(1)
          else ! Not going to do this in heating, so do it now.
            !clean up heating rate
            solnData(PHHE_VAR,i,j,k) = 0d0
#ifdef PE_HEAT
            ! Clean up fluxes too. - JW
            solnData(PEFL_VAR,i,j,k) = 0d0
#endif
          end if
          
          solnData(PHIO_VAR,i,j,k) = 0d0

!           end if


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
!    call Eos_wrapped(MODE_DENS_EI, blkLimits, tmpID)
  
  enddo block ! block

#ifdef timing
  t_stop = MPI_Wtime()
  ion_time = ion_time + t_stop - t_start
#endif

! Now subcycle heating and cooling. - JW
#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  write(*,'(A,ES12.3E3)') "[RadTrans]: Calling RadHeat with dt =", rt_dt
#endif

#ifdef timing
  t_start = MPI_Wtime()
#endif
if (.not. first_loop) &
  call RadHeat(nblk, blklst, rt_dt, total_timestep)

#ifdef timing
  t_stop = MPI_Wtime()
  heat_time = heat_time + t_stop - t_start
#endif

#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  write(*,*) "[RadTrans]: Leaving RadHeat."
#endif


! If early exit, we don't count the time since we are going to redo it
! at a smaller timestep.

!  if (early_exit) &
!    total_timestep = total_timestep - rt_dt
if (.not. first_loop) &
  total_timestep = total_timestep + rt_dt

  
  ! Unfortunately we can't get around a blocking allreduce here. Its
  ! important that if the ionization changes enough to allow radiation
  ! to travel from one processor to another that it is caught, and therefore
  ! all processors are going to have to take the same sized steps. Honestly
  ! though, this is still much cheaper than calling all of Driver_evolveFlash
  ! on small steps. - JW

#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  write(*,*) "[RadTrans]: Reducing dt_save."
#endif

  call MPI_ALLREDUCE(MPI_IN_PLACE, dt_save, 1, FLASH_REAL, &
                     MPI_MIN, dr_globalcomm, ierr)

#ifdef DEBUG_RADTRANS

  call MPI_allReduce(MPI_IN_PLACE, largest_xion, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

  call MPI_allReduce(MPI_IN_PLACE, xion_save, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

    write(*,'(A,ES12.3E3)') 'rt_dt =', rt_dt
    write(*,'(A,ES12.3E3)') 'dt_save =', dt_save
    write(*,'(A,ES12.3E3)') 'total_timestep =', total_timestep
    write(*,'(A,ES12.3E3)') 'hydro timestep =', dt
    write(*,'(A,ES12.3E3)') 'highest xion =', largest_xion
    write(*,'(A,ES12.3E3)') 'highest xion change =', xion_save
    write(*,'(A,ES12.3E3)') "nsteps =", real(nsteps)
#endif
                     
  rt_dt = min(dt-total_timestep, dt_save)

#ifdef DEBUG_RADTRANS
if (rt_meshMe .eq. MASTER_PE) &
  write(*,*) "[RadTrans]: After dt_save. Reducing xion stuff"
#endif

  call Timers_stop("solving_ionization")

!!!!!!!!!!!!!!!!!!!
!!! Debugging stop
!!!!!!!!!!!!!!!!!!!
  
!    stop

  if (first_loop) then
    first_loop = .false.
    !continue
  end if

end do !subcycle ! do while (dt-rt_dt .lt. 1e-6)

#ifdef timing
call MPI_allReduce(MPI_IN_PLACE, ion_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

call MPI_allReduce(MPI_IN_PLACE, heat_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

call MPI_allReduce(MPI_IN_PLACE, raytrace_time, 1, FLASH_REAL, &
                     MPI_MAX, dr_globalcomm, ierr)

if (rt_meshMe .eq. MASTER_PE) then
    write(*, '(A,ES12.3E3, X, A)') "Ionization solver took", ion_time, 'seconds.'
    write(*, '(A,ES12.3E3, X, A)') "Heating and cooling solver took", heat_time, 'seconds.'
    write(*, '(A,ES12.3E3, X, A)') "Ray tracing took", raytrace_time, 'seconds.'
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
end subroutine RadTrans
