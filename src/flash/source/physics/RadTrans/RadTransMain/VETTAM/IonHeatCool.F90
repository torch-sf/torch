!!****if* source/physics/RadTrans/RadTransMain/Vettam/IonHeatCool
!!
!!  NAME 
!!
!!  IonHeatCool
!!
!!  SYNOPSIS
!!
!!  call IonHeatCool( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!
!!      Modified from original Torch RadTrans.F90 to only calculate 
!!      heating and cooling and ionization.
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
!#define DEBUG_RADTRANS
#ifdef ONE_CELL_TESTING
#define DEBUG_RADTRANS
#endif
#define timing

subroutine IonHeatCool(nblk, blklst, dt, pass)

#include "Flash.h"

! extra unit specific data
  use rt_data, only : rt_protonMass, abu_c, rt_abar, rt_idealgas, &
                      rt_dt, rt_dt_pos, rt_maxHchange, &
                      rt_rayTrace, rt_heatInRad, rt_neutral_min, &
                      rt_ion_min, rt_ion_threshold

! general radiation transport data
  use RadTrans_data,  ONLY : rt_meshMe, rt_globalComm
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getDeltas

  use Driver_interface, ONLY : Driver_abortFlash
  use Diffuse_interface, ONLY: Diffuse_solveScalar, Diffuse_fluxLimiter

  use RadTrans_data, ONLY: rt_useRadTrans, use_uv_bkgd, uv_bkgd_ion_rate, uv_bkgd_heat_rate
  use Eos_interface, ONLY: Eos_wrapped
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Driver_data, ONLY: dr_simTime, dr_globalcomm, dr_globalMe, dr_nStep

  use Logfile_interface, ONLY: Logfile_stamp ! Lets start recording some info.

! for subcycling internal energy
  use Heat_interface, ONLY : RadHeat
! for background CR ionization
  use Heat_data !, ONLY : he_crIonRate, he_crIonNH, he_crIonExp, he_use_cr_heating, he_subfactor, he_stratifyHeat, he_h_UV
  use heatCool !, ONLY : approx_column_dens, dei_dt
  use cool_vars
  use calc_ion


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
  real, allocatable,dimension(:)          :: dx, dy, dz
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
           NumH0, phoio_old, phoio_new, del_phoio, phoio_save
  real, save :: dt_save=1d99
  real  :: t_start, t_stop, ion_time, heat_time, raytrace_time
  logical :: first_loop, converged, early_exit, all_done, global_all_done, fully_ionized, fileexists
  integer :: nsteps, global_nsteps, step, heatcool_nsteps
  integer :: p, ierr, stat(MPI_STATUS_SIZE), sendproc, dt_request, done_request
  character(len=MAX_STRING_LENGTH*2) :: strbuff
  
  real, parameter :: ion_threshold=1d-10, ion_frac_change_accept=1d-1

  !=========================================================================


  if(.not. rt_useRadTrans) return
! Don't ever trust a function written by someone else to ensure guardcells
! are filled. Always fill on entry to a routine.

  call Grid_fillGuardCells(CENTER, ALLDIR)

!  reset radiation timestep

  rt_dt = min(0.1*rt_maxHchange*dt,dt_save) 
! rt_dt = dt

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


call Logfile_stamp("Entering IonHeatCool.", "[IonHeatCool]")


  call Timers_start("IonHeatCool")

! Subcycle the timestep over the ionization calculation until we hit the
! hydro timestep or we hit convergence everywhere on the grid.

! This loop needs to cycle over the heating and cooling solver on timesteps
! based on the change in internal energy.
do while ((dt-total_timestep) .gt. (1d-6*dt))
  xion_save = 1.0d99

! Lets just run til convergence, or til the cows come home, to solve ionization fraction.  
  do while (xion_save .gt. rt_ion_threshold) !ion_frac_change_accept) ! .or. phoio_save > 1d-3 .or. nsteps < 2)
  
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
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

! not sure if I need the blk pointers or just a solnData array
! temp and density calls for multiple sources here
! additional state variables
          dens = solnData(DENS_VAR,i,j,k)
          temp = solnData(TEMP_VAR,i,j,k)

! neutral and ionised hydrogen fraction
          xH0  = solnData(IHA_SPEC,i,j,k)
          xH1  = solnData(IHP_SPEC,i,j,k)

! Define mean molecular weight from temperature and ionization fraction
! since we don't actually follow the elements of the gas. - JW

      !    if (xHp .gt. 0.5) then ! More than half ionized.
      !       mu_mol = 0.61  ! Ionized
      !    else
            !if (temp .gt. 100.0) then ! Probably not molecular
               mu_mol = 1.3  ! Atomic
        
            !else
          !    mu_mol = 24.0d0/11.0d0   ! Molecular
            !end if
          !end if

! Convert mass density to number density
          !numdens = dens /(rt_abar* rt_protonMass)
          numdens = dens/ (mu_mol * rt_protonMass)
          xh = (/xH0,xH1/)

! phot ionization rate. 
          phih   = solnData(PHIO_VAR,i,j,k)

! Set vars to approx rate of change in energy.
          rho    = dens
          ndens  = numdens
          ei     = solnData(EINT_VAR,i,j,k)
          tdust  = solnData(TDUS_VAR,i,j,k)
          xHp    = xH1
          ephen  = solnData(PHHE_VAR,i,j,k)
          TtoEI   = ei / temp
          conf    = ndens*ndens/rho
          
!          ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
!                    & solnData(VELY_VAR,i,j,k)**2 + &
!                    & solnData(VELZ_VAR,i,j,k)**2)
                    
! convert to erg /(g s)
          ephen  = ephen/rho
! Get the photoelectric flux from radiation sources. -JW
          Gflux = solnData(PEFL_VAR,i,j,k)

! You cannot photoionize more hydrogen than neutral hydrogen exists in the cell!
          !fully_ionized = (phih*rt_dt .ge. NumH0) .and. (xH1 .ge. 1.0d0) .and. (xH0 .le. 0.0d0)
          fully_ionized = ((xH1 .ge. 1.0d0) .or. (xH0 .le. 0.0d0))
! Add background cosmic ray ionization rate to the photoionization rate.
          if (he_use_cr_heating .and. .not. fully_ionized) then
            if (he_crIonNH == 1.0 .and. he_crIonExp == 1.0) then
              phih = phih + he_crIonRate ! Uniform background ionization from CRs
            else ! Use eqn 27 from Padovani et. al. 2009 
              phih = phih + he_crIonRate*(max(he_crIonNH,approx_column_dens(dens, numdens, temp)) / he_crIonNH)**(-he_crIonExp)
            end if
          end if

! Add background UV ionization and heating rate
          if (use_uv_bkgd .and. .not. fully_ionized) then
                  phih = phih + uv_bkgd_ion_rate
                  ephen = ephen + uv_bkgd_heat_rate/(dens*del(IAXIS)*del(JAXIS)*del(KAXIS)) ! convert to erg/(g s)
          endif

! Use this to estimate the radiation transport timestep. Figure 10 % change in ionization is the floor.
          ! (# of neutral H - # of ionizations / s * dt)  / # of neutral H 
          ! is the fractional change in ionization due strictly to radiation transport, 
          ! which is the dt we are looking for.

          !frac_change = abs(NumH0 - phih*rt_dt)/NumH0

! Initialize new ionization fractions and new temperature to initial ones.
! Needed for convergence and flip-flop check below.
! Hopefully no more flip-flopping with implicit solver! - JW

          hFracNew = xh

! abu_c for stabiliy (non zero electron density)
!          eldens = ndens * hFracNew(1) + abu_c
! Calculate the new and mean ionization state and the new electron
! density.
! hFracNew is changed by calc_ionization, xh is the initial value 

! Here we call the new implicit solver for ionization. It will return
! the subcycling timestep if it is smaller than the hydro timestep,
! otherwise it returns the hydro timestep. - JW
!          call calc_ionization(dt, temp, eldens, ndens, hFracNew, xh, phih )
          
          call calc_ionization(rt_dt, ion_dt, temp, ndens, hFracNew, xh, phih)
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
          

! Don't let neutral fraction ever be less than some given value.
! 1. This should be an input runtime parameter and
! 2. With proper CR heating and ionization this would not be needed. - JW

! used to use ion_threshold, now its rt_neutral_min. - JW

          if ((hFracNew(1) .gt. 1.0d0-rt_neutral_min) .or. (hFracNew(0) .lt. rt_neutral_min)) then
            !print*, "[RadTrans]: Warning! Both xHp and xH0 wrong!"
            hFracNew(1) = 1.0d0 -rt_neutral_min !ion_threshold
            hFracNew(0) = rt_neutral_min !ion_threshold
          end if

          if (hFracNew(1) .lt. rt_ion_min) then
            hFracNew(1) = rt_ion_min
            hFracNew(0) = 1.0d0 - rt_ion_min
          end if

          ! NOTE: Having this check for a fractional change in the ionization
          !       fraction is absolutely critical to getting the right
          !       heating and cooling. Failing to control the change in
          !       ionization to be somewhat slow (~10% a step or less) can
          !       lead to large errors in the heating and cooling that are
          !       hard to decipher. - JW
          !frac_change = min(1.0,max(frac_change,abs(hFracNew(1) - xh(1))))
          frac_change = abs(hFracNew(0) - xh(0)) !/max(1d-50,max(hFracNew(0),xh(0)))
          
          !phoio_new = solnData(PHHE_VAR,i,j,k)
          !phoio_old = solnData(OPHH_VAR,i,j,k)
          
          !if (phoio_new .gt. 0.0) sub_dt = min(sub_dt,he_subfactor*phoio_new/max(1d-50,abs(phoio_new - phoio_old)))

! Look I only care about the dangerous instability due to rapidly varying radiation fields
! at very high ionization. Therefore, I'm going to try and ignore timestepping criterion for
! everyone else and I think this is safe (and will be much faster).
        if (ephen .gt. 0d0) then
! Now get the estimated heating and cooling timestep from ei /  d(ei)/dt.
          hc_dt   = he_subfactor * ei / max(abs(dei_dt(0.0,ei)), 1d-50)          
          frac_dt = rt_dt * rt_maxHchange / max(1e-50,frac_change)
        
#ifdef ONE_CELL_TESTING
          write(*,'(A,ES18.6E3,I4)') 'ei / de_dt   =', hc_dt, dr_globalMe
          write(*,'(A,ES18.6E3,I4)') 'ion_dt       =', ion_dt, dr_globalMe
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

#ifdef ONE_CELL_TESTING
          write(*,'(A,ES18.6E3,I4)') 'second sub_dt =', sub_dt, dr_globalMe
#endif

          largest_xion = max(largest_xion, hFracNew(1))
          xion_save = max(xion_save, frac_change)
!          phoio_save = max(phoio_save, del_phoio)
          
! update the ionisation state
! Unless your the first loop, in which just calculate the
! timestep that's proper for the next loop. - JW
          if (.not. first_loop) then
            solnData(IHA_SPEC,i,j,k) = hFracNew(0)
            solnData(IHP_SPEC,i,j,k) = hFracNew(1)
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
!    call Eos_wrapped(MODE_DENS_EI, blkLimits, tmpID)

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

  call MPI_allReduce(MPI_IN_PLACE, xion_save, 1, FLASH_REAL, &
                   MPI_MAX, dr_globalcomm, ierr)
!  call MPI_allReduce(MPI_IN_PLACE, phoio_save, 1, FLASH_REAL, &
!                   MPI_MAX, dr_globalcomm, ierr)

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

  ! On the first loop we saved no
  ! data, so we need to make sure that
  ! xion_save is reset back to a large
  ! number to guarantee this runs for
  ! one more loop.
  if (first_loop) then
    xion_save  = 1d99
    first_loop = .false.
  end if

  end do ! End ionization convergence loop. 

  call MPI_allReduce(MPI_IN_PLACE, sub_dt, 1, FLASH_REAL, &
                   MPI_MIN, dr_globalcomm, ierr)
                   
! Now subcycle heating and cooling. - JW
#ifdef DEBUG_RADTRANS
!if (rt_meshMe .eq. MASTER_PE) &
  write(*,'(A,ES12.3E3,I4)') "[IonHeatCool]: Calling RadHeat with dt =", rt_dt, dr_globalMe
#endif

if (rt_heatInRad) then
#ifdef timing
  t_start = MPI_Wtime()
#endif
  call RadHeat(nblk, blklst, rt_dt, total_timestep)
#ifdef timing
  t_stop = MPI_Wtime()
  heat_time = heat_time + t_stop - t_start
#endif
end if

#ifdef DEBUG_RADTRANS
!if (rt_meshMe .eq. MASTER_PE) &
  write(*,*) "[IonHeatCool]: Leaving RadHeat.", dr_globalMe
#endif

  total_timestep = total_timestep + rt_dt
  rt_dt = min(dt-total_timestep, sub_dt) !, rt_dt*1.5)

end do ! subcycle over heating and cooling ! do while (dt-rt_dt .lt. 1e-6)

dt_save = sub_dt

#ifdef DEBUG_RADTRANS
!if (rt_meshMe .eq. MASTER_PE) then
  write(*,*) "[IonHeatCool]: After subcycle loop.", dr_globalMe
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

  call Timers_stop("IonHeatCool")

  call Logfile_stamp("Leaving IonHeatCool.", "[IonHeatCool]")

  return
end subroutine IonHeatCool
