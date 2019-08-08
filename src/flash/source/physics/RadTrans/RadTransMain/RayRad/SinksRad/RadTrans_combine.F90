!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans
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
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step using multigroup diffusion theory.
!!
!!      Modified to properly cycle over the ray tracing and ionization
!!      calculation here, instead of running the entirety of Flash on
!!      the very short fractional ionization change timesteps. Also,
!!      the ionization solver was made implicit so as to take larger
!!      timesteps overall. Finally, coupling to sinks was added. - Josh Wall
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
subroutine RadTrans(nblk, blklst, dt, pass)

#include "Flash.h"

! extra unit specific data
  use rt_data, only : rt_protonMass, abu_c, rt_abar, rt_idealgas, rt_dt, rt_dt_pos, rt_maxHchange, &
                      rt_rayTrace

! general radiation transport data
  use RadTrans_data,  ONLY : rt_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getDeltas

  use Driver_interface, ONLY : Driver_abortFlash
  use Diffuse_interface, ONLY: Diffuse_solveScalar, Diffuse_fluxLimiter

  use RadTrans_data, ONLY: rt_useRadTrans
  use Eos_interface, ONLY: Eos_wrapped
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_rayAdvance
  use Driver_data, ONLY: dr_simTime, dr_globalcomm

  use Logfile_interface, ONLY: Logfile_stamp ! Lets start recording some info.

! For sink particles - JW
#ifdef SINK_PART_TYPE
  use Particles_sinkData, only : localnpf, localnp, particles_local
  use Grid_interface
#endif

! for subcycling internal energy
  !use Heat_interface, ONLY : Heat, radloss
  use Heat_data

  implicit none

#include "Flash_mpi.h"

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
  real  :: send_buff(2), rec_buff(2), pos_save(3)
  real  :: dt_old, frac_old_save(2), frac_new_save(2)
  real  :: total_timestep, xion_save
  real  :: TtoEI, eplus, emin, ei, dei, dt_dei, dt1, conf
  real, save :: dt_save=1e99
  logical :: first_loop, converged, early_exit
  integer :: p, ierr, stat(MPI_STATUS_SIZE), sendproc
  character(len=MAX_STRING_LENGTH*2) :: strbuff

  !=========================================================================

!  reset radiation timestep
  rt_dt = min(dt,dt_save)
  total_timestep = 0.0
  pos_save = 0.0
!  fac_old = 1.0
  dt_old = 0.0
  dt_save = 1e99
  xion_save = 0.0
  first_loop = .true.
  
  

  if(.not. rt_useRadTrans) return

call Logfile_stamp("Entering RadRay.", "[RadTrans]")

#ifdef SINK_PART_TYPE
! If there are no sink particles or massive particles, return
! immediately. - JW
! Note this assumes if you wanted sinks, that they are the only
! sources of ionization. If that's not true, you'll want to take
! this out. - JW
#ifndef ACTIVE_PART_TYPE
  if (localnpf == 0) then
    !print*, "No particles so no rays generated this step."
    !print*, "Also don't mess with ionization fraction either, its too damn slow."
    return
  end if
#endif
#endif

! if single source solving stop here

#ifdef DEBUG_RADTRANS
  print*,'entering rad. solver'
#endif

  call Timers_start("RadTrans")

! Subcycle the timestep over the ray tracing method
! and the ionization calculation until we hit the
! hydro timestep or we hit convergence everywhere on the
! grid.

!  dt_save = 1e99

  do while ((dt-total_timestep) .gt. 1e-6)
  
! First step we try to go the whole dt.
    converged = .false. !.true. ! If you converge, you get to go home!
    early_exit = .false. ! If you're bad, you have to start over!

    xion_save = 0.0
    total_timestep = total_timestep + rt_dt

! Try to run the whole hydro loop first, since
! we need x_ion_new to calculate the subcycle dt anyways.
!    if (first_loop) then
!      rt_dt = dt 
!    end if
    
#ifdef DEBUG_RADTRANS
  print*,'entering raytracing'
#endif

! raytracing, should go in Particles_advance.F90 but then
! order in Driver is screwed up, as source term is calculated after call to hydro solver
! as long as the raytracing manipulates the data structure orderly then nothing should happen

! Actually we'd like to do this here, lets keep this together with
! solving dx_ion. - JW

! Actual actual radiation transport. -JW
  if (rt_rayTrace) then
#ifdef DEBUG_RADTRANS    
    if (rt_Meshme == MASTER_PE) print*, "Calling ray tracing with dt=.", rt_dt
#endif
    call Timers_start("raytracing")
    call Particles_rayAdvance(rt_dt)
    call Timers_stop("raytracing")
  endif

#ifdef DEBUG_RADTRANS
  print*,'leaving raytracing'
#endif

!  call Heat(nblk, blklst, rt_dt, dt)

  call Timers_start("solving_ionization") 
!===========================================
! actual radiation transport (ionization calc... -JW)
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
          ei   = solnData(EINT_VAR,i,j,k)
! neutral and ionised hydrogen fraction
          xH0  = solnData(IHA_SPEC,i,j,k)
          xH1  = solnData(IHP_SPEC,i,j,k)

! phot ionization/ heating rates
          phih   = solnData(PHIO_VAR,i,j,k)

! Convert mass density to number density
          ndens = dens /(rt_abar* rt_protonMass)
          xh = (/xH0,xH1/)

! ei to T conversion, if eos was called both should be proportional
          TtoEI  = ei/temp
          
! conversion factor for cooling rate
! instead of n^2 it should be n_ion*n_elec 
! but we do not calculate the electron fraction explicitly
! also implies this is all recombination cooling which I am not sure of
          conf   = ndens*ndens/dens
          
          if (temp .lt. he_absTmin) then
! fix internal energy to 10 K equivalent
            ei  = he_absTmin*TtoEI
            temp = he_absTmin
          endif

! get heating from feedback processes, i.e. radiation, SN
! this is in [erg/(s cm^3)], as radiation heating only sees atomic Hydrogen
          eplus  = solnData(PHHE_VAR,i,j,k)

! convert to erg /(g s)
          eplus  = eplus/dens

! add UV heating
           if ( (temp < he_theatmax) .and. (temp > he_theatmin) ) then
! add any additional heating terms
! photoelectric heating
! stratify ?
             if(he_stratifyHeat) then
! converts to [erg/(g s)] 
! he_peheat is in [erg/s]
               eplus = eplus +  he_peheat * exp(-abs(zz)/(he_h_UV)) / (dens/ndens)
             else
               eplus = eplus +  he_peheat / (dens/ndens)
             endif
           endif

           if (  (temp   <= he_tradmax .AND. temp   >= he_tradmin) & 
             .AND. (ndens >= he_dradmin .AND. ndens <= he_dradmax) ) then

             ! get da cooling rate [erg cm^3/s] (volumetric)
             if(he_coolOff) then
               emin = 0d0
             else
             ! emin is in erg/s(?), also output for given temperature
               call radloss(temp,emin)
             endif

             ! convert to erg/(g s)
             emin = emin*conf
! calculate energy change rate, at 1e4 K rapid change might be unstable if temperature is just above
! with a much shallower slope
             dei    = max(1e-50,abs(eplus-emin))

             ! get timestep, large change -> small time step
             dt_dei = ei/abs(dei)

            ! subcycle timestep size 
            dt1    = he_subfactor*dt_dei

             ! change internal energy by a fraction
             ei     = ei + rt_dt*(eplus-emin)
             ! change temperature 
             temp = ei/TtoEI


             ! over cooled, fix	
             if (temp .lt. he_absTmin) then
             ! fix internal energy to 10 K equivalent
               ei  = he_absTmin*TtoEI
               temp = he_absTmin
             endif
           end if
! Initialize new ionization fractions and new temperature to initial ones.
! Needed for convergence and flip-flop check below.
! Hopefully no more flip-flopping with implicit solver! - JW
!    			store = xh(0)

          hFracNew = 0.0 !xh

! abu_c for stabiliy (non zero electron density)
!          eldens = ndens * hFracNew(1) + abu_c
! Calculate the new and mean ionization state and the new electron
! density.
! hFracNew is changed by calc_ionization, xh is the initial value 

! Here we call the new implicit solver for ionization. It will return
! the subcycling timestep if it is smaller than the hydro timestep,
! otherwise it returns the hydro timestep. - JW
!          call calc_ionization(dt, temp, eldens, ndens, hFracNew, xh, phih )          
          call calc_ionization(rt_dt, sub_dt, temp, ndens, hFracNew, xh, phih)
          
          ! Now set rt_dt for comparison in the while loop. If  current 
          ! rt_dt+sub_dt > dt-rt_dt, we'll just finish the last loop on
          ! dt-rt_dt. Note that if the first subcycle sub_dt > dt, then
          ! rt_dt = dt which ends the while loop. Note dt_old allows
          ! for global comparison of cell timesteps. - JW


          !rt_dt  = min(rt_dt+sub_dt, dt-rt_dt, dt_old)
          ! Set for global comparision on this processor. - JW
          dt_save = min(sub_dt, dt1) !dt_save, 

          
          if(hFracNew(0) + hFracNew(1) .gt. 1.01) then
            print*,phih,dens,temp,hFracNew,xh
            print*,'ion wrong'
            stop
          endif
          
! If the step was too large (because dt_save < rt_dt or because we
! overshoot the fraction and xion > 1.0, exit early.

!          if (rt_dt/dt_save .gt. 10) then !.or. (hFracNew(1) .gt. 1.0)) &
!!               .and. (frac_change .gt. 1e-4)) then
!!!#ifdef DEBUG_RADTRANS
!!            print*, "Exiting early."
!!            print*, "Xion =", hFracNew(1)
!!            print*, "dt_save, rt_dt =", dt_save, rt_dt
!!!#endif
!             early_exit = .true.
!!            if (hFracNew(1) .gt. 1.0) dt_save = rt_dt*(xh(1)/hFracNew(1))
!! Clean up on our way out.
            
!              solnData(PHIO_VAR,i,j,k) = 0d0
            
!              call Grid_releaseBlkPtr(tmpID,solnData)
!              deallocate(xCoord)
!              deallocate(yCoord)
!              deallocate(zCoord)            
            
!              exit block
!            else
!              early_exit = .false.
!          end if
          
            if (hFracNew(1) .gt. 1.0) then
              hFracNew(1)=1.0
              hFracNew(0)=0.0
            end if
          xion_save = max(xion_save, hFracNew(1))
          ! Is this cell converged to a solution?
!          frac_change = abs(hFracNew(1)-xh(1))/xh(1)
!          if (frac_change .gt. 1e-4) converged = .false.

! update the ionisation state
          solnData(IHA_SPEC,i,j,k) = hFracNew(0)
          solnData(IHP_SPEC,i,j,k) = hFracNew(1)
          solnData(PHIO_VAR,i,j,k) = 0d0
          solnData(PHHE_VAR,i,j,k) = 0d0
  
!           end if
! reset/treated in phen heating
! 					solnData(PHHE_VAR,i,j,k) = 0d0
			
! time step criterion change of ionized hydrogen fraction in timestep
! this is after changing the ionisation fraction, corrective method
! if change in neutral fraction in current timestep was greater than 0.1 than 
! slow down timestep

! Unless the cell is occupied by a source. -JW


!          if(rt_rayTrace ) then
!            frac_change = abs(xH0 - hFracNew(0))

!            if(frac_change .gt. rt_maxHchange) then

!              fac = rt_maxHchange/frac_change ! Factor is max fractional change over
                                              ! actual fractional change. - JW
!              if (dt .lt. dt_old) then ! Are we smaller than previous
                                         ! factors? - JW
! Is there a source in this cell? If so no limiter applied. - JW
! Note here it'd be easy to add another check for active/massive sources. -JW
!#ifdef SINK_PART_TYPE
!                x_dis = 1e99
!                y_dis = 1e99
!                z_dis = 1e99
!                ! Find the closest particle on this processor.
!                ! Tacit assumption that no particle is too close on neighboring
!                ! proc here, which is okay since we are upping the timestep
!                ! if a particle is too close, not lowering it. -JW
!                do p=1, localnp
!
!                  x_dis = min(x_dis,abs(xx-particles_local(POSX_PART_PROP, p)))
!                  y_dis = min(y_dis,abs(yy-particles_local(POSY_PART_PROP, p)))
!                  z_dis = min(z_dis,abs(zz-particles_local(POSZ_PART_PROP, p)))
!
!                  
!                end do
!
!                if ((x_dis .le. 0.5*del(1)) .and. (y_dis .le. 0.5*del(2)) &
!                  .and. (z_dis .le. 0.5*del(3))) then
!
!                  if (frac_change .gt. 0.5) then
!                    fac = 0.5 / frac_change
!                  end if ! Is frac_change > 0.5?
!                  ! Is the timestep factor larger than any previous? - JW
!                  if (fac .ge. fac_old) then
!                    cycle ! Then go to the next cell. - JW
!                  end if ! Is fac < fac_old?  
!#endif
! change in current timstep, next timestep should be smaller 
!                      rt_dt = dt !fac*
!                      rt_dt_pos(1) = 1
!                      rt_dt_pos(2) = 1
!                      rt_dt_pos(3) = 1
!                      rt_dt_pos(4) = tmpID
!                      rt_dt_pos(5) = rt_meshMe
!                  ! Update the old fac.
!                  !    fac_old = fac
!                      pos_save(1) = xx
!                      pos_save(2) = yy
!                      pos_save(3) = zz
!                      frac_old_save = xh
!                      frac_new_save = hFracNew
!                      dt_old = dt

!#ifdef SINK_PART_TYPE
!                endif ! Source check. - JW
!#endif
!              endif ! If the frac is smaller than previous fracs. -JW
!            endif ! If fractional change is > 0.1.
!          endif ! If ray tracing is on.
        enddo ! coord loops
      enddo
    enddo

!  clean up memory 
    call Grid_releaseBlkPtr(tmpID,solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)
  enddo block ! block

!===========================================
! call to EOS to set zone
!===========================================
! no call to eos hydrogen does not partake in hydrodynamics it is too good for that.
    call Eos_wrapped(MODE_DENS_EI, blkLimits, tmpID)



! If early exit, we don't count the time since we are going to redo it
! at a smaller timestep.

  if (early_exit) &
    total_timestep = total_timestep - rt_dt

  rt_dt = min(dt-total_timestep, dt_save)
  !rt_dt = max(rt_dt, 1e-8)
  
  ! Moved up to top of the loop.
  !total_timestep = total_timestep + rt_dt

! If the first step was too large, don't count it
! since we are going to start over at smaller steps.
!  if (first_loop .and. (early_exit)) then
!    total_timestep = 0.0
!    first_loop = .false.
!  end if
!#ifdef DEBUG_RADTRANS
  if (rt_meshMe == MASTER_PE) then
    write(*,'(A,ES12.3E3)') 'rt_dt =', rt_dt
    write(*,'(A,ES12.3E3)') 'sub_dt =', dt_save
    write(*,'(A,ES12.3E3)') 'total_timestep =', total_timestep
    write(*,'(A,ES12.3E3)') 'highest xion =', xion_save
  end if
!#endif
  call Timers_stop("solving_ionization")
  
  if (converged) then
!#ifdef DEBUG_RADTRANS
    print*, "We converged! Exiting the loop."
!#endif
    exit
  end if
  
end do !subcycle ! do while (dt-rt_dt .lt. 1e-6)

!call Grid_fillGuardCells(CENTER, ALLDIR)
! Set rt_dt for comparison with other timesteps (hydro, etc).
rt_dt = 1e99 !dt_save

#ifdef DEBUG_RADTRANS
if (rt_meshMe == MASTER_PE) print*, "Leaving RadTrans!"
#endif

! Lets see where the minimum is. -JW

!send_buff = [rt_dt, real(rt_meshMe)]

!call MPI_ALLREDUCE(send_buff, rec_buff, 2, MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0, dr_globalcomm , ierr)

!sendproc = int(rec_buff(2))

!if (rt_meshMe == sendproc) then
!!   print*, "Proc ", sendproc, rt_meshMe, " sending x, y, z ", pos_save
!!   print*, "Proc ", sendproc, " sending neu, ion, diff", frac_old_save, frac_new_save, abs(frac_old_save - frac_new_save)
!   call MPI_Send(pos_save, 3, MPI_DOUBLE_PRECISION, 0, 0, dr_globalcomm, ierr)
!   call MPI_Send(frac_old_save, 2, MPI_DOUBLE_PRECISION, 0, 1, dr_globalcomm, ierr)
!   call MPI_Send(frac_new_save, 2, MPI_DOUBLE_PRECISION, 0, 2, dr_globalcomm, ierr)
!end if

!if (rt_meshMe == 0) then
!!   print*, "Proc ", rt_meshMe, " recieving."
!   call MPI_Recv(pos_save, 3, MPI_DOUBLE_PRECISION, sendproc, 0, dr_globalcomm, stat, ierr)
!   call MPI_Recv(frac_old_save, 2, MPI_DOUBLE_PRECISION, sendproc, 1, dr_globalcomm, stat, ierr)
!   call MPI_Recv(frac_new_save, 2, MPI_DOUBLE_PRECISION, sendproc, 2, dr_globalcomm, stat, ierr)
   
!   print*, "In rad trans, smallest timestep is ", rec_buff(1), " on proc ", &
!           rec_buff(2), " at ", pos_save
!   print*, "Old and new neu and ion fracs, diff ", frac_old_save, frac_new_save, abs(frac_old_save-frac_new_save)

!   write(strbuff,'(a,e10.3)') "Smallest timestep = ", rec_buff(1)
!   call Logfile_stamp(trim(strbuff), "[RadTrans]")
!   write(strbuff,'(a,3e10.3)') "At location x, y, z = ", pos_save
!   call Logfile_stamp(trim(strbuff), "[RadTrans]")

!end if


#ifdef DEBUG_RADTRANS
	print*,'leaving rad. solver'
#endif


  call Timers_stop("RadTrans")

  call Logfile_stamp("Leaving RadRay.", "[RadTrans]")

  return
end subroutine RadTrans
