!!  ++ / 
!!  + / -   Heating/Cooling
!!   / --
!!written by C. Baczynski, 2012-2013

!! Modified by J. Wall 05/2016 to no longer be called by
!! Driver_sourceTerms but instead to be subcycled by
!! RadTrans along with ray tracing and ionization frac
!! calculation. - JW

!! Further modified by J. Wall 06/2016 to include:
!!
!! molecular cooling from: Neufield et. al. 1995
!!  H2 dissociation from Lepp & Shull ApJ 270, 1983
!!                  or Shapiro & Kang, ApJ 318, 1987
!! and the dust cooling from Goldsmith ApJ 557
!!
!! that was implemented in Flash by Robi Banerjee
!! and Daniel Seifried, who graciously shared their work with me.

!! Note we want to implement these together since we are subcycling
!! heating and cooling to a solution at the same time.

!! Description:
!!   uses heating rates from SN and ionising radiation to 
!!   calculate new temperature
!!   cooling according to Dalgarno&McCray 1972 is applied
!!
!! Input: 
!!   dt: current simulation timestep
!!   blockcount: number of local blocks
!!   blockList: list of local block IDs
!!   time: global simulation time
!!

!#define debug
!#define DEBUG_OVERHEAT
!#define DEBUG_RADHEAT
!#define DEBUG_DUST
!#define DEBUG_STEP
!#define timing
!#define ONE_CELL_TESTING
!#define verbose

subroutine RadHeat (blockCount,blockList,dt,time)
  use Heat_data ! Use Christian's heating data file.
  use Cool_data ! Use Robi's cooling data file.
  use cool_vars ! Use my variable storage.
  
  
  use Grid_interface,  ONLY : Grid_fillGuardCells, Grid_getCellCoords, &
      Grid_releaseBlkPtr,Grid_getBlkPtr, Grid_getDeltas, &
      Grid_getBlkIndexLimits
  use Grid_data, ONLY : gr_meshComm  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Eos_interface, ONLY : Eos_wrapped
  use Driver_data,   ONLY : dr_simTime, dr_nStep, dr_globalMe, dr_globalcomm
  use Driver_interface, ONLY : Driver_abortFlash

  use heatCool  ! Use my heating and cooling stuff.
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Flash_mpi.h"

! arguments  
  integer,intent(IN)  :: blockCount
  integer,dimension(blockCount),intent(IN):: blockList
  real,intent(IN)  :: dt,time

! block data
  integer  :: blockID, thisBlock
  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:)  :: xCoord, yCoord, zCoord
  real, allocatable,dimension(:)  :: dx, dy, dz
  integer :: xSizeCoord, ySizeCoord, zSizeCoord
  logical :: getGuardCells = .true.
  integer, dimension(2,MDIM)  :: blkLimits, blkLimitsGC
  real, dimension(MDIM)  :: del

! iterators
! TODO see if they can be reused
  integer :: i, j, k, l, m
  integer :: Ncycles
  real    :: dtsub, eldens, neutr
  logical :: fileexists

! variables for the SN explosion
  real    :: xx, yy, zz

! heating and cooling variables
  !real :: tranheat, sdot, sheat,tdepheat
  real :: ei, ek, temp !, timestep
  !real :: dt0, dt1, dt_dei

  real :: radia, tcool, theat, ttherm, scdot, dei, xH0 !, xHp
  real ::  convf, tmpnew, eiold, tmpold
  integer :: nstep

! Code timing
#ifdef timing
  real :: t_start, t_stop, t1, t2
  real :: heat_time, heat_time_max, ttempstart, ttempend
  real, save :: maxtempever = 0d0
  integer :: tnstep
  logical :: got_hotter
#endif

#if defined(timing) || defined(DEBUG_OVERHEAT)
  integer :: ierr
#endif

#ifdef DEBUG_OVERHEAT
  logical :: overheated
#endif

!  logical :: implicit_on = .true.
!  real    :: emin_impl

!interface
!    function dei_dt(del_t, ei_in)
!    implicit none
!    integer, parameter    :: dp=KIND(1d0)
!    real(dp), intent(in) :: del_t, ei_in
!    real(dp)             :: dei_dt
!    end function dei_dt
!end interface

  if (.not. he_useHeat) return

#ifdef timing
  heat_time = 0.0
  tnstep    = 0
  ttempstart = 0.0
  ttempend   = 0.0
  got_hotter = .false.
#endif

#ifdef DEBUG_OVERHEAT
  overheated = .false.
#endif

! Don't ever trust a function written by someone else to ensure guardcells
! are filled. Always fill on entry to a routine.

!  call Grid_fillGuardCells(CENTER, ALLDIR)

!  if ( dt .lt.  he_dtThres) then
!    call Timers_start("heat")
!     if (he_meshMe == MASTER_PE ) then
!       print*,'turning off cooling for this timestep, dt < dtThres'
!     endif
!     do thisBlock = 1, blockCount
!       blockID = blockList(thisBlock)
!! Get a pointer to solution data 
!       call Grid_getBlkPtr(blockID,solnData)
!       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!         do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!             solnData(PHHE_VAR,i,j,k) = 0d0
!           enddo ! coord loops
!         enddo
!       enddo
!      call Grid_releaseBlkPtr(blockID,solnData)
!    enddo ! block loop
!    call Timers_stop("heat")
!    return
!  endif

  ! start the timer ticking
  call Timers_start("heat")


! loop over local blocks in domain, so the blocklist and apply heating and cooling as well as stellar winds
! this loop should go into its own subroutine 
  do thisBlock = 1, blockCount
    blockID = blockList(thisBlock)

! Get a pointer to solution data 
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) !indices for the interior zones and all zones including guard zones
    call Grid_getDeltas(blockID, del) !grid spacing dx dz dy

! Get a pointer to solution data 
    xSizeCoord = blkLimitsGC(HIGH,IAXIS)
    ySizeCoord = blkLimitsGC(HIGH,JAXIS)
    zSizeCoord = blkLimitsGC(HIGH,KAXIS)

! allocate space for dimensions
    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))
    !allocate(dx(xSizeCoord))
    !allocate(dy(ySizeCoord))
    !allocate(dz(zSizeCoord))

! actually only one dx would be needed, -SizeCoord might not be right size, for cubic zones irrelephant
    !dx(:) = del(IAXIS)
    !dy(:) = del(JAXIS)
    !dz(:) = del(KAXIS)

    call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)     

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

!  get zone variables
          temp   = solnData(TEMP_VAR,i,j,k)
          
         
          rho    = solnData(DENS_VAR,i,j,k)
          ei     = solnData(EINT_VAR,i,j,k)
          tdust  = solnData(TDUS_VAR,i,j,k)
#ifdef IHP_SPEC
          xH0    = solnData(IHA_SPEC,i,j,k)
          xHp    = solnData(IHP_SPEC,i,j,k)
#endif
          
! Go away you nasty NaN's. - JW
!          if ( (xHp .le. 0.0) .or. (xH0 .ge. 1.0) ) then
!              xHp = 1e-4
!          end if
          
          eiold  = ei
          tmpold = temp

! not sure ek is needed
          ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
                  & solnData(VELY_VAR,i,j,k)**2 + &
                  & solnData(VELZ_VAR,i,j,k)**2)

! get heating from feedback processes, i.e. radiation, SN
! this is in [erg/(s cm^3)], as radiation heating only sees atomic Hydrogen
          ephen  = solnData(PHHE_VAR,i,j,k)
          
          !if (ephen > 0.0d0) print*, "[RadHeat]: ephen =", ephen
#ifdef DEBUG
          if (ephen < 0.0) then
            print*, "Heating rate less than zero in RadHeat!", ephen
            print*, "blk, i, j, k", blockID, i, j, k, he_meshMe
            call flush(6)
            stop
          end if
#endif
! convert to erg /(g s)
          ephen  = ephen/rho

#ifdef PE_HEAT
! Get the photoelectric flux from radiation sources. -JW
          Gflux = solnData(PEFL_VAR,i,j,k)
#elif defined(VETTAM)
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
          !       mu_mol = 2.3   ! Molecular
          !end if

#ifdef DEBUG_OVERHEAT
! Debugging. -JW
          if (temp .gt. 1e9) then
#ifdef IHP_SPEC
             write(*,'(A, 5ES12.3E3, X, A, 5I4, X, 3ES14.5E3)') &
                "[RadHeat] Too hot at loop begin.  temp, ei, ephen, xH0, xHp:", &
                temp, ei, ephen, xH0, xHp, &
                "  Proc, blk, i, j, k, x, y, z:", &
                he_meshMe, blockID, i, j, k, xx, yy, zz
#else
             write(*,'(A, 3ES12.3E3, X, A, 5I4, X, 3ES14.5E3)') &
                "[RadHeat] Too hot at loop begin.  temp, ei, ephen:", &
                temp, ei, ephen, &
                "  Proc, blk, i, j, k, x, y, z:", &
                he_meshMe, blockID, i, j, k, xx, yy, zz
#endif
             overheated = .true.
          end if
#endif 
#ifdef DEBUG_RADHEAT
               !if (nstep .gt. 100) then
               write(*,'(A,ES12.3E3,X,4I3)') "Heat dt =", dt, blockID, i, j, k
			   print*, "temp, ei, ephen, xH0, xHp:"
               write(*,'(5ES12.3E3)') temp, ei, ephen, xH0, xHp
               call flush()
               !stop
               !end if
#endif

#ifdef ONE_CELL_TESTING
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
                if (temp > maxtempever) then
                !if (t_stop-t_start > heat_time) then
                !if (nstep > tnstep) then
                    heat_time  = t_stop-t_start
                    ttempstart = tmpold
                    ttempend   = temp
                    tnstep     = nstep
                    got_hotter = .true.
                    maxtempever = temp
                end if
#endif

! Debugging. -JW

#ifdef DEBUG_RADHEAT
               !if (nstep .gt. 100) then
               write(*,'(A,ES12.3E3,X,4I3)') "Heat dt =", dt, blockID, i, j, k
			   print*, "temp, ei:"
               write(*,'(2ES12.3E3)') temp, ei !, emin, eplus
               call flush()
               call Driver_abortFlash("[RadHeat]: Solve a single cell for debugging.")
               !end if
#endif
#ifdef DEBUG_OVERHEAT          
           if (temp .gt. 1e9) then
#ifdef IHP_SPEC
             write(*,'(A, 5ES12.3E3, X, A, 5I4, X, 3ES14.5E3, X, A, ES12.3E3)') &
               "[RadHeat] Too hot at loop end.    temp, ei, ephen, xH0, xHp:", &
                temp, ei, ephen, xH0, xHp, &
                "  Proc, blk, i, j, k, x, y, z:", &
                he_meshMe, blockID, i, j, k, xx, yy, zz, &
                "  initial temp:", tmpold
#else
             write(*,'(A, 3ES12.3E3, X, A, 5I4, X, 3ES14.5E3, X, A, ES12.3E3)') &
               "[RadHeat] Too hot at loop end.    temp, ei, ephen:", &
                temp, ei, ephen, &
                "  Proc, blk, i, j, k, x, y, z:", &
                he_meshMe, blockID, i, j, k, xx, yy, zz, &
                "  initial temp:", tmpold
#endif
             overheated = .true.
           end if

#endif
! print warning if temperature too high, maybe fix to upper limit
           if (temp > he_absTmax) then
             inquire(file=trim(he_outfile),exist=fileexists)
             if(fileexists .eqv. .false.) then
               open(he_funit_log, file=trim(he_outfile), position='APPEND')
               write(he_funit_log,'(9(1X,A16))') '[00]time','[01]tempnow' ,'[02]tempcorr', '[03]step',  &
                                     '[04]i', '[05]j', '[06]k', '[07]blockID', '[08]proc'
               close(he_funit_log)
             endif

! write to file some info
             open(he_funit_log, file=trim(he_outfile), position='APPEND')
             write(he_funit_log,'(3(1X,ES16.9),6(1X,I16))') &
             dr_simTime, temp, he_absTmax, dr_nStep, i, j, k, blockID, he_meshMe
             close(he_funit_log)

! fix it and call eos with density temperature mode, afterwards gracefully exit
             ei  = he_absTmax*TtoEI
             temp = he_absTmax
! adjust internal energy

#ifdef verbose
             print*,'\\Correction to temp =', he_absTmax,' in Heat'
             print*,'\\purely from shock heating'
             print*,'\\core ID, block ID, zone ID xyz',he_meshMe,blockID,i,j,k
#endif
           endif

! floor internal energy
           if(ei .lt. he_smallpres/solnData(DENS_VAR,i,j,k)) then
             ei =  he_smallpres/solnData(DENS_VAR,i,j,k)
           endif

! update the global thermodynamic quantities due to the phen. heating
            solnData(ENER_VAR,i,j,k) = ei + ek
            solnData(EINT_VAR,i,j,k) = ei
  					solnData(TEMP_VAR,i,j,k) = temp
            solnData(TDUS_VAR,i,j,k) = tdust

!clean up heating rate
#ifdef FERVENT
            solnData(PHHE_VAR,i,j,k) = 0d0
#endif
#ifdef PE_HEAT
! Clean up fluxes too. - JW
            solnData(PEFL_VAR,i,j,k) = 0d0
#endif

#ifdef DEBUG_STEP
  write(*,'(A,I4)') "[RadHeat]: Number of steps ", nstep
#endif


#ifdef ONE_CELL_TESTING
            write(*,'(A,ES13.3E3)') "[RadHeat]: End energy =", ei
            write(*,'(A,ES13.3E3)') "[RadHeat]: End temp =", temp
            write(*,'(A,ES13.3E3)') "[RadHeat]: End dust temp =", tdust
            write(*,'(A,ES13.3E3)') "[RadHeat]: End PEFL =", solnData(PEFL_VAR,i,j,k)
            write(*,'(A,ES13.3E3)') "[RadHeat]: End PHHE =", solnData(PHHE_VAR,i,j,k)
            call flush(6)
        !stop
#endif

          enddo ! coord loops
        enddo
      enddo

!  clean up memory 
     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
     !deallocate(dx)
     !deallocate(dy)
     !deallocate(dz) 

!  crank changed state variables through EOS
     call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
     !call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)
     
   enddo ! block loop

   call Timers_stop("heat")
   
   !print*, "[RadHeat]: Leaving RadHeat." 

#ifdef DEBUG_OVERHEAT
   call MPI_allReduce(MPI_IN_PLACE, overheated, 1, FLASH_LOGICAL, &
                      MPI_LOR, dr_globalcomm, ierr)
   if (overheated) then
      if (dr_globalMe .eq. MASTER_PE) then
        print *, "Writing checkpoint and plot files, hope it helps find what broke!"
      end if
      call IO_writeCheckpoint()
      call IO_writePlotfile(.true.)  ! mark forced
      ! WEIRD BEHAVIOR: if Driver_abortFlash(...) is wrapped within
      !
      !    if (dr_globalMe .eq. MASTER_PE) then
      !      call Driver_abortFlash("Temp too high at beginning or end of RadHeat loop.")
      !    end if
      !
      ! the code stalls in IO_writeCheckpoint(); therefore make sure
      ! all procs call Driver_abortFlash(...).  ! -AT 2019 June 06
      call Driver_abortFlash("Temp too high at beginning or end of RadHeat loop.")
   endif
#endif

#ifdef timing
  call MPI_Allreduce(heat_time, heat_time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, maxtempever, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, ttempstart, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, ttempend, 1, MPI_DOUBLE_PRECISION, MPI_MAX, gr_meshComm, ierr)

  call MPI_Allreduce(MPI_IN_PLACE, got_hotter, 1, MPI_LOGICAL, MPI_LOR, gr_meshComm, ierr)
  !if (heat_time_max == heat_time) then
  if (got_hotter .and. dr_globalMe == 0) then
  write(*,'(A,ES12.3E3)') "Max time in loop =", heat_time
  write(*,'(A,ES12.3E3)') "Begin loop temp  =", ttempstart
  write(*,'(A,ES12.3E3)') "End loop temp    =", ttempend
  write(*,'(A,ES12.3E3)') "max temp    =", maxtempever
  write(*,'(A,I10)') "End loop nstep   =", tnstep
  call flush(6)
  end if
#endif
   return
end subroutine RadHeat
