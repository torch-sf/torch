!!****f* source/physics/RadTrans/RadTransMain/VETTAM/rt_dustTemperature
!!
!! NAME
!!  
!!  rt_dustTemperature
!!
!!
!! SYNOPSIS
!! 
!!  rt_dustTemperature()
!!  
!! DESCRIPTION
!!
!!  This subroutine updates the dust/gas temperature coupled to the radiation field
!!  
!!
!! AUTHOR
!!  Shyam Harimohan Menon (2020)
!! 
!!
!!***

#undef DEBUG_PICARD

#include "constants.h"
#include "Flash.h"
SUBROUTINE rt_dustTemperature(dt)

  use RadTrans_data
  use Hydro_data, ONLY: hy_useHydro
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_fillGuardCells
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_stampVarMask
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_data, ONLY: dr_globalMe, dr_nStep
#ifdef UEUV_VAR
  use rt_ionisedata, ONLY: useEUVIonize
#endif
  use Eos_interface, ONLY:Eos_wrapped
#ifndef IHA_SPEC
  use Eos_data, ONLY: eos_gamma,eos_singleSpeciesA
#endif


  implicit none

  real, intent(in) :: dt
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: b,blockID,i,j,k,l, iter,error
  integer, parameter :: communication_count = 1
  real :: tdust, tnew, dEr,dTemp,dFrx,dFry,dFrz,rel,del_abs, dele
  real :: max_change_local, max_change_global
  logical :: convergeStatus, dust_addsrcterms
  integer :: picard_substeps, substep_counter
  logical :: gcmask(NUNK_VARS), gcMaskLogged = .false.
  real :: dt_substep, min_dt_substep
  real,dimension(NDIM,blockCount,GRID_ILO_GC:GRID_IHI_GC,&
    GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) :: fl_guess1,fl_guess2,fl_guess3
  real,dimension(blockCount,GRID_ILO_GC:GRID_IHI_GC,&
    GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC) :: local_picard_theta
#ifdef DEBUG_PICARD  
  integer, parameter :: funit = 22
  character(len=80)  :: outfile
  character(len=10)  :: file_ID
  !Tag to differentiate b/w the 2 sweeps in the split solver
  character(len=1), save :: splitSolverTag = 'a'
#endif

  SELECT CASE(rt_hydro_type)
  CASE(1)
    call RadiativeEq() !Radiative equilibrium; i..e Tdust=Trad assumed at all space and time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Picard's fixed-point iteration
  CASE(2)

    picard_substeps = 1
    dt_substep = dt
    !TODO: min subcycling dt hardcoded here for now, make this a user parameter maybe
    min_dt_substep = dt_substep/1000.
    do 
      do substep_counter = 1, picard_substeps
        if(picard_substeps .eq. 1) then 
          call Picard_Scheme(dt_substep,.false.)
        else
          if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,A,I0)") 'Picard Subcycle: ', substep_counter,'/', picard_substeps
          call Picard_Scheme(dt_substep,.true.)
        endif
        !Exit and recycle subcycling at any point scheme doesn't converge
        if(.not. convergeStatus) exit
      end do
      if(convergeStatus) then 
        !Subcycling worked. Move on.
        exit
      else
        !Subcycling didn't work. Try again with smaller timestep.
        picard_substeps = picard_substeps*2
        dt_substep = dt/picard_substeps
        if(dt_substep .lt. min_dt_substep) call Driver_abortFlash("[VET]: Picard's iteration subcycling limit reached.")
      endif

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  CASE DEFAULT 
    call Driver_abortFlash("[VET]: Invalid hydro update type.")
  END SELECT
  !Update eint and pres by calling EOS
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Eos_wrapped(MODE_DENS_TEMP,blkLimits, blockID)
  end do
  
CONTAINS

  SUBROUTINE RadiativeEq()
    do b = 1, blockCount
      blockID = blockList(b)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      call Grid_getBlkPtr(blockID,solnData)
      do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
          do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
            tdust = (solnData(ERAD_VAR,i,j,k)/rt_radconst)**0.25
            if(tdust<rt_T_min) tdust = rt_T_min

            !if KDUS_VAR defined (i.e. Krome Compiled), don't update this temperature
#if defined(TEMP_VAR) && !defined(KDUS_VAR)
            solnData(TEMP_VAR,i,j,k) = tdust !Tdust = Tgas = Trad
#endif
            !TODO: Implement other versions of dust-gas coupling here (e.g. tdust=trad, but tgas and trad interact via collision term)
          end do
        end do
      end do
      call Grid_releaseBlkPtr(blockID,solnData)
    end do
  END SUBROUTINE RadiativeEq
  
  SUBROUTINE Picard_Scheme(dtstep,newSubstep)
    real, intent(in) :: dtstep
    logical, intent(in) :: newSubstep
    iter = 0    
    if(newSubstep) then 
      !Reset solution vectors to previous time solution.
      do b = 1, blockCount
        
        blockID = blockList(b)
        call Grid_getBlkPtr(blockID,solnData)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        !Set old time quantities
#ifdef TEMP_VAR        
        solnData(TEMP_VAR,:,:,:) = solnData(TMPO_VAR,:,:,:)
#endif        
        solnData(ERAD_VAR,:,:,:) = solnData(ERAO_VAR,:,:,:)
        solnData(MOHX_VAR,:,:,:) = solnData(MHXO_VAR,:,:,:)
#if NDIM>1
        solnData(MOHY_VAR,:,:,:) = solnData(MHYO_VAR,:,:,:)
#if NDIM>2
        solnData(MOHZ_VAR,:,:,:) = solnData(MHZO_VAR,:,:,:)
#endif
#endif      
        call Grid_releaseBlkPtr(blockID,solnData)
      end do
      call Timers_start("vet_implicit")
      !Perform first Radiation subsystem step if subcycling is being done (otherwise Cool.F90 does this)
      call Radiation_step(dtstep)
      call Timers_stop("vet_implicit")
    endif

    ! Set the previous iteration value to the old time solution 
    do b = 1, blockCount
      blockID = blockList(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      !Set old time quantities
      solnData(TMPP_VAR,:,:,:) = solnData(TMPO_VAR,:,:,:)
      solnData(ERAP_VAR,:,:,:) = solnData(ERAO_VAR,:,:,:)
      solnData(MHXP_VAR,:,:,:) = solnData(MHXO_VAR,:,:,:)
      fl_guess1(1,b,:,:,:) = solnData(MHXO_VAR,:,:,:)
#if NDIM>1
      solnData(MHYP_VAR,:,:,:) = solnData(MHYO_VAR,:,:,:)
      fl_guess1(2,b,:,:,:) = solnData(MHYO_VAR,:,:,:)
#if NDIM>2
      solnData(MHZP_VAR,:,:,:) = solnData(MHZO_VAR,:,:,:)
      fl_guess1(3,b,:,:,:) = solnData(MHZO_VAR,:,:,:)
#endif
#endif
      local_picard_theta(b,:,:,:) = rt_picard_theta
      call Grid_releaseBlkPtr(blockID,solnData)
    end do

    !Save convergence pattern if debug switched on
#ifdef DEBUG_PICARD
    write(file_ID,'(I0.4)') (dr_nStep)
    outfile = 'PicardIterations_' // trim(adjustl(file_ID)) // trim(adjustl(splitSolverTag)) // '.dat'
    if (dr_globalMe .eq. MASTER_PE) then
      write(*,"(A,A)") 'DEBUG_PICARD:Writing Picard iteration residuals to ',trim(outfile)
      open(funit, file=trim(outfile), status='replace')
      write(funit,'(A14,3X,A8)') '#[00]Iteration', '[01]dMax'
      close(funit)
    endif
    !Split solver tag alternated for 2 files per simulation time step
    if(splitSolverTag .eq. 'a') then 
      splitSolverTag = 'b'
    else 
      splitSolverTag = 'a'
    endif
#endif

    !Start Picard iteration scheme
    iter = 0
    do
      ! Update temperature over all blocks
#ifdef TEMP_VAR
      call Timers_start("vet_nrtempupdate")
      do b = 1, blockCount
        blockID = blockList(b)
        call Temperature_Update(dtstep,blockID)
      end do
      call Timers_stop("vet_nrtempupdate")
#endif
      !Perform check if Picard scheme is converged
      call Timers_start("vet_picardcheck")
      call Picard_CheckConvergence()
      call Timers_stop("vet_picardcheck")

      !Check for convergence
      if(max_change_global<rt_picard_rtol) then
        convergeStatus = .true.
        exit

      else if(iter>rt_picard_maxits) then
        if (dr_globalMe .eq. MASTER_PE) then
          write(*,"(A,F5.2)") 'Picard scheme did not converge, residual change: ',max_change_global
        endif
        convergeStatus = .false.
        return
        !call Driver_abortFlash("[VET]: Picard's iteration has not converged after reaching max iterations")

      endif
      

      iter = iter + 1

      if(rt_picard_correction) then
        do b = 1, blockCount
          blockID = blockList(b)
          call Grid_getBlkPtr(blockID,solnData)
          call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

          do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
            do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                dele = 0.0
#ifdef TEMP_VAR
                dele = solnData(TAUP_VAR,i,j,k)*rt_speedlt*rt_radconst *&
                (solnData(TEMP_VAR,i,j,k)**4 - solnData(TMPP_VAR,i,j,k)**4)
#endif
                solnData(ERAD_VAR,i,j,k) = max(solnData(ERAD_VAR,i,j,k) + dele*dtstep,rt_smalle)

              end do
            end do
          end do
          call Grid_releaseBlkPtr(blockID,solnData)
        end do
      endif

      !Set vector x_k = x_k+1, i.e. as new guess 
      !Also update guesses for fluxes
      do b = 1, blockCount
          
        blockID = blockList(b)
        call Grid_getBlkPtr(blockID,solnData)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
#ifdef TEMP_VAR
        solnData(TMPP_VAR,:,:,:) = solnData(TEMP_VAR,:,:,:)
#endif
        solnData(ERAP_VAR,:,:,:) = solnData(ERAD_VAR,:,:,:)
        solnData(MHXP_VAR,:,:,:) = solnData(MOHX_VAR,:,:,:)
#if NDIM>1
        solnData(MHYP_VAR,:,:,:) = solnData(MOHY_VAR,:,:,:)
#if NDIM>2
        solnData(MHZP_VAR,:,:,:) = solnData(MOHZ_VAR,:,:,:)
#endif
#endif
        call Grid_releaseBlkPtr(blockID,solnData)
      end do

      !Update radiation quantities with new temperature
      call Timers_start("vet_implicit")
      call Radiation_Step(dtstep)
      call Timers_stop("vet_implicit")
    end do
    ! Picard iteration complete
    

    ! Do one last energy synchronisation step by explicitly adding energy locally to all cells
    if(rt_picard_correction) then
      do b = 1, blockCount
        blockID = blockList(b)
        call Grid_getBlkPtr(blockID,solnData)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
          do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
            do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              dele = 0.0
#ifdef TEMP_VAR
              dele = solnData(TAUP_VAR,i,j,k)*rt_speedlt*rt_radconst *&
              (solnData(TEMP_VAR,i,j,k)**4 - solnData(TMPP_VAR,i,j,k)**4)
#endif
              solnData(ERAD_VAR,i,j,k) = max(solnData(ERAD_VAR,i,j,k) + dele*dtstep,rt_smalle)
            end do
          end do
        end do
        call Grid_releaseBlkPtr(blockID,solnData)
      end do
    endif

  END SUBROUTINE Picard_Scheme

  SUBROUTINE Picard_CheckConvergence()
    real :: dErmax, dFrxmax, dFrymax, dFrzmax, dTempmax
    max_change_global = 0.0
    max_change_local = 0.0
    dEr = 0.0
    dFrx = 0.0
    dFry = 0.0
    dFrz = 0.0
    dTemp = 0.0

    do b = 1, blockCount
      
      blockID = blockList(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

      ! Find max change in all 5 quantities in processor
      do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
          do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        
            ! Damp new guess if Picard guesses oscillating about 0 for any of the components
            ! of radiation flux
            if(iter>5) then 
              if(((fl_guess1(1,b,i,j,k)*fl_guess2(1,b,i,j,k) .lt. 0) .and. &
                (fl_guess2(1,b,i,j,k)*fl_guess3(1,b,i,j,k) .lt. 0)) &
#if NDIM>1
                .or. ((fl_guess1(2,b,i,j,k)*fl_guess2(2,b,i,j,k) .lt. 0) .and. &
                (fl_guess2(2,b,i,j,k)*fl_guess3(2,b,i,j,k) .lt. 0)) &
#if NDIM>2
                .or. ((fl_guess1(3,b,i,j,k)*fl_guess2(3,b,i,j,k) .lt. 0) .and. &
                (fl_guess2(3,b,i,j,k)*fl_guess3(3,b,i,j,k) .lt. 0)) & 
#endif
#endif
                ) then

                del_abs = abs((solnData(MOHX_VAR,i,j,k)-solnData(MHXP_VAR,i,j,k)))
                if(fl_guess1(1,b,i,j,k)>rt_picard_rtol .and. del_abs>rt_picard_abstol) then
                  local_picard_theta(b,i,j,k) = local_picard_theta(b,i,j,k)/2.
                endif

              endif  

            endif

            !theta scheme on guess 
            if(iter>1) then 
#ifdef TEMP_VAR
              solnData(TEMP_VAR,i,j,k) = local_picard_theta(b,i,j,k)*solnData(TEMP_VAR,i,j,k) + &
                (1.-local_picard_theta(b,i,j,k))*solnData(TMPP_VAR,i,j,k)
#endif               
              solnData(ERAD_VAR,i,j,k) = local_picard_theta(b,i,j,k)*solnData(ERAD_VAR,i,j,k) + &
                (1.-local_picard_theta(b,i,j,k))*solnData(ERAP_VAR,i,j,k)
              solnData(MOHX_VAR,i,j,k) = local_picard_theta(b,i,j,k)*solnData(MOHX_VAR,i,j,k) + &
                (1.-local_picard_theta(b,i,j,k))*solnData(MHXP_VAR,i,j,k)
#if NDIM>1
              solnData(MOHY_VAR,i,j,k) = local_picard_theta(b,i,j,k)*solnData(MOHY_VAR,i,j,k) + &
                (1.-local_picard_theta(b,i,j,k))*solnData(MHYP_VAR,i,j,k)
#if NDIM>2
              solnData(MOHZ_VAR,i,j,k) = local_picard_theta(b,i,j,k)*solnData(MOHZ_VAR,i,j,k) + &
                (1.-local_picard_theta(b,i,j,k))*solnData(MHZP_VAR,i,j,k)
#endif 
#endif 
            endif        
            
            rel = abs((solnData(ERAD_VAR,i,j,k)-solnData(ERAP_VAR,i,j,k))/solnData(ERAP_VAR,i,j,k))
            del_abs = abs(solnData(ERAD_VAR,i,j,k)-solnData(ERAP_VAR,i,j,k))
            !Set reltol = 0.0 if abstol condition satisfied
            if (del_abs < rt_picard_abstol) rel = 0.0
            dEr = max(rel,dEr)

            
            
#ifdef TEMP_VAR
            rel = abs((solnData(TEMP_VAR,i,j,k)-solnData(TMPP_VAR,i,j,k))/solnData(TMPP_VAR,i,j,k)) 
            del_abs = abs((solnData(TEMP_VAR,i,j,k)-solnData(TMPP_VAR,i,j,k)))
            if (del_abs < rt_picard_abstol) rel = 0.0
            dTemp = max(dTemp,rel)
#endif              
            
           
            rel = abs((solnData(MOHX_VAR,i,j,k)-solnData(MHXP_VAR,i,j,k))/solnData(MHXP_VAR,i,j,k))
            del_abs = abs((solnData(MOHX_VAR,i,j,k)-solnData(MHXP_VAR,i,j,k)))
            if (del_abs < rt_picard_abstol) rel = 0.0
            !If flux low enough allow Picard iterations to converge
            if (abs(solnData(MOHX_VAR,i,j,k)) < rt_smallf) rel = 0.0
            dFrx = max(dFrx,rel)
            fl_guess3(1,b,i,j,k) = fl_guess2(1,b,i,j,k)
            fl_guess2(1,b,i,j,k) = fl_guess1(1,b,i,j,k)
            fl_guess1(1,b,i,j,k) = (solnData(MOHX_VAR,i,j,k)-solnData(MHXP_VAR,i,j,k))/solnData(MHXP_VAR,i,j,k)
            
#if NDIM>1
            rel = abs((solnData(MOHY_VAR,i,j,k)-solnData(MHYP_VAR,i,j,k))/solnData(MHYP_VAR,i,j,k))
            del_abs = abs((solnData(MOHY_VAR,i,j,k)-solnData(MHYP_VAR,i,j,k)))
            if (del_abs < rt_picard_abstol) rel = 0.0
            if (abs(solnData(MOHY_VAR,i,j,k)) < rt_smallf) rel = 0.0
            dFry = max(dFry,rel)
            fl_guess3(2,b,i,j,k) = fl_guess2(2,b,i,j,k)
            fl_guess2(2,b,i,j,k) = fl_guess1(2,b,i,j,k)
            fl_guess1(2,b,i,j,k) = (solnData(MOHY_VAR,i,j,k)-solnData(MHYP_VAR,i,j,k))/solnData(MHYP_VAR,i,j,k)
            
            
#if NDIM>2
            rel = abs((solnData(MOHZ_VAR,i,j,k)-solnData(MHZP_VAR,i,j,k))/solnData(MHZP_VAR,i,j,k))
            del_abs = abs((solnData(MOHZ_VAR,i,j,k)-solnData(MHZP_VAR,i,j,k)))
            if (del_abs < rt_picard_abstol) rel = 0.0
            if (abs(solnData(MOHZ_VAR,i,j,k)) < rt_smallf) rel = 0.0
            dFrz = max(dFrz,rel)
            fl_guess3(3,b,i,j,k) = fl_guess2(3,b,i,j,k)
            fl_guess2(3,b,i,j,k) = fl_guess1(3,b,i,j,k)
            fl_guess1(3,b,i,j,k) = (solnData(MOHZ_VAR,i,j,k)-solnData(MHZP_VAR,i,j,k))/solnData(MHZP_VAR,i,j,k)
                                      
#endif
#endif
  
          end do
        end do
      end do
      call Grid_releaseBlkPtr(blockID,solnData)
    end do
    ! Convergence check finished
    max_change_local = max(dEr,dFrx,dFry,dFrz,dTemp)
    
    ! This will give us max change in the entire domain
    call MPI_AllReduce(max_change_local,max_change_global,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)

!Save convergence pattern if debug switched on
#ifdef DEBUG_PICARD
    !Get global max of each quantity
    call MPI_AllReduce(dEr,dErmax,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call MPI_AllReduce(dFrx,dFrxmax,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call MPI_AllReduce(dFry,dFrymax,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call MPI_AllReduce(dFrz,dFrzmax,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)
    call MPI_AllReduce(dTemp,dTempmax,communication_count,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,error)

    if ((dr_globalMe .eq. MASTER_PE) .and. (iter>0)) then
      open(funit, file=trim(outfile), position='APPEND')
      write(funit,'(I0,1X,ES16.9,1X,ES16.9,1X,ES16.9,1X,ES16.9,1X,ES16.9,1X,ES16.9)') iter+1, &
        dErmax, dFrxmax, dFrymax, dFrzmax, dTempmax, max_change_global
      close(funit)
    endif
#endif
  END SUBROUTINE Picard_CheckConvergence

END SUBROUTINE rt_dustTemperature

! Update temperature for all cells
! Can update temperature through : 1. Newton-Rhapson solve, 2. Linearised Commercon method
! 3. LTE temperature update, 4. Explicit temperature update
#ifdef TEMP_VAR
SUBROUTINE Temperature_Update(dt,blockID)

  use RadTrans_data
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
#ifdef RAYTRACE_3DRT  
  use RadTrans_data, ONLY : rt_useRadTrans
#endif

#ifdef SINK_PART_TYPE  
  use Particles_sinkData, ONLY : useSinkParticles
#endif  

#ifdef IHA_SPEC
  use Eos_interface, ONLY:Eos_wrapped, Eos_getAbarZbar
  use Eos_data, ONLY: eos_gasConstant
#else
  use Eos_interface, ONLY:Eos_wrapped
  use Eos_data, ONLY: eos_gamma,eos_singleSpeciesA, eos_gasConstant
#ifdef TEMP_VAR
  use eos_idealGammaData, ONLY: eos_gammam1
#endif
#endif
  implicit none


#include "Flash.h"
#include "constants.h"


  real, intent(in) :: dt
  integer,intent(in) :: blockID
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(NXB,NYB,NZB):: tnew
  real, dimension(:,:,:,:), pointer :: solnData

  integer :: i,j,k,l, iter,i_b,j_b,k_b
  real :: kB, mH, abar
  tnew(:,:,:) = 0.0

  !Get some constants
  call PhysicalConstants_get("Boltzmann",kB)
  call PhysicalConstants_get("proton mass",mH)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData)
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
    do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
      do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        i_b = i-blkLimits(LOW,IAXIS)+1
        j_b = j-blkLimits(LOW,JAXIS)+1
        k_b = k-blkLimits(LOW,KAXIS)+1
        if(rt_temp_type .eq. 1) then
          call NR_tempupdate(tnew(i_b,j_b,k_b))
        else if(rt_temp_type .eq. 2) then
          call Linearised_tempupdate()
        else if(rt_temp_type .eq. 3) then
          call LTE_tempupdate()
        else if(rt_temp_type .eq. 4) then
          call ExplicitTempUpdate(tnew(i_b,j_b,k_b))
        endif
        if (tnew(i_b,j_b,k_b) < rt_T_min) tnew(i_b,j_b,k_b) = rt_T_min
        ! Safety check
        if(isnan(tnew(i_b,j_b,k_b))) then 
          call Driver_abortFlash("[VET]: Unphysical temperature after NR update")
        endif
#ifdef TEMP_VAR
        solnData(TEMP_VAR,i,j,k) = tnew(i_b,j_b,k_b)
#endif
      end do
    end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData)
  

CONTAINS

  !Subroutine to perform Newton-Rhapson for a single cells
  SUBROUTINE NR_tempupdate(tnew)
    real, intent(out) :: tnew
    real :: told,rho,opac_planck,opac_rosseland,d2,alpha,a,b,c,fx,fdash_x, Er,rel_t,abs_t
    real :: dE, vx, vy, vz, ggprod, toldsave
    integer :: iter

    !Initial temperature before update
#ifdef TEMP_VAR
    told = solnData(TEMP_VAR,i,j,k)
#endif
    rho = solnData(DENS_VAR,i,j,k)
    Er = solnData(ERAD_VAR,i,j,k)
    opac_planck = solnData(TAUP_VAR,i,j,k)
    opac_rosseland = solnData(TAUR_VAR,i,j,k)

#ifdef VELX_VAR 
    vx = solnData(VELX_VAR,i,j,k)
#endif          
#ifdef VELY_VAR
    vy = solnData(VELY_VAR,i,j,k)
#endif
#ifdef VELZ_VAR
    vz = solnData(VELZ_VAR,i,j,k)
#endif


    !Compute coefficients of f(x) = ax**4 + bx + c
    d2 = dt*rt_speedlt

#ifdef IHA_SPEC
    ggprod = (1.0 / (solnData(GAMC_VAR,i,j,k) - 1.0)) * eos_gasConstant
    call Eos_getAbarZbar(solnData(:,i,j,k), abar)
    alpha = (rho * ggprod)/abar
#else
    ggprod = eos_gammam1 * eos_gasConstant
    alpha = (rho * ggprod)/eos_singleSpeciesA
#endif

    ! Velocity dependent terms 
    !Computing explicit terms to add to the energy to gas
    dE = 0.0
    if(rt_ovcterms) then 
    
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vx*solnData(MOHX_VAR,i,j,k)/(rt_speedlt**2) + &
        d2*(opac_planck - opac_rosseland)*vx*(vx*(1+solnData(XXED_VAR,i,j,k)))*&
        solnData(ERAD_VAR,i,j,k)/(rt_speedlt**2)

#if NDIM>1
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vy*solnData(MOHY_VAR,i,j,k)/(rt_speedlt**2) + &
         d2*(opac_planck - opac_rosseland)* solnData(ERAD_VAR,i,j,k)*&
         vy*(vy*(1+solnData(YYED_VAR,i,j,k)) + vx*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2) &
         + d2*(opac_planck - opac_rosseland)*solnData(ERAD_VAR,i,j,k)*vx* &
         (vy*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2)

#if NDIM>2
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vz*solnData(MOHZ_VAR,i,j,k)/(rt_speedlt**2) + &
         d2*(opac_planck - opac_rosseland)* vz* solnData(ERAD_VAR,i,j,k)*(vz*(1+solnData(ZZED_VAR,i,j,k)) + & 
          vx*solnData(XZED_VAR,i,j,k) + vy* solnData(YZED_VAR,i,j,k))/(rt_speedlt**2) + &
          d2*(opac_planck - opac_rosseland)*vx*solnData(ERAD_VAR,i,j,k)* &
          (vz*solnData(XZED_VAR,i,j,k))/(rt_speedlt**2) + d2*(opac_planck - opac_rosseland)* &
          vy*solnData(ERAD_VAR,i,j,k)*(vz*solnData(YZED_VAR,i,j,k))/(rt_speedlt**2)
#endif
#endif
    endif


    a = d2*solnData(TAUP_VAR,i,j,k)*rt_radconst
    b = alpha
    c = -1.*(alpha*solnData(TMPO_VAR,i,j,k) + d2*opac_planck*Er + dE)

!Add contributions from point sources to heating
#ifdef FSPT_VAR
    if(rt_sinkheat .and. rt_sink_implicit) then 
      c = c - solnData(FSPT_VAR,i,j,k)*dt
    endif
#endif

    ! Newton Rhapson iterations start
    iter = 0
    do

      fx = a*told**4 + b*told + c 
      fdash_x = 4*a*told**3 + b

      tnew = told - fx/fdash_x
      rel_t = abs((tnew-told)/told)
      abs_t = abs(tnew-told)

      !TODO: Remove this later
      if (told .eq. 0.0) rel_t = 0.0
      !Check exit conditions
      
      ! Set told to tnew now
      toldsave = told 
      told = tnew
      iter = iter+1

      ! Reltol reached
      if(abs_t .lt. 1.e-10) then 
        exit
      else if(rel_t .lt. rt_nr_rtol) then
        exit 

      else if(rel_t .gt. rt_nr_dtol) then
        call Driver_abortFlash("[VET]: Newton Rhapson temperature update relative change too high.")

      else if(iter .eq. rt_nr_maxits) then
        write(*,"(A,ES9.2,A,ES9.2,ES9.2)") 'Newton scheme did not converge, residual change: ',rel_t,' old,New T',toldsave,tnew
        call Driver_abortFlash("[VET]: Newton Rhapson temperature update does not converge after maxits.")        

      endif

    end do


  END SUBROUTINE NR_tempupdate

  SUBROUTINE ExplicitTempUpdate(tguess_n)
    real, intent(out) :: tguess_n

    real :: told, rho, opac_planck, opac_rosseland, d2, alpha, &
    Er, tguess_o, vx, vy, vz, dE

#ifdef TEMP_VAR
    !Temperature at previous timestep
    told = solnData(TMPO_VAR,i,j,k)
    !Previous guess for temperature
    tguess_o = solnData(TEMP_VAR,i,j,k)
#endif
    rho = solnData(DENS_VAR,i,j,k)
    !Newest guess for Er
    Er = solnData(ERAD_VAR,i,j,k)
    opac_planck = solnData(TAUP_VAR,i,j,k)
    opac_rosseland = solnData(TAUR_VAR,i,j,k)

#ifdef VELX_VAR 
    vx = solnData(VELX_VAR,i,j,k)
#endif
#ifdef VELY_VAR
    vy = solnData(VELY_VAR,i,j,k)
#endif
#ifdef VELZ_VAR
    vz = solnData(VELZ_VAR,i,j,k)
#endif
    d2 = dt*rt_speedlt
#ifdef IHA_SPEC
    call Eos_getAbarZbar(solnData(:,i,j,k), abar)
    alpha = (rho * kB)/((solnData(GAMC_VAR,i,j,k)-1.0)*mH*abar)
#else
    abar = eos_singleSpeciesA
    alpha = (rho * kB)/((eos_gamma-1.0)*mH*abar)
#endif

    tguess_n = (1./alpha)*(alpha*told - opac_planck*d2*(rt_radconst*tguess_o**4 - Er))

    ! Velocity dependent terms 
    !Computing explicit terms to add to the energy to gas
    dE = 0.0 
    if(rt_ovcterms) then
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vx*solnData(MOHX_VAR,i,j,k)/(rt_speedlt**2) + &
        d2*(opac_planck - opac_rosseland)*vx*(vx*(1+solnData(XXED_VAR,i,j,k)))*&
        solnData(ERAD_VAR,i,j,k)/(rt_speedlt**2)

#if NDIM>1
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vy*solnData(MOHY_VAR,i,j,k)/(rt_speedlt**2) + &
          d2*(opac_planck - opac_rosseland)* solnData(ERAD_VAR,i,j,k)*&
          vy*(vy*(1+solnData(YYED_VAR,i,j,k)) + vx*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2) &
          + d2*(opac_planck - opac_rosseland)*solnData(ERAD_VAR,i,j,k)*vx* &
          (vy*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2)

#if NDIM>2
      dE = dE - d2*(2*opac_planck - opac_rosseland)*vz*solnData(MOHZ_VAR,i,j,k)/(rt_speedlt**2) + &
          d2*(opac_planck - opac_rosseland)* vz* solnData(ERAD_VAR,i,j,k)*(vz*(1+solnData(ZZED_VAR,i,j,k)) + & 
          vx*solnData(XZED_VAR,i,j,k) + vy* solnData(YZED_VAR,i,j,k))/(rt_speedlt**2) + &
          d2*(opac_planck - opac_rosseland)*vx*solnData(ERAD_VAR,i,j,k)* &
          (vz*solnData(XZED_VAR,i,j,k))/(rt_speedlt**2) + d2*(opac_planck - opac_rosseland)* &
          vy*solnData(ERAD_VAR,i,j,k)*(vz*solnData(YZED_VAR,i,j,k))/(rt_speedlt**2)
#endif
#endif
    endif

    tguess_n = tguess_n + (1./alpha)*(dE)

    !Add contributions from point sources to heating
    !TODO: This may actually not conserve energy, unless I add this term in the RHD eqns as well
#ifdef FSPT_VAR
    if(rt_useRadTrans .and. useSinkParticles .and. rt_sink_implicit) then 
      tguess_n = tguess_n + solnData(FSPT_VAR,i,j,k)*dt
    endif
#endif

  END SUBROUTINE ExplicitTempUpdate

  SUBROUTINE Linearised_tempupdate()

    real :: told, mu,rho,opac_planck,d2,alpha,Er

    !Initial temperature before update
    told = solnData(TMPO_VAR,i,j,k)
    rho = solnData(DENS_VAR,i,j,k)
    Er = solnData(ERAD_VAR,i,j,k)
    opac_planck = solnData(TAUP_VAR,i,j,k)


    d2 = dt*rt_speedlt
#ifdef IHA_SPEC
    call Eos_getAbarZbar(solnData(:,i,j,k), abar)
    alpha = (solnData(GAMC_VAR,i,j,k)-1.0)*mH*abar/(rho * kB)
#else
    abar = eos_singleSpeciesA
    alpha = (eos_gamma-1.0)*mH*abar/(rho * kB)
#endif

    mu = d2*opac_planck*alpha

    ! Commerçon et al 2011 linearisation

    tnew(i_b,j_b,k_b) = (told + 3*mu*rt_radconst*(told**4) + mu*Er) / (1 + 4*mu*rt_radconst*(told**3))

  END SUBROUTINE Linearised_tempupdate

  SUBROUTINE LTE_tempupdate()

    tnew(i_b,j_b,k_b) = (solnData(ERAD_VAR,i,j,k)/rt_radconst)**0.25

  END SUBROUTINE LTE_tempupdate


END SUBROUTINE Temperature_Update
#endif
