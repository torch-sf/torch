!!****f* source/physics/RadTrans/RadTransMain/VETTAM/rt_petsc
!!
!! NAME
!!  
!!  rt_petsc
!!
!!
!! SYNOPSIS
!! 
!!  rt_petsc()
!!  
!! DESCRIPTION
!!
!!  This contains the subroutines for the allocation, deallocation, initialising and solving of the matrices/vectors
!!  of the linearised system of equations of RHD solved with VETTAM
!!***

#include "Flash.h"
#include "constants.h"
#include "petsc/finclude/petscksp.h"
#undef VET_DEBUG
! Initial PETSc data structures. Called by RadTrans_init or when the AMR Structure changes due to refinement. 
SUBROUTINE Petsc_init(blockCount_,blockList_)
  use RadTrans_data
  use tree, only : grid_changed
  use Timers_interface, ONLY: Timers_start,Timers_stop
  implicit none
  integer, INTENT(IN) :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  type(PetscErrorCode) :: ierr

  !blockCount and blockList are the global block info in the VET unit. 
  !blockCount_ and blockList_ are the temporary ones passed. 
  !TODO: grid_changed does not get updated for 2 steps. Maybe I can limit the grid_changed
  !      condition to only act if nSteps>2
  IF(grid_changed.NE.0.or.(.not.alloced)) THEN
    blockCount = blockCount_
    blockList(1:blockCount) = blockList_(1:blockCount)
    ! Make sure surr_blks is updated.
    call Timers_start("vet_neghinit")
    call gr_ensureValidNeighborInfo(10)
    call InitNeghLevels()
    call Timers_stop("vet_neghinit")
 
    call Petsc_alloc()
  END IF
END SUBROUTINE Petsc_init

! Allocate memory for PETSc vectors and matrices.
SUBROUTINE Petsc_alloc()
  use RadTrans_data
  use Driver_data, ONLY:dr_globalMe
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
  integer :: ierr, band_count
  type(tPC) :: pc
  character(len=MAX_STRING_LENGTH) :: current_band_temp
  !Preallocating memory based on maximum possible neighbours in diagonal/off-diagonal
  !Max neghs = NDIM+NDIM*2**(NDIM-1)) considering cell at 3 simultaneous block boundaries
  !Own cell stencil columns has NDIM+1 non-zero stencil elements, neghs have 2 each 
  !TODO: This does max possible mem allocation. Can optimise if required by reading the grid info.
  type(PetscInt), parameter :: diagonal = (NDIM+1)+2*(NDIM+NDIM*2**(NDIM-1)),&
                               offdiagonal= 2*(NDIM*2**(NDIM-1))
  EXTERNAL :: MyConvergenceTest

  IF(alloced) THEN
    call Petsc_dealloc()
    if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Petsc reallocating memory and setup as grid changed.", "[VET] ")
  ELSE 
    if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Petsc allocating memory and setup on first call.", "[VET] ")
  END IF 

  ! Create the KSP instance and set solver/PC info
  call KSPCreate(comm,ksp,ierr) 
  call KSPSetType(ksp, rt_linearsolver_str, ierr)
  call KSPGetPC(ksp,pc,ierr)
  call PCSetType(pc,rt_preconditioner_str,ierr)

  ! Set user defined convergence test. Not sure if we need this. 
  call KSPConvergedDefaultCreate(ctx,ierr)
  call KSPSetConvergenceTest(ksp, MyConvergenceTest, ctx, PETSC_NULL_FUNCTION, ierr)

  ! Create the solution vector, RHS and the coefficient matrix
  ! Number of cells in solution vector/no of equations = number of elements (for Erad and Frads)
  ncells = NXB*NYB*NZB*blockCount*(NDIM+1)

  !Create a vector for each band
  do band_count=1, rt_freqbands
    SELECT CASE(band_count)
    !Note: the initial guess is set to non-zero for the IR band only, and is set to the previous time soln.
    !For the UV bands, the previous soln is not available; so the guess is set to zero. Subsequent steps should have previous solutions in the vector as guesses.
    !This might lead to a higher number of iterations for the UV bands when allocation happens (start of sims or refinement), but it is transitory
    CASE(1)
        current_band_temp = rt_band1
        call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Xvec_Band1, ierr)
        if(current_band_temp .eq. 'IR') then
          call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
          call CopyGuess(Xvec_Band1)
        endif
      CASE(2)
        current_band_temp = rt_band2
        call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Xvec_Band2, ierr)
        if(current_band_temp .eq. 'IR') then
          call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
          call CopyGuess(Xvec_Band2)
        endif
      CASE(3)
        current_band_temp = rt_band3
        call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Xvec_Band3, ierr)
        if(current_band_temp .eq. 'IR') then
          call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
          call CopyGuess(Xvec_Band3)
        endif
      CASE(4)
        current_band_temp = rt_band4
        call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Xvec_Band4, ierr)
        if(current_band_temp .eq. 'IR') then
          call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
          call CopyGuess(Xvec_Band4)
        endif
      CASE(5)
        current_band_temp = rt_band5
        call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Xvec_Band5, ierr)
        if(current_band_temp .eq. 'IR') then
          call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
          call CopyGuess(Xvec_Band5)
        endif
      CASE DEFAULT
        if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,I0)") 'Band Count, No of bands: ', band_count, rt_freqbands
        call Driver_abortFlash("[VETTAM]: Something is wrong; band_count is greater than 5, should not be happening!")
    END SELECT
  end do
  !Create RHS Vector
  call VecCreateMPI(comm, ncells, PETSC_DETERMINE, Yvec, ierr)

  !Allocate matrix
  ! Let Petsc take care of all initialisation. Not fast.
  call MatCreateAIJ(comm,ncells,ncells,PETSC_DETERMINE,PETSC_DETERMINE,&
    diagonal,PETSC_NULL_INTEGER,offdiagonal,PETSC_NULL_INTEGER,Amat,ierr)
  
  
  call KSPSetTolerances(ksp, rt_rtol, rt_abstol, rt_dtol, rt_maxits, ierr)
  call KSPSetFromOptions(ksp,ierr)
  
#ifdef VET_DEBUG
  call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
#endif

  alloced = .TRUE.

END SUBROUTINE Petsc_alloc



SUBROUTINE Petsc_dealloc()
  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe
  implicit none
  integer :: ierr, band_count
  character(len=MAX_STRING_LENGTH) :: current_band_temp

  call MatDestroy(Amat,ierr)

  !Destroy all solution vectors
  do band_count=1, rt_freqbands
    SELECT CASE(band_count)
    CASE(1)
        current_band_temp = rt_band1
        call VecDestroy(Xvec_Band1,ierr)
      CASE(2)
        current_band_temp = rt_band2
        call VecDestroy(Xvec_Band2,ierr)
      CASE(3)
        current_band_temp = rt_band3
        call VecDestroy(Xvec_Band3,ierr)
      CASE(4)
        current_band_temp = rt_band4
        call VecDestroy(Xvec_Band4,ierr)
      CASE(5)
        current_band_temp = rt_band5
        call VecDestroy(Xvec_Band5,ierr)
      CASE DEFAULT
        if (dr_globalMe .eq. MASTER_PE) write(*,"(A,I0,I0)") 'Band Count, No of bands: ', band_count, rt_freqbands
        call Driver_abortFlash("[VETTAM]: Something is wrong; band_count is greater than 5, should not be happening!")
    END SELECT
  end do

  call VecDestroy(Yvec,ierr)

  call KSPDestroy(ksp,ierr)

  ! At this point there should be no allocated memory from Petsc be left.
  ! This will (for debugging purposes) show left over memory.
#ifdef VET_DEBUG
  !call PetscMallocDump(ierr)
#endif

  alloced = .FALSE.

END SUBROUTINE Petsc_dealloc


SUBROUTINE Petsc_step(dt,reason)
  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe
  use Timers_interface, ONLY: Timers_start,Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
  real, intent(IN) :: dt
  type(PetscInt), intent(OUT) :: reason
  type(PetscErrorCode) :: ierr
  type(PetscInt) :: iterations
  call MatZeroEntries(Amat,ierr)

#ifdef VET_DEBUG
  call Logfile_stamp("Matrix setting starts", "[VET_DEBUG] ")
#endif

  call Timers_start("vet_setcoeff")
  call rt_fillMatrix(dt)
  call Timers_stop("vet_setcoeff")
#ifdef VET_DEBUG
  call Logfile_stamp("Matrix setting done", "[VET_DEBUG] ")
#endif

  call Timers_start("vet_vecassemble")
  call VecAssemblyBegin(Yvec,ierr)
  call VecAssemblyEnd(Yvec,ierr)
  call Timers_stop("vet_vecassemble")
  IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while assembling Petsc Vector.")

#ifdef VET_DEBUG
  call Logfile_stamp("Vector assembly done", "[VET_DEBUG] ")
#endif

  call Timers_start("vet_matassemble")
  call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
  call Timers_stop("vet_matassemble")
  IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while assembling Petsc Matrix.")
#ifdef VET_DEBUG
  call Logfile_stamp("Matrix assembly done", "[VET_DEBUG] ")
#endif

  if(rt_debug) then 
    print *, 'YVECTOR......'
    call VecView(YVec,PETSC_VIEWER_STDOUT_WORLD,ierr)     
    print *, 'Amat......'
    call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr)
  endif

  call Timers_start("vet_kspsetup")
  call KSPSetOperators(ksp,Amat,Amat,ierr)
  call KSPSetUp(ksp, ierr)
  call Timers_stop("vet_kspsetup")
  call Timers_start("vet_kspsolve")
  if(current_band .eq. rt_band1) then
    call KSPSolve(ksp,Yvec,Xvec_Band1,ierr)
    call KSPGetSolution(ksp,Xvec_Band1,ierr)
  else if(current_band .eq. rt_band2) then
    call KSPSolve(ksp,Yvec,Xvec_Band2,ierr)
    call KSPGetSolution(ksp,Xvec_Band2,ierr)
  else if(current_band .eq. rt_band3) then
    call KSPSolve(ksp,Yvec,Xvec_Band3,ierr)
    call KSPGetSolution(ksp,Xvec_Band3,ierr)
  else if(current_band .eq. rt_band4) then
    call KSPSolve(ksp,Yvec,Xvec_Band4,ierr)
    call KSPGetSolution(ksp,Xvec_Band4,ierr)
  else if(current_band .eq. rt_band5) then
    call KSPSolve(ksp,Yvec,Xvec_Band5,ierr)
    call KSPGetSolution(ksp,Xvec_Band5,ierr)
  endif
  call Timers_stop("vet_kspsolve")
  IF(ierr.NE.0) call Driver_abortFlash("[VET]: Error in KSPSolve.")
#ifdef VET_DEBUG
  call Logfile_stamp("Operator assembly done", "[VET_DEBUG] ")
#endif

  call KSPGetIterationNumber(ksp, iterations, ierr)
  call KSPGetConvergedReason(ksp, reason, ierr)

#ifdef VET_DEBUG
  call Logfile_stamp("Matrix operations done.", "[VET_DEBUG] ")
#endif

  ! The iterations returned by above is number_iters - 1
  iterations = iterations + 1

  !Iteration log stuff
  if(current_band .eq. rt_band1) then 
    cumulative_iter_Band1 = cumulative_iter_Band1 + iterations
    if(min_iter_Band1>iterations) min_iter_Band1 = iterations
    if(max_iter_Band1<iterations) max_iter_Band1 = iterations
  else if(current_band .eq. rt_band2) then 
    cumulative_iter_Band2 = cumulative_iter_Band2 + iterations
    if(min_iter_Band2>iterations) min_iter_Band2 = iterations
    if(max_iter_Band2<iterations) max_iter_Band2 = iterations
  else if(current_band .eq. rt_band3) then 
    cumulative_iter_Band3 = cumulative_iter_Band3 + iterations
    if(min_iter_Band3>iterations) min_iter_Band3 = iterations
    if(max_iter_Band3<iterations) max_iter_Band3 = iterations
  else if(current_band .eq. rt_band4) then 
    cumulative_iter_Band4 = cumulative_iter_Band4 + iterations
    if(min_iter_Band4>iterations) min_iter_Band4 = iterations
    if(max_iter_Band4<iterations) max_iter_Band4 = iterations
  else if(current_band .eq. rt_band5) then 
    cumulative_iter_Band5 = cumulative_iter_Band5 + iterations
    if(min_iter_Band5>iterations) min_iter_Band5 = iterations
    if(max_iter_Band5<iterations) max_iter_Band5 = iterations
  endif
  

  !Maximum iterations reached and not converged
  IF(reason .eq. -3) then 
    if(dr_globalMe .eq. MASTER_PE) print *, 'PETSC did not converge, attempting subcycling...'
    return
  ENDIF

  IF(dr_globalMe.eq.MASTER_PE) THEN

    IF(reason .lt. 0) THEN
      print *,"No. of PETSC iterations: ", iterations
      print *,"PETSC Converged reason: ", reason


      SELECT CASE(reason)

      CASE(-2)
        call Driver_abortFlash("[VET]: KSP diverged due to Null.")
      CASE(-3)
       call Driver_abortFlash("[VET]: Maximum Iterations reached, and no convergence.")
      CASE(-4)
        call Driver_abortFlash("[VET]: Solution is diverging with residual>dtol.")
      CASE(-5)
        call Driver_abortFlash("[VET]: Breakdown detected, could be due to singular matrix.")
      CASE(-6)
        call Driver_abortFlash("[VET]: Breakdown detected in BICG method.")
      CASE(-7)
        call Driver_abortFlash("[VET]: Breakdown detected in BICG method.")
      CASE(-8)
        call Driver_abortFlash("[VET]: PC is indefinite, and method requires positive definite.")
      CASE(-9)
        call Driver_abortFlash("[VET]: NAN or infinite values detected in matrix.")
      CASE(-10)
        call Driver_abortFlash("[VET]: Matrix is indefinite, and method requires positive definite.")
      CASE(-11)
        call Driver_abortFlash("[VET]: Unable to build requested PC.")
      CASE DEFAULT
        call Driver_abortFlash("[VET] : Unrecognised reason returned by Petsc.")

      END SELECT 
  
    ELSE
#ifdef VET_DEBUG      
      print *,"No. of PETSC iterations: ", iterations
      print *,"PETSC Converged reason: ", reason
#endif      

    ENDIF
  END IF  

  if(rt_debug) then   
    print *, 'XVECTOR......'
    call VecView(Xvec_Band1,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if(current_band .eq. rt_band1) then
      call VecView(Xvec_Band1,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else if(current_band .eq. rt_band2) then
      call VecView(Xvec_Band2,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else if(current_band .eq. rt_band3) then
      call VecView(Xvec_Band3,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else if(current_band .eq. rt_band4) then
      call VecView(Xvec_Band4,PETSC_VIEWER_STDOUT_WORLD,ierr)
    else if(current_band .eq. rt_band5) then
      call VecView(Xvec_Band5,PETSC_VIEWER_STDOUT_WORLD,ierr)
    endif
  endif
#ifdef VET_DEBUG
  call Logfile_stamp("Result being copied", "[VET_DEBUG] ")
#endif
  if(current_band .eq. rt_band1) then
    call CopyResult(Xvec_Band1)
  else if(current_band .eq. rt_band2) then
    call CopyResult(Xvec_Band2)
  else if(current_band .eq. rt_band3) then
    call CopyResult(Xvec_Band3)
  else if(current_band .eq. rt_band4) then
    call CopyResult(Xvec_Band4)
  else if(current_band .eq. rt_band5) then
    call CopyResult(Xvec_Band5)
  endif
#ifdef VET_DEBUG
  call Logfile_stamp("Result copying done", "[VET_DEBUG] ")
#endif

END SUBROUTINE Petsc_step


! Set Radiation quantities to the solution vector values. 
SUBROUTINE CopyResult(X)
  use RadTrans_data
  use Grid_interface
#ifdef UEUV_VAR
  use rt_ionisedata, ONLY: multiple_ionbands
#endif
  implicit none
  type(tVec), INTENT(IN) :: X
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b,blockID,i,j,k,l
  real, pointer :: solnData(:,:,:,:)
  type(PetscErrorCode) ::  ierr
  type(PetscScalar), dimension(:), POINTER :: arr
  integer, dimension(MDIM) :: cornerID, stride, Ncell, globalIndexLimits
  INTEGER, EXTERNAL :: localblockindex
  
  call VecGetArrayReadF90(X, arr, ierr)

  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkCornerID(blockID,cornerID,stride)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    cornerID = (cornerID-1)/(stride) + 1
    ! Number of cells assuming grid fully refined.

    !For getting global row index
    call Grid_getGlobalIndexLimits(globalIndexLimits)
    Ncell(1) = globalIndexLimits(IAXIS)/(stride(IAXIS))
    Ncell(2) = globalIndexLimits(JAXIS)/(stride(JAXIS))
    Ncell(3) = globalIndexLimits(KAXIS)/(stride(KAXIS))
    do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
      do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
          ! flatten the index
          l = localblockindex(i-blkLimits(LOW,IAXIS),&
                      j-blkLimits(LOW,JAXIS),&
                      k-blkLimits(LOW,KAXIS),&
                      b-1)
          solnData(ERAD_VAR,i,j,k) = arr(l+1)
          solnData(ERAD_VAR,i,j,k) = max(solnData(ERAD_VAR,i,j,k),rt_smalle)
          solnData(MOHX_VAR,i,j,k) = arr(l+2)
          

#if NDIM>1          
          solnData(MOHY_VAR,i,j,k) = arr(l+3)
          
#if NDIM>2
          solnData(MOHZ_VAR,i,j,k) = arr(l+4)
          
#endif
#endif
          !If current band is EUV, store the volumetric rate of absorbed EUV photons
#ifdef UEUV_VAR
          !Stores energy absorbed per unit area per unit volume by the gas -- used to compute ionisation rate
          if(current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
          & current_band .eq. 'EUV_15P2_INFTY') then
            if(current_band .eq. 'EUV_13P6_15P2' .and. multiple_ionbands) then
              solnData(HEUV_VAR,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(ERAD_VAR,i,j,k)
            else
              solnData(UEUV_VAR,i,j,k) = solnData(TAUH_VAR,i,j,k) * rt_speedlt * solnData(ERAD_VAR,i,j,k)
            endif
          endif
#endif 

          !Store solutions in different containers if they are defined:
          !This is relevant if you have multiple bands and need the solution, as ERAD and MOHX/Y/Z are overwritten
          !Note: These need to be defined in the config of the simulation by the user -- they are not present by default

          if(current_band .eq. rt_band1) then
#ifdef B1ER_VAR
            solnData(B1ER_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
#endif
#ifdef B1FX_VAR
            solnData(B1FX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k)
#endif
#ifdef B1FY_VAR
            solnData(B1FY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k)
#endif
#ifdef B1FZ_VAR
            solnData(B1FZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k)
#endif

          else if(current_band .eq. rt_band2) then
#ifdef B2ER_VAR
            solnData(B2ER_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
#endif
#ifdef B2FX_VAR
            solnData(B2FX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k)
#endif
#ifdef B2FY_VAR
            solnData(B2FY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k)
#endif
#ifdef B2FZ_VAR
            solnData(B2FZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k)
#endif
          else if(current_band .eq. rt_band3) then
#ifdef B3ER_VAR
            solnData(B3ER_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
#endif
#ifdef B3FX_VAR
            solnData(B3FX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k)
#endif
#ifdef B3FY_VAR
            solnData(B3FY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k)
#endif
#ifdef B3FZ_VAR
            solnData(B3FZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k)
#endif

          else if(current_band .eq. rt_band4) then
#ifdef B4ER_VAR
            solnData(B4ER_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
#endif
#ifdef B4FX_VAR
            solnData(B4FX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k)
#endif
#ifdef B4FY_VAR
            solnData(B4FY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k)
#endif
#ifdef B4FZ_VAR
            solnData(B4FZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k)
#endif

          else if(current_band .eq. rt_band5) then
#ifdef B5ER_VAR
            solnData(B5ER_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
#endif
#ifdef B5FX_VAR
            solnData(B5FX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k)
#endif
#ifdef B5FY_VAR
            solnData(B5FY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k)
#endif
#ifdef B5FZ_VAR
            solnData(B5FZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k)
#endif
          else
            call Driver_abortFlash("[rt_petsc]: current_band >3; Maximum bands is set to 3.")
          endif

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

  call VecRestoreArrayReadF90(X, arr, ierr)

END SUBROUTINE CopyResult

! Set solution vector guess to radiation quantities in initial state 
SUBROUTINE CopyGuess(X)
  use RadTrans_data
  use Grid_interface
  use Driver_data, ONLY: dr_globalMe
  implicit none
  type(tVec), INTENT(IN) :: X
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: b,blockID,i,j,k,l
  real, pointer :: solnData(:,:,:,:)
  type(PetscErrorCode) ::  ierr
  type(PetscScalar), dimension(:), POINTER :: arr
  integer, dimension(MDIM) :: cornerID, stride, Ncell, globalIndexLimits
  INTEGER, EXTERNAL :: localblockindex

  call VecGetArrayReadF90(X, arr, ierr)
  CHKERRA(ierr)
  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkCornerID(blockID,cornerID,stride)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    cornerID = (cornerID-1)/(stride) + 1
    ! Number of cells assuming grid fully refined.

    !For getting global row index
    call Grid_getGlobalIndexLimits(globalIndexLimits)
    Ncell(1) = globalIndexLimits(IAXIS)/(stride(IAXIS))
    Ncell(2) = globalIndexLimits(JAXIS)/(stride(JAXIS))
    Ncell(3) = globalIndexLimits(KAXIS)/(stride(KAXIS))
    do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
      do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
          ! flatten the index
           l = localblockindex(i-blkLimits(LOW,IAXIS),&
                      j-blkLimits(LOW,JAXIS),&
                      k-blkLimits(LOW,KAXIS),&
                      b-1)
          arr(l+1) = solnData(ERAD_VAR,i,j,k) 
          arr(l+2) = solnData(MOHX_VAR,i,j,k) 
#if NDIM>1          
          arr(l+3) = solnData(MOHY_VAR,i,j,k) 
#if NDIM>2
          arr(l+4) = solnData(MOHZ_VAR,i,j,k)
#endif
#endif          

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

  call VecRestoreArrayReadF90(X, arr, ierr)
END SUBROUTINE CopyGuess

SUBROUTINE InitNeghLevels()
  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  IMPLICIT NONE
  integer :: b, s, p, blockID, ierr

  IF(ALLOCATED(NeghLevels)) DEALLOCATE(NeghLevels)
  IF(ALLOCATED(SurrBlkSum)) DEALLOCATE(SurrBlkSum)
  !! Stores the refinement level of neighbor along each direction.
  allocate (NeghLevels(blockCount, 1:RIGHT_EDGE, 1:RIGHT_EDGE, 1:RIGHT_EDGE))
  allocate (SurrBlkSum(blockCount))

  ! Given a Blkno and a Procno rblockListAll (reverse block list all)
  ! will return the global contigous number of that block
  IF(.NOT.ALLOCATED(rblockListAll)) &
    ALLOCATE(rblockListAll(MAXBLOCKS,0:dr_globalNumProcs-1))
  rblockListAll(:,:) = -1
  DO b = 1, blockCount
    blockID = blockList(b)
    rblockListAll(blockID, dr_globalMe) = 1 ! This is a leaf block
  END DO
  call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, &
                     rblockListAll, MAXBLOCKS, MPI_INTEGER, &
                     MPI_COMM_WORLD,ierr)
  ! All leaf blocks are marked with a "1", now make the cummulative Sum.
  ! In principle any other unique numbering would also be fine, it just has to
  ! be contingues and the same on all procs.
  s = 0
  DO p = 0, dr_globalNumProcs-1
    DO b = 1, MAXBLOCKS
      IF(rblockListAll(b,p).EQ.1) THEN
        rblockListAll(b,p) = s
!        print *,"b -> s: ", b, " -> ", s
        s = s + 1
      END IF
    END DO
  END DO

  call getNeghLevels()

END SUBROUTINE InitNeghLevels

SUBROUTINE getNeghLevels()
  use RadTrans_data
  use gr_interface,     ONLY : gr_findAllNeghID, gr_getBlkHandle
  use tree,             ONLY:  surr_blks
  use Grid_data,        ONLY : gr_meshMe
  IMPLICIT NONE
  INTEGER :: lb,k,j,i,numNegh,eachNegh,neghBlk,allCenters,level,blockID,lrefine
  integer, dimension(BLKNO:PROCNO) :: neghBlkProc
  integer, dimension(BLKNO:TYPENO) :: negh_prop
! This is adapted from gr_hypreSetupGrid, l. 359ff 
  !!-----------------------------------------------------------------------
  !!     19. This entire step is devoted to computing the refinement
  !!         levels (i.e HYPRE part number) of our neighbors. The process is
  !!         a bit warped as we need off processor information.
  !!         NOTE: i)  A dependency on PARAMESH is built into the system, as a
  !!                    tradeoff we avoid off processor communication.
  !!               ii) Some block information might be cached by PARAMESH depending
  !!                   on last performed operation.
  !!               iii)We store this information locally.
  !!----------------------------------------------------------------------

  allCenters = 2**NDIM
  do lb=1, blockCount
     blockID = blockList(lb)
     call Grid_getBlkRefineLevel(blockList(lb),lrefine)
     call gr_findAllNeghID(blockID, SurrBlkSum(lb))
     do k = LEFT_EDGE , LEFT_EDGE+K3D*(RIGHT_EDGE-1)
        do j = LEFT_EDGE , LEFT_EDGE+K2D*(RIGHT_EDGE-1)
           do i = LEFT_EDGE , RIGHT_EDGE
              if (allCenters /= i*j*k) then
                 numNegh = SurrBlkSum(lb) % regionInfo(i,j,k) % numNegh
                 do eachNegh = 1, numNegh
                    neghBlkProc(BLKNO:PROCNO) = SurrBlkSum(lb) % regionInfo(i,j,k) % details(BLKNO:PROCNO,eachNegh)
                    if (neghBlkProc (PROCNO) .EQ. gr_meshMe) then
                       call Grid_getBlkRefineLevel(neghBlkProc(BLKNO),level)
                    else
                       !! Search for cached blocks
                       neghBlk = -1
                       call gr_getBlkHandle(neghBlkProc(BLKNO), neghBlkProc(PROCNO), neghBlk)

                       if (neghBlk == -1) then !! information is not locally cached.
                          negh_prop(:) = surr_blks(:,i,j,k,blockID)
                          !! piece of code extracted from gr_findAllNeghID.F90
                          if (negh_prop(BLKNO) == NONEXISTENT) then
                             !(Assumes NONEXISTENT is equal to -1).
                             !! If this query returns NONEXISTENT, that means that the
                             !! the neighbor is at a lower resolution.
                             level = lrefine - 1
                          else if (negh_prop(TYPENO) == PARENT_BLK) then
                             !! PARENT_BLK indicates that there may be more than one
                             !! block in the neighborhood here, and they will be
                             !! at a higher resolution.
                             level = lrefine + 1
                          else if (negh_prop(TYPENO) == LEAF) then
                             !! This is the situation when the neighbor is at the
                             !! same level of resulution. There is only one neighbor
                             !! very simply found.
                             level = lrefine
                          else
                             call Driver_abortFlash("unable to find negh")
                          end if
                       else
                          !! information is locally cached.
                          call Grid_getBlkRefineLevel(neghBlk,level)
                       end if
                    end if
                    NeghLevels(lb,i,j,k) = level
                 end do
              end if
           end do
        end do
     end do
  end do
END SUBROUTINE getNeghLevels

SUBROUTINE MyConvergenceTest(ksp_, it, rnorm, reason, ctx_, ierr)
  use RadTrans_data
  implicit none
  type(tKSP) :: ksp_
  PetscInt it
  PetscFortranAddr ctx_
  PetscErrorCode ierr
  PetscReal rnorm
  KSPConvergedReason :: reason
  call KSPConvergedDefault(ksp_,it,rnorm,reason,ctx_,ierr)
  IF(it.lt.rt_minits) THEN
    ierr = 0
    reason = 0
  END IF 
END SUBROUTINE MyConvergenceTest

SUBROUTINE Print_PETSCConvergence()
  use RadTrans_data
  use Driver_data, ONLY: dr_globalMe, dr_nstep, dr_simTime, dr_dt
#ifdef UEUV_VAR
  use rt_ionisedata, ONLY: multiple_ionbands
#endif

  implicit none
  integer :: min_iter, max_iter, cumulative_iter
  character(len=MAX_STRING_LENGTH), save :: PetscPrint
  !Variables for convergence file logging
  logical, save :: firstCall = .true.
  integer, save :: firstcallband_count = 0
  integer, parameter :: funit = 25
  character(len=80) :: outfile
  logical :: exist

  if(firstCall .and. Petsc_Log) then
    if(current_band .eq. rt_band1) then 
      outfile = 'PetscConvergence_' // trim(rt_band1) // '.log'
    else if(current_band .eq. rt_band2) then 
      outfile = 'PetscConvergence_' // trim(rt_band2) // '.log'
    else if(current_band .eq. rt_band3) then 
      outfile = 'PetscConvergence_' // trim(rt_band3) // '.log'
    else if(current_band .eq. rt_band4) then
      outfile = 'PetscConvergence_' // trim(rt_band4) // '.log'
    else if(current_band .eq. rt_band5) then
      outfile = 'PetscConvergence_' // trim(rt_band5) // '.log'
    endif
    !Open file
    
    if(dr_globalMe .eq. MASTER_PE) then
      write(*,"(A,A)") 'Writing PETSC convergence info to ',trim(outfile)
      inquire(file=outfile,exist=exist)
      !Create the file only if it does not exist
      if(.not. exist) then
        open(funit, file=trim(outfile), status='new')
        write(funit,'(A10,3X,A5,3X,A6,3X,A9,3X,A12,3X,A12)') '#[00]nStep', '[01]t','[02]dt',&
          '[03]nIter','[04]min_Iter','[05]max_Iter'
        close(funit)
      endif
    endif
    firstcallband_count = firstcallband_count + 1
    !Switch off firstCall only if looped through all the bands once, since each band has a separate log file
    if(firstcallband_count .eq. rt_freqbands) &
      & firstCall = .false.
  endif

  if(current_band .eq. rt_band1) then
    outfile = 'PetscConvergence_' // trim(rt_band1) // '.log'
    cumulative_iter = cumulative_iter_band1
    max_iter = max_iter_band1
    min_iter = min_iter_band1
  else if(current_band .eq. rt_band2) then
    outfile = 'PetscConvergence_' // trim(rt_band2) // '.log'
    cumulative_iter = cumulative_iter_band2
    max_iter = max_iter_band2
    min_iter = min_iter_band2
  else if(current_band .eq. rt_band3) then
    outfile = 'PetscConvergence_' // trim(rt_band3) // '.log'
    cumulative_iter = cumulative_iter_band3
    max_iter = max_iter_band3
    min_iter = min_iter_band3
  else if(current_band .eq. rt_band4) then
    outfile = 'PetscConvergence_' // trim(rt_band4) // '.log'
    cumulative_iter = cumulative_iter_band4
    max_iter = max_iter_band4
    min_iter = min_iter_band4
  else if(current_band .eq. rt_band5) then
    outfile = 'PetscConvergence_' // trim(rt_band5) // '.log'
    cumulative_iter = cumulative_iter_band5
    max_iter = max_iter_band5
    min_iter = min_iter_band5
  endif

  PetscPrint = 'PETSC iterations Used (' // trim(current_band) // '): '

  !Log number of iterations
  if(Petsc_Log) then 
    if(dr_globalMe .eq. MASTER_PE) then
      open(funit, file=trim(outfile), position='APPEND')
      !TODO: Get time spent in setting matrix and obtaining soln
      write(funit,'(I0,1X,ES16.9,1X,ES16.9,1X,I0,1X,I0,1X,I0)') dr_nstep, &
        dr_dt, dr_simTime,cumulative_iter,min_iter,max_iter
      close(funit)
    endif
  endif

  if(dr_globalMe .eq. MASTER_PE) &
    & write(*,"(A,I0,1X,I0,1X,I0)") trim(PetscPrint), cumulative_iter

END SUBROUTINE Print_PETSCConvergence



