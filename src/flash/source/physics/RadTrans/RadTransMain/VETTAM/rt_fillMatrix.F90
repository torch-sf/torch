!!****f* source/physics/RadTrans/RadTransMain/VETTAM/rt_fillMatrix
!!
!! NAME
!!  
!!  rt_fillMatrix
!!
!!
!! SYNOPSIS
!! 
!!  rt_fillMatrix(dt)
!!  
!! DESCRIPTION
!!
!!  This set of subroutine calculates and allocates the coefficient matrix 
!!  values of the linearised radiation moment equation subsystem, i.e. the 
!!  matrix stencil is prepared here. 
!!  The matrix is then passed to the PETSc library for solving the linear system.
!!
!!***
!! NOTE: Currently the matrix values for each row among the NDIM+1 equations are
!!       set seperately. Might be more efficient to pass it in one call to MatSet.
!!       For this however, I would need to set a large coeff/ID array with all possible
!!       P&C of indices - i.e. 28 indices, and set unwanted indices to negative values.
!!       So far setting the matrix does not affect performance, but can be revisited in the future.
!!***


#include "Flash.h"
#include "constants.h"
#include "petsc/finclude/petscksp.h"

SUBROUTINE rt_fillMatrix(dt)

  use RadTrans_data
  use Grid_interface
  use SemenovOpacities, only: getOpacity_planck, getOpacity_rosseland
#ifdef IONY_MSCALAR
  use Eos_data, only: eos_singleSpeciesA
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use rt_ionisedata, ONLY: alpha_ground_constant, alpha_type, ion_ots, hnu, alpha_B
  use rt_ionisemodule, ONLY: get_recombination_ground, get_recombination_coefficient
#elif defined(IHA_SPEC) && defined(UEUV_VAR)
  use rt_ionisedata, ONLY : hpA, elecA, alpha_ground_constant, alpha_type, ion_ots, hnu, alpha_B
  use rt_ionisemodule, ONLY: get_recombination_ground, get_recombination_coefficient
#endif
  implicit none

  real, intent(in) :: dt
  type(PetscErrorCode) ::  ierr
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, faces_BC,onBoundary
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: i, j, k, b, blockID, mylevel, axis, cell_offset, myPE
  real :: cellsize(MDIM), coords(MDIM)
  real,dimension(5*NDIM+1) :: er_coeff, frx_coeff,fry_coeff,frz_coeff
  type(PetscInt) :: coeff_ID(NDIM+1,5*NDIM+1),vecID(NDIM+1)
  type(PetscInt) :: no_rows = NDIM +1, no_columns = 5*NDIM+1
  real :: coeff(NDIM+1,5*NDIM+1),vec(NDIM+1)
  integer, dimension(MDIM) :: cornerID, stride, Ncell
  real :: d1x, d1y, d1z, d2, Ci0, Ci1, Cj0, Cj1, Ck0, Ck1, huge_val, tau, tau_avg
  real :: opac_planck, opac_rosseland, vx, vy, vz, tau_hlle, wavespeed_cells
  real, dimension(2,NDIM) :: wavespeed, wavespeed_negh
  integer, dimension(MDIM) :: globalIndexLimits
  logical,DIMENSION(2,NDIM) :: boundary_inblock, FCInterface
  INTEGER, EXTERNAL :: globalmatrixindex
  integer, dimension(5*NDIM+1) :: er_columnID,frx_columnID,fry_columnID,frz_columnID
  ! Stores Eddington Tensor components for each direction neighbours i.e. i/j/k-1,and i/j/k+1
  real, dimension(2,NDIM) :: xxed, xyed, xzed, yyed, zzed, yzed, tau_hlle_neigh
  real :: mH, iony, nH, ne, nHplus, diffuse_rec_term, temp, alpha_cell
  

  ! Set opacity everywhere. This makes the opacities implicit for Picard iterations
  if(current_band .eq. 'IR' .and. rt_hydro_type .eq. 2) call rt_setOpacity()

  !Reset values in matrix to zero: this is done as coefficients are given an ADD_VALUES below
  call MatZeroEntries(Amat,ierr)
  CHKERRA(ierr)
  do b = 1, blockCount

    blockID = blockList(b)
    ! get block limits
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    ! Get pointer to solution data
    call Grid_getBlkPtr(blockID,solnData)
    ! get block cell size and refine level.
    call Grid_getDeltas(blockID,cellsize)
    call Grid_getBlkRefineLevel(blockID,mylevel)
    
    ! Flag to determine whether physical boundary at any of block faces.
    call Grid_getBlkBC(blockID,faces_BC,onBoundary) 

    call Grid_getGlobalIndexLimits(globalIndexLimits)

    !Get processor ID
    call Driver_getMype(GLOBAL_COMM,myPE)

    ! Locate the corner cell of the block in the matrix subsystem.
    call Grid_getBlkCornerID(blockID,cornerID,stride)
    ! Appropriate corner ID when stride !=1
    cornerID = (cornerID-1)/(stride) + 1

    ! Number of cells assuming grid fully refined.
    Ncell(1) = globalIndexLimits(IAXIS)/(stride(IAXIS))
    Ncell(2) = globalIndexLimits(JAXIS)/(stride(JAXIS))
    Ncell(3) = globalIndexLimits(KAXIS)/(stride(KAXIS))

    !call locate_globalindex()
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          ! By default assume it is not a boundary
          boundary_inblock(LOW:HIGH,IAXIS:NDIM) = .false.
          
          !Check for possible boundaries
          if(i .eq. blkLimits(LOW, IAXIS)) call CheckBoundary(LOW,IAXIS,boundary_inblock)
          if(i .eq. blkLimits(HIGH, IAXIS)) call CheckBoundary(HIGH,IAXIS,boundary_inblock)
#if NDIM>1
          if(j .eq. blkLimits(LOW, JAXIS)) call CheckBoundary(LOW,JAXIS,boundary_inblock)
          if(j .eq. blkLimits(HIGH, JAXIS)) call CheckBoundary(HIGH,JAXIS,boundary_inblock)
#if NDIM>2
          if(k .eq. blkLimits(LOW, KAXIS)) call CheckBoundary(LOW,KAXIS,boundary_inblock)
          if(k .eq. blkLimits(HIGH, KAXIS)) call CheckBoundary(HIGH,KAXIS,boundary_inblock)
#endif
#endif

          call getrowcolumnIDs(er_columnID,frx_columnID,fry_columnID,frz_columnID)

          d1x = dt/(2*cellsize(IAXIS))
          d1y = dt/(2*cellsize(JAXIS))
          d1z = dt/(2*cellsize(KAXIS))
          d2 = dt * rt_speedlt

          opac_planck = solnData(TAUP_VAR,i,j,k)
          opac_rosseland = solnData(TAUR_VAR,i,j,k)

          do axis =1, NDIM
            do cell_offset=LOW, HIGH
              if(axis .eq. 1) then
                xxed(cell_offset,axis) = solnData(XXED_VAR,i+2*cell_offset-3,j,k)
                xyed(cell_offset,axis) = solnData(XYED_VAR,i+2*cell_offset-3,j,k)
                xzed(cell_offset,axis) = solnData(XZED_VAR,i+2*cell_offset-3,j,k)
                yyed(cell_offset,axis) = solnData(YYED_VAR,i+2*cell_offset-3,j,k)
                yzed(cell_offset,axis) = solnData(YZED_VAR,i+2*cell_offset-3,j,k)
                zzed(cell_offset,axis) = solnData(ZZED_VAR,i+2*cell_offset-3,j,k)

                tau_hlle_neigh(cell_offset,axis) = (10.0*cellsize(axis)*solnData(TAUR_VAR,i+2*cell_offset-3,j,k))**2 /&
                                              (2*xxed(cell_offset,axis)) 
#if NDIM>1
              else if(axis .eq. 2) then
                xxed(cell_offset,axis) = solnData(XXED_VAR,i,j+2*cell_offset-3,k)
                xyed(cell_offset,axis) = solnData(XYED_VAR,i,j+2*cell_offset-3,k)
                xzed(cell_offset,axis) = solnData(XZED_VAR,i,j+2*cell_offset-3,k)
                yyed(cell_offset,axis) = solnData(YYED_VAR,i,j+2*cell_offset-3,k)
                yzed(cell_offset,axis) = solnData(YZED_VAR,i,j+2*cell_offset-3,k)
                zzed(cell_offset,axis) = solnData(ZZED_VAR,i,j+2*cell_offset-3,k)

                tau_hlle_neigh(cell_offset,axis) = (10.0*cellsize(axis)*solnData(TAUR_VAR,i,j+2*cell_offset-3,k))**2 /&
                                              (2*yyed(cell_offset,axis))
#if NDIM>2
              else if(axis .eq. 3) then
                xxed(cell_offset,axis) = solnData(XXED_VAR,i,j,k+2*cell_offset-3)
                xyed(cell_offset,axis) = solnData(XYED_VAR,i,j,k+2*cell_offset-3)
                xzed(cell_offset,axis) = solnData(XZED_VAR,i,j,k+2*cell_offset-3)
                yyed(cell_offset,axis) = solnData(YYED_VAR,i,j,k+2*cell_offset-3)
                yzed(cell_offset,axis) = solnData(YZED_VAR,i,j,k+2*cell_offset-3)
                zzed(cell_offset,axis) = solnData(ZZED_VAR,i,j,k+2*cell_offset-3)

                tau_hlle_neigh(cell_offset,axis) = (10.0*cellsize(axis)*solnData(TAUR_VAR,i,j,k+2*cell_offset-3))**2 /&
                                              (2*zzed(cell_offset,axis))
#endif
#endif
              endif

            end do
          end do


          if(i .eq. blkLimits(LOW,IAXIS)) call SetEddBoundary(LOW,IAXIS)
          
          if(i .eq. blkLimits(HIGH,IAXIS)) call SetEddBoundary(HIGH,IAXIS)
#if NDIM>1
          if(j .eq. blkLimits(LOW,JAXIS)) call SetEddBoundary(LOW,JAXIS)

          if(j .eq. blkLimits(HIGH,JAXIS)) call SetEddBoundary(HIGH,JAXIS)
          
#if NDIM>2
          if(k .eq. blkLimits(LOW,KAXIS)) call SetEddBoundary(LOW,KAXIS)

          if(k .eq. blkLimits(HIGH,KAXIS)) call SetEddBoundary(HIGH,KAXIS)
          
#endif
#endif
          !Set wavespeeds at Left and Right interfaces for this cell
          call Set_Wavespeeds()
          !Set the coefficients for the equation governing this cell
          call Set_Coeffs()

          !Combine coefficient ID's into a single container
          coeff_ID(1,:) = er_columnID(:)
          coeff_ID(2,:) = frx_columnID(:)
#if NDIM>1  
          coeff_ID(3,:) = fry_columnID(:)
#if NDIM>2
          coeff_ID(4,:) = frz_columnID(:)
#endif
#endif

          !Combine coefficients into a single container
          coeff(1,:) = er_coeff(:)
          coeff(2,:) = frx_coeff(:)
#if NDIM>1  
          coeff(3,:) = fry_coeff(:)
#if NDIM>2
          coeff(4,:) = frz_coeff(:)
#endif
#endif
          ! Vector values

          if(current_band .eq. 'IR') then

            !Terms coming from the (first-order) discretisation of dE_r/dt, dF/dt. Only present for IR band; other bands assumed instantaneous 
            vec(1) = solnData(ERAO_VAR,i,j,k)
            vec(2) = solnData(MHXO_VAR,i,j,k)
#if NDIM>1
            vec(3) = solnData(MHYO_VAR,i,j,k)
#if NDIM>2
            vec(4) = solnData(MHZO_VAR,i,j,k)
#endif
#endif
            !Gas/dust emission term; i.e. the a_RT^4 term in the radiation energy density equation
#ifdef TEMP_VAR
            !Only add gas-radiation exchange term if hydro switched on; also do not do this if radiative equilibrium assumed (hydro_type .eq. 1)
            if(rt_update_hydro .and. rt_hydro_type .ne. 1) then
              vec(1) = vec(1) + opac_planck*d2*rt_radconst*(solnData(TEMP_VAR,i,j,k)**4)
            endif
#endif
          ! Contribution from recombinations:
          ! i) Recombinations to ground state for EUV, ii) Recombinations to all other levels (B-type) for FUV (TODO: This is not implemented now)
#ifdef UEUV_VAR                   
          else if((current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
      & current_band .eq. 'EUV_15P2_INFTY') .and. .not. ion_ots) then
            vec(1) = 0.0
            !Get the coefficient for recombinations to ground state
#ifdef TEMP_VAR
            temp = solnData(TEMP_VAR,i,j,k)
#elif defined(TGAS_VAR)
            temp = solnData(TGAS_VAR,i,j,k)
#endif
            if(alpha_type .eq. 'default') then
              alpha_cell = get_recombination_ground(temp,ion_ots) !alpha_A
            else if(alpha_type .eq. 'constant') then
              alpha_cell = alpha_ground_constant !alpha_A
            else
              print *, alpha_type
              call Driver_abortFlash("[rt_fillMatrix]:alpha_type should be 'constant' or 'default'; check!")
            endif
#ifdef IONY_MSCALAR
            call PhysicalConstants_get("proton mass",mH)
            iony = solnData(IONY_MSCALAR,i,j,k)
            nH = solnData(DENS_VAR,i,j,k)/(eos_singleSpeciesA*mH)
            nHplus = nH*(1.-iony)
            ne = nHplus
            diffuse_rec_term = alpha_cell * nHplus* ne * hnu * dt
            
            !KROME Version
#elif defined(IHA_SPEC) && defined(UEUV_VAR)
            nHplus = solnData(DENS_VAR,i,j,k)*solnData(IHP_SPEC,i,j,k)/hpA
            ne     = solnData(DENS_VAR,i,j,k)*solnData(ELEC_SPEC,i,j,k)/(elecA)
            diffuse_rec_term = alpha_cell * nHplus* ne * hnu * dt
#else
            diffuse_rec_term = 0.0
#endif
            vec(1) = vec(1) + diffuse_rec_term
            vec(2) = 0.0
#if NDIM>1
            vec(3) = 0.0
#if NDIM>2
            vec(4) = 0.0
#endif
#endif

#endif !End UEUV_VAR IF condition
          else
            vec(1) = 0.0
            vec(2) = 0.0
#if NDIM>1
            vec(3) = 0.0
#if NDIM>2
            vec(4) = 0.0
#endif
#endif
          endif



          !Contributions from sinks (UV bands) or reprocessed higher frequency bands (for IR)
#ifdef JSTR_VAR
          vec(1) = vec(1) + solnData(JSTR_VAR,i,j,k)*dt
#endif
          !REIR is already multiplied by time, so no dt required here
#ifdef REIR_VAR
          if(current_band .eq. 'IR') vec(1) = vec(1) + solnData(REIR_VAR,i,j,k)
#endif
          !Set VecIDs
          vecID(1) = coeff_ID(1,3)
          vecID(2) = vecID(1)+1
#if NDIM>1
          vecID(3) = vecID(2)+1
#if NDIM>2
          vecID(4) = vecID(3)+1
#endif
#endif

          ! Some safety checks
          if(solnData(ERAO_VAR,i,j,k) .gt. huge_val .or. solnData(ERAO_VAR,i,j,k) .ne. solnData(ERAO_VAR,i,j,k)) then
            call Grid_getSingleCellCoords((/i,j,k/),blockID,CENTER,EXTERIOR,coords)
            print *,'ERAO infinite or nan. ERAO_VAR = ' ,solnData(ERAO_VAR,i,j,k)
            print *, 'My coordinates are ', coords
            call Driver_abortFlash("ERAO is nan or infinite")
          endif

          if(solnData(MHXO_VAR,i,j,k) .gt. huge_val .or. solnData(MHXO_VAR,i,j,k) .ne. solnData(MHXO_VAR,i,j,k)) then
            call Grid_getSingleCellCoords((/i,j,k/),blockID,CENTER,EXTERIOR,coords)
            print *,'MHXO infinite or nan. MHXO_VAR = ' ,solnData(MHXO_VAR,i,j,k)
            print *, 'My coordinates are ', coords
            call Driver_abortFlash("MHXO is nan or infinite")
          endif

          if(solnData(MHYO_VAR,i,j,k) .gt. huge_val .or. solnData(MHYO_VAR,i,j,k) .ne. solnData(MHYO_VAR,i,j,k)) then
            call Grid_getSingleCellCoords((/i,j,k/),blockID,CENTER,EXTERIOR,coords)
            print *,'MHYO infinite or nan. MHYO_VAR = ' ,solnData(MHYO_VAR,i,j,k)
            print *, 'My coordinates are ', coords
            call Driver_abortFlash("MHYO is nan or infinite")
          endif

          if(solnData(MHZO_VAR,i,j,k) .gt. huge_val .or. solnData(MHZO_VAR,i,j,k) .ne. solnData(MHZO_VAR,i,j,k)) then
            call Grid_getSingleCellCoords((/i,j,k/),blockID,CENTER,EXTERIOR,coords)
            print *,'MHZO infinite or nan. MHZO_VAR = ' ,solnData(MHZO_VAR,i,j,k)
            print *, 'My coordinates are ', coords
            call Driver_abortFlash("MHZO is nan or infinite")
          endif
          ! Safety checks complete



           !!!!! Boundary Conditions

           !Check for possible boundaries
          if(i .eq. blkLimits(LOW, IAXIS)) then
            if(boundary_inblock(LOW,IAXIS)) then
              call SetPhysicalBoundary(LOW,IAXIS)
            else
              call SetBlockBoundary(LOW,IAXIS)
            endif
          endif

          if(i .eq. blkLimits(HIGH, IAXIS)) then
            if(boundary_inblock(HIGH,IAXIS)) then
              call SetPhysicalBoundary(HIGH,IAXIS)
            else
              call SetBlockBoundary(HIGH,IAXIS)
            endif
          endif
#if NDIM>1
          if(j .eq. blkLimits(LOW, JAXIS)) then
            if(boundary_inblock(LOW,JAXIS)) then
              call SetPhysicalBoundary(LOW,JAXIS)
            else
              call SetBlockBoundary(LOW,JAXIS)
            endif
          endif

          if(j .eq. blkLimits(HIGH, JAXIS)) then
            if(boundary_inblock(HIGH,JAXIS)) then
              call SetPhysicalBoundary(HIGH,JAXIS)
            else
              call SetBlockBoundary(HIGH,JAXIS)
            endif
          endif

#if NDIM>2
          if(k .eq. blkLimits(LOW, KAXIS)) then
            if(boundary_inblock(LOW,KAXIS)) then
              call SetPhysicalBoundary(LOW,KAXIS)
            else
              call SetBlockBoundary(LOW,KAXIS)
            endif
          endif

          if(k .eq. blkLimits(HIGH, KAXIS)) then
            if(boundary_inblock(HIGH,KAXIS)) then
              call SetPhysicalBoundary(HIGH,KAXIS)
            else
              call SetBlockBoundary(HIGH,KAXIS)
            endif
          endif
#endif
#endif

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Consider the case where 1D problem initiated by setting NxB and NyB =1
          ! Force Petsc to ignore the contributions in the boundaries then, and add it
          ! to the boundary point.
          ! For this case set the contributions of boundary cells to j=k=1
#if NDIM>1

          if(Ncell(1) .eq. 1) then
            coeff_ID(1,1) = -1
            coeff_ID(1,2) = -1
            coeff_ID(1,5) = -1
            coeff_ID(1,6) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,1) + coeff(1,5)
            coeff(1,4) = coeff(1,4) + coeff(1,2) + coeff(1,6)

            coeff_ID(2:NDIM+1,1) = -1
            coeff_ID(2:NDIM+1,2) = -1
            coeff_ID(2:NDIM+1,5) = -1
            coeff_ID(2:NDIM+1,6) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,1) + coeff(2:NDIM+1,5)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,2) + coeff(2:NDIM+1,6) 
          endif


          if(Ncell(2) .eq. 1) then
            coeff_ID(1,7) = -1
            coeff_ID(1,8) = -1
            coeff_ID(1,10) = -1
            coeff_ID(1,11) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,7) + coeff(1,10)
            coeff(1,9) = coeff(1,9) + coeff(1,8) + coeff(1,11)

            coeff_ID(2:NDIM+1,7) = -1
            coeff_ID(2:NDIM+1,8) = -1
            coeff_ID(2:NDIM+1,9) = -1
            coeff_ID(2:NDIM+1,10) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,7) + coeff(2:NDIM+1,9)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,8) + coeff(2:NDIM+1,10) 
          endif
#if NDIM>2
          if(Ncell(3) .eq. 1) then 
            coeff_ID(1,12) = -1
            coeff_ID(1,13) = -1
            coeff_ID(1,15) = -1
            coeff_ID(1,16) = -1
            coeff(1,3) = coeff(1,3)+ coeff(1,12) + coeff(1,15)
            coeff(1,14) = coeff(1,14) + coeff(1,13) + coeff(1,16)

            coeff_ID(2:NDIM+1,11) = -1
            coeff_ID(2:NDIM+1,12) = -1
            coeff_ID(2:NDIM+1,13) = -1
            coeff_ID(2:NDIM+1,14) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,11) + coeff(2:NDIM+1,13)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,12) + coeff(2:NDIM+1,14) 
          endif
#endif
#endif
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !Put values in Vector and Matrix
          !TODO: Possibly can optimise filling the matrix faster?
          call VecSetValues(YVec,NDIM+1,vecID,vec,INSERT_VALUES,ierr)
          IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Vector.")
          no_rows = 1
          call MatSetValues(Amat,no_rows,vecID(1),no_columns,coeff_ID(1,:),coeff(1,:),ADD_VALUES,ierr)
          IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix.")
          call MatSetValues(Amat,no_rows,vecID(2),no_columns,coeff_ID(2,:),coeff(2,:),ADD_VALUES,ierr)
          IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix.")
#if NDIM>1
          call MatSetValues(Amat,no_rows,vecID(3),no_columns,coeff_ID(3,:),coeff(3,:),ADD_VALUES,ierr)
          IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix.")
#if NDIM>2          
          call MatSetValues(Amat,no_rows,vecID(4),no_columns,coeff_ID(4,:),coeff(4,:),ADD_VALUES,ierr)
          IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix.")
#endif
#endif
        end do
      end do
    end do
  call Grid_releaseBlkPtr(blockID,solnData)
  end do




CONTAINS 
  !Computes the L/R wavespeeds for the local cell
  SUBROUTINE Set_Wavespeeds()
    !Reset FCInterface
    FCInterface(LOW:HIGH,IAXIS:NDIM) = .false.
    !Check for fine-coarse interfaces
    if(i .eq. blkLimits(LOW, IAXIS) .and. .not. boundary_inblock(LOW,IAXIS)) call CheckFCInterface(LOW,IAXIS,FCInterface(LOW,IAXIS))
    if(i .eq. blkLimits(HIGH, IAXIS) .and. .not. boundary_inblock(HIGH,IAXIS)) call CheckFCInterface(HIGH,IAXIS,FCInterface(HIGH,IAXIS))
#if NDIM>1
    if(j .eq. blkLimits(LOW,JAXIS) .and. .not. boundary_inblock(LOW,JAXIS)) call CheckFCInterface(LOW,JAXIS,FCInterface(LOW,JAXIS))

    if(j .eq. blkLimits(HIGH,JAXIS) .and. .not. boundary_inblock(HIGH,JAXIS)) call CheckFCInterface(HIGH,JAXIS,FCInterface(HIGH,JAXIS))
    
#if NDIM>2
    if(k .eq. blkLimits(LOW,KAXIS) .and. .not. boundary_inblock(LOW,KAXIS)) call CheckFCInterface(LOW,KAXIS,FCInterface(LOW,KAXIS))

    if(k .eq. blkLimits(HIGH,KAXIS) .and. .not. boundary_inblock(HIGH,KAXIS)) call CheckFCInterface(HIGH,KAXIS,FCInterface(HIGH,KAXIS))
    
#endif
#endif

    wavespeed(LOW:HIGH,IAXIS) = SQRT(solnData(XXED_VAR,i,j,k))
    wavespeed_negh(LOW:HIGH,IAXIS) = SQRT(xxed(LOW:HIGH,IAXIS))
#if NDIM>1
    wavespeed(LOW:HIGH,JAXIS) = SQRT(solnData(YYED_VAR,i,j,k))
    wavespeed_negh(LOW:HIGH,JAXIS) = SQRT(yyed(LOW:HIGH,JAXIS))

#if NDIM>2
    wavespeed(LOW:HIGH,KAXIS) = SQRT(solnData(ZZED_VAR,i,j,k))
    wavespeed_negh(LOW:HIGH,KAXIS) = SQRT(zzed(LOW:HIGH,KAXIS))
#endif
#endif

    ! Number of cells to resolve the propogation speed on (see Jiang et al 2014)
    wavespeed_cells = 10
    !Only IR band requires this fix
    if(use_hlle_correction .and. current_band .eq. 'IR') then 
      tau_hlle = (wavespeed_cells*cellsize(IAXIS)*solnData(TAUR_VAR,i,j,k))**2/(2*solnData(XXED_VAR,i,j,k))
      !-X interface
      call CorrectedWaveSpeed(tau_hlle_neigh(LOW,IAXIS),tau_hlle,FCInterface(LOW,IAXIS),&
        wavespeed_negh(LOW,IAXIS),wavespeed(LOW,IAXIS)) 
      !+X interface
      call CorrectedWaveSpeed(tau_hlle,tau_hlle_neigh(HIGH,IAXIS),FCInterface(HIGH,IAXIS),&
        wavespeed(HIGH,IAXIS),wavespeed_negh(HIGH,IAXIS))
#if NDIM>1
      tau_hlle = (wavespeed_cells*cellsize(JAXIS)*solnData(TAUR_VAR,i,j,k))**2/(2*solnData(YYED_VAR,i,j,k))
      !-Y interface
      call CorrectedWaveSpeed(tau_hlle_neigh(LOW,JAXIS),tau_hlle,FCInterface(LOW,JAXIS),&
        wavespeed_negh(LOW,JAXIS),wavespeed(LOW,JAXIS)) 
      !+Y interface
      call CorrectedWaveSpeed(tau_hlle,tau_hlle_neigh(HIGH,JAXIS),FCInterface(HIGH,JAXIS),&
        wavespeed(HIGH,JAXIS),wavespeed_negh(HIGH,JAXIS))
#if NDIM>2
      tau_hlle = (wavespeed_cells*cellsize(KAXIS)*solnData(TAUR_VAR,i,j,k))**2/(2*solnData(ZZED_VAR,i,j,k))
      !-Z interface
      call CorrectedWaveSpeed(tau_hlle_neigh(LOW,KAXIS),tau_hlle,FCInterface(LOW,KAXIS),&
        wavespeed_negh(LOW,KAXIS),wavespeed(LOW,KAXIS))
       !+Z interface
      call CorrectedWaveSpeed(tau_hlle,tau_hlle_neigh(HIGH,KAXIS),FCInterface(HIGH,KAXIS),&
        wavespeed(HIGH,KAXIS),wavespeed_negh(HIGH,KAXIS))
#endif
#endif
    endif
    
    !Compute some factors
    if(solnData(XXED_VAR,i,j,k) .gt. 0.0 .or. xxed(LOW,IAXIS) .gt. 0.0) then 
      Ci0 = (wavespeed(LOW,IAXIS) - wavespeed_negh(LOW,IAXIS))/ &
            (wavespeed(LOW,IAXIS) + wavespeed_negh(LOW,IAXIS))
    else
      Ci0 = 0.0
    endif

    if(solnData(XXED_VAR,i,j,k) .gt. 0.0 .or. xxed(HIGH,IAXIS) .gt. 0.0) then
      Ci1 = (wavespeed_negh(HIGH,IAXIS) - wavespeed(HIGH,IAXIS))/ &
            (wavespeed_negh(HIGH,IAXIS) + wavespeed(HIGH,IAXIS))
    else 
      Ci1 = 0.0
    endif

#if NDIM>1

    if(solnData(YYED_VAR,i,j,k) .gt. 0.0 .or. yyed(LOW,JAXIS) .gt. 0.0) then
      Cj0 = (wavespeed(LOW,JAXIS) - wavespeed_negh(LOW,JAXIS))/ &
            (wavespeed(LOW,JAXIS) + wavespeed_negh(LOW,JAXIS))
    else 
      Cj0 = 0.0
    endif

    if(solnData(YYED_VAR,i,j,k) .gt. 0.0 .or. yyed(HIGH,JAXIS) .gt. 0.0) then
      Cj1 = (wavespeed_negh(HIGH,JAXIS) - wavespeed(HIGH,JAXIS))/ &
            (wavespeed_negh(HIGH,JAXIS) + wavespeed(HIGH,JAXIS))
    else 
      Cj1 = 0.0
    endif

#if NDIM>2

    if(solnData(ZZED_VAR,i,j,k) .gt. 0.0 .or. zzed(LOW,KAXIS) .gt. 0.0) then
      Ck0 = (wavespeed(LOW,KAXIS) - wavespeed_negh(LOW,KAXIS))/ &
            (wavespeed(LOW,KAXIS) + wavespeed_negh(LOW,KAXIS))
    else 
      Ck0 = 0.0
    endif

    if(solnData(ZZED_VAR,i,j,k) .gt. 0.0 .or. zzed(HIGH,KAXIS) .gt. 0.0) then
       Ck1 = (wavespeed_negh(HIGH,KAXIS) - wavespeed(HIGH,KAXIS))/ &
            (wavespeed_negh(HIGH,KAXIS) + wavespeed(HIGH,KAXIS))
    else 
      Ck1 = 0.0
    endif

#endif
#endif

    ! Some safety checks for infinites or nans in wavespeeds
    huge_val = huge(Ci0)
    if(Ci0 .gt. huge_val .or. Ci1 .gt. huge_val & 
#if NDIM>1
      .or. Cj0 .gt. huge_val .or. Cj1 .gt. huge_val &
#if NDIM>2
      .or. Ck0 .gt. huge_val .or. Ck1 .gt. huge_val &
#endif
#endif
      ) then
      print *, 'Coefficients', Ci0, Ci1, Cj0, Cj1, Ck0, Ck1
      print *, 'Local R/L interface wavespeeds'
      do axis =1, NDIM
        do cell_offset=LOW, HIGH
          print *, wavespeed(cell_offset,axis)
        end do 
      end do 
      print *, 'wavespeeds of neighbours'
      do axis =1, NDIM
        do cell_offset=LOW, HIGH
          print *, wavespeed_negh(cell_offset,axis)
        end do 
      end do 
      call Driver_abortFlash("Wavespeed differences infinite. This might mean eddington tensors zero.")
    endif

    if(Ci0 .ne. Ci0 .or. Ci1 .ne. Ci1 &
#if NDIM>1
      .or. Cj0 .ne. Cj0 .or. Cj1 .ne. Cj1 &
#if NDIM>2
      .or. Ck0 .ne. Ck0 .or. Ck1 .ne. Ck1 &
#endif
#endif
      ) then
      print *, 'Coefficients', Ci0, Ci1, Cj0, Cj1, Ck0, Ck1
      if(use_hlle_correction) then 
        print *, 'tau factors',tau_hlle,tau_hlle_neigh, solnData(TAUR_VAR,i-1,j,k),solnData(TAUR_VAR,i+1,j,k)
      endif
      print *, 'Local R/L interface wavespeeds'
      do axis =1, NDIM
        do cell_offset=LOW, HIGH
          print *, wavespeed(cell_offset,axis)
        end do 
      end do 
      print *, 'wavespeeds of neighbours'
      do axis =1, NDIM
        do cell_offset=LOW, HIGH
          print *, wavespeed_negh(cell_offset,axis)
        end do 
      end do 
      call Driver_abortFlash("Wavespeed differences nans. This might mean eddington tensors/mean intensity is nan.")
    endif
    ! Safety checks complete

  END SUBROUTINE Set_Wavespeeds

  ! Sets the coefficients in the matrix for a given cell
  SUBROUTINE Set_Coeffs()
#ifdef VELX_VAR
    vx = solnData(VELX_VAR,i,j,k)
#endif
#ifdef VELY_VAR
    vy = solnData(VELY_VAR,i,j,k)
#endif
#ifdef VELZ_VAR
    vz = solnData(VELZ_VAR,i,j,k)
#endif

    !Matrix coefficients and corresponding global index for non-boundary pts

    !er_coeff order : Er,Fr1 all i, Er,Fr2 all j, Er, Fr3 all k - 5*NDIM+1 components
!================================================================================================
    !Er equation

    !Time discretisation term (d/dt). Only present for IR band
    if(current_band .eq. 'IR') then 
      er_coeff(3) = 1
    else
      er_coeff(3) = 0
    endif

    !Transport terms
    er_coeff(1) = -1 * d1x * (1+Ci0) * rt_speedlt * wavespeed_negh(LOW,IAXIS)
    er_coeff(2) = -1 * d1x * (1+Ci0)
    er_coeff(3) = er_coeff(3) + rt_speedlt * (d1x*(1+Ci1)*wavespeed(HIGH,IAXIS) &
                  + d1x*(1-Ci0)*wavespeed(LOW,IAXIS))
    er_coeff(4) = d1x*(1+Ci1) - d1x*(1-Ci0)
    er_coeff(5) = -1 * d1x * (1-Ci1) * rt_speedlt * wavespeed_negh(HIGH,IAXIS)
    er_coeff(6) = d1x * (1-Ci1)
#if NDIM>1
    er_coeff(7) = -1 * d1y * (1+Cj0) * rt_speedlt * wavespeed_negh(LOW,JAXIS)
    er_coeff(8) = -1 * d1y * (1+Cj0)
    er_coeff(3) = er_coeff(3) + rt_speedlt * (d1y*(1+Cj1)*wavespeed(HIGH,JAXIS) &
                  + d1y*(1-Cj0)*wavespeed(LOW,JAXIS))
    er_coeff(9) = d1y*(1+Cj1) - d1y*(1-Cj0)
    er_coeff(10) = -1 * d1y * (1-Cj1) * rt_speedlt * wavespeed_negh(HIGH,JAXIS)
    er_coeff(11) = d1y * (1-Cj1)
#if NDIM>2
    er_coeff(12) = -1 * d1z * (1+Ck0) * rt_speedlt * wavespeed_negh(LOW,KAXIS)
    er_coeff(13) = -1 * d1z * (1+Ck0)
    er_coeff(3) = er_coeff(3) + rt_speedlt * (d1z * (1+Ck1)*wavespeed(HIGH,KAXIS) + & 
                    d1z*(1-Ck0)*wavespeed(LOW,KAXIS))
    er_coeff(14) = d1z*(1+Ck1) - d1z*(1-Ck0)
    er_coeff(15) = -1 * d1z * (1-Ck1) * rt_speedlt * wavespeed_negh(HIGH,KAXIS)
    er_coeff(16) = d1z*(1-Ck1)
#endif
#endif

    !Source terms
    ! Radiation absorption term. Include for UV bands. Also include for IR band if radiative equilibrium is not assumed.
    if(rt_update_hydro .and. ((current_band .eq. 'IR' .and. rt_hydro_type .ne. 1) .or. &
      & current_band .eq. 'FUV' .or. current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
      & current_band .eq. 'EUV_15P2_INFTY' .or. current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW' .or. &
      current_band .eq. 'PE')) then
      er_coeff(3) = er_coeff(3) + d2*opac_planck
    endif
    if(rt_ovcterms) then
    !Advection of radiation enthalpy
      er_coeff(3) = er_coeff(3) + d2*(opac_planck - opac_rosseland)*vx*(vx*(1+solnData(XXED_VAR,i,j,k)))/(rt_speedlt**2)
#if NDIM>1
      er_coeff(3) = er_coeff(3) + d2*(opac_planck - opac_rosseland)* &
                    vy*(vy*(1+solnData(YYED_VAR,i,j,k)) + vx*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2) &
                    + d2*(opac_planck - opac_rosseland)*vx*(vy*solnData(XYED_VAR,i,j,k))/(rt_speedlt**2)
#if NDIM>2
      er_coeff(3) = er_coeff(3) + d2*(opac_planck - opac_rosseland)* &
                      vz*(vz*(1+solnData(ZZED_VAR,i,j,k)) + vx*solnData(XZED_VAR,i,j,k) + &
                      vy* solnData(YZED_VAR,i,j,k))/(rt_speedlt**2) + d2*(opac_planck - opac_rosseland)* &
                      vx*(vz*solnData(XZED_VAR,i,j,k))/(rt_speedlt**2) + d2*(opac_planck - opac_rosseland)* &
                      vy*(vz*solnData(YZED_VAR,i,j,k))/(rt_speedlt**2)
#endif
#endif
    !Work done by radiation on gas
      er_coeff(4) = er_coeff(4) - d2*(2*opac_planck - opac_rosseland)*vx/(rt_speedlt**2)
#if NDIM>1
      er_coeff(9) = er_coeff(9) - d2*(2*opac_planck - opac_rosseland)*vy/(rt_speedlt**2)
#if NDIM>2
      er_coeff(14) = er_coeff(14) - d2*(2*opac_planck - opac_rosseland)*vz/(rt_speedlt**2)
#endif
#endif
    endif

!================================================================================================
    !Frx
    !Time discretisation term (d/dt). Only present for IR band
    if(current_band .eq. 'IR') then 
      frx_coeff(4) = 1
    else
      frx_coeff(4) = 0
    endif
    !Transport terms
    frx_coeff(1) = -d1x*(1+Ci0)*(rt_speedlt**2)*xxed(LOW,IAXIS)
    frx_coeff(2) = -d1x*(1+Ci0)*rt_speedlt*wavespeed_negh(LOW,IAXIS)
    frx_coeff(3) = rt_speedlt**2 * (d1x*(Ci1+Ci0)*solnData(XXED_VAR,i,j,k))
    frx_coeff(4) = frx_coeff(4) + rt_speedlt * d1x*((1+Ci1)*wavespeed(HIGH,IAXIS)+(1-Ci0)*wavespeed(LOW,IAXIS))
    frx_coeff(5) = d1x*(1-Ci1)*(rt_speedlt**2)*xxed(HIGH,IAXIS)
    frx_coeff(6) = -1*d1x*(1-Ci1)*rt_speedlt*wavespeed_negh(HIGH,IAXIS)

#if NDIM>1
    frx_coeff(7) = -d1y*(1+Cj0)*(rt_speedlt**2)*xyed(LOW,JAXIS)
    frx_coeff(8) = -d1y*(1+Cj0)*rt_speedlt * wavespeed_negh(LOW,JAXIS)
    frx_coeff(3) = frx_coeff(3)+ rt_speedlt**2 * (d1y*(Cj1+Cj0)*solnData(XYED_VAR,i,j,k))
    frx_coeff(4) = frx_coeff(4)+ rt_speedlt*d1y*((1+Cj1)*wavespeed(HIGH,JAXIS)+(1-Cj0)*wavespeed(LOW,JAXIS))
    frx_coeff(9) = d1y*(1-Cj1)*(rt_speedlt**2)*xyed(HIGH,JAXIS)
    frx_coeff(10) = -1*d1y*(1-Cj1)*rt_speedlt*wavespeed_negh(HIGH,JAXIS)
    frx_coeff(11) = 0
#if NDIM>2
    frx_coeff(11) = -1*d1z*(1+Ck0)*(rt_speedlt**2)*xzed(LOW,KAXIS)
    frx_coeff(12) = -1*d1z*(1+Ck0)*rt_speedlt*wavespeed_negh(LOW,KAXIS)
    frx_coeff(3) = frx_coeff(3)+ rt_speedlt**2 * (d1z*(Ck1+Ck0)*solnData(XZED_VAR,i,j,k))
    frx_coeff(4) = frx_coeff(4)+ rt_speedlt*d1z*((1+Ck1)*wavespeed(HIGH,KAXIS)+(1-Ck0)*wavespeed(LOW,KAXIS))
    frx_coeff(13) = d1z*(1-Ck1)*(rt_speedlt**2)*xzed(HIGH,KAXIS)
    frx_coeff(14) = -1*d1z*(1-Ck1)*rt_speedlt*wavespeed_negh(HIGH,KAXIS)
#endif
#endif

    !Source terms
    !Gas-Radiation momentum exchange term
    frx_coeff(4) = frx_coeff(4) + d2*opac_rosseland
    if(rt_ovcterms) then
      !Radiation enthalpy related term
      frx_coeff(3) = frx_coeff(3) - d2*opac_rosseland *vx*(1+solnData(XXED_VAR,i,j,k))
#if NDIM>1
      frx_coeff(3) = frx_coeff(3) - d2*opac_rosseland*vy*solnData(XYED_VAR,i,j,k)
#if NDIM>2
      frx_coeff(3) = frx_coeff(3) - d2*opac_rosseland*vz*solnData(XZED_VAR,i,j,k)
#endif
#endif
    endif
!================================================================================================
#if NDIM>1
    !Fry equation
    !Time discretisation term (d/dt). Only present for IR band
    if(current_band .eq. 'IR') then 
      fry_coeff(4) = 1
    else
      fry_coeff(4) = 0
    endif
    !Transport Terms
    fry_coeff(1) = -d1x*(1+Ci0)*(rt_speedlt**2)*xyed(LOW,IAXIS)
    fry_coeff(2) = -d1x*(1+Ci0)*rt_speedlt*wavespeed_negh(LOW,IAXIS)
    fry_coeff(3) = rt_speedlt**2 * (d1x*(Ci1+Ci0)*solnData(XYED_VAR,i,j,k) + & 
                   d1y*(Cj1+Cj0)*solnData(YYED_VAR,i,j,k))
    fry_coeff(4) = fry_coeff(4) + rt_speedlt * d1x*((1+Ci1)*wavespeed(HIGH,IAXIS)+(1-Ci0)*wavespeed(LOW,IAXIS)) + &
                   rt_speedlt * d1y * ((1+Cj1)*wavespeed(HIGH,JAXIS)+(1-Cj0)*wavespeed(LOW,JAXIS))
    fry_coeff(5) = d1x*(1-Ci1)*(rt_speedlt**2)*xyed(HIGH,IAXIS)
    fry_coeff(6) = -1*d1x*(1-Ci1)*rt_speedlt*wavespeed_negh(HIGH,IAXIS)
    fry_coeff(7) = -d1y*(1+Cj0)*(rt_speedlt**2)*yyed(LOW,JAXIS)
    fry_coeff(8) = -d1y*(1+Cj0)*rt_speedlt * wavespeed_negh(LOW,JAXIS)
    fry_coeff(9) = d1y*(1-Cj1)*(rt_speedlt**2)*yyed(HIGH,JAXIS)
    fry_coeff(10) = -1*d1y*(1-Cj1)*rt_speedlt*wavespeed_negh(HIGH,JAXIS)
    fry_coeff(11) = 0
#if NDIM>2
    fry_coeff(11) = -1*d1z*(1+Ck0)*(rt_speedlt**2)*yzed(LOW,KAXIS)
    fry_coeff(12) = -1*d1z*(1+Ck0)*rt_speedlt*wavespeed_negh(LOW,KAXIS)
    fry_coeff(3) = fry_coeff(3) + rt_speedlt**2 * (d1z*(Ck1+Ck0)*solnData(YZED_VAR,i,j,k))
    fry_coeff(4) = fry_coeff(4) + rt_speedlt*d1z*((1+Ck1)*wavespeed(HIGH,KAXIS)+(1-Ck0)*wavespeed(LOW,KAXIS))
    fry_coeff(13) = d1z*(1-Ck1)*(rt_speedlt**2)*yzed(HIGH,KAXIS)
    fry_coeff(14) = -1*d1z*(1-Ck1)*rt_speedlt*wavespeed_negh(HIGH,KAXIS)
#endif
!if NDIM>2 end

    !Source Terms
    !Gas-Radiation momentum exchange term
    fry_coeff(4) = fry_coeff(4) + d2*opac_rosseland
    if(rt_ovcterms) then
    !Radiation enthalpy related term
      fry_coeff(3) = fry_coeff(3) - d2*opac_rosseland*vx*solnData(XYED_VAR,i,j,k)
      fry_coeff(3) = fry_coeff(3) - d2*opac_rosseland *vy*(1+solnData(YYED_VAR,i,j,k))
#if NDIM>2
      fry_coeff(3) = fry_coeff(3) - d2*opac_rosseland*vz*solnData(YZED_VAR,i,j,k)
#endif
    endif
!if NDIM>2 end
#endif
!if NDIM>1 end
!================================================================================================
#if NDIM>2
    !Frz equation
    !Time discretisation term (d/dt). Only present for IR band
    if(current_band .eq. 'IR') then 
      frz_coeff(4) = 1
    else
      frz_coeff(4) = 0
    endif
    !Transport terms
    frz_coeff(1) = -d1x*(1+Ci0)*(rt_speedlt**2)*xzed(LOW,IAXIS)
    frz_coeff(2) = -d1x*(1+Ci0)*rt_speedlt*wavespeed_negh(LOW,IAXIS)
    frz_coeff(3) = rt_speedlt**2 * (d1x*(Ci1+Ci0)*solnData(XZED_VAR,i,j,k) + & 
                   d1y*(Cj1+Cj0)*solnData(YZED_VAR,i,j,k)+d1z*(Ck1+Ck0)*solnData(ZZED_VAR,i,j,k))
    frz_coeff(4) = frz_coeff(4) + rt_speedlt * d1x*((1+Ci1)*wavespeed(HIGH,IAXIS)+(1-Ci0)*wavespeed(LOW,IAXIS)) + &
                   rt_speedlt * d1y*((1+Cj1)*wavespeed(HIGH,JAXIS)+(1-Cj0)*wavespeed(LOW,JAXIS)) + &
                   rt_speedlt *d1z*((1+Ck1)*wavespeed(HIGH,KAXIS)+(1-Ck0)*wavespeed(LOW,KAXIS))
    frz_coeff(5) = d1x*(1-Ci1)*(rt_speedlt**2)*xzed(HIGH,IAXIS)
    frz_coeff(6) = -1*d1x*(1-Ci1)*rt_speedlt*wavespeed_negh(HIGH,IAXIS)
    frz_coeff(7) = -d1y*(1+Cj0)*(rt_speedlt**2)*yzed(LOW,JAXIS)
    frz_coeff(8) = -d1y*(1+Cj0)*rt_speedlt * wavespeed_negh(LOW,JAXIS)
    frz_coeff(9) = d1y*(1-Cj1)*(rt_speedlt**2)*yzed(HIGH,JAXIS)
    frz_coeff(10) = -1*d1y*(1-Cj1)*rt_speedlt*wavespeed_negh(HIGH,JAXIS)
    frz_coeff(11) = -1*d1z*(1+Ck0)*(rt_speedlt**2)*zzed(LOW,KAXIS)
    frz_coeff(12) = -1*d1z*(1+Ck0)*rt_speedlt*wavespeed_negh(LOW,KAXIS)
    frz_coeff(13) = d1z*(1-Ck1)*(rt_speedlt**2)*zzed(HIGH,KAXIS)
    frz_coeff(14) = -1*d1z*(1-Ck1)*rt_speedlt*wavespeed_negh(HIGH,KAXIS)

    frx_coeff(15) = 0
    fry_coeff(15) = 0
    frz_coeff(15) = 0
    frx_coeff(16) = 0
    fry_coeff(16) = 0
    frz_coeff(16) = 0

    !Source terms
    !Gas-Radiation momentum exchange term
    frz_coeff(4) = frz_coeff(4) + d2*opac_rosseland
    if(rt_ovcterms) then
    !Radiation enthalpy related term
      frz_coeff(3) = frz_coeff(3) - d2*opac_rosseland * (vz*(1+solnData(ZZED_VAR,i,j,k)) &
                     + vx*solnData(XZED_VAR,i,j,k) + vy*solnData(YZED_VAR,i,j,k))
    endif
#endif
!if NDIM>2 end
!================================================================================================

  END SUBROUTINE Set_Coeffs
  
  SUBROUTINE getrowcolumnIDs(er_columnID,frx_columnID,fry_columnID,frz_columnID)
    ! implicit none
    ! integer, intent(in) :: blockID, i ,j ,k 
    ! type(PetscInt), dimension(NDIM+1,5*NDIM+1), intent(out) :: coeff_ID
    integer, dimension(5*NDIM+1), intent(out) :: er_columnID,frx_columnID,fry_columnID,frz_columnID
    integer :: i_b, j_b, k_b
    
    ! call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    i_b = i-blkLimits(LOW, IAXIS)
    j_b = j-blkLimits(LOW, JAXIS)
    k_b = k-blkLimits(LOW, KAXIS)

    
    er_columnID(1) = globalmatrixindex(i_b-1,j_b,k_b,blockID,myPE) 
    er_columnID(2) = er_columnID(1) + 1
    er_columnID(3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
    er_columnID(4) = er_columnID(3) + 1
    er_columnID(5) = globalmatrixindex(i_b+1,j_b,k_b,blockID,myPE)
    er_columnID(6) = er_columnID(5) + 1

    frx_columnID(1:6) = er_columnID(1:6) 

#if NDIM>1
    er_columnID(7) = globalmatrixindex(i_b,j_b-1,k_b,blockID,myPE)
    er_columnID(8) = er_columnID(7)+ 2
    er_columnID(9) = er_columnID(3) + 2
    er_columnID(10) = globalmatrixindex(i_b,j_b+1,k_b,blockID,myPE)
    er_columnID(11) = er_columnID(10) + 2   

    frx_columnID(7) =  er_columnID(7)
    frx_columnID(8) =  er_columnID(7)+1
    frx_columnID(9) =  er_columnID(10)
    frx_columnID(10) =  er_columnID(10)+1

    fry_columnID(1) = er_columnID(1)
    fry_columnID(2) = er_columnID(1)+2
    fry_columnID(3) = er_columnID(3)
    fry_columnID(4) = er_columnID(3)+2
    fry_columnID(5) = er_columnID(5)
    fry_columnID(6) = er_columnID(5)+2
    fry_columnID(7) = er_columnID(7)
    fry_columnID(8) = er_columnID(7)+2
    fry_columnID(9) = er_columnID(10)
    fry_columnID(10) = er_columnID(10)+2

    !Set to negative to ignore this while filling matrix
    frx_columnID(11) = -1
    fry_columnID(11) = -1
#if NDIM>2
    er_columnID(12) = globalmatrixindex(i_b,j_b,k_b-1,blockID,myPE)
    er_columnID(13) = er_columnID(12) + 3
    er_columnID(14) = er_columnID(3) + 3
    er_columnID(15) = globalmatrixindex(i_b,j_b,k_b+1,blockID,myPE)
    er_columnID(16) = er_columnID(15) + 3

    frx_columnID(11) = er_columnID(12)
    frx_columnID(12) = er_columnID(12)+1
    frx_columnID(13) = er_columnID(15)
    frx_columnID(14) = er_columnID(15)+1

    fry_columnID(11) = er_columnID(12)
    fry_columnID(12) = er_columnID(12)+2  
    fry_columnID(13) = er_columnID(15)
    fry_columnID(14) = er_columnID(15)+2

    frz_columnID(1) = er_columnID(1)
    frz_columnID(2) = er_columnID(1)+3
    frz_columnID(3) = er_columnID(3)
    frz_columnID(4) = er_columnID(3) + 3
    frz_columnID(5) = er_columnID(5)
    frz_columnID(6) = er_columnID(5) + 3
    frz_columnID(7) = er_columnID(7)
    frz_columnID(8) = er_columnID(7) +3
    frz_columnID(9) = er_columnID(10)
    frz_columnID(10) = er_columnID(10) +3
    frz_columnID(11) = er_columnID(12)
    frz_columnID(12) = er_columnID(12)+3
    frz_columnID(13) = er_columnID(15)
    frz_columnID(14) = er_columnID(15) + 3

    frx_columnID(15) = -1
    fry_columnID(15) = -1
    frz_columnID(15) = -1

    frx_columnID(16) = -1
    fry_columnID(16) = -1
    frz_columnID(16) = -1
#endif
#endif
  END SUBROUTINE getrowcolumnIDs

  SUBROUTINE SetBlockBoundary(which_boundary,which_axis)
    use RadTrans_data
    use Grid_interface, ONLY: Grid_getSingleCellCoords, Grid_getBlkIDFromPos
    INTEGER, INTENT(IN) :: which_axis, which_boundary
    type(PetscErrorCode) :: ierr
    INTEGER, EXTERNAL :: globalmatrixindex
    INTEGER :: i_b,j_b,k_b, baseindex
    real :: negh_values(NDIM+1,4), geometrical_factor, cell_length_factor
    type(PetscInt) :: negh_cID(NDIM+1,4),negh_rID(NDIM+1),negh_no_rows,negh_no_columns
    INTEGER, DIMENSION(IAXIS:KAXIS), parameter :: NX = (/NXB,NYB,NZB/)
    INTEGER, DIMENSION(MDIM) :: edge, negh_cell
    integer, dimension(IAXIS:KAXIS) :: child_pos
    INTEGER, DIMENSION(BLKNO:PROCNO) :: neghBlkProc
    INTEGER :: numNegh,n

    i_b = i-blkLimits(LOW, IAXIS)
    j_b = j-blkLimits(LOW, JAXIS)
    k_b = k-blkLimits(LOW, KAXIS)


    edge(:) = (/ 1+K1D, 1+K2D, 1+K3D /)
    negh_cell(:) = (/ i_b, j_b, k_b /)

    if(which_boundary .eq. LOW) then
      edge(which_axis) = LEFT_EDGE
      negh_cell(which_axis) = NX(which_axis) - 1
    else 
      edge(which_axis) = RIGHT_EDGE
      negh_cell(which_axis) = 0
    endif

    numNegh = SurrBlkSum(b) % regionInfo(edge(1),edge(2),edge(3)) % numNegh

    SELECT CASE(NeghLevels(b,edge(1),edge(2),edge(3)) - mylevel)

    CASE(-1) ! Neighbour lower refinement than current block (fine-coarse boundary)
      if(numNegh .ne. 1) then 
        call Driver_abortFlash("[VET] : ERROR, Neighbour lower refinement level - more than one neighbour found.")
      endif

      neghBlkProc(BLKNO:PROCNO) = SurrBlkSum(b) &
            % regionInfo(edge(1),edge(2),edge(3)) % details(BLKNO:PROCNO,numNegh)

#if NDIM > 1
        ! In case of multidimensional simulations account for a
        ! shift of the smaller/finer block to the larger/coarse block in the
        ! directions along the edge surface.
        call whichChild(blockID,child_pos)
        DO n=1,MDIM
          IF(n.NE. which_axis) THEN
            IF(child_pos(n).eq.HIGH) THEN
              negh_cell(n) = negh_cell(n)/2 + NX(n)/2
            ELSE ! LOW
              negh_cell(n) = negh_cell(n)/2
            END IF
          END IF
        END DO
#endif

      !Get flattened index of neighbouring active cell
      baseindex = globalmatrixindex(negh_cell(1),negh_cell(2), &
        negh_cell(3),neghBlkProc(BLKNO),neghBlkProc(PROCNO))

      !Set neighbour cell coefficients to this cell (as for same level case)
      ! i.e.: set coefficients in corresponding columns for the current cell row
      if(which_boundary .eq. LOW) then
        if(which_axis .eq. IAXIS) then

          coeff_ID(1,1) = baseindex
          coeff_ID(1,2) = coeff_ID(1,1) + 1
          coeff_ID(2,1) = coeff_ID(1,1)
          coeff_ID(2,2) = coeff_ID(1,2)

#if NDIM>1
          coeff_ID(3,1) = coeff_ID(1,1)
          coeff_ID(3,2) = coeff_ID(1,1) + 2
#if NDIM>2
          coeff_ID(4,1) = coeff_ID(1,1)
          coeff_ID(4,2) = coeff_ID(1,1) + 3
#endif
#endif
        endif 

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,7) = baseindex
          coeff_ID(1,8) = coeff_ID(1,7) + 2
          coeff_ID(2,7) = coeff_ID(1,7)
          coeff_ID(2,8) = coeff_ID(1,7) + 1

#if NDIM>1
          coeff_ID(3,7) = coeff_ID(1,7)
          coeff_ID(3,8) = coeff_ID(1,7) + 2
#if NDIM>2
          coeff_ID(4,7) = coeff_ID(1,7)
          coeff_ID(4,8) = coeff_ID(1,7) + 3
#endif
#endif
      
        endif 

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,12) = baseindex
          coeff_ID(1,13) = coeff_ID(1,12) + 3
          coeff_ID(2,11) = coeff_ID(1,12)
          coeff_ID(2,12) = coeff_ID(1,12) + 1

#if NDIM>1
          coeff_ID(3,11) = coeff_ID(1,12)
          coeff_ID(3,12) = coeff_ID(1,12) + 2
#if NDIM>2
          coeff_ID(4,11) = coeff_ID(1,12)
          coeff_ID(4,12) = coeff_ID(1,12) + 3
#endif
#endif
      
        endif

      endif

      if(which_boundary .eq. HIGH) then
        if(which_axis .eq. IAXIS) then
          
          coeff_ID(1,5) = baseindex
          coeff_ID(1,6) = coeff_ID(1,5) + 1
          coeff_ID(2,5) = coeff_ID(1,5)
          coeff_ID(2,6) = coeff_ID(1,6)

#if NDIM>1
          coeff_ID(3,5) = coeff_ID(1,5)
          coeff_ID(3,6) = coeff_ID(1,5) + 2
#if NDIM>2
          coeff_ID(4,5) = coeff_ID(1,5)
          coeff_ID(4,6) = coeff_ID(1,5) + 3
#endif
#endif

        endif

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,10) = baseindex
          coeff_ID(1,11) = coeff_ID(1,10) + 2
          coeff_ID(2,9) = coeff_ID(1,10)
          coeff_ID(2,10) = coeff_ID(1,10) + 1

#if NDIM>1
          coeff_ID(3,9) = coeff_ID(1,10)
          coeff_ID(3,10) = coeff_ID(1,10) + 2
#if NDIM>2
          coeff_ID(4,9) = coeff_ID(1,10)
          coeff_ID(4,10) = coeff_ID(1,10) + 3
#endif
#endif
      
        endif

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,15) = baseindex
          coeff_ID(1,16) = coeff_ID(1,15) + 3
          coeff_ID(2,13) = coeff_ID(1,15)
          coeff_ID(2,14) = coeff_ID(1,15) + 1

#if NDIM>1
          coeff_ID(3,13) = coeff_ID(1,15)
          coeff_ID(3,14) = coeff_ID(1,15) + 2
#if NDIM>2
          coeff_ID(4,13) = coeff_ID(1,15)
          coeff_ID(4,14) = coeff_ID(1,15) + 3
#endif
#endif
      
        endif

      endif


      !Now fill coefficients for coarse-fine interface here for simplicity
      !i.e. fill in for row baseindex (i.e. for the eqn concerning the neighbour cell) - 
      !       1. The coefficient concerned with the flux across the face with a -ve sign 
      !          on this cell column. This will include NDIM +1 terms to fill - E_r + fluxes
      !          i.e. (baseindex,currentcell_index) -  i.e. off-diagonal term of neighbour
      !       2. Add the corresponding term in the coefficient for the neighbour cell column, 
      !          i.e. (baseindex, baseindex)- i.e. diagonal term of neighbour

      ! negh_values contains the values to be filled in different columns in the neighbour cell row
      ! In order: coarse cell, Er/Fr of each neighbour cell

      ! Now we are filling coefficients of neighbours, for the neighbour cell

      if(which_boundary .eq. HIGH) then
        if(which_axis .eq. IAXIS) then
          ! this cell as negh column index
          !E_r eqn
          negh_values(1,1) = -1*er_coeff(5)
          negh_values(1,2) = -1.*er_coeff(6)
          negh_values(1,3) = -1.*d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = -1.*d1x*(1+Ci1)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(5)
          negh_values(2,2) = -1.*frx_coeff(6)
          negh_values(2,3) = -1.*d1x*(1+Ci1)*solnData(XXED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt


          ! Column IDs
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 1
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+1
          negh_cID(2,1:4) = negh_cID(1,1:4)

          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1
#if NDIM>1

          !Fry eqn
          negh_values(3,1) = -1*fry_coeff(5)
          negh_values(3,2) = -1*fry_coeff(6)
          negh_values(3,3) = -1.*d1x*(1+Ci1)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Column IDs
          negh_cID(3,1:4) = negh_cID(1,1:4)
          negh_cID(3,2) = negh_cID(3,1)+2
          negh_cID(3,4) = negh_cID(3,3)+2

          !Row ids
          negh_rID(3) = baseindex+2
          
#if NDIM>2

          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(5)
          negh_values(4,2) = -1*frz_coeff(6)
          negh_values(4,3) = -1.*d1x*(1+Ci1)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt


          !Column IDs
          negh_cID(4,1:4) = negh_cID(1,1:4)
          negh_cID(4,2) = negh_cID(4,1)+3
          negh_cID(4,4) = negh_cID(4,3)+3

          !Row ids
          negh_rID(4) = baseindex+3
        
#endif
#endif
        endif

        if(which_axis .eq. JAXIS) then

          !E_r eqn
          negh_values(1,1) = -1.*er_coeff(10)
          negh_values(1,2) = -1.*er_coeff(11)
          negh_values(1,3) = -1.*d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = -1.*d1y*(1+Cj1)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(9)
          negh_values(2,2) = -1.*frx_coeff(10)
          negh_values(2,3) = -1.*d1y*(1+Cj1)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Fry eqn 
          negh_values(3,1) = -1*fry_coeff(9)
          negh_values(3,2) = -1.*fry_coeff(10)
          negh_values(3,3) = -1.*d1y*(1+Cj1)*solnData(YYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          ! Column IDs
          !Er eqn
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 2
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+2

          !Frx eqn
          negh_cID(2,1:4) = negh_cID(1,1:4)
          negh_cID(2,2) = negh_cID(2,1)+1
          negh_cID(2,4) = negh_cID(2,3)+1

          !Fry eqn
          negh_cID(3,1:4) = negh_cID(1,1:4)
          
          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1
          negh_rID(3) = baseindex+2
#if NDIM>2
          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(9)
          negh_values(4,2) = -1.*frz_coeff(10)
          negh_values(4,3) = -1.*d1y*(1+Cj1)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt


          !Column Ids
          !Frz eqn
          negh_cID(4,1:4) = negh_cID(1,1:4)
          negh_cID(4,2) = negh_cID(4,1) + 3
          negh_cID(4,4) = negh_cID(4,3) + 3

          !Row ids       
          negh_rID(4) = baseindex+3

#endif
      
        endif

        if(which_axis .eq. KAXIS) then
          !E_r eqn
          negh_values(1,1) = -1.*er_coeff(15)
          negh_values(1,2) = -1.*er_coeff(16)
          negh_values(1,3) = -1.*d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = -1.*d1z*(1+Ck1)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(13)
          negh_values(2,2) = -1.*frx_coeff(14)
          negh_values(2,3) = -1.*d1z*(1+Ck1)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Fry eqn 
          negh_values(3,1) = -1*fry_coeff(13)
          negh_values(3,2) = -1.*fry_coeff(14)
          negh_values(3,3) = -1.*d1z*(1+Ck1)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(13)
          negh_values(4,2) = -1.*frz_coeff(14)
          negh_values(4,3) = -1.*d1z*(1+Ck1)*solnData(ZZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt


          !Column Ids
          !Er eqn
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 3
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+3

          !Frx eqn
          negh_cID(2,1:4) = negh_cID(1,1:4)
          negh_cID(2,2) = negh_cID(2,1)+1
          negh_cID(2,4) = negh_cID(2,3)+1

          !Fry eqn
          negh_cID(3,1:4) = negh_cID(1,1:4)
          negh_cID(3,2) = negh_cID(3,1) + 2
          negh_cID(3,4) = negh_cID(3,3) + 2

          !Frz eqn
          negh_cID(4,1:4) = negh_cID(1,1:4)
          
          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1
          negh_rID(3) = baseindex+2
          negh_rID(4) = baseindex+3

        endif 

      endif

      if(which_boundary .eq. LOW) then
        if(which_axis .eq. IAXIS) then
        !E_r eqn
          negh_values(1,1) = -1*er_coeff(1)
          negh_values(1,2) = -1*er_coeff(2)
          negh_values(1,3) = -1.*d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = d1x*(1-Ci0)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(1)
          negh_values(2,2) = -1*frx_coeff(2)
          negh_values(2,3) = d1x*(1-Ci0)*solnData(XXED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt


          ! Column IDs
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 1
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+1
          negh_cID(2,1:4) = negh_cID(1,1:4)

          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1


          
#if NDIM>1

          !Fry eqn
          negh_values(3,1) = -1*fry_coeff(1)
          negh_values(3,2) = -1*fry_coeff(2)
          negh_values(3,3) = d1x*(1-Ci0)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt


          !Column IDs
          negh_cID(3,1:4) = negh_cID(1,1:4)
          negh_cID(3,2) = negh_cID(3,1)+2
          negh_cID(3,4) = negh_cID(3,3)+2

          !Row ids
          negh_rID(3) = baseindex+2
          
#if NDIM>2

          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(1)
          negh_values(4,2) = -1*frz_coeff(2)
          negh_values(4,3) = d1x*(1-Ci0)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt


          !Column IDs
          negh_cID(4,1:4) = negh_cID(1,1:4)
          negh_cID(4,2) = negh_cID(4,1)+3
          negh_cID(4,4) = negh_cID(4,3)+3

          !Row ids
          negh_rID(4) = baseindex+3
        
#endif
#endif
        endif

        if(which_axis .eq. JAXIS) then

          !E_r eqn
          negh_values(1,1) = -1.*er_coeff(7)
          negh_values(1,2) = -1.*er_coeff(8)
          negh_values(1,3) = -1.*d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = d1y*(1-Cj0)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(7)
          negh_values(2,2) = -1*frx_coeff(8)
          negh_values(2,3) = d1y*(1-Cj0)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Fry eqn 
          negh_values(3,1) = -1*fry_coeff(7)
          negh_values(3,2) = -1.*fry_coeff(8)
          negh_values(3,3) = d1y*(1-Cj0)*solnData(YYED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt


          ! Column IDs
          !Er eqn
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 2
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+2

          !Frx eqn
          negh_cID(2,1:4) = negh_cID(1,1:4)
          negh_cID(2,2) = negh_cID(2,1)+1
          negh_cID(2,4) = negh_cID(2,3)+1

          !Fry eqn
          negh_cID(3,1:4) = negh_cID(1,1:4)
          
          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1
          negh_rID(3) = baseindex+2
#if NDIM>2
          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(7)
          negh_values(4,2) = -1.*frz_coeff(8)
          negh_values(4,3) = d1y*(1-Cj0)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Column Ids
          !Frz eqn
          negh_cID(4,1:4) = negh_cID(1,1:4)
          negh_cID(4,2) = negh_cID(4,1) + 3
          negh_cID(4,4) = negh_cID(4,3) + 3

          !Row ids       
          negh_rID(4) = baseindex+3

#endif
        endif 

        if(which_axis .eq. KAXIS) then

          !E_r eqn
          negh_values(1,1) = -1.*er_coeff(12)
          negh_values(1,2) = -1.*er_coeff(13)
          negh_values(1,3) = -1.*d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          negh_values(1,4) = d1z*(1-Ck0)

          !Frx eqn
          negh_values(2,1) = -1*frx_coeff(11)
          negh_values(2,2) = -1.*frx_coeff(12)
          negh_values(2,3) = d1z*(1-Ck0)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(2,4) = -1.*d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Fry eqn 
          negh_values(3,1) = -1*fry_coeff(11)
          negh_values(3,2) = -1.*fry_coeff(12)
          negh_values(3,3) = d1z*(1-Ck0)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(3,4) = -1.*d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          !Frz eqn
          negh_values(4,1) = -1*frz_coeff(11)
          negh_values(4,2) = -1.*frz_coeff(12)
          negh_values(4,3) = d1z*(1-Ck0)*solnData(ZZED_VAR,i,j,k)*rt_speedlt**2
          negh_values(4,4) = -1.*d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          !Column Ids
          !Er eqn
          negh_cID(1,1) = baseindex
          negh_cID(1,2) = negh_cID(1,1) + 3
          negh_cID(1,3) = globalmatrixindex(i_b,j_b,k_b,blockID,myPE)
          negh_cID(1,4) = negh_cID(1,3)+3

          !Frx eqn
          negh_cID(2,1:4) = negh_cID(1,1:4)
          negh_cID(2,2) = negh_cID(2,1)+1
          negh_cID(2,4) = negh_cID(2,3)+1

          !Fry eqn
          negh_cID(3,1:4) = negh_cID(1,1:4)
          negh_cID(3,2) = negh_cID(3,1) + 2
          negh_cID(3,4) = negh_cID(3,3) + 2

          !Frz eqn
          negh_cID(4,1:4) = negh_cID(1,1:4)
          
          !Row ids
          negh_rID(1) = baseindex
          negh_rID(2) = baseindex+1
          negh_rID(3) = baseindex+2
          negh_rID(4) = baseindex+3 
      
        endif

      endif

      ! Correct for factors related to fine-coarse cell differences
      ! 1/2 factor because denominator of d1x contains fine-cell width = 2*coarse one
      cell_length_factor = 1/2.
      ! 1/2^(ndim-1) factor due to the fact that the face area is divided into this many parts
      geometrical_factor = (1./2**(NDIM-1))
      negh_values = negh_values*cell_length_factor*geometrical_factor

      ! Set the coefficients for neighbour cell
      negh_no_rows = NDIM+1
      negh_no_columns = 4
      call MatSetValues(Amat,1,negh_rID(1),4,negh_cID(1,:),negh_values(1,:),ADD_VALUES,ierr)
      IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix at Bboundary.")
      call MatSetValues(Amat,1,negh_rID(2),4,negh_cID(2,:),negh_values(2,:),ADD_VALUES,ierr)
      IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix at Bboundary.")
#if NDIM>1      
      call MatSetValues(Amat,1,negh_rID(3),4,negh_cID(3,:),negh_values(3,:),ADD_VALUES,ierr)
      IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix at Bboundary.")
#if NDIM>2      
      call MatSetValues(Amat,1,negh_rID(4),4,negh_cID(4,:),negh_values(4,:),ADD_VALUES,ierr)
      IF(ierr.NE.0) call Driver_abortFlash("[VET]:Error while initializing Petsc Matrix at Bboundary.")
#endif
#endif

    CASE(1) ! Neighbour at higher refinement

      ! In multidimensions neighbours should be greater than 1
      if(numNegh .eq. 1 .and. NDIM .gt. 1) then 
        call Driver_abortFlash("[VET] : ERROR, Neighbour higher refinement level - only one neighbour found with NDIM>1.")
      endif

      ! Now subtract the terms that are already by the fine cell for this coarse cell
      ! 1. Contribution from boundary related flux terms removed
      ! 2. HLLE flux terms arising in current cell coefficients removed for appropriate boundary
      if(which_boundary .eq. LOW) then 
        if(which_axis .eq. IAXIS) then
        
          coeff_ID(1,1) = -1
          coeff_ID(1,2) = -1
          coeff(1,3) = coeff(1,3) - d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,4) = coeff(1,4) + d1x*(1-Ci0)

          coeff_ID(2,1) = -1
          coeff_ID(2,2) = -1
          coeff(2,3) = coeff(2,3) + d1x*(1-Ci0)*solnData(XXED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt

#if NDIM>1
          coeff_ID(3,1) = -1
          coeff_ID(3,2) = -1
          coeff(3,3) = coeff(3,3) + d1x*(1-Ci0)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt
#if NDIM>2          
          coeff_ID(4,1) = -1
          coeff_ID(4,2) = -1
          coeff(4,3) = coeff(4,3) + d1x*(1-Ci0)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1x*(1-Ci0)*wavespeed(which_boundary,which_axis)*rt_speedlt
#endif
#endif


        endif

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,7) = -1
          coeff_ID(1,8) = -1
          coeff(1,3) = coeff(1,3) - d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,9) = coeff(1,9) + d1y*(1-Cj0)

          coeff_ID(2,7) = -1
          coeff_ID(2,8) = -1
          coeff(2,3) = coeff(2,3) + d1y*(1-Cj0)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          coeff_ID(3,7) = -1
          coeff_ID(3,8) = -1
          coeff(3,3) = coeff(3,3) + d1y*(1-Cj0)*solnData(YYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt
#if NDIM>2          
          coeff_ID(4,7) = -1
          coeff_ID(4,8) = -1
          coeff(4,3) = coeff(4,3) + d1y*(1-Cj0)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1y*(1-Cj0)*wavespeed(which_boundary,which_axis)*rt_speedlt
#endif


        endif

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,12) = -1
          coeff_ID(1,13) = -1
          coeff(1,3) = coeff(1,3) - d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,14) = coeff(1,14) + d1z*(1-Ck0)

          coeff_ID(2,11) = -1
          coeff_ID(2,12) = -1
          coeff(2,3) = coeff(2,3) + d1z*(1-Ck0)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt

          coeff_ID(3,11) = -1
          coeff_ID(3,12) = -1
          coeff(3,3) = coeff(3,3) + d1z*(1-Ck0)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt
    
          coeff_ID(4,11) = -1
          coeff_ID(4,12) = -1
          coeff(4,3) = coeff(4,3) + d1z*(1-Ck0)*solnData(ZZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1z*(1-Ck0)*wavespeed(which_boundary,which_axis)*rt_speedlt

        endif
      endif


      if(which_boundary .eq. HIGH) then 
        if(which_axis .eq. IAXIS) then

          coeff_ID(1,5) = -1
          coeff_ID(1,6) = -1
          coeff(1,3) = coeff(1,3) - d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,4) = coeff(1,4) - d1x*(1+Ci1)

          coeff_ID(2,5) = -1
          coeff_ID(2,6) = -1
          coeff(2,3) = coeff(2,3) - d1x*(1+Ci1)*solnData(XXED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt

#if NDIM>1
          coeff_ID(3,5) = -1
          coeff_ID(3,6) = -1
          coeff(3,3) = coeff(3,3) - d1x*(1+Ci1)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt
#if NDIM>2          
          coeff_ID(4,5) = -1
          coeff_ID(4,6) = -1
          coeff(4,3) = coeff(4,3) - d1x*(1+Ci1)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1x*(1+Ci1)*wavespeed(which_boundary,which_axis)*rt_speedlt
#endif
#endif


        endif

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,10) = -1
          coeff_ID(1,11) = -1
          coeff(1,3) = coeff(1,3) - d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,9) = coeff(1,9) - d1y*(1+Cj1)

          coeff_ID(2,9) = -1
          coeff_ID(2,10) = -1
          coeff(2,3) = coeff(2,3) - d1y*(1+Cj1)*solnData(XYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          coeff_ID(3,9) = -1
          coeff_ID(3,10) = -1
          coeff(3,3) = coeff(3,3) - d1y*(1+Cj1)*solnData(YYED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt
#if NDIM>2          
          coeff_ID(4,9) = -1
          coeff_ID(4,10) = -1
          coeff(4,3) = coeff(4,3) - d1y*(1+Cj1)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1y*(1+Cj1)*wavespeed(which_boundary,which_axis)*rt_speedlt
#endif


        endif

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,15) = -1
          coeff_ID(1,16) = -1
          coeff(1,3) = coeff(1,3) - d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt
          coeff(1,14) = coeff(1,14) - d1z*(1+Ck1)

          coeff_ID(2,13) = -1
          coeff_ID(2,14) = -1
          coeff(2,3) = coeff(2,3) - d1z*(1+Ck1)*solnData(XZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(2,4) = coeff(2,4) - d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt

          coeff_ID(3,13) = -1
          coeff_ID(3,14) = -1
          coeff(3,3) = coeff(3,3) - d1z*(1+Ck1)*solnData(YZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(3,4) = coeff(3,4) - d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt
         
          coeff_ID(4,13) = -1
          coeff_ID(4,14) = -1
          coeff(4,3) = coeff(4,3) - d1z*(1+Ck1)*solnData(ZZED_VAR,i,j,k)*rt_speedlt**2 
          coeff(4,4) = coeff(4,4) - d1z*(1+Ck1)*wavespeed(which_boundary,which_axis)*rt_speedlt


        endif
      endif
    CASE(0) ! Neighbour at same refinement level as current block  

      if(numNegh .ne. 1) then 
        call Driver_abortFlash("[VET] : ERROR, Neighbour same refinement level - more than one neighbour found.")
      endif

      !This gives the block no and processor number of the neighbouring block appropriate for the direction
      neghBlkProc(BLKNO:PROCNO) = SurrBlkSum(b) &
            % regionInfo(edge(1),edge(2),edge(3)) % details(BLKNO:PROCNO,numNegh)

      baseindex = globalmatrixindex(negh_cell(1),negh_cell(2),negh_cell(3),neghBlkProc(BLKNO),neghBlkProc(PROCNO))

      if(which_boundary .eq. LOW) then
        if(which_axis .eq. IAXIS) then

          coeff_ID(1,1) = baseindex
          coeff_ID(1,2) = coeff_ID(1,1) + 1
          coeff_ID(2,1) = coeff_ID(1,1)
          coeff_ID(2,2) = coeff_ID(1,2)

#if NDIM>1
          coeff_ID(3,1) = coeff_ID(1,1)
          coeff_ID(3,2) = coeff_ID(1,1) + 2
#if NDIM>2
          coeff_ID(4,1) = coeff_ID(1,1)
          coeff_ID(4,2) = coeff_ID(1,1) + 3
#endif
#endif
        endif

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,7) = baseindex
          coeff_ID(1,8) = coeff_ID(1,7) + 2
          coeff_ID(2,7) = coeff_ID(1,7)
          coeff_ID(2,8) = coeff_ID(1,7) + 1

#if NDIM>1
          coeff_ID(3,7) = coeff_ID(1,7)
          coeff_ID(3,8) = coeff_ID(1,7) + 2
#if NDIM>2
          coeff_ID(4,7) = coeff_ID(1,7)
          coeff_ID(4,8) = coeff_ID(1,7) + 3
#endif
#endif
      
        endif

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,12) = baseindex
          coeff_ID(1,13) = coeff_ID(1,12) + 3
          coeff_ID(2,11) = coeff_ID(1,12)
          coeff_ID(2,12) = coeff_ID(1,12) + 1

#if NDIM>1
          coeff_ID(3,11) = coeff_ID(1,12)
          coeff_ID(3,12) = coeff_ID(1,12) + 2
#if NDIM>2
          coeff_ID(4,11) = coeff_ID(1,12)
          coeff_ID(4,12) = coeff_ID(1,12) + 3
#endif
#endif
      
        endif

      endif        

      if(which_boundary .eq. HIGH) then
        if(which_axis .eq. IAXIS) then
          
          coeff_ID(1,5) = baseindex
          coeff_ID(1,6) = coeff_ID(1,5) + 1
          coeff_ID(2,5) = coeff_ID(1,5)
          coeff_ID(2,6) = coeff_ID(1,6)

#if NDIM>1
          coeff_ID(3,5) = coeff_ID(1,5)
          coeff_ID(3,6) = coeff_ID(1,5) + 2
#if NDIM>2
          coeff_ID(4,5) = coeff_ID(1,5)
          coeff_ID(4,6) = coeff_ID(1,5) + 3
#endif
#endif

        endif

        if(which_axis .eq. JAXIS) then

          coeff_ID(1,10) = baseindex
          coeff_ID(1,11) = coeff_ID(1,10) + 2
          coeff_ID(2,9) = coeff_ID(1,10)
          coeff_ID(2,10) = coeff_ID(1,10) + 1

#if NDIM>1
          coeff_ID(3,9) = coeff_ID(1,10)
          coeff_ID(3,10) = coeff_ID(1,10) + 2
#if NDIM>2
          coeff_ID(4,9) = coeff_ID(1,10)
          coeff_ID(4,10) = coeff_ID(1,10) + 3
#endif
#endif
      
        endif

        if(which_axis .eq. KAXIS) then

          coeff_ID(1,15) = baseindex
          coeff_ID(1,16) = coeff_ID(1,15) + 3
          coeff_ID(2,13) = coeff_ID(1,15)
          coeff_ID(2,14) = coeff_ID(1,15) + 1

#if NDIM>1
          coeff_ID(3,13) = coeff_ID(1,15)
          coeff_ID(3,14) = coeff_ID(1,15) + 2
#if NDIM>2
          coeff_ID(4,13) = coeff_ID(1,15)
          coeff_ID(4,14) = coeff_ID(1,15) + 3
#endif
#endif

        endif
      endif


    CASE DEFAULT 
        call Driver_abortFlash("[VET] : Unrecognised level boundary condition.")

    END SELECT 




  END SUBROUTINE SetBlockBoundary



  SUBROUTINE SetPhysicalBoundary(which_boundary,which_axis)

    use RadTrans_data
    use Grid_interface, ONLY: Grid_getSingleCellCoords, Grid_getBlkIDFromPos
    INTEGER, INTENT(IN) :: which_axis, which_boundary
    type(PetscErrorCode) :: ierr
    type(PetscScalar), parameter :: one = 1
    INTEGER, EXTERNAL :: globalmatrixindex
    integer :: i_b, j_b, k_b,baseindex
    INTEGER, DIMENSION(IAXIS:KAXIS), parameter :: NX = (/NXB,NYB,NZB/)
    INTEGER, DIMENSION(MDIM) :: edge, negh_cell
    INTEGER, DIMENSION(BLKNO:PROCNO) :: neghBlkProc
    INTEGER :: numNegh
    REAL :: delx, taur, F_inc, rho, a1, b1, kE, kF

      SELECT CASE(rad_BC(which_boundary,which_axis))
      ! Reflecting BC
      CASE(1)
        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,1) = -1
            coeff_ID(1,2) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,1)
            coeff(1,4) = coeff(1,4) - coeff(1,2)

            coeff_ID(2:NDIM+1,1) = -1
            coeff_ID(2:NDIM+1,2) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,1)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,2)
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,7) = -1
            coeff_ID(1,8) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,7)
            coeff(1,9) = coeff(1,9) - coeff(1,8)

            coeff_ID(2:NDIM+1,7) = -1
            coeff_ID(2:NDIM+1,8) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,7)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,8)
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,12) = -1
            coeff_ID(1,13) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,12)
            coeff(1,14) = coeff(1,14) - coeff(1,13)

            coeff_ID(2:NDIM+1,11) = -1
            coeff_ID(2:NDIM+1,12) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,11)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,12)
          endif

        endif

        if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,5) = -1
            coeff_ID(1,6) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,5)
            coeff(1,4) = coeff(1,4) - coeff(1,6)

            coeff_ID(2:NDIM+1,5) = -1
            coeff_ID(2:NDIM+1,6) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,5)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,6)
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,10) = -1
            coeff_ID(1,11) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,10)
            coeff(1,9) = coeff(1,9) - coeff(1,11)

            coeff_ID(2:NDIM+1,9) = -1
            coeff_ID(2:NDIM+1,10) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,9)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,10)
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,15) = -1
            coeff_ID(1,16) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,15)
            coeff(1,14) = coeff(1,14) - coeff(1,16)

            coeff_ID(2:NDIM+1,13) = -1
            coeff_ID(2:NDIM+1,14) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,13)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) - coeff(2:NDIM+1,14)
          endif

        endif
      !Outflow BC
      CASE(2)   
        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,1) = -1
            coeff_ID(1,2) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,1)
            coeff(1,4) = coeff(1,4) + coeff(1,2)

            coeff_ID(2:NDIM+1,1) = -1
            coeff_ID(2:NDIM+1,2) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,1)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,2)
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,7) = -1
            coeff_ID(1,8) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,7)
            coeff(1,9) = coeff(1,9) + coeff(1,8)

            coeff_ID(2:NDIM+1,7) = -1
            coeff_ID(2:NDIM+1,8) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,7)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,8)
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,12) = -1
            coeff_ID(1,13) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,12)
            coeff(1,14) = coeff(1,14) + coeff(1,13)

            coeff_ID(2:NDIM+1,11) = -1
            coeff_ID(2:NDIM+1,12) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,11)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,12)
          endif

        endif

        if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,5) = -1
            coeff_ID(1,6) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,5)
            coeff(1,4) = coeff(1,4) + coeff(1,6)

            coeff_ID(2:NDIM+1,5) = -1
            coeff_ID(2:NDIM+1,6) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,5)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,6)
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,10) = -1
            coeff_ID(1,11) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,10)
            coeff(1,9) = coeff(1,9) + coeff(1,11)

            coeff_ID(2:NDIM+1,9) = -1
            coeff_ID(2:NDIM+1,10) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,9)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,10)
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,15) = -1
            coeff_ID(1,16) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,15)
            coeff(1,14) = coeff(1,14) + coeff(1,16)

            coeff_ID(2:NDIM+1,13) = -1
            coeff_ID(2:NDIM+1,14) = -1
            coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,13)
            coeff(2:NDIM+1,4) = coeff(2:NDIM+1,4) + coeff(2:NDIM+1,14)
          endif

        endif
      ! Dirichlet BC
      CASE(3) 
        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1:NDIM+1,1) = -1
            coeff_ID(1:NDIM+1,2) = -1
            vec(:) = vec(:) - coeff(:,1)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,2)*fradx_boundary_value
            vec(2) = vec(2) - coeff(2,2)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,2)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,2)*fradz_boundary_value
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1:NDIM+1,7) = -1
            coeff_ID(1:NDIM+1,8) = -1
            coeff_ID(1,8) = -1
            coeff_ID(2:NDIM+1,8) = -1

            vec(:) = vec(:) - coeff(:,7)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,8)*frady_boundary_value
            vec(2) = vec(2) - coeff(2,8)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,8)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,8)*fradz_boundary_value
#endif
#endif
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,12:13) = -1
            coeff_ID(2:NDIM+1,11) = -1
            coeff_ID(2:NDIM+1,12) = -1

            vec(1) = vec(1) - coeff(1,12)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,13)*fradz_boundary_value

            vec(2:NDIM+1) = vec(2:NDIM+1) - coeff(2:NDIM+1,11)*rt_radconst*rt_boundary_T**4
            vec(2) = vec(2) - coeff(2,12)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,12)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,12)*fradz_boundary_value
#endif
#endif
          endif
        endif

        if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1:NDIM+1,5) = -1
            coeff_ID(1:NDIM+1,6) = -1
            vec(:) = vec(:) - coeff(:,5)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,6)*fradx_boundary_value
            vec(2) = vec(2) - coeff(2,6)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,6)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,6)*fradz_boundary_value
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,10:11) = -1
            coeff_ID(2:NDIM+1,9) = -1
            coeff_ID(2:NDIM+1,10) = -1

            vec(1) = vec(1) - coeff(1,10)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,11)*frady_boundary_value

            vec(2:NDIM+1) = vec(2:NDIM+1) - coeff(2:NDIM+1,9)*rt_radconst*rt_boundary_T**4
            vec(2) = vec(2) - coeff(2,10)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,10)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,10)*fradz_boundary_value
#endif
#endif
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,15:16) = -1
            coeff_ID(2:NDIM+1,13) = -1
            coeff_ID(2:NDIM+1,14) = -1

            vec(1) = vec(1) - coeff(1,15)*rt_radconst*rt_boundary_T**4
            vec(1) = vec(1) - coeff(1,16)*fradz_boundary_value

            vec(2:NDIM+1) = vec(2:NDIM+1) - coeff(2:NDIM+1,13)*rt_radconst*rt_boundary_T**4
            vec(2) = vec(2) - coeff(2,14)*fradx_boundary_value
#if NDIM>1
            vec(3) = vec(3) - coeff(3,14)*frady_boundary_value
#if NDIM>2
            vec(4) = vec(4) - coeff(4,14)*fradz_boundary_value
#endif
#endif
          endif
        endif
      !Periodic BC
      CASE(4)
        
        i_b = i-blkLimits(LOW, IAXIS)
        j_b = j-blkLimits(LOW, JAXIS)
        k_b = k-blkLimits(LOW, KAXIS)

        edge(:) = (/ 1+K1D, 1+K2D, 1+K3D /)
        negh_cell(:) = (/ i_b, j_b, k_b /)

        if(which_boundary .eq. LOW) then
          edge(which_axis) = LEFT_EDGE
          negh_cell(which_axis) = NX(which_axis) - 1
        else 
          edge(which_axis) = RIGHT_EDGE
          negh_cell(which_axis) = 0
        endif

        if((NeghLevels(b,edge(1),edge(2),edge(3)) - mylevel) .ne. 0) then 
          call Driver_abortFlash("[VET] : Neighbouring block at different level. Not implemented yet.")
        endif

        numNegh = SurrBlkSum(b) % regionInfo(edge(1),edge(2),edge(3)) % numNegh

        if(numNegh .ne. 1) then 
          call Driver_abortFlash("[VET] : More than one neighbour found here.")
        endif

        !This gives the block no and processor number of the neighbouring block appropriate for the direction
        neghBlkProc(BLKNO:PROCNO) = SurrBlkSum(b) &
              % regionInfo(edge(1),edge(2),edge(3)) % details(BLKNO:PROCNO,numNegh)

        baseindex = globalmatrixindex(negh_cell(1),negh_cell(2),negh_cell(3),neghBlkProc(BLKNO),neghBlkProc(PROCNO))

        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then

            coeff_ID(1,1) = baseindex
            coeff_ID(1,2) = coeff_ID(1,1) + 1 !Frx
            coeff_ID(2,1) = coeff_ID(1,1) !Er 
            coeff_ID(2,2) = coeff_ID(1,1) + 1 !Frx
#if NDIM>1
            coeff_ID(3,1) = coeff_ID(1,1)!Er 
            coeff_ID(3,2) = coeff_ID(1,1) + 2 !Fry
#if NDIM>2
            coeff_ID(4,1) = coeff_ID(1,1) !Er 
            coeff_ID(4,2) = coeff_ID(1,1) + 3 !Frz
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then
          
            coeff_ID(1,7) = baseindex
            coeff_ID(1,8) = coeff_ID(1,7) + 2 !Fry
            coeff_ID(2,7) = coeff_ID(1,7) !Er 
            coeff_ID(2,8) = coeff_ID(1,7) + 1 !Frx
#if NDIM>1
            coeff_ID(3,7) = coeff_ID(1,7) !Er 
            coeff_ID(3,8) = coeff_ID(1,7) + 2 !Fry
#if NDIM>2
            coeff_ID(4,7) = coeff_ID(1,7) !Er 
            coeff_ID(4,8) = coeff_ID(1,7) + 3 !Frz
#endif
#endif
          endif

          if(which_axis .eq. KAXIS) then

            coeff_ID(1,12) = baseindex
            coeff_ID(1,13) = coeff_ID(1,12) + 3 !Frz
            coeff_ID(2,11) = coeff_ID(1,12) !Er 
            coeff_ID(2,12) = coeff_ID(1,12) + 1 !Frx
#if NDIM>1
            coeff_ID(3,11) = coeff_ID(1,12) !Er 
            coeff_ID(3,12) = coeff_ID(1,12) + 2 !Fry
#if NDIM>2
            coeff_ID(4,11) = coeff_ID(1,12) !Er 
            coeff_ID(4,12) = coeff_ID(1,12) + 3 !Frz
#endif
#endif
          endif
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then

            coeff_ID(1,5) = baseindex
            coeff_ID(1,6) = coeff_ID(1,5) + 1 !Frx
            coeff_ID(2,5) = coeff_ID(1,5) !Er 
            coeff_ID(2,6) = coeff_ID(1,5) + 1 !Frx
#if NDIM>1
            coeff_ID(3,5) = coeff_ID(1,5)!Er 
            coeff_ID(3,6) = coeff_ID(1,5) + 2 !Fry
#if NDIM>2
            coeff_ID(4,5) = coeff_ID(1,5) !Er 
            coeff_ID(4,6) = coeff_ID(1,5) + 3 !Frz
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then

            coeff_ID(1,10) = baseindex
            coeff_ID(1,11) = coeff_ID(1,10) + 2 !Fry
            coeff_ID(2,9) = coeff_ID(1,10) !Er 
            coeff_ID(2,10) = coeff_ID(1,10) + 1 !Frx
#if NDIM>1
            coeff_ID(3,9) = coeff_ID(1,10) !Er 
            coeff_ID(3,10) = coeff_ID(1,10) + 2 !Fry
#if NDIM>2
            coeff_ID(4,9) = coeff_ID(1,10) !Er 
            coeff_ID(4,10) = coeff_ID(1,10) + 3 !Frz
#endif
#endif
          endif

          if(which_axis .eq. KAXIS) then

            coeff_ID(1,15) = baseindex
            coeff_ID(1,16) = coeff_ID(1,15) + 3 !Frx
            coeff_ID(2,13) = coeff_ID(1,15) !Er 
            coeff_ID(2,14) = coeff_ID(1,15) + 1 !Frx
#if NDIM>1
            coeff_ID(3,13) = coeff_ID(1,15) !Er 
            coeff_ID(3,14) = coeff_ID(1,15) + 2 !Fry
#if NDIM>2
            coeff_ID(4,13) = coeff_ID(1,15) !Er 
            coeff_ID(4,14) = coeff_ID(1,15) + 3 !Frz
#endif
#endif 
          endif
        endif

      ! Vacuum BCs
      CASE(5)
        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,1) = -1
            coeff_ID(1,2) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,1) + coeff(1,2)*wavespeed_negh(LOW,IAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,1) = -1
            coeff_ID(2:NDIM+1,2) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,1) + coeff(2,2)*wavespeed_negh(LOW,IAXIS)*rt_speedlt
#if NDIM>1
            coeff(3,3) = coeff(3,3) + coeff(3,1)
            coeff(3,4) = coeff(3,4) + coeff(3,2)
#if NDIM>2
            coeff(4,3) = coeff(4,3) + coeff(4,1)
            coeff(4,4) = coeff(4,4) + coeff(4,2)
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,7) = -1
            coeff_ID(1,8) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,7) + coeff(1,8)*wavespeed_negh(LOW,JAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,7) = -1
            coeff_ID(2:NDIM+1,8) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,7)
            coeff(2,4) = coeff(2,4) + coeff(2,8)
            coeff(3,3) = coeff(3,3) + coeff(3,7) + coeff(3,8)*wavespeed_negh(LOW,JAXIS)*rt_speedlt
#if NDIM>2
            coeff(4,3) = coeff(4,3) + coeff(4,7)
            coeff(4,4) = coeff(4,4) + coeff(4,8)
#endif
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,12) = -1
            coeff_ID(1,13) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,12) + coeff(1,13)*wavespeed_negh(LOW,KAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,11) = -1
            coeff_ID(2:NDIM+1,12) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,11)
            coeff(2,4) = coeff(2,4) + coeff(2,12)
            coeff(3,3) = coeff(3,3) + coeff(3,11)
            coeff(3,4) = coeff(3,4) + coeff(3,12)
            coeff(4,3) = coeff(4,3) + coeff(4,11) + coeff(4,12)*wavespeed_negh(LOW,KAXIS)*rt_speedlt
          endif

        endif

        if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,5) = -1
            coeff_ID(1,6) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,5) + coeff(1,6)*wavespeed_negh(HIGH,IAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,5) = -1
            coeff_ID(2:NDIM+1,6) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,5) + coeff(2,6)*wavespeed_negh(HIGH,IAXIS)*rt_speedlt
#if NDIM>1
            coeff(3,3) = coeff(3,3) + coeff(3,5)
            coeff(3,4) = coeff(3,4) + coeff(3,6)
#if NDIM>2
            coeff(4,3) = coeff(4,3) + coeff(4,5)
            coeff(4,4) = coeff(4,4) + coeff(4,6)
#endif
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,10) = -1
            coeff_ID(1,11) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,10) + coeff(1,11)*wavespeed_negh(HIGH,JAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,9) = -1
            coeff_ID(2:NDIM+1,10) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,9) 
            coeff(2,4) = coeff(2,4) + coeff(2,10)
            coeff(3,3) = coeff(3,3) + coeff(3,9) + coeff(3,10)*wavespeed_negh(HIGH,JAXIS)*rt_speedlt
#if NDIM>2
            coeff(4,3) = coeff(4,3) + coeff(4,9) 
            coeff(4,4) = coeff(4,4) + coeff(4,10)
#endif
          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,15) = -1
            coeff_ID(1,16) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,15) + coeff(1,16)*wavespeed_negh(HIGH,KAXIS)*rt_speedlt

            coeff_ID(2:NDIM+1,13) = -1
            coeff_ID(2:NDIM+1,14) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,13)
            coeff(2,4) = coeff(2,4) + coeff(2,14)
            coeff(3,3) = coeff(3,3) + coeff(3,13)
            coeff(3,4) = coeff(3,4) + coeff(3,14)
            coeff(4,3) = coeff(4,3) + coeff(4,13) + coeff(4,14)*wavespeed_negh(HIGH,KAXIS)*rt_speedlt
          endif

        endif
      ! Marshak Boundary Conditions.
      CASE(6)
        taur = solnData(TAUR_VAR,i,j,k)
        rho = solnData(DENS_VAR,i,j,k)
        delx = cellsize(which_axis)
        !Set boundary radiation field to zero for UV bands; for IR user-set boundary condition
        if(current_band .eq. 'IR') then
          F_inc = rt_radconst*rt_speedlt/4 *rt_boundary_T**4
        else
          !TODO: This can be something else if a backround UV field is required.
          F_inc = 0.0
        endif
        a1 = (2*rt_speedlt/(3*taur*delx) - rt_speedlt/2.)/(2*rt_speedlt/(3*taur*delx) + rt_speedlt/2.)
        b1 = (4*rt_speedlt)/(3*taur*delx + 4)
        kE = (4*F_inc)/(2*rt_speedlt/(3*taur*delx) + rt_speedlt/2.)
        kF = (16*F_inc)/(3*taur*delx + 4.)
        if(which_boundary .eq. LOW) then 
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,1) = -1
            coeff_ID(1,2) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,1)* a1 - coeff(1,2)* b1
            coeff(1,4) = coeff(1,4) - coeff(1,2)
            vec(1) = vec(1) - coeff(1,1)*kE - coeff(1,2) * kF

            coeff_ID(2,1) = -1
            coeff_ID(2,2) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,1)* a1 - coeff(2,2)*b1
            coeff(2,4) = coeff(2,4) - coeff(2,2)
            vec(2) = vec(2) - coeff(2,1)*kE - coeff(2,2) * kF

            !Reflecting BC for Fry and Frz, and same condition as above for E_r
#if NDIM>1
            coeff_ID(3:NDIM+1,1) = -1
            coeff_ID(3:NDIM+1,2) = -1
            !Same condition for E_r as above
            coeff(3:NDIM+1,3) = coeff(3:NDIM+1,3) + coeff(3:NDIM+1,1)*a1
            !Reflecting for the flux
            coeff(3:NDIM+1,4) = coeff(3:NDIM+1,4) - coeff(3:NDIM+1,2)
            vec(3:NDIM+1) = vec(3:NDIM+1) - coeff(3:NDIM+1,1)*kE
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,7) = -1
            coeff_ID(1,8) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,7)* a1 - coeff(1,8)*b1
            coeff(1,9) = coeff(1,9) - coeff(1,8)
            vec(1) = vec(1) - coeff(1,7)* kE - coeff(1,8) * kF

            coeff_ID(2,7) = -1
            coeff_ID(2,8) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,7) * a1
            coeff(2,4) = coeff(2,4) - coeff(2,8)
            vec(2) = vec(2) - coeff(2,7) * kE
            
            !Marshak for Fy
            coeff_ID(3,7) = -1
            coeff_ID(3,8) = -1
            coeff(3,3) = coeff(3,3) + coeff(3,7)* a1 - coeff(3,8)* b1
            coeff(3,4) = coeff(3,4) - coeff(3,8)
            vec(3) = vec(3) - coeff(3,7)*kE - coeff(3,8) * kF

#if NDIM>2
            coeff_ID(4,7) = -1
            coeff_ID(4,8) = -1
            coeff(4,3) = coeff(4,3) + coeff(4,7) * a1
            coeff(4,4) = coeff(4,4) - coeff(4,8)
            vec(4) = vec(4) - coeff(4,7) * kE
#endif

          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,12) = -1
            coeff_ID(1,13) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,12)* a1 - coeff(1,13)*b1
            coeff(1,14) = coeff(1,14) - coeff(1,13)
            vec(1) = vec(1) - coeff(1,12)*kE - coeff(1,13) * kF

            coeff_ID(2:3,11) = -1
            coeff_ID(2:3,12) = -1
            coeff(2:3,3) = coeff(2:3,3) + coeff(2:3,11) * a1
            coeff(2:3,4) = coeff(2:3,4) - coeff(2:3,12)
            vec(2:3) = vec(2:3) - coeff(2:3,11)*kE
            
            !Marshak for Fz
            coeff_ID(4,11) = -1
            coeff_ID(4,12) = -1
            coeff(4,3) = coeff(4,3) + coeff(4,11)*a1 - coeff(4,12)*b1
            coeff(4,4) = coeff(4,4) - coeff(4,12)
            vec(4) = vec(4) - coeff(4,11)* kE - coeff(4,12) * kF
          endif
        endif

        if(which_boundary .eq. HIGH) then 
          if(which_axis .eq. IAXIS) then
            coeff_ID(1,5) = -1
            coeff_ID(1,6) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,5)* a1 + coeff(1,6)*b1
            coeff(1,4) = coeff(1,4) - coeff(1,6)
            vec(1) = vec(1) - coeff(1,5)*kE + coeff(1,6) * kF

            coeff_ID(2,5) = -1
            coeff_ID(2,6) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,5)*a1 + coeff(2,6)*b1
            coeff(2,4) = coeff(2,4) - coeff(2,6)
            vec(2) = vec(2) - coeff(2,5)*kE + coeff(2,6) * kF

#if NDIM>1
            coeff_ID(3:NDIM+1,5) = -1
            coeff_ID(3:NDIM+1,6) = -1
            coeff(3:NDIM+1,3) = coeff(3:NDIM+1,3) + coeff(3:NDIM+1,5) * a1
            coeff(3:NDIM+1,4) = coeff(3:NDIM+1,4) - coeff(3:NDIM+1,6)
            vec(3:NDIM+1) = vec(3:NDIM+1) - coeff(3:NDIM+1,5) * kE
#endif
          endif

          if(which_axis .eq. JAXIS) then
            coeff_ID(1,10) = -1
            coeff_ID(1,11) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,10)* a1 + coeff(1,11)*b1
            coeff(1,9) = coeff(1,9) - coeff(1,11)
            vec(1) = vec(1) - coeff(1,10)* kE + coeff(1,11) * kF

            ! Outflow for Fx
            coeff_ID(2,9) = -1
            coeff_ID(2,10) = -1
            coeff(2,3) = coeff(2,3) + coeff(2,9) * a1
            coeff(2,4) = coeff(2,4) - coeff(2,10)
            vec(2) = vec(2) - coeff(2,9) * kE
            
            !Marshak for Fy
            coeff_ID(3,9) = -1
            coeff_ID(3,10) = -1
            coeff(3,3) = coeff(3,3) + coeff(3,9)*a1 + coeff(3,10)*b1
            coeff(3,4) = coeff(3,4) - coeff(3,10)
            vec(3) = vec(3) - coeff(3,9)*kE + coeff(3,10) * kF

#if NDIM>2
            coeff_ID(4,9) = -1
            coeff_ID(4,10) = -1
            coeff(4,3) = coeff(4,3) + coeff(4,9) * a1
            coeff(4,4) = coeff(4,4) - coeff(4,10)
            vec(4) = vec(4) - coeff(4,9) * kE
#endif

          endif

          if(which_axis .eq. KAXIS) then
            coeff_ID(1,15) = -1
            coeff_ID(1,16) = -1
            coeff(1,3) = coeff(1,3) + coeff(1,15)* a1 + coeff(1,16)*b1
            coeff(1,14) = coeff(1,14) - coeff(1,16)
            vec(1) = vec(1) - coeff(1,15)*kE + coeff(1,16) * kF

            coeff_ID(2:3,13) = -1
            coeff_ID(2:3,14) = -1
            coeff(2:3,3) = coeff(2:3,3) + coeff(2:3,13) * a1
            coeff(2:3,4) = coeff(2:3,4) - coeff(2:3,14)
            vec(2:3) = vec(2:3) - coeff(2:3,13) * kE
            
            !Marshak for Fz
            coeff_ID(4,13) = -1
            coeff_ID(4,14) = -1
            coeff(4,3) = coeff(4,3) + coeff(4,13)* a1 + coeff(4,14)* b1
            coeff(4,4) = coeff(4,4) - coeff(4,14)
            vec(4) = vec(4) - coeff(4,13)*kE + coeff(4,14) * kF
          endif
        endif
      !User-defined BC for problem: streaming radiation inwards with given flux, zero-gradient BCs for radiation energy
      CASE(7) 
        if(which_boundary .eq. LOW .and. which_axis .eq. IAXIS) then
          
          coeff_ID(1:NDIM+1,1) = -1
          coeff_ID(1:NDIM+1,2) = -1
          coeff(1,3) = coeff(1,3) + coeff(1,1)
          coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,1)
          vec(1) = vec(1) - coeff(1,2)*fradx_boundary_value
          vec(2) = vec(2) - coeff(2,2)*fradx_boundary_value
#if NDIM>1
          vec(3) = vec(3) - coeff(3,2)*frady_boundary_value
#if NDIM>2
          vec(4) = vec(4) - coeff(4,2)*fradz_boundary_value
#endif
#endif
        else if(which_boundary .eq. LOW .and. which_axis .eq. JAXIS) then
                  
          coeff_ID(1:NDIM+1,7) = -1
          coeff_ID(1:NDIM+1,8) = -1
          coeff(1,3) = coeff(1,3) + coeff(1,7)
          coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,7)
          vec(1) = vec(1) - coeff(1,8)*frady_boundary_value
          vec(2) = vec(2) - coeff(2,8)*fradx_boundary_value
#if NDIM>1
          vec(3) = vec(3) - coeff(3,8)*frady_boundary_value
#if NDIM>2
          vec(4) = vec(4) - coeff(4,8)*fradz_boundary_value
#endif
#endif
        

        else if(which_boundary .eq. LOW .and. which_axis .eq. KAXIS) then
          
          coeff_ID(1,12) = -1
          coeff_ID(1,13) = -1
          coeff_ID(2:NDIM+1,11) = -1
          coeff_ID(2:NDIM+1,12) = -1
          coeff(1,3) = coeff(1,3) + coeff(1,12)
          coeff(2:NDIM+1,3) = coeff(2:NDIM+1,3) + coeff(2:NDIM+1,11)
          vec(1) = vec(1) - coeff(1,13)*fradz_boundary_value
          vec(2) = vec(2) - coeff(2,12)*fradx_boundary_value
#if NDIM>1
          vec(3) = vec(3) - coeff(3,12)*frady_boundary_value
#if NDIM>2
          vec(4) = vec(4) - coeff(4,12)*fradz_boundary_value
#endif
#endif
      
        else
          call Driver_abortFlash("[VET] : User-defined free-flowing BC only available for -ve x direction BC.")
        endif


      CASE DEFAULT 
        print *, rad_BC(which_boundary,which_axis)
        call Driver_abortFlash("[VET] : Check Radiation BC defined in flash.par.")

      END SELECT 

  END SUBROUTINE SetPhysicalBoundary

  SUBROUTINE CheckBoundary(which_boundary,which_axis,boundary_inblock)

    use RadTrans_data
    INTEGER, INTENT(IN) :: which_axis, which_boundary
    LOGICAL, DIMENSION(2,NDIM), INTENT(OUT) :: boundary_inblock

    IF(faces_BC(which_boundary,which_axis) .ne. NOT_BOUNDARY & 
      .or. onBoundary(which_boundary,which_axis) &
      .ne. NOT_BOUNDARY) THEN

      boundary_inblock(which_boundary,which_axis) = .true.

    ELSE
      boundary_inblock(which_boundary,which_axis) = .false.


    ENDIF

  END SUBROUTINE CheckBoundary
   
  SUBROUTINE SetEddBoundary(which_boundary,which_axis)
  use RadTrans_data
  INTEGER, INTENT(IN) :: which_axis, which_boundary
  INTEGER, DIMENSION(MDIM) :: edge
  INTEGER :: boundary_cell, i_b, j_b, k_b
  REAL  :: del_negh

  if(boundary_inblock(which_boundary,which_axis)) then

    SELECT CASE(rad_BC(which_boundary,which_axis))
      !Reflecting 
      CASE(1)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)

        !delx would change here in AMR


        if(which_axis .eq. IAXIS) then
          xyed(which_boundary,which_axis) = -solnData(XYED_VAR,i,j,k)
          xzed(which_boundary,which_axis) = -solnData(XZED_VAR,i,j,k)
          yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))

        else if(which_axis .eq. JAXIS) then
          
          xyed(which_boundary,which_axis) = -solnData(XYED_VAR,i,j,k)
          xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
          yzed(which_boundary,which_axis) = -solnData(YZED_VAR,i,j,k)
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))

        else if(which_axis .eq. KAXIS) then
          xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
          xzed(which_boundary,which_axis) = -solnData(XZED_VAR,i,j,k)
          yzed(which_boundary,which_axis) = -solnData(YZED_VAR,i,j,k)
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif
      !Outflow
      CASE(2)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)
        xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
        xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
        yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)

        if(which_axis .eq. IAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))
        else if(which_axis .eq. JAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))
        else if(which_axis .eq. KAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif
      !Dirichlet
      !Conditions for edd tensor kept identical to outflow BCs
      !TODO: Confirm if this is alright
      CASE(3)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)
        xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
        xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
        yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)

        if(which_axis .eq. IAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))
        else if(which_axis .eq. JAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))
        else if(which_axis .eq. KAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif

      !Periodic
      CASE(4)

        
        if(which_boundary .eq. LOW) then
          if(which_axis .eq. IAXIS) then

            if(Ncell(1) .gt. 1) then
              boundary_cell = i-1
            else
              boundary_cell = i
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,boundary_cell,j,k)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,boundary_cell,j,k)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,boundary_cell,j,k)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,boundary_cell,j,k)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,boundary_cell,j,k)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,boundary_cell,j,k)
            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,boundary_cell,j,k))**2 /&
                                              (2.0*xxed(which_boundary,which_axis))

          else if(which_axis .eq. JAXIS) then
            if(Ncell(2) .gt. 1) then
              boundary_cell = j-1
            else
              boundary_cell = j
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,boundary_cell,k)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,boundary_cell,k)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,boundary_cell,k)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,boundary_cell,k)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,boundary_cell,k)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,boundary_cell,k)

            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,boundary_cell,k))**2 /&
                                              (2.0*yyed(which_boundary,which_axis))

          else if(which_axis .eq. KAXIS) then
            if(Ncell(3) .gt. 1) then
              boundary_cell = k-1
            else
              boundary_cell = k
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,boundary_cell)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,boundary_cell)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,boundary_cell)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,boundary_cell)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,boundary_cell)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,boundary_cell)
            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,boundary_cell))**2 /&
                                              (2.0*zzed(which_boundary,which_axis))
          endif

        else if(which_boundary .eq. HIGH) then
          if(which_axis .eq. IAXIS) then

            if(Ncell(1) .gt. 1) then
              boundary_cell = i+1
            else
              boundary_cell = i
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,boundary_cell,j,k)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,boundary_cell,j,k)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,boundary_cell,j,k)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,boundary_cell,j,k)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,boundary_cell,j,k)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,boundary_cell,j,k) 

            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,boundary_cell,j,k))**2 /&
                                              (2.0*xxed(which_boundary,which_axis))

          else if(which_axis .eq. JAXIS) then
            if(Ncell(2) .gt. 1) then
              boundary_cell = j+1
            else
              boundary_cell = j
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,boundary_cell,k)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,boundary_cell,k)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,boundary_cell,k)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,boundary_cell,k)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,boundary_cell,k)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,boundary_cell,k)

            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,boundary_cell,k))**2 /&
                                              (2.0*yyed(which_boundary,which_axis))

          else if(which_axis .eq. KAXIS) then
            if(Ncell(3) .gt. 1) then
              boundary_cell = k+1
            else
              boundary_cell = k
            endif
            xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,boundary_cell)
            xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,boundary_cell)
            xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,boundary_cell)
            yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,boundary_cell)
            zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,boundary_cell)
            yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,boundary_cell) 
            tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,boundary_cell))**2 /&
                                              (2.0*zzed(which_boundary,which_axis))
          endif
        endif
      ! Vacuum BCs
      CASE(5)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)
        xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
        xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
        yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)

        if(which_axis .eq. IAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))
        else if(which_axis .eq. JAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))
        else if(which_axis .eq. KAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif

      CASE(6)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)
        xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
        xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
        yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)

        if(which_axis .eq. IAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))
        else if(which_axis .eq. JAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))
        else if(which_axis .eq. KAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif

      CASE(7)
        xxed(which_boundary,which_axis) = solnData(XXED_VAR,i,j,k)
        yyed(which_boundary,which_axis) = solnData(YYED_VAR,i,j,k)
        zzed(which_boundary,which_axis) = solnData(ZZED_VAR,i,j,k)
        xyed(which_boundary,which_axis) = solnData(XYED_VAR,i,j,k)
        xzed(which_boundary,which_axis) = solnData(XZED_VAR,i,j,k)
        yzed(which_boundary,which_axis) = solnData(YZED_VAR,i,j,k)

        if(which_axis .eq. IAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(XXED_VAR,i,j,k))
        else if(which_axis .eq. JAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(YYED_VAR,i,j,k))
        else if(which_axis .eq. KAXIS) then 
          tau_hlle_neigh(which_boundary,which_axis) = (10.0*cellsize(which_axis)*solnData(TAUR_VAR,i,j,k))**2 /&
                                              (2.0*solnData(ZZED_VAR,i,j,k))
        endif


      CASE DEFAULT
        call Driver_abortFlash("[VET] : Check Radiation BC defined in flash.par.")
      END SELECT

    else
      i_b = i-blkLimits(LOW, IAXIS)
      j_b = j-blkLimits(LOW, JAXIS)
      k_b = k-blkLimits(LOW, KAXIS)
      edge(:) = (/ 1+K1D, 1+K2D, 1+K3D /)
      if(which_boundary .eq. LOW) then
        edge(which_axis) = LEFT_EDGE
      else
        edge(which_axis) = RIGHT_EDGE
      endif

      SELECT CASE(NeghLevels(b,edge(1),edge(2),edge(3)) - mylevel)
        CASE(-1) ! Neighbour lower refinement than current block (fine-coarse boundary)
          del_negh = cellsize(which_axis)*2.0

        CASE(1)
          del_negh = cellsize(which_axis)/2.0

        CASE(0)
          del_negh = cellsize(which_axis)

        CASE DEFAULT 
          call Driver_abortFlash("[VET] : Unrecognised level boundary condition.")

      END SELECT

      if(which_axis .eq. IAXIS) then
        tau_hlle_neigh(which_boundary,which_axis) = (10.0*del_negh*solnData(TAUR_VAR,i+2*which_boundary-3,j,k))**2 /&
                                              (2*xxed(which_boundary,which_axis))
      else if(which_axis .eq. JAXIS) then 
        tau_hlle_neigh(which_boundary,which_axis) = (10.0*del_negh*solnData(TAUR_VAR,i,j+2*which_boundary-3,k))**2 /&
                                              (2*yyed(which_boundary,which_axis))
      else if(which_axis .eq. KAXIS) then
        tau_hlle_neigh(which_boundary,which_axis) = (10.0*del_negh*solnData(TAUR_VAR,i,j,k+2*which_boundary-3))**2 /&
                                              (2*zzed(which_boundary,which_axis)) 
      endif

    endif



  END SUBROUTINE SetEddBoundary

  SUBROUTINE CheckFCInterface(which_boundary,which_axis,FCInterface)
    use RadTrans_data
    INTEGER, INTENT(IN) :: which_axis, which_boundary
    LOGICAL, INTENT(INOUT) :: FCInterface
    INTEGER, DIMENSION(MDIM) :: edge

    edge(:) = (/ 1+K1D, 1+K2D, 1+K3D /)
    if(which_boundary .eq. LOW) then
      edge(which_axis) = LEFT_EDGE
    else
      edge(which_axis) = RIGHT_EDGE
    endif
    !Check if fine-cell interfaces present anywhere
    SELECT CASE(NeghLevels(b,edge(1),edge(2),edge(3)) - mylevel)
      CASE(-1) ! Neighbour lower refinement than current block (fine-coarse boundary)
        FCInterface = .true.

      CASE(1)
        FCInterface = .false.

      CASE(0)
        FCInterface = .false.

      CASE DEFAULT 
        call Driver_abortFlash("[VET] : Unrecognised level boundary condition.")

    END SELECT
  END SUBROUTINE CheckFCInterface

END SUBROUTINE rt_fillMatrix

!Returns Corrected Left and Right wave speed at an interface 
SUBROUTINE CorrectedWaveSpeed(taul,taur,fcflag,wavespeedl,wavespeedr)
  use RadTrans_data
  implicit none
  REAL, INTENT(IN) :: taul, taur
  LOGICAL, INTENT(IN) :: fcflag
  REAL, INTENT(INOUT) :: wavespeedl, wavespeedr
  REAL :: tauavg, taulwave, taurwave

  !Average value of tau of L and R cell centres at interface
  tauavg = (taul + taur)/2.
  
  !Left going wavespeed at interface
  if(fcflag) then 
    taulwave = taur !Upstream wavespeed
  else
    taulwave = tauavg
  endif
  
  !Safety Check: Return to uncorrected wavespeed to asymptotic optically thin limit
  if(exp(-taulwave) .eq. 1) then 
    wavespeedl = wavespeedl
  !If all good, correct wavespeed
  else
    wavespeedl = wavespeedl*SQRT((1.0-exp(-taulwave))/taulwave)
  endif

  !Right going wavespeed at interface
  if(fcflag) then 
    taurwave = taul !Upstream wavespeed
  else
    taurwave = tauavg
  endif
  
  !Return to uncorrected wavespeed to asymptotic optically thin limit
  if(exp(-taurwave) .eq. 1) then 
    wavespeedr = wavespeedr
  !If all good, correct wavespeed
  else
    wavespeedr = wavespeedr*SQRT((1.0-exp(-taurwave))/taurwave)
  endif

END SUBROUTINE CorrectedWaveSpeed

! Copied & modified from gr_findAllNeghID: subroutine gr_getPartensGuardCellView
! Find out which position this child block has in its parten block, e.g.
! LOW/HIGH in each dimension.
subroutine whichChild(srcBlkID, blkChild)

  use tree, only : which_child

  implicit none
  integer, intent(IN) :: srcBlkID
!  integer, dimension(MDIM), intent(IN) :: srcGuardCellID
!  integer, dimension(MDIM), intent(OUT) :: parentGuardCellID
  integer, dimension(MDIM), intent(OUT) :: blkChild

  integer, parameter :: TOTAL_CHILDREN = 2**MDIM
  integer, parameter, dimension(MDIM,TOTAL_CHILDREN) :: blockChild = &
       reshape (source = (/&
       LOW,  LOW,  LOW , &
       HIGH, LOW,  LOW , &
       LOW,  HIGH, LOW , &
       HIGH, HIGH, LOW , &
       LOW,  LOW,  HIGH, &
       HIGH, LOW,  HIGH, &
       LOW,  HIGH, HIGH, &
       HIGH, HIGH, HIGH  &
       /), shape = (/MDIM,TOTAL_CHILDREN/))
  integer :: eachAxis, child

  child = which_child(srcBlkID)

  !Obtain the source blocks position in the parent block.
  !e.g. child=5 is at location: LOW,LOW,HIGH.
  !We use this location to figure out how the source block's 
  !neighbor appears with respect to the parent.

  blkChild(1:MDIM) = blockChild(1:MDIM,child)

end subroutine whichChild


FUNCTION globalmatrixindex(i,j,k,blockID,myPE) RESULT(l)
  use RadTrans_data
  implicit none
  INTEGER, INTENT(IN) :: i,j,k,blockID,myPE
  INTEGER :: l
  l = (NDIM+1)*(i + NXB*(j + NYB*(k + NZB*rblockListAll(blockID,myPE))))
END FUNCTION globalmatrixindex

FUNCTION localblockindex(i,j,k,blockID) RESULT(l)
  implicit none
  INTEGER, INTENT(IN) :: i,j,k,blockID
  INTEGER :: l
  ! Corner ID's subtracted by one to transform to C-type matrix indexing 
  l = (NDIM+1)*(i + NXB*(j + NYB*(k + NZB*blockID)))
END FUNCTION localblockindex
