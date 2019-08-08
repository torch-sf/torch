!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_dataReconstOneStep
!!
!! NAME
!!
!!  hy_uhd_dataReconstOneStep
!!
!! SYNOPSIS
!!
!!  hy_uhd_dataReconstOneStep(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(:,:),
!!                            integer(IN) :: order,
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: dt,
!!                            real(IN)    :: del(MDIM),
!!                            real(IN)    :: ogravX(:,:,:),
!!                            real(IN)    :: ogravY(:,:,:),
!!                            real(IN)    :: ogravZ(:,:,:),
!!                            real(IN)    :: DivV(:,:,:),
!!                            real(IN)    :: FlatCoeff(:,:,:,:),
!!                            logical(IN) :: TransX_updateOnly,
!!                            logical(IN) :: TransY_updateOnly,
!!                            logical(IN) :: TransZ_updateOnly,
!!                            real(OUT)   :: Wp(HY_VARINUMMAX,NDIM),
!!                            real(OUT)   :: Wn(HY_VARINUMMAX,NDIM),
!!                            real(OUT)   :: sig(HY_VARINUMMAX,NDIM),
!!                            real(OUT)   :: lambda(HY_WAVENUM,NDIM),
!!                            real(OUT)   :: leftEig(HY_VARINUM,HY_VARINUM,NDIM),
!!                            real(OUT)   :: rghtEig(HY_VARINUM,HY_WAVENUM,NDIM),
!!                   OPTIONAL,real,pointer:: hy_SpcR(:,:,:,:,:),
!!                   OPTIONAL,real,pointer:: hy_SpcL(:,:,:,:,:),
!!                   OPTIONAL,real,pointer:: hy_SpcSig(:,:,:,:,:))
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimitsGC - block limits including guardcells
!!  order       - order of reconstruction method
!!  ix,iy,iz    - local indices
!!  dt          - timestep
!!  del         - deltas in each {x,y,z} direction
!!  ogravX,ogravY,ogravZ - gravity components at n-step in each direction
!!  DivV        - divergence of velocity fields used for hybrid order method
!!  FlatCoeff   - flattening parameters
!!  TransX_updateOnly - a switch for a selective transverse update in x direction
!!  TransY_updateOnly - a switch for a selective transverse update in y direction
!!  TransZ_updateOnly - a switch for a selective transverse update in z direction
!!  Wp,Wn       - right (plus) and left(negative) Riemann states
!!  sig         - a transverse flux term
!!  lambda      - eigenvalue
!!  leftEig     - left eigenvector
!!  rghtEig     - right eigenvector
!!  hy_SpcR,hy_SpcL,hy_SpcSig - pointers of Riemann states for Species and mass scalar
!!
!! DESCRIPTION
!!
!!  This onestep data reconstruction routine evolves the cell centered 
!!  values by dt/2 time step at each cell interface using characteristic 
!!  tracing method based on 2nd order MUSCL-Hancock, 3rd order PPM,
!!  or 5th order WENO schemes.
!!  For MHD, the cell interface values are reconstructed in multidimensional
!!  fashion with inclusions of MHD multidimensional terms that are
!!  proportional to divB for the unsplit MHD implementation.
!!
!! REFERENCES
!!
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella and Woodward, 54, 174 (1984), JCP
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***
#include "Flash.h"
#ifdef FLASH_USM_MHD
!! REORDER(4): B[xyz]
#endif

#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif
#ifdef DEBUG_HYDRO
#define DEBUG_HYDRO_POSITIVITY
#endif


Subroutine hy_uhd_dataReconstOneStep(blockID,blkLimitsGC,order,ix,iy,iz, &
                                     dt,del,ogravX,ogravY,ogravZ,&
                                     DivU, FlatCoeff,  &
                                     TransX_updateOnly,&
                                     TransY_updateOnly,&
                                     TransZ_updateOnly,&
                                     Wp,Wn,sig,&
                                     lambda,leftEig,rghtEig,&
                                     cellCfl,&
                                     hy_SpcR,hy_SpcL,hy_SpcSig)

  use Hydro_data,        ONLY : hy_cfl,hy_cfl_original,hy_wenoMethod,&
                                hy_cflFallbackFactor,&
                                hy_fallbackLowerCFL, &
                                hy_meshMe, &
                                hy_fullSpecMsFluxHandling,hy_numXN,&
                                hy_useHybridOrder,hy_shockDetectOn,&
                                hy_hybridOrderKappa,hy_eswitch, &
                                hy_useAuxEintEqn, hy_radiusGP
#ifdef FLASH_USM_MHD
  use Hydro_data,        ONLY : hy_killdivb
#endif
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  use Hydro_data,        ONLY : hy_forceHydroLimit
#endif

  use Timers_interface,  ONLY : Timers_start, Timers_stop

  use Grid_interface,    ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                                Grid_getCellCoords

  !! Use short nicknames for the reconstruction subroutines.
  !! Some compilers (e.g., Absoft) can be confused with similar subroutine names.
  use hy_uhd_interface,  ONLY : hy_uhd_checkRHjumpCond,&
                                MH   => hy_uhd_DataReconstructNormalDir_MH,&
                                PPM  => hy_uhd_DataReconstructNormalDir_PPM,&
                                WENO => hy_uhd_DataReconstructNormalDir_WENO,&
                                GP   => hy_uhd_DataReconstructNormalDir_GP

  implicit none

#include "constants.h"
#include "UHD.h"


  !!-----Arguments---------------------------------------------------------
  integer, intent(IN) :: blockID
  integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
  integer, intent(IN) :: order,ix,iy,iz
  real,    intent(IN) :: dt
  real,    intent(IN), dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
  real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN), target :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: FlatCoeff
  real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: DivU
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(HIGH,KAXIS)),intent(IN), target :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(HIGH,KAXIS)),intent(IN) :: FlatCoeff
  real, dimension(blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(HIGH,KAXIS)), intent(IN) :: DivU
#endif
  logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly

  real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: Wp, Wn
  real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: sig
  real, dimension(              HY_WAVENUM,NDIM),intent(OUT) :: lambda
  real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: leftEig
  real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: rghtEig
  real, optional, intent(INOUT) :: cellCfl ! We may pass back a lowered CFL factor here - KW
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  !!------------------------------------------------------------------------


  real, pointer, dimension(:)       :: DataGrav1D  
  real, pointer, dimension(:,:)     :: Data1D, LambdaArr
  real, pointer, dimension(:,:,:)   :: LeigArr, ReigArr
  real, pointer, dimension(:,:,:,:) :: U

  real    :: Cs,dx,dy,dz,hdt,hdtx,hdty,hdtz,idx,idy,idz
  real    :: epsilon ! = 1.e-8, kappa
  real    :: cflIn, maxCfl
  integer :: k2, k3, k4, kOrder, kOrder_orig, kLambda, iDim,kGrav
  integer :: hyEndVar,hyEndFlux,hyEndPrimVar
  logical :: doFullSpecMs
  logical, dimension(NDIM) :: SWp,SWn
  logical, dimension(MDIM) :: Trans_updateOnly
  real, pointer :: SpcSigPtr(:)
  real, target  :: SpcSigIgnore(1:0)
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /* More local variables for MHD */
  real    :: uB
  real, dimension(MDIM) :: dnBn
  real, dimension(HY_VARINUM,NDIM) :: aBn
#ifdef FLASH_UGLM_MHD /* extra source term in GLM */
  real, dimension(NDIM) :: dnGLM
  real, dimension(HY_VARINUM,NDIM) :: aGLM
#endif
#if NFACE_VARS > 0 /* only needed for USM */
#if NDIM > 1 /* pointers for face-centered magnetic fields */
  real, dimension(:,:,:,:), pointer :: Bx, By, Bz
#endif /* NDIM > 1      */
#endif /* NFACE_VARS > 0 */
#endif /* FLASH_USM_MHD || FLASH_UGLM_MHD*/

 
  !****************************************************************!
  ! ADDITIONAL ARGUMENTS NEEDED FOR GAUSSIAN PROCESS RECONSTRUCTION!
  ! GP stuff initialized here. DEV NOTE: For the                   !
  ! first implementation we are hardwiring 2D                      !
  ! ***************************************************************!
  real    :: shockSwitchGP
  integer :: radius, iSizeGC, jSizeGC,kSizeGC
  real, pointer, dimension(:)      :: x1, x2, x3
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits

#if NDIM < 3 
  ! Although we define this array for 1D here for the sake of compilations,
  ! GP does not support 1D calculations.
  real, pointer, dimension(:,:,:)   :: DataMultiD,DataGravMultiD
#elif NDIM == 3
  real, pointer, dimension(:,:,:,:) :: DataMultiD,DataGravMultiD
#endif

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC), target :: xCenter
  real, dimension(GRID_JHI_GC), target :: yCenter
  real, dimension(GRID_KHI_GC), target :: zCenter
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS)), target :: xCenter
  real, dimension(blkLimitsGC(HIGH,JAXIS)), target :: yCenter
  real, dimension(blkLimitsGC(HIGH,KAXIS)), target :: zCenter
#endif
  !****************** End of GP declarations **********************!

#ifdef FLASH_USM_MHD
#if NDIM == 2
  !!NAG complains if this isn't here, since Bz is uninitialized.
  nullify(Bz)
#endif /* NDIM < 3 */
#endif

  ! Set index range depending on hydro or MHD
  ! default for for hydro
  hyEndVar  = HY_ENER
  hyEndFlux = F05ENER_FLUX
#ifdef FLASH_USM_MHD /* for USM-MHD */
  hyEndVar  = HY_MAGZ
  hyEndFlux = F08MAGZ_FLUX
#elif defined(FLASH_UGLM_MHD)
  hyEndVar  = HY_GLMP
  hyEndFlux = F09GLMP_FLUX
#endif

  hyEndPrimVar = hyEndVar
  if (hy_useAuxEintEqn) hyEndVar = HY_EINT
#ifdef FLASH_UHD_3T
  hyEndVar = HY_END_VARS
#endif
#ifdef GRAVITY
  hyEndVar = HY_GRAV
#endif

  kGrav = 0
#ifdef GRAVITY
  kGrav = 1
#endif

  !!-------------------------------------------!!
  !! DEV: Gravity with HEDP does not work.     !!
  !!-------------------------------------------!!
  ! Switch for transverse update
  Trans_updateOnly(DIR_X) = TransX_updateOnly
  Trans_updateOnly(DIR_Y) = TransY_updateOnly
  Trans_updateOnly(DIR_Z) = TransZ_updateOnly

  ! dimensionality tuning parameters
  k2 = 0
  k3 = 0
  k4 = order - 1
  if (order > 3) k4 = 2 !(i.e., assume order = 3)
#ifdef FLASH_USM_MHD
  if (.NOT. hy_forceHydroLimit) k4 = min(NGUARD-3,k4)
#endif

  ! half delta t
  hdt=.5*dt
  
  ! grid deltas
  dx  =del(DIR_X)
  hdtx=0.5*dt/dx
  idx = 1./dx
#if NDIM > 1
  dy  =del(DIR_Y)
  hdty=0.5*dt/dy
  idy = 1./dy
  k2 = 1
#if NDIM == 3
  dz  =del(DIR_Z)
  hdtz=0.5*dt/dz
  idz = 1./dz
  k3 = 1
#endif
#endif

  !! Set indices for 1D stencil arrays
  !! Stencil size for reconstruction: [i-kOrder:i+kOrder]
  if (NGUARD == 4) then 
     if (order <= 2) then
        ! Second order MH method
        kOrder = 1
     elseif (order > 2) then
        ! PPM with 4 guardcells
        kOrder = 2
     endif
  elseif (NGUARD == 6) then
     if (order <=2) then
        kOrder = 1
     elseif (order > 2) then !if (order > 2) then
        ! PPM, WENO, GP with 6 guardcells
        kOrder = 3
     endif
  endif
  kOrder_orig = kOrder

  !! Array bound check for hybrid order in checking DivU
  if (ix-k4    < blkLimitsGC( LOW,IAXIS) .or. &
      ix+k4    > blkLimitsGC(HIGH,IAXIS) .or. &
      iy-k4*k2 < blkLimitsGC( LOW,JAXIS) .or. &
      iy+k4*k2 > blkLimitsGC(HIGH,JAXIS) .or. &
      iz-k4*k3 < blkLimitsGC( LOW,KAXIS) .or. &
      iz+k4*k3 > blkLimitsGC(HIGH,KAXIS)) then
     k4 = k4-1
  endif

  if (present(cellCfl)) then
     cflIn = cellCfl
  else
     cflIn = hy_cfl
  end if
  maxCfl = cflIn

  doFullSpecMs = (hy_fullSpecMsFluxHandling .AND. hy_numXN > 0 .AND. present(hy_SpcR)) 
  if (doFullSpecMs) doFullSpecMs = associated(hy_SpcR)

  !!*************************************************************!
  !! First, get block pointer for cell-centered variables        !
  !!*************************************************************!
  call Grid_getBlkPtr(blockID,U,CENTER)


  !!*************************************************************!
  !!  BEGIN ADDITIONAL IMPLEMENTATIONS FOR USM-MHD               !
  !!*************************************************************!
#ifdef FLASH_USM_MHD /* extra definitions for MHD */
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
#if NDIM == 3
  call Grid_getBlkPtr(blockID,Bz,FACEZ)
#endif /* NDIM == 3      */
#endif /* NDUM > 1       */
#endif /* NFACE_VARS > 0 */

  !! Dot product of U & B fields
  uB = dot_product(U(VELX_VAR:VELZ_VAR,ix,iy,iz),&
                   U(MAGX_VAR:MAGZ_VAR,ix,iy,iz))

  !! Multidimensional MHD terms
  if (.not. hy_forceHydroLimit) then
     aBn(:,DIR_X)= (/0. , &
          -U(MAGX_VAR,ix,iy,iz)/U(DENS_VAR,ix,iy,iz), &
          -U(MAGY_VAR,ix,iy,iz)/U(DENS_VAR,ix,iy,iz), &
          -U(MAGZ_VAR,ix,iy,iz)/U(DENS_VAR,ix,iy,iz), &
          (U(GAME_VAR,ix,iy,iz)-1.)*uB , &     
           0.                          , &
          -U(VELY_VAR,ix,iy,iz)        , &
          -U(VELZ_VAR,ix,iy,iz)          &
#ifdef FLASH_UGLM_MHD /* for GLM-MHD */
          ,0.                            &
#endif
          /)
#ifdef FLASH_UGLM_MHD /* for GLM-MHD */
     aGLM(:,DIR_X) = 0.
!!$     aGLM(HY_GLMP,DIR_X) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGX_VAR,ix,iy,iz)
     aGLM(HY_PRES,DIR_X) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGX_VAR,ix,iy,iz)
#endif
#if NDIM > 1
     aBn(:,DIR_Y)=aBn(:,DIR_X)
     aBn(HY_MAGX,DIR_Y)=-U(VELX_VAR,ix,iy,iz)
     aBn(HY_MAGY,DIR_Y)=0.
#ifdef FLASH_UGLM_MHD /* for GLM-MHD */
     aGLM(:,DIR_Y) = 0.
!!$     aGLM(HY_GLMP,DIR_Y) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGY_VAR,ix,iy,iz)
     aGLM(HY_PRES,DIR_Y) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGY_VAR,ix,iy,iz)
#endif

#if NDIM == 3
     aBn(:,DIR_Z)=aBn(:,DIR_X)
     aBn(HY_MAGX,DIR_Z)=-U(VELX_VAR,ix,iy,iz)
     aBn(HY_MAGZ,DIR_Z)=0.
#ifdef FLASH_UGLM_MHD /* for GLM-MHD */
     aGLM(:,DIR_Z) = 0.
!!$     aGLM(HY_GLMP,DIR_Z) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGZ_VAR,ix,iy,iz)
     aGLM(HY_PRES,DIR_Z) = (1.-U(GAME_VAR,ix,iy,iz))*U(MAGZ_VAR,ix,iy,iz)
#endif
#endif /* NDIM == 3 */
#endif /* NDIM > 1 */
  else
     aBn = 0.
  endif


  !! Define the derivative of the normal B field
  dnBn(:) = 0.

#if NFACE_VARS > 0
#if NDIM > 1
  if (hy_killdivb) then
     dnBn(DIR_X) = Bx(MAG_FACE_VAR,ix+1,iy,  iz  )-Bx(MAG_FACE_VAR,ix,iy,iz)
     if (NDIM >= 2) &
     dnBn(DIR_Y) = By(MAG_FACE_VAR,ix,  iy+1,iz  )-By(MAG_FACE_VAR,ix,iy,iz)
     if (NDIM == 3) &
     dnBn(DIR_Z) = Bz(MAG_FACE_VAR,ix,  iy,  iz+1)-Bz(MAG_FACE_VAR,ix,iy,iz)
  else
     dnBn(DIR_X) = .5*(U(MAGX_VAR,ix+1,iy,  iz  )-U(MAGX_VAR,ix-1,iy,iz))
     if (NDIM >= 2) &
     dnBn(DIR_Y) = .5*(U(MAGY_VAR,ix,  iy+1,iz  )-U(MAGY_VAR,ix,iy-1,iz))
     if (NDIM == 3) &
     dnBn(DIR_Z) = .5*(U(MAGZ_VAR,ix,  iy,  iz+1)-U(MAGZ_VAR,ix,iy,iz-1))
  endif
#endif/* if NDIM > 1 */
#else /* if NFACE_VARS = 0, that is, 1D MHD USM or GLM */
#ifdef FLASH_UGLM_MHD
  dnGLM(DIR_X) = .5*(U(GLMP_VAR,ix+1,iy,iz)-U(GLMP_VAR,ix-1,iy,iz))
  dnBn (DIR_X) = .5*(U(MAGX_VAR,ix+1,iy,iz)-U(MAGX_VAR,ix-1,iy,iz))
  if (NDIM >= 2) then
  dnGLM(DIR_Y) = .5*(U(GLMP_VAR,ix,iy+1,iz)-U(GLMP_VAR,ix,iy-1,iz))
  dnBn (DIR_Y) = .5*(U(MAGY_VAR,ix,iy+1,iz)-U(MAGY_VAR,ix,iy-1,iz))
  if (NDIM == 3) then
  dnGLM(DIR_Z) = .5*(U(GLMP_VAR,ix,iy,iz+1)-U(GLMP_VAR,ix,iy,iz-1))
  dnBn (DIR_Z) = .5*(U(MAGZ_VAR,ix,iy,iz+1)-U(MAGZ_VAR,ix,iy,iz-1))
  endif
  endif
#endif
#endif /* if NFACE_VARS > 0 */
#endif /* FLASH_USM_MHD || FLASH_UGLM_MHD */
  !!*************************************************************!
  !! DONE WITH ADDITIONAL IMPLEMENTATIONS FOR USM-MHD            !
  !!*************************************************************!


  !!*************************************************************!
  !! (I) Perform data reconstruction in each normal direction    !
  !!*************************************************************!

  !! PREPARE FOR GP RECONSTRUCTION
  if (order == 6) then
     ! GP radius
     radius = int(hy_radiusGP)

#if NDIM == 2
     DataMultiD => U(:,ix-radius:ix+radius,&
                       iy-radius:iy+radius,&
                       iz)
#elif NDIM == 3
     DataMultiD => U(:,ix-radius:ix+radius,&
                       iy-radius:iy+radius,&
                       iz-radius:iz+radius)
#endif

     iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     if (NDIM > 1) &
     jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     if (NDIM == 3) &
     kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
    
     call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
    if (NDIM > 1) &
     call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)
     if (NDIM == 3) &
     call Grid_getCellCoords(KAXIS,blockID, CENTER, .true.,zCenter,kSizeGC)

     x1 => xCenter(ix-radius:ix+radius)
     if (NDIM > 1) &
     x2 => yCenter(iy-radius:iy+radius)
     if (NDIM == 3) &
     x3 => zCenter(iz-radius:iz+radius)

#ifdef SHOK_VAR
#if NDIM == 2
     shockSwitchGP = maxval(U(SHOK_VAR,ix-radius:ix+radius,&
                                       iy-radius:iy+radius,&
                                       iz))
#elif NDIM == 3
     shockSwitchGP = maxval(U(SHOK_VAR,ix-radius:ix+radius,&
                                       iy-radius:iy+radius,&
                                       iz-radius:iz+radius))
#endif
#else
     shockSwitchGP = 0.0        !fallback if there is no SHOK_VAR
#endif
  endif !END OF GP PREPARATION


  !! Begin reconstruction
  DO iDim=DIR_X,NDIM
     If (Trans_updateOnly(iDim)) Then
        kOrder = 1
     Else
        kOrder = kOrder_orig
     Endif

     !! Set 1D arrays based on each direction
     select case (iDim)
     ! X-direction
     case (DIR_X)
        Data1D => U(:,ix-kOrder:ix+kOrder,iy,iz)
#ifdef GRAVITY
        DataGrav1D => ogravX(ix-kOrder:ix+kOrder,iy,iz)
#endif

     ! Y-direction
     case (DIR_Y)
        Data1D => U(:,ix,iy-kOrder:iy+kOrder,iz)
#ifdef GRAVITY
        DataGrav1D => ogravY(ix,iy-kOrder:iy+kOrder,iz)
#endif

     ! Z-direction
     case (DIR_Z)
        Data1D => U(:,ix,iy,iz-kOrder:iz+kOrder)
#ifdef GRAVITY
        DataGrav1D => ogravZ(ix,iy,iz-kOrder:iz+kOrder)
#endif
     end select

     if (doFullSpecMs) then
        SpcSigPtr => hy_SpcSig  (:,ix,iy,iz,iDim)
     else
        SpcSigPtr => SpcSigIgnore
     end if

#ifdef FLASH_USM_MHD
     !! CASE 1: MHD calls the reconstruction subroutines WITH all optional arguments for MHD 
     !!         BUT WITHOUT those for species.
     !! THE CONVENTIONAL INTERPOLATION-BASED RECONSTRUCTION, AS WELL AS GP:
     !! 2nd order MUSCL-Hancock, 3rd order PPM, 5th order WENO, and GP
     if (doFullSpecMs .AND. (.NOT. Trans_updateOnly(iDim))) then
        if (order == 2) then
           !! For second-order MH
           call MH&
                (iDim,dt,del(iDim),Data1D,DataGrav1D,&
                 FlatCoeff(iDim,ix,iy,iz), &
                 Trans_updateOnly(iDim),&
                 lambda (  1,iDim),&
                 leftEig(1,1,iDim),&
                 rghtEig(1,1,iDim),&
                 Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
                 dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
                 Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
                 Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
                 SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        else if (order == 3) then
           call PPM&
                (iDim,dt,del(iDim),Data1D,DataGrav1D,&
                 FlatCoeff(iDim,ix,iy,iz), &
                 Trans_updateOnly(iDim),&
                 lambda (  1,iDim),&
                 leftEig(1,1,iDim),&
                 rghtEig(1,1,iDim),&
                 Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
                 dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
                 Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
                 Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
                 SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        else if (order == 5) then
           !! wenoMethod = "WENO5" (or "weno5") will call weno5 by Jiang and Shu.
           !! wenoMethod = "WENOZ" (or "wenoz") will call weno-Z by Borges et al.
           call WENO&
                (hy_wenoMethod,iDim,dt,del(iDim),Data1D,DataGrav1D,&
                 FlatCoeff(iDim,ix,iy,iz), &
                 Trans_updateOnly(iDim),&
                 lambda (  1,iDim),&
                 leftEig(1,1,iDim),&
                 rghtEig(1,1,iDim),&
                 Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
                 dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
                 Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
                 Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
                 SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        else if (order == 6) then ! GP
           ! check if we are near shocks:
           ! (1) use GP if not near shocks when shockDetect = .true., or
           ! (2) use GP everywhere when shockDetect = .false.
           IF (.not. hy_shockDetectOn .or. &
                (hy_shockDetectOn .and. shockSwitchGP == 0.)) THEN
#ifdef GPRO_VAR
              U(GPRO_VAR,ix,iy,iz) = 1.
#endif
              call GP&
                (iDim,dt,del(iDim),DataMultiD,DataGravMultiD,&
                 x1,x2,x3,radius,&
                 FlatCoeff(iDim,ix,iy,iz), &
                 Trans_updateOnly(iDim),&
                 lambda (  1,iDim),&
                 leftEig(1,1,iDim),&
                 rghtEig(1,1,iDim),&
                 Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
                 dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
                 Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
                 Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
                 SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))
           ELSE !use PPM instead near shocks
              call PPM&
                (iDim,dt,del(iDim),Data1D,DataGrav1D,&
                 FlatCoeff(iDim,ix,iy,iz), &
                 Trans_updateOnly(iDim),&
                 lambda (  1,iDim),&
                 leftEig(1,1,iDim),&
                 rghtEig(1,1,iDim),&
                 Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
                 dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
                 Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
                 Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
                 SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))
           ENDIF
        endif ! if (order = 2, 3, 5, or 6)
     elseif (order == 2) then
        !! For second-order MH
        call MH&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
              SpcSig=SpcSigPtr)

     elseif (order == 3) then
        call PPM&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
              SpcSig=SpcSigPtr)

     elseif (order == 5) then
        !! wenoMethod = "WENO5" (or "weno5") will call weno5 by Jiang and Shu
        !! wenoMethod = "WENOZ" (or "wenoz") will call weno-Z by Borges et al.
        call WENO&
             (hy_wenoMethod,iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
              SpcSig=SpcSigPtr)

     elseif (order == 6) then ! GP
        ! check if we are near shocks:
        ! (1) use GP if not near shocks when shockDetect = .true., or
        ! (2) use GP everywhere when shockDetect = .false.
        IF (.not. hy_shockDetectOn .or. &
            (hy_shockDetectOn .and. shockSwitchGP == 0.)) THEN
#ifdef GPRO_VAR
           U(GPRO_VAR,ix,iy,iz) = 1.
#endif
           call GP&
             (iDim,dt,del(iDim),DataMultiD,DataGravMultiD,&
              x1,x2,x3,radius,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
              SpcSig=SpcSigPtr)
        ELSE !use PPM instead near shocks
           call PPM&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              dnBnormal=dnBn(iDim),aBnormal=aBn(1,iDim),&
              SpcSig=SpcSigPtr)
        ENDIF
     endif ! if (order = 2, 3, 5, or 6)
#endif

#ifdef FLASH_UHD_HYDRO
     if (doFullSpecMs .AND. (.NOT. Trans_updateOnly(iDim))) then
     !! CASE 2: Hydro calls the reconstruction subroutines WITHOUT 
     !!         those optional arguments for MHD, and WITH those for species (e.g., hy_Spc*).
     !!         This is a *default* setup without requiring to use scratch variables.
     !! THE CONVENTIONAL INTERPOLATION-BASED RECONSTRUCTION, AS WELL AS GP:
     !! 2nd order MUSCL-Hancock, 3rd order PPM, 5th order WENO, and GP
        if (order == 2) then
        !! For second-order MH
           call MH&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
              Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
              SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        elseif (order == 3) then
           call PPM&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
              Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
              SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        elseif (order == 5) then
        !! wenoMethod = "WENO5" (or "weno5") will call weno5 by Jiang and Shu
        !! wenoMethod = "WENOZ" (or "wenoz") will call weno-Z by Borges et al.
           call WENO&
             (hy_wenoMethod,iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
              Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
              SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))

        elseif (order == 6) then ! GP
        ! check if we are near shocks:
        ! (1) use GP if not near shocks when shockDetect = .true., or
        ! (2) use GP everywhere when shockDetect = .false.
           IF (.not. hy_shockDetectOn .or. &
            (hy_shockDetectOn .and. shockSwitchGP == 0.)) THEN
#ifdef GPRO_VAR
              U(GPRO_VAR,ix,iy,iz) = 1.
#endif
              call GP&
                   (iDim,dt,del(iDim),DataMultiD,DataGravMultiD,&
              x1,x2,x3,radius,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
              Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
              SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))
           ELSE !use PPM instead near shocks
              call PPM&
                   (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              Sr    =hy_SpcR  (:,ix,iy,iz,iDim),&
              Sl    =hy_SpcL  (:,ix,iy,iz,iDim),&
              SpcSig=hy_SpcSig(:,ix,iy,iz,iDim))
           ENDIF

        endif ! if (order = 2, 3, 5, or 6)

     else
     !! CASE 3: Hydro calls the reconstruction subroutines WITHOUT any of those optional arguments.
     !!         This is a *non-default* setup requiring to use scratch variables.
     !! THE CONVENTIONAL INTERPOLATION-BASED RECONSTRUCTION, AS WELL AS GP:
     !! 2nd order MUSCL-Hancock, 3rd order PPM, 5th order WENO, and GP
        if (order == 2) then
        !! For second-order MH
           call MH&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              SpcSig=SpcSigPtr)

        elseif (order == 3) then
           call PPM&
             (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              SpcSig=SpcSigPtr)

        elseif (order == 5) then
        !! wenoMethod = "WENO5" (or "weno5") will call weno5 by Jiang and Shu
        !! wenoMethod = "WENOZ" (or "wenoz") will call weno-Z by Borges et al.
           call WENO&
             (hy_wenoMethod,iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              SpcSig=SpcSigPtr)

        elseif (order == 6) then ! GP
        ! check if we are near shocks:
        ! (1) use GP if not near shocks when shockDetect = .true., or
        ! (2) use GP everywhere when shockDetect = .false.
           IF (.not. hy_shockDetectOn .or. &
            (hy_shockDetectOn .and. shockSwitchGP == 0.)) THEN
#ifdef GPRO_VAR
              U(GPRO_VAR,ix,iy,iz) = 1.
#endif
              call GP&
                   (iDim,dt,del(iDim),DataMultiD,DataGravMultiD,&
              x1,x2,x3,radius,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              SpcSig=SpcSigPtr)
           ELSE !use PPM instead near shocks
              call PPM&
                   (iDim,dt,del(iDim),Data1D,DataGrav1D,&
              FlatCoeff(iDim,ix,iy,iz), &
              Trans_updateOnly(iDim),&
              lambda (  1,iDim),&
              leftEig(1,1,iDim),&
              rghtEig(1,1,iDim),&
              Wp(1,iDim),Wn(1,iDim),sig(1,iDim),&
              SpcSig=SpcSigPtr)
           ENDIF

        endif ! if (order = 2, 3, 5, or 6)
     end if
#endif

  ENDDO ! do iDim = DIR_X, NDIM



  !! *************************************************************
  !! (Ia) Avoid any possible negative states:                    *
  !!   - Use first-order Godunov scheme if they occur            *
  !! *************************************************************
  if (hy_useHybridOrder) then
     ! kappa = 0. is default; kappa = 0.4 is what Balsara uses.
     ! kappa = -0.1 for mh gives almost same bw shock profile without hybrid order
     ! kappa = -0.1 for ppm more overshoots than without it; -0.05 is bit more overshooting but almost identical
     !! epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)
!!$     epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)
     Cs = U(GAMC_VAR,ix,iy,iz)*U(PRES_VAR,ix,iy,iz)
#ifdef FLASH_USM_MHD
     Cs = Cs + dot_product(U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),U(MAGX_VAR:MAGZ_VAR,ix,iy,iz))
#endif
     Cs = sqrt(Cs/U(DENS_VAR,ix,iy,iz))
     epsilon = hy_hybridOrderKappa*Cs
          

     if (minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < epsilon ) then
        DO iDim = DIR_X,NDIM
           call hy_uhd_checkRHjumpCond(iDim,U(DENS_VAR,ix,iy,iz),&
                                            U(VELX_VAR:VELZ_VAR,ix,iy,iz),&
                                            U(PRES_VAR,ix,iy,iz),&
!!$#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
!!$                                            U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),&
!!$#ifdef FLASH_UGLM_MHD
!!$                                            U(GLMP_VAR,ix,iy,iz),&
!!$#endif
!!$#endif
                                            U(GAMC_VAR,ix,iy,iz),&
                                            Wp(HY_DENS:HY_EINT,iDim),&
                                            Wn(HY_DENS:HY_EINT,iDim),&
                                            SWp(iDim),SWn(iDim))
           if (SWp(iDim) .or. SWn(iDim)) then

              if (.NOT. Trans_updateOnly(iDim)) then
#ifdef DEBUG_HYDRO
                 print*,'[dR1St]h fallback to order 1 for iDim=',iDim,'at ix,iy=',ix,iy,' in Block',blockID,'@',hy_meshMe
                 print*,'[dR1St]h',Wn(HY_DENS,iDim),Wp(HY_DENS,iDim), &
                      U(DENS_VAR,ix,iy,iz), U(PRES_VAR,ix,iy,iz), &
                      Cs, U(GAMC_VAR,ix,iy,iz),SWn(iDim),SWp(iDim)
#endif
                 if (hy_fallbackLowerCFL) maxCfl = min(maxCfl,hy_cflFallbackFactor/real(NDIM))
#ifdef SHOK_VAR
                 ! Add 3 to the shock variable if RH conditions are violated
                 ! to ensure that HLL is called if RiemannSolver=HYBRID - JW
                 U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) = &
                   U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) + 3.0
#endif
              end if

              ! First-order for the right state
              Wp(HY_DENS:HY_EINT,iDim) = (/U(DENS_VAR,ix,iy,iz),&
                                           U(VELX_VAR:VELZ_VAR,ix,iy,iz),&
                                           U(PRES_VAR,ix,iy,iz),&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                                           U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),&
#ifdef FLASH_UGLM_MHD
                                           U(GLMP_VAR,ix,iy,iz),&
#endif
#endif
                                           U(GAMC_VAR:GAME_VAR,ix,iy,iz),&
                                           U(EINT_VAR,ix,iy,iz)/)
              ! Same for the left state
              Wn(HY_DENS:HY_EINT,iDim) = Wp(HY_DENS:HY_EINT,iDim)
           endif
        ENDDO !do iDim

     endif !(minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < epsilon )
  endif !if (hy_useHybridOrder) then






#if NDIM > 1
  !! *************************************************************
  !! (II) Perform to make contributions from transverse fluxes   *
  !! *************************************************************
  !! This is needed only for 2D & 3D
  !! For 2D
#ifdef FLASH_UGLM_MHD
  sig(HY_GLMP,:) = 0.
#endif

  Wn(HY_DENS:HY_END_VARS,DIR_X) = Wn(HY_DENS:HY_END_VARS,DIR_X) - hdty*sig(HY_DENS:HY_END_VARS,DIR_Y)
  Wp(HY_DENS:HY_END_VARS,DIR_X) = Wp(HY_DENS:HY_END_VARS,DIR_X) - hdty*sig(HY_DENS:HY_END_VARS,DIR_Y)

  Wn(HY_DENS:HY_END_VARS,DIR_Y) = Wn(HY_DENS:HY_END_VARS,DIR_Y) - hdtx*sig(HY_DENS:HY_END_VARS,DIR_X)
  Wp(HY_DENS:HY_END_VARS,DIR_Y) = Wp(HY_DENS:HY_END_VARS,DIR_Y) - hdtx*sig(HY_DENS:HY_END_VARS,DIR_X)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) 
  ! Need to add multidimensional MHD terms in MHD calculations
  Wn(HY_DENS:hyEndPrimVar,DIR_X) = Wn(HY_DENS:hyEndPrimVar,DIR_X) - hdty*aBn(HY_DENS:hyEndPrimVar,DIR_Y)*dnBn(DIR_Y)
  Wp(HY_DENS:hyEndPrimVar,DIR_X) = Wp(HY_DENS:hyEndPrimVar,DIR_X) - hdty*aBn(HY_DENS:hyEndPrimVar,DIR_Y)*dnBn(DIR_Y)
  Wn(HY_DENS:hyEndPrimVar,DIR_Y) = Wn(HY_DENS:hyEndPrimVar,DIR_Y) - hdtx*aBn(HY_DENS:hyEndPrimVar,DIR_X)*dnBn(DIR_X)
  Wp(HY_DENS:hyEndPrimVar,DIR_Y) = Wp(HY_DENS:hyEndPrimVar,DIR_Y) - hdtx*aBn(HY_DENS:hyEndPrimVar,DIR_X)*dnBn(DIR_X)
#ifdef FLASH_UGLM_MHD
  Wn(HY_DENS:hyEndPrimVar,DIR_X) = Wn(HY_DENS:hyEndPrimVar,DIR_X) - hdty*aGLM(HY_DENS:hyEndPrimVar,DIR_Y)*dnGLM(DIR_Y)
  Wp(HY_DENS:hyEndPrimVar,DIR_X) = Wp(HY_DENS:hyEndPrimVar,DIR_X) - hdty*aGLM(HY_DENS:hyEndPrimVar,DIR_Y)*dnGLM(DIR_Y)
  Wn(HY_DENS:hyEndPrimVar,DIR_Y) = Wn(HY_DENS:hyEndPrimVar,DIR_Y) - hdtx*aGLM(HY_DENS:hyEndPrimVar,DIR_X)*dnGLM(DIR_X)
  Wp(HY_DENS:hyEndPrimVar,DIR_Y) = Wp(HY_DENS:hyEndPrimVar,DIR_Y) - hdtx*aGLM(HY_DENS:hyEndPrimVar,DIR_X)*dnGLM(DIR_X)
#endif
#endif

#if NDIM > 2
  Wn(HY_DENS:HY_END_VARS,DIR_X) = Wn(HY_DENS:HY_END_VARS,DIR_X) - hdtz*sig(HY_DENS:HY_END_VARS,DIR_Z)
  Wp(HY_DENS:HY_END_VARS,DIR_X) = Wp(HY_DENS:HY_END_VARS,DIR_X) - hdtz*sig(HY_DENS:HY_END_VARS,DIR_Z)

  Wn(HY_DENS:HY_END_VARS,DIR_Y) = Wn(HY_DENS:HY_END_VARS,DIR_Y) - hdtz*sig(HY_DENS:HY_END_VARS,DIR_Z)
  Wp(HY_DENS:HY_END_VARS,DIR_Y) = Wp(HY_DENS:HY_END_VARS,DIR_Y) - hdtz*sig(HY_DENS:HY_END_VARS,DIR_Z)

  Wn(HY_DENS:HY_END_VARS,DIR_Z) = Wn(HY_DENS:HY_END_VARS,DIR_Z) - hdtx*sig(HY_DENS:HY_END_VARS,DIR_X)&
                                                                - hdty*sig(HY_DENS:HY_END_VARS,DIR_Y)

  Wp(HY_DENS:HY_END_VARS,DIR_Z) = Wp(HY_DENS:HY_END_VARS,DIR_Z) - hdtx*sig(HY_DENS:HY_END_VARS,DIR_X)&
                                                                - hdty*sig(HY_DENS:HY_END_VARS,DIR_Y)

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  ! Need to add multidimensional MHD terms in MHD calculations
  Wn(HY_DENS:hyEndPrimVar,DIR_X) = Wn(HY_DENS:hyEndPrimVar,DIR_X) - hdtz*aBn(HY_DENS:hyEndPrimVar,DIR_Z)*dnBn(DIR_Z)
  Wp(HY_DENS:hyEndPrimVar,DIR_X) = Wp(HY_DENS:hyEndPrimVar,DIR_X) - hdtz*aBn(HY_DENS:hyEndPrimVar,DIR_Z)*dnBn(DIR_Z)

  Wn(HY_DENS:hyEndPrimVar,DIR_Y) = Wn(HY_DENS:hyEndPrimVar,DIR_Y) - hdtz*aBn(HY_DENS:hyEndPrimVar,DIR_Z)*dnBn(DIR_Z)
  Wp(HY_DENS:hyEndPrimVar,DIR_Y) = Wp(HY_DENS:hyEndPrimVar,DIR_Y) - hdtz*aBn(HY_DENS:hyEndPrimVar,DIR_Z)*dnBn(DIR_Z)

  Wn(HY_DENS:hyEndPrimVar,DIR_Z) = Wn(HY_DENS:hyEndPrimVar,DIR_Z) - hdtx*aBn(HY_DENS:hyEndPrimVar,DIR_X)*dnBn(DIR_X)&
                                                                  - hdty*aBn(HY_DENS:hyEndPrimVar,DIR_Y)*dnBn(DIR_Y)
  Wp(HY_DENS:hyEndPrimVar,DIR_Z) = Wp(HY_DENS:hyEndPrimVar,DIR_Z) - hdtx*aBn(HY_DENS:hyEndPrimVar,DIR_X)*dnBn(DIR_X)&
                                                                  - hdty*aBn(HY_DENS:hyEndPrimVar,DIR_Y)*dnBn(DIR_Y)

  ! corrections for 3D MHD source terms
  sig(HY_DENS:hyEndPrimVar,DIR_X) = sig(HY_DENS:hyEndPrimVar,DIR_X) + aBn(HY_DENS:hyEndPrimVar,DIR_X)*dnBn(DIR_X)
  sig(HY_DENS:hyEndPrimVar,DIR_Y) = sig(HY_DENS:hyEndPrimVar,DIR_Y) + aBn(HY_DENS:hyEndPrimVar,DIR_Y)*dnBn(DIR_Y)
  sig(HY_DENS:hyEndPrimVar,DIR_Z) = sig(HY_DENS:hyEndPrimVar,DIR_Z) + aBn(HY_DENS:hyEndPrimVar,DIR_Z)*dnBn(DIR_Z)


#ifdef FLASH_UGLM_MHD
  Wn(HY_DENS:hyEndPrimVar,DIR_X) = Wn(HY_DENS:hyEndPrimVar,DIR_X) - hdtz*aGLM(HY_DENS:hyEndPrimVar,DIR_Z)*dnGLM(DIR_Z)
  Wp(HY_DENS:hyEndPrimVar,DIR_X) = Wp(HY_DENS:hyEndPrimVar,DIR_X) - hdtz*aGLM(HY_DENS:hyEndPrimVar,DIR_Z)*dnGLM(DIR_Z)

  Wn(HY_DENS:hyEndPrimVar,DIR_Y) = Wn(HY_DENS:hyEndPrimVar,DIR_Y) - hdtz*aGLM(HY_DENS:hyEndPrimVar,DIR_Z)*dnGLM(DIR_Z)
  Wp(HY_DENS:hyEndPrimVar,DIR_Y) = Wp(HY_DENS:hyEndPrimVar,DIR_Y) - hdtz*aGLM(HY_DENS:hyEndPrimVar,DIR_Z)*dnGLM(DIR_Z)

  Wn(HY_DENS:hyEndPrimVar,DIR_Z) = Wn(HY_DENS:hyEndPrimVar,DIR_Z) - hdtx*aGLM(HY_DENS:hyEndPrimVar,DIR_X)*dnGLM(DIR_X)
                                                                  - hdty*aGLM(HY_DENS:hyEndPrimVar,DIR_Y)*dnGLM(DIR_Y)
  Wp(HY_DENS:hyEndPrimVar,DIR_Z) = Wp(HY_DENS:hyEndPrimVar,DIR_Z) - hdtx*aGLM(HY_DENS:hyEndPrimVar,DIR_X)*dnGLM(DIR_X)
                                                                  - hdty*aGLM(HY_DENS:hyEndPrimVar,DIR_Y)*dnGLM(DIR_Y)

  ! corrections for 3D MHD source terms
  sig(HY_DENS:hyEndPrimVar,DIR_X) = sig(HY_DENS:hyEndPrimVar,DIR_X) + aGLM(HY_DENS:hyEndPrimVar,DIR_X)*dnGLM(DIR_X)
  sig(HY_DENS:hyEndPrimVar,DIR_Y) = sig(HY_DENS:hyEndPrimVar,DIR_Y) + aGLM(HY_DENS:hyEndPrimVar,DIR_Y)*dnGLM(DIR_Y)
  sig(HY_DENS:hyEndPrimVar,DIR_Z) = sig(HY_DENS:hyEndPrimVar,DIR_Z) + aGLM(HY_DENS:hyEndPrimVar,DIR_Z)*dnGLM(DIR_Z)
#endif

#endif /* FLASH_USM_MHD || FLASH_UGLM_MHD*/

#endif /* NDIM > 2 */
#endif /* NDIM > 1 */


#if (NSPECIES+NMASS_SCALARS) > 0
#if NDIM > 1
  if (hy_fullSpecMsFluxHandling) then
     if (.NOT. Trans_updateOnly(DIR_X)) then
     hy_SpcL(:,ix,iy,iz,DIR_X) = hy_SpcL(:,ix,iy,iz,DIR_X) - hdty*hy_SpcSig(:,ix,iy,iz,DIR_Y)
     hy_SpcR(:,ix,iy,iz,DIR_X) = hy_SpcR(:,ix,iy,iz,DIR_X) - hdty*hy_SpcSig(:,ix,iy,iz,DIR_Y)
     end if

     if (.NOT. Trans_updateOnly(DIR_Y)) then
     hy_SpcL(:,ix,iy,iz,DIR_Y) = hy_SpcL(:,ix,iy,iz,DIR_Y) - hdtx*hy_SpcSig(:,ix,iy,iz,DIR_X)
     hy_SpcR(:,ix,iy,iz,DIR_Y) = hy_SpcR(:,ix,iy,iz,DIR_Y) - hdtx*hy_SpcSig(:,ix,iy,iz,DIR_X)
     end if

#if NDIM > 2
     if (.NOT. Trans_updateOnly(DIR_X)) then
     hy_SpcL(:,ix,iy,iz,DIR_X) = hy_SpcL(:,ix,iy,iz,DIR_X) - hdtz*hy_SpcSig(:,ix,iy,iz,DIR_Z)
     hy_SpcR(:,ix,iy,iz,DIR_X) = hy_SpcR(:,ix,iy,iz,DIR_X) - hdtz*hy_SpcSig(:,ix,iy,iz,DIR_Z)
     end if

     if (.NOT. Trans_updateOnly(DIR_Y)) then
     hy_SpcL(:,ix,iy,iz,DIR_Y) = hy_SpcL(:,ix,iy,iz,DIR_Y) - hdtz*hy_SpcSig(:,ix,iy,iz,DIR_Z)
     hy_SpcR(:,ix,iy,iz,DIR_Y) = hy_SpcR(:,ix,iy,iz,DIR_Y) - hdtz*hy_SpcSig(:,ix,iy,iz,DIR_Z)
     end if

     if (.NOT. Trans_updateOnly(DIR_Z)) then
     hy_SpcL(:,ix,iy,iz,DIR_Z) = hy_SpcL(:,ix,iy,iz,DIR_Z) - hdtx*hy_SpcSig(:,ix,iy,iz,DIR_X)&
                                                           - hdty*hy_SpcSig(:,ix,iy,iz,DIR_Y)

     hy_SpcR(:,ix,iy,iz,DIR_Z) = hy_SpcR(:,ix,iy,iz,DIR_Z) - hdtx*hy_SpcSig(:,ix,iy,iz,DIR_X)&
                                                           - hdty*hy_SpcSig(:,ix,iy,iz,DIR_Y)
     end if
#endif /* NDIM > 2 */
  end if
#endif /* NDIM > 1 */
#endif /* (NSPECIES+NMASS_SCALARS) > 0 */



#ifdef FLASH_USM_MHD
#if NFACE_VARS == 0 /* if NFACE_VAR = 0 */
  ! 1D MHD or pure hydro mode
  Wn(HY_MAGX,DIR_X) = U(MAGX_VAR,ix,iy,iz)
  Wp(HY_MAGX,DIR_X) = U(MAGX_VAR,ix,iy,iz)
#if NDIM > 1
  Wn(HY_MAGY,DIR_Y) = U(MAGY_VAR,ix,iy,iz)
  Wp(HY_MAGY,DIR_Y) = U(MAGY_VAR,ix,iy,iz)
#if NDIM == 3
  Wn(HY_MAGZ,DIR_Z) = U(MAGZ_VAR,ix,iy,iz)
  Wp(HY_MAGZ,DIR_Z) = U(MAGZ_VAR,ix,iy,iz)
#endif /* NDIM = 3 */
#endif /* NDIM > 1 */
#endif /* NFACE_VARS == 0 */
#endif /* FLASH_USM_MHD */

  !! *************************************************************
  !! (Ia) Avoid any possible negative states:                    *
  !!   - Use first-order Godunov scheme if they occur            *
  !! *************************************************************
#if NDIM > 1
  if (hy_useHybridOrder) then
     ! kappa = 0. is default; kappa = 0.4 is what Balsara uses.
     ! kappa = -0.1 for mh gives almost same bw shock profile without hybrid order
     ! kappa = -0.1 for ppm more overshoots than without it;
     ! kappa = -0.05 is bit more overshooting but almost identical
     !! epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)
     Cs = U(GAMC_VAR,ix,iy,iz)*U(PRES_VAR,ix,iy,iz)
#ifdef FLASH_USM_MHD
     Cs = Cs + dot_product(U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),U(MAGX_VAR:MAGZ_VAR,ix,iy,iz))
#endif
     Cs = sqrt(Cs/U(DENS_VAR,ix,iy,iz))
     epsilon = hy_hybridOrderKappa*Cs

     if (minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < epsilon ) then
        DO iDim = DIR_X,NDIM
           call hy_uhd_checkRHjumpCond(iDim,U(DENS_VAR,ix,iy,iz),&
                                            U(VELX_VAR:VELZ_VAR,ix,iy,iz),&
                                            U(PRES_VAR,ix,iy,iz),&
                                            U(GAMC_VAR,ix,iy,iz),&
                                            Wp(HY_DENS:HY_EINT,iDim),&
                                            Wn(HY_DENS:HY_EINT,iDim),&
                                            SWp(iDim),SWn(iDim))
           if (SWp(iDim) .or. SWn(iDim)) then
              ! First-order for the right state
#ifdef DEBUG_HYDRO
              if (.NOT. Trans_updateOnly(iDim)) then
                 print*,'[dR1St]H fallback to order 1 for iDim=',iDim,'at ix,iy=',ix,iy,' in Block',blockID,'@',hy_meshMe
                 print*,'[dR1St]H',Wn(HY_DENS,iDim),Wp(HY_DENS,iDim), &
                      U(DENS_VAR,ix,iy,iz), U(PRES_VAR,ix,iy,iz), &
                      Cs, U(GAMC_VAR,ix,iy,iz),SWn(iDim),SWp(iDim)
#ifdef SHOK_VAR
                 ! Add 3 to the shock variable if RH conditions are violated
                 ! to ensure that HLL is called if RiemannSolver=HYBRID - JW
                 U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) = &
                   U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) + 3.0
#endif
              end if
#endif

              Wp(HY_DENS:HY_EINT,iDim) = (/U(DENS_VAR,ix,iy,iz),&
                                           U(VELX_VAR:VELZ_VAR,ix,iy,iz),&
                                           U(PRES_VAR,ix,iy,iz),&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                                           U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),&
#ifdef FLASH_UGLM_MHD
                                           U(GLMP_VAR,ix,iy,iz),&
#endif
#endif
                                           U(GAMC_VAR:GAME_VAR,ix,iy,iz),&
                                           U(EINT_VAR,ix,iy,iz)/)
              ! Same for the left state
              Wn(HY_DENS:HY_EINT,iDim) = Wp(HY_DENS:HY_EINT,iDim)
           endif
        ENDDO !do iDim


        !! Lower CFL if needed to reduce to first-order (i.e., donor-cell)
        if ((SWp(DIR_X) .or. SWn(DIR_X)  .or. &
             SWp(DIR_Y) .or. SWn(DIR_Y)) .and. maxCfl > 0.45) &
             maxCfl = 0.45
#if NDIM > 2
        if ((SWp(DIR_X) .or. SWn(DIR_X)  .or. &
             SWp(DIR_Y) .or. SWn(DIR_Y)  .or. &
             SWp(DIR_Z) .or. SWn(DIR_Z)) .and. maxCfl > 0.25) &
             maxCfl = 0.25
#endif
     endif !if maxval(divv) < epsilon
  endif !if hy_useHybridOrder
#endif /*if NDIM > 1 */

!!$  if (maxCFL < cflIn) then
!!$     print*,'Would set hy_cfl = maxCfl:',hy_cfl,maxCfl
!!$  end if
  if (present(cellCfl)) cellCfl = maxCfl


  !! *************************************************************************
  !! (IV) The last check for negative pressure and density in reconstruction *
  !! *************************************************************************
  do iDim = DIR_X,NDIM
     if (Wn(HY_DENS,iDim) < 0. .or. Wp(HY_DENS,iDim) < 0. .or. &
         Wn(HY_PRES,iDim) < 0. .or. Wp(HY_PRES,iDim) < 0. ) then

        if (.NOT. Trans_updateOnly(iDim)) then
#ifdef DEBUG_HYDRO_POSITIVITY
           print*,'[dR1St] fallback to order 1 for iDim=',iDim,'at ix,iy=',ix,iy,' in Block',blockID,'@',hy_meshMe
           print*,'[dR1St]',Wn(HY_DENS,iDim),Wp(HY_DENS,iDim), &
                Wn(HY_PRES,iDim), Wp(HY_PRES,iDim)
#endif
           if (present(cellCfl)) cellCfl = min(cellCfl,hy_cflFallbackFactor/real(NDIM))
#ifdef SHOK_VAR
           ! Add 2 to the shock variable if RH conditions are violated
           ! to ensure that HLL is called if RiemannSolver=HYBRID - JW
           U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) = &
              U(SHOK_VAR,ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3) + 2.0
#endif
        end if

        ! First-order for the right state
        Wn(HY_DENS:HY_END_VARS-kGrav,iDim) = (/U(DENS_VAR,ix,iy,iz)&
                                              ,U(VELX_VAR:VELZ_VAR,ix,iy,iz)&
                                              ,U(PRES_VAR,ix,iy,iz)&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                                              ,U(MAGX_VAR:MAGZ_VAR,ix,iy,iz)&
#ifdef FLASH_UGLM_MHD
                                              ,U(GLMP_VAR,ix,iy,iz)&
#endif
#endif
                                              ,U(GAMC_VAR:GAME_VAR,ix,iy,iz)&
                                              ,U(EINT_VAR,ix,iy,iz)&
#ifdef FLASH_UHD_3T
                                              ,U(EELE_VAR,ix,iy,iz)&
                                              ,U(EION_VAR,ix,iy,iz)&
                                              ,U(ERAD_VAR,ix,iy,iz)&
#endif
                                              /)

        ! Same for the left state
        Wp(HY_DENS:HY_END_VARS-kGrav,iDim) = Wn(HY_DENS:HY_END_VARS-kGrav,iDim)

     endif
  enddo

#ifdef FLASH_UHD_3T
  !! *************************************************************************
  !! (V) ADDED LAST check for negative energy components !!!
  !! *************************************************************************
  do iDim = DIR_X,NDIM
     if (Wn(HY_EION,iDim) < 0. .or. Wp(HY_EION,iDim) < 0. .or. &
         Wn(HY_EELE,iDim) < 0. .or. Wp(HY_EELE,iDim) < 0. .or. &
         Wn(HY_ERAD,iDim) < 0. .or. Wp(HY_ERAD,iDim) < 0. ) then

        if (.NOT. Trans_updateOnly(iDim)) then
           print*,'[dR1St]+ fallback to order 1 for iDim=',iDim,'at ix,iy=',ix,iy,' in Block',blockID,'@',hy_meshMe
           print*,'[dR1St]+',Wn(HY_EION,iDim),Wp(HY_EION,iDim), &
                Wn(HY_EELE,iDim), Wp(HY_EELE,iDim), &
                Wn(HY_ERAD,iDim), Wp(HY_ERAD,iDim)
           if (present(cellCfl)) then
              if (hy_fallbackLowerCFL) cellCfl = min(cellCfl,hy_cflFallbackFactor/real(NDIM))
           end if
        end if

        ! First-order for the right state
        Wn(HY_DENS:HY_END_VARS-kGrav,iDim) = (/U(DENS_VAR,ix,iy,iz)&
                                              ,U(VELX_VAR:VELZ_VAR,ix,iy,iz)&
                                              ,U(PRES_VAR,ix,iy,iz)&
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
                                              ,U(MAGX_VAR:MAGZ_VAR,ix,iy,iz)&
#ifdef FLASH_UGLM_MHD
                                              ,U(GLMP_VAR,ix,iy,iz)&
#endif
#endif
                                              ,U(GAMC_VAR:GAME_VAR,ix,iy,iz)&
                                              ,U(EINT_VAR,ix,iy,iz)&
                                              ,U(EELE_VAR,ix,iy,iz)&
                                              ,U(EION_VAR,ix,iy,iz)&
                                              ,U(ERAD_VAR,ix,iy,iz)&
                                              /)

        ! Same for the left state
        Wp(HY_DENS:HY_END_VARS-kGrav,iDim) = Wn(HY_DENS:HY_END_VARS-kGrav,iDim)

     endif
  enddo
#endif

  call Grid_releaseBlkPtr(blockID,U,CENTER)

#ifdef FLASH_USM_MHD /* extra definitions for MHD */
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
#endif /* NDIM == 3      */
#endif /* NDUM > 1       */
#endif /* NFACE_VARS > 0 */
#endif

End Subroutine hy_uhd_dataReconstOneStep
