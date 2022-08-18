!!****ih* source/physics/Hydro/HydroMain/unsplit/hy_uhd_interface
!!
!! NAME
!!  hy_uhd_interface
!!
!! SYNOPSIS
!!  use hy_uhd_interface
!!
!! DESCRIPTION
!!  This is an interface specific for the two unsplit solvers:
!!   1. MHD_StaggeredMesh
!!   2. Hydro_Unsplit
!!  This interface defines its public interfaces.
!!
!!***

!!REORDER(4): [xyz]flux, faceFlux, fl[xyz], Sp, U

Module hy_uhd_interface

  implicit none


#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  interface
     subroutine hy_uhd_avgState(sweepDir,VL,VR,Vavg)
       implicit none
       integer, intent(IN) :: sweepDir
       real, dimension(HY_VARINUM3), intent(IN)  :: VL,VR
       real, dimension(HY_VARINUM2), intent(OUT) :: Vavg
     end subroutine hy_uhd_avgState
  end interface



  interface
     subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                       ogravX,ogravY,ogravZ,&
                                       scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                       hy_SpcR,&
                                       hy_SpcL,hy_SpcSig,normalFieldUpdate)

       implicit none
       integer, intent(IN)   :: blockID
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN)   :: dt
       real,    intent(IN),dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) &
            :: ogravX,ogravY,ogravZ
#else
       real, dimension(blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: ogravX,ogravY,ogravZ
#endif
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, intent(IN), optional :: normalFieldUpdate
     end subroutine hy_uhd_getRiemannState
 end interface



  interface
     subroutine hy_uhd_entropyFix(lambda,lambdaL,lambdaR)
       implicit none
       real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
       real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
     end subroutine hy_uhd_entropyFix
  end interface



  interface
     subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC,datasize,del,&
                                     xflux, yflux,zflux,&
                                     scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                                     hy_SpcR,hy_SpcL,hy_SpcSig,lastCall)
       implicit none
       integer, intent(IN)  :: blockID
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       integer, dimension(MDIM), intent(IN)         :: datasize
       real,    dimension(MDIM), intent(IN)         :: del
#ifdef FIXEDBLOCKSIZE
       real, intent(OUT) :: xflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)
       real, intent(OUT) :: yflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)
       real, intent(OUT) :: zflux (NFLUXES, &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC)                
#else
       real, intent(OUT) :: xflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real, intent(OUT) :: yflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
       real, intent(OUT) :: zflux(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
#endif
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, optional, intent(IN) :: lastCall
     end subroutine hy_uhd_getFaceFlux
  end interface



  interface
     Subroutine hy_uhd_dataReconstOneStep(blockID,blkLimitsGC,order,ix,iy,iz, &
                                          dt,del,ogravX,ogravY,ogravZ,&
                                          DivU,FlatCoeff,   &
                                          TransX_updateOnly,&
                                          TransY_updateOnly,&
                                          TransZ_updateOnly,&
                                          Wp,Wn,sig,&
                                          lambda,leftEig,rghtEig,&
                                          cellCfl,&
                                          hy_SpcR,hy_SpcL,hy_SpcSig)

       implicit none
       integer,intent(IN) :: blockID
       integer,intent(IN),dimension(LOW:HIGH,MDIM):: blkLimitsGC
       integer,intent(IN) :: order,ix,iy,iz
       real,   intent(IN) :: dt
       real,   intent(IN), dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
       real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN), target :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: FlatCoeff
       real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: DivU
#else
       real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN), target :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: FlatCoeff
       real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: DivU
#endif
       logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: Wp, Wn
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: sig
       real, dimension(              HY_WAVENUM,NDIM),intent(OUT) :: lambda
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: leftEig
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: rghtEig
       real, optional, intent(INOUT) :: cellCfl ! We may pass back a lowered CFL factor here - KW
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
     end subroutine hy_uhd_dataReconstOnestep
  end interface



  interface
     Subroutine hy_uhd_DataReconstructNormalDir_MH&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_uhd_DataReconstructNormalDir_MH
  end interface




  interface
     Subroutine hy_uhd_DataReconstructNormalDir_PPM&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_uhd_DataReconstructNormalDir_PPM
  end interface



    interface
     Subroutine hy_uhd_DataReconstructNormalDir_PCM&
          (wenoMethod,dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: wenoMethod,dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_uhd_DataReconstructNormalDir_PCM
  end interface

  

  interface
     Subroutine hy_uhd_DataReconstructNormalDir_WENO&
          (wenoMethod,dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: wenoMethod,dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_uhd_DataReconstructNormalDir_WENO
  end interface



  interface
     Subroutine hy_uhd_DataReconstructNormalDir_GP&
          (radius,dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: radius,dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     end Subroutine hy_uhd_DataReconstructNormalDir_GP
    end interface



  interface
     subroutine hy_uhd_Roe(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm,Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_Roe
  end interface



  interface
     subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1),  intent(OUT) :: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLL
  end interface



  interface
     subroutine hy_uhd_HLLC(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLLC
  end interface



  interface
     subroutine hy_uhd_Marquina(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_Marquina
  end interface


  interface
     subroutine hy_uhd_MarquinaModified(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_MarquinaModified
  end interface


  interface
     subroutine hy_uhd_LLF(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_LLF
  end interface



  interface
     subroutine hy_uhd_setMinTimeStep(blockID,i,j,k,delta,speed)
       implicit none
       integer, INTENT(in) :: blockID,i,j,k
       real, INTENT(in) :: delta,speed
     end subroutine hy_uhd_setMinTimeStep
  end interface




!!$  interface
!!$     subroutine hy_uhd_upwindTransverseFlux(order,vector,lambda,leig,reig,sig)
!!$       implicit none
!!$       integer,intent(IN) :: order
!!$       real,intent(IN),dimension(HY_VARINUMMAX,-2:2)   :: vector
!!$       real,intent(IN),dimension(HY_WAVENUM) :: lambda
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,intent(OUT), dimension(HY_VARINUMMAX) :: sig
!!$     end subroutine hy_uhd_upwindTransverseFlux
!!$  end interface


  interface
     subroutine hy_uhd_upwindTransverseFlux&
          (dir,order,vm1,vc0,vp1,lambda,leig,reig,sigSize,sig,speciesScalar)
       implicit none
       integer,intent(IN) :: dir,order
       real,intent(IN),dimension(:)  :: vm1,vc0,vp1
       real,intent(IN),dimension(HY_WAVENUM) :: lambda
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,pointer,dimension(:) :: lambda
!!$       real,pointer,dimension(:,:) :: leig,reig
       integer,intent(IN) :: sigSize
       real,intent(OUT),dimension(sigSize) :: sig
       logical, intent(IN),optional :: speciesScalar
     end subroutine hy_uhd_upwindTransverseFlux
  end interface



  interface
     subroutine hy_uhd_TVDslope(dir,VL,V0,VR,lambda,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX),intent(IN)  :: VL,V0,VR
       real, dimension(HY_WAVENUM),   intent(IN)  :: lambda
       real, dimension(HY_VARINUM,HY_WAVENUM), intent(IN) :: leig
       real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_uhd_TVDslope
  end interface


  interface
     subroutine hy_uhd_TVDslopeUpwind(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real,dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
       real,dimension(HY_WAVENUM),  intent(IN)  :: lambdaL,lambda,lambdaR
       real,dimension(HY_VARINUM,HY_WAVENUM),intent(IN) :: leig
       real,dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_uhd_TVDslopeUpwind
  end interface



  interface
     subroutine hy_uhd_unsplit( blockCount, blockList, dt, dtOld )
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt, dtOld
     end subroutine hy_uhd_unsplit
  end interface



  interface
     subroutine hy_uhd_unsplitUpdate(blockID,rangeSwitch,dt,del,dataSize,blkLimits,&    
                                     blkLimitsGC,xflux,yflux,zflux,gravX,gravY,gravZ,&
                                     scrch_Ptr)
       implicit none
       integer,intent(IN) :: blockID
       integer,intent(IN) :: rangeSwitch
       real, intent(IN)   :: dt
       real, intent(IN)   :: del(MDIM)
       integer,dimension(MDIM),intent(IN) :: dataSize
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
       integer,intent(IN) :: blkLimitsGC(LOW:HIGH,MDIM)
#ifdef FIXEDBLOCKSIZE
  real, intent(in) ::xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) ::yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, intent(in) ::zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: gravX,gravY,gravZ

#else
  real, intent(in) :: xflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))  
  real, intent(in) :: yflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))  
  real, intent(in) :: zflux(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS))
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)),& 
                  intent(IN) :: gravX,gravY,gravZ
#endif
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
     end subroutine hy_uhd_unsplitUpdate
  end interface





    interface
       subroutine hy_uhd_unsplitUpdateMultiTemp&
            (blockID,rangeSwitch,blkLimits,datasize,dt,del,xflux,yflux,zflux, scrch_Ptr)
         implicit none
         integer,intent(IN) :: blockID, rangeSwitch
         integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
         integer,intent(IN) :: datasize(MDIM)      
         real, intent(IN) :: dt
         real, intent(IN) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
         real, intent(in) :: xflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
         real, intent(in) :: yflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
         real, intent(in) :: zflux(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
         real, intent(in) :: xflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
         real, intent(in) :: yflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
         real, intent(in) :: zflux(NFLUXES,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
         real, pointer, dimension(:,:,:,:) :: scrch_Ptr
       end subroutine hy_uhd_unsplitUpdateMultiTemp
    end interface

    interface
       subroutine hy_uhd_multiTempAfter(blockCount, blockList, dt)
         implicit none
         integer, intent(in) :: blockCount
         integer, intent(in) :: blockList(blockCount)
         real,    intent(in) :: dt
       end subroutine hy_uhd_multiTempAfter
    end interface

    interface
       subroutine hy_uhd_ragelike(uextra, soln, dens_old, &
            duele_adv, duion_adv, durad_adv, &
            xc, yc, zc, &
            uele_new, uion_new, urad_new, &
            stepFactor)
         implicit none
         real, intent(in)  :: uextra
         real, intent(in)  :: soln(NUNK_VARS)
         real, intent(in)  :: dens_old
         real, intent(in)  :: duele_adv
         real, intent(in)  :: duion_adv
         real, intent(in)  :: durad_adv
         real, intent(in)  :: xc 
         real, intent(in)  :: yc 
         real, intent(in)  :: zc
         real, intent(out) :: uele_new
         real, intent(out) :: uion_new
         real, intent(out) :: urad_new
         real, intent(in),OPTIONAL  :: stepFactor
       end subroutine hy_uhd_ragelike
    end interface

  interface
     subroutine hy_uhd_updateSpeciesMassScalar&
          (order,densNew,Sp,U,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt,SpNew)
  implicit none
  integer, intent(IN) :: order
  real, intent(IN) :: densNew
#if NDIM == 1
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,1,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,1,1) :: U
#elif NDIM == 2
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,1) :: U
#elif NDIM == 3
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,-3:3)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,-3:3)   :: U
#endif
  real, intent(IN)  :: FL,FR,GL,GR,HL,HR,dx,dy,dz,dt
  real, intent(OUT), dimension(NSPECIES+NMASS_SCALARS) :: SpNew
     end subroutine hy_uhd_updateSpeciesMassScalar
  end interface


  interface
     subroutine hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,eosMode,rangeSwitch,noMHD)
       implicit none
       integer, intent(IN) :: blockID
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       real, intent(IN) :: dt,dtOld
       real, dimension(MDIM), intent(IN) :: del
       integer, intent(IN) :: eosMode
       integer,OPTIONAL,intent(IN) :: rangeSwitch
       logical,OPTIONAL,intent(IN) :: noMHD
     end subroutine hy_uhd_energyFix
  end interface



  interface 
     subroutine hy_uhd_unitConvert(blockID,convertDir)
       implicit none
       integer, intent(IN) :: blockID, convertDir
     end subroutine hy_uhd_unitConvert
  end interface
  


  interface
     subroutine hy_uhd_addViscousFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,mu,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: mu
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: mu
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_uhd_addViscousFluxes
    end interface



  interface
     subroutine hy_uhd_addThermalFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,kappa,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: kappa
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: kappa
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_uhd_addThermalFluxes
    end interface



    interface
       subroutine hy_uhd_eigenParameters&
            (V,dir,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM2), intent(IN)  :: V
         integer, intent(IN)  :: dir
         real, intent(OUT) :: U_normal,C_fast
         real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(OUT),optional :: B_beta
         real, intent(OUT), optional :: C_hyp
       end subroutine hy_uhd_eigenParameters
    end interface



    interface
       subroutine hy_uhd_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow,C_hyp)
         implicit none
         real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
         real,intent(IN) :: U_normal,C_fast
         real,intent(IN), optional :: C_alfn,C_slow,C_hyp
       end subroutine hy_uhd_eigenValue
    end interface



    interface
       subroutine hy_uhd_eigenVector&
            (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: LeftEigvec
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
         real, dimension(HY_VARINUM2), intent(IN) :: V
         integer, intent(IN) :: dir
         logical, intent(IN) :: cons
         real,    intent(IN) :: C_fast
         real,    intent(IN), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(IN), optional  :: B_beta
         real,    intent(IN), optional :: C_hyp
       end subroutine hy_uhd_eigenVector
    end interface


    interface
       subroutine hy_uhd_prepareNewGravityAccel(blockCount,blockList,gcMaskLogged)
         integer, INTENT(IN) ::  blockCount
         integer, INTENT(IN), dimension(blockCount) :: blockList
         logical,intent(in) :: gcMaskLogged
       end subroutine hy_uhd_prepareNewGravityAccel
    end interface


    interface
       subroutine hy_uhd_putGravityUnsplit&
            (blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ,potentialIndex,&
             lastCall)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
         integer, dimension(MDIM), intent(IN) :: dataSize
         real,    intent(IN) :: dt, dtOld
#ifdef FIXEDBLOCKSIZE
         real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(OUT) :: &
              gravX,gravY,gravZ
#else
         real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(OUT) :: &
              gravX,gravY,gravZ
#endif
         integer, intent(IN), OPTIONAL :: potentialIndex
         logical, intent(IN), OPTIONAL :: lastCall
       end subroutine hy_uhd_putGravityUnsplit
    end interface



    interface
       Subroutine hy_uhd_addGravityUnsplit&
            (blockID,blkLimitsGC,dataSize,dt,gravX,gravY,gravZ)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
         integer, dimension(MDIM), intent(IN) :: dataSize
         real,    intent(IN) :: dt
#ifdef FIXEDBLOCKSIZE
         real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
              gravX,gravY,gravZ
#else
         real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) :: &
              gravX,gravY,gravZ
#endif
       end Subroutine hy_uhd_addGravityUnsplit
    end interface



    interface
       subroutine hy_uhd_shockDetect(blockID)
         implicit none
         integer, INTENT(IN) :: blockID
       end subroutine hy_uhd_shockDetect
    end interface



    interface
       subroutine hy_uhd_prim2con(V,CU)
         implicit none
         real ,dimension(HY_VARINUM2), intent(IN) :: V
         real ,dimension(HY_VARINUM),  intent(OUT) :: CU
       end subroutine hy_uhd_prim2con
    end interface



    interface
       subroutine hy_uhd_con2prim(CU,game,V)
         implicit none
         real ,dimension(HY_VARINUM), intent(IN)  :: CU
         real, intent(IN) :: game
         real ,dimension(HY_VARINUM), intent(OUT) :: V
       end subroutine hy_uhd_con2prim
    end interface



    interface
       subroutine hy_uhd_prim2flx(dir,V,F)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(*),            intent(IN) :: V
         real, dimension(HY_VARINUM1),  intent(OUT) :: F
       end subroutine hy_uhd_prim2flx
    end interface



    interface
       Subroutine hy_uhd_checkRHjumpCond(dir,dens,velocity,pres,gamc,Wp,Wn,SWp,SWn)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(MDIM), intent(IN) :: velocity
         real, intent(IN) :: dens, pres, gamc
         real, dimension(HY_DENS:HY_EINT), intent(IN) :: Wp,Wn
         logical, intent(OUT) :: SWp,SWn
       end Subroutine hy_uhd_checkRHjumpCond
    end interface


#ifdef FLASH_UGLM_MHD
    interface
       Subroutine hy_uhd_updateSourceGLM(blockID,dt,del,blkLimits,blkLimitsGC,C_hyp,C_par)
         implicit none
         integer, intent(IN) :: blockID
         real,    intent(IN) :: dt
         real,    dimension(MDIM),   intent(IN) :: del
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
         real,    intent(IN) :: C_hyp,C_par
       end Subroutine hy_uhd_updateSourceGLM
    end interface
#endif


  ! Dongwook thinks that the following two interfaces should be taken care of for
  ! general implementation later.
!!$  interface
!!$     Subroutine hy_uhd_counterGP_init(radiusGP, counterGP, blkLimitsGC)
!!$       implicit none
!!$       real,    intent(IN)  :: radiusGP
!!$       integer, intent(OUT) :: counterGP
!!$       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
!!$     end subroutine hy_uhd_counterGP_init
!!$  end interface
!!$
  interface
     Subroutine hy_uhd_initGP(Zk, Wk, Pk)
       implicit none
       real, intent(INOUT), dimension(:,:,:) :: Zk
       real, intent(INOUT), dimension(:,:) :: Wk
       real, intent(INOUT), dimension(:,:) :: Pk
     end Subroutine hy_uhd_initGP
  end interface


!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  interface
     subroutine hy_uhd_addResistiveFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,eta,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE 
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, & 
                  GRID_JLO_GC:GRID_JHI_GC, & 
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: eta
#else 
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), & 
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), & 
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),& 
                  intent(IN) :: eta
#endif 
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_uhd_addResistiveFluxes
    end interface



  interface
     subroutine hy_uhd_addOhmicHeating&
          (blockID,blkLimits,ix,iy,iz,Qohm,eta)
     implicit none
     integer, INTENT(IN) :: blockID,ix,iy,iz
     integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
     real, intent(INOUT) :: Qohm
     real, intent(IN)    :: eta  
     end subroutine hy_uhd_addOhmicHeating
  end interface


  interface
     subroutine hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
       implicit none
       integer, intent(IN) :: blockID
       real,    intent(IN) :: dt
       real,    dimension(MDIM),   intent(IN) :: del
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
       logical, intent(IN) :: halfTimeAdvance
     end subroutine hy_uhd_staggeredDivb
  end interface



  interface
     subroutine hy_uhd_HLLD(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_uhd_HLLD
  end interface



  interface
     subroutine hy_uhd_getElectricFields( blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
       implicit none
       integer, intent(IN)  :: blockID
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN), dimension(MDIM)  :: del  
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,               &
                       GRID_ILO_GC:GRID_IHI_GC,  &
                       GRID_JLO_GC:GRID_JHI_GC,  &
                       GRID_KLO_GC:GRID_KHI_GC), &
                       intent(IN) :: flx,fly,flz
#else
       real, DIMENSION(NFLUXES,               &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: flx,fly,flz
#endif
     end subroutine hy_uhd_getElectricFields
  end interface



  interface
     subroutine hy_uhd_getCurrents( &
        blockID, rangeSwitch, blkLimits,datasize, del, Jp, Jm, mode_switch,&
        scrch_Ptr,&
        ix,iy,iz)
        implicit none
        ! Arguments:
        integer,intent(in) :: blockID, rangeSwitch
        integer,intent(in) :: blkLimits  (LOW:HIGH,MDIM)
        integer,intent(in) :: datasize(MDIM), mode_switch
        real,   intent(in) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
        real, intent(inout) :: Jp(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
        real, intent(inout) :: Jm(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
        real, intent(inout) :: Jp(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
        real, intent(inout) :: Jm(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
        real, pointer, dimension(:,:,:,:) :: scrch_Ptr
        integer, intent(in), OPTIONAL :: ix,iy,iz 
     end subroutine hy_uhd_getCurrents
  end interface



  interface
     subroutine hy_uhd_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
                                     fluxType,DerivDir,   &
                                     faceFlux,            &
                                     Flux1Deriv,          &
                                     Flux2Deriv           )
       implicit none
       integer, intent(IN) :: ix,iy,iz
       integer, intent(IN) :: fluxType,DerivDir
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real,    dimension(NFLUXES,                  &
                          GRID_ILO_GC:GRID_IHI_GC,  &
                          GRID_JLO_GC:GRID_JHI_GC,  &
                          GRID_KLO_GC:GRID_KHI_GC), &
                          intent(IN) :: faceFlux
#else
       real,    dimension(NFLUXES,                  &
                          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                          intent(IN) :: faceFlux
#endif
       real,    intent(OUT):: Flux1Deriv,Flux2Deriv
     end subroutine hy_uhd_getFluxDeriv
  end interface
  
  interface
     subroutine hy_uhd_biermannSource(blockCount, blockList, dt)
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt
     end subroutine hy_uhd_biermannSource
  end interface

  interface
    subroutine hy_uhd_addHallFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
       integer, INTENT(IN) :: sweepDir       
     end subroutine hy_uhd_addHallFluxes
  end interface
  
  interface
     Subroutine hy_uhd_addBiermannFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir,dt)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
       real, dimension(NFLUXES), intent(INOUT) :: Flux
       integer, INTENT(IN) :: sweepDir
       real, optional, intent(IN) :: dt
     end Subroutine hy_uhd_addBiermannFluxes
  end interface
 
#endif 
!endif #ifdef FLASH_USM_MHD

End Module hy_uhd_interface
