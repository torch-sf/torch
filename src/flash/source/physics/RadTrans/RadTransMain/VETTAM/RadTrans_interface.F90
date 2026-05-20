!!****h* source/physics/RadTrans/RadTrans_interface
!!
!! NAME
!!   RadTrans_interface
!!
!! SYNOPSIS
!!   use RadTrans_interface : ONLY [...]
!!
!! DESCRIPTION
!!   Public interface for the RadTrans unit
!!***
module RadTrans_interfaceTypeDecl
  implicit none
  type RadTrans_dbgContext_t
     integer :: step
     integer :: group
     integer :: component
     integer :: libErrCode
     integer :: flashErrCode
     integer :: retriable       ! 0 for NO, 1 for YES
     logical :: willingToRetry
  end type RadTrans_dbgContext_t
end module RadTrans_interfaceTypeDecl


module RadTrans_interface
  use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
  implicit none

#include "constants.h"

  interface
     subroutine RadTrans_getDbgContext(context)
       use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
!!$       import RadTrans_dbgContext_t ! a modern alternative
       implicit none
       type(RadTrans_dbgContext_t),intent(OUT) :: context
     end subroutine RadTrans_getDbgContext
     subroutine RadTrans_getDbgContextPtr(context)
       use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
!!$       import RadTrans_dbgContext_t ! a modern alternative
       implicit none
       type(RadTrans_dbgContext_t),pointer :: context
     end subroutine RadTrans_getDbgContextPtr
  end interface

  interface
     subroutine RadTrans(nblk, blklst, dt, pass)
       implicit none
       integer, intent(in) :: nblk
       integer, intent(in) :: blklst(nblk)
       real,    intent(in) :: dt
       integer, intent(in), optional :: pass
     end subroutine RadTrans
  end interface

  interface
     subroutine IonHeatCool(nblk, blklst, dt, pass)
       implicit none
       integer, intent(in) :: nblk
       integer, intent(in) :: blklst(nblk)
       real,    intent(in) :: dt
       integer, intent(in), optional :: pass
     end subroutine IonHeatCool
  end interface

  interface
     subroutine RadTrans_computeDt(blockID,  blkLimits,blkLimitsGC, &
          solnData, dt_radtrans, dt_minloc)
       implicit none
       integer, intent(IN) :: blockID
       integer, intent(IN) :: blkLimits(2,MDIM)
       integer, intent(IN) :: blkLimitsGC(2,MDIM)
       real,  pointer :: solnData(:,:,:,:) 
       real, intent(INOUT) :: dt_radtrans
       integer, intent(INOUT)  :: dt_minloc(5)
     end subroutine RadTrans_computeDt
  end interface
  
  interface
     subroutine RadTrans_computeFluxLimiter(ifl, iflOut, ieddi3, solnData, blockID, gcLayers)
       implicit none
       integer, intent(in) :: ifl
       integer, intent(in) :: iflOut
       integer, intent(in) :: ieddi3
       real,    intent(INOUT) :: solnData(:,1:,1:,1:)
       integer, intent(IN) :: blockID
       integer, intent(IN),OPTIONAL :: gcLayers
     end subroutine RadTrans_computeFluxLimiter
  end interface

  interface
     subroutine RadTrans_init()
       implicit none
     end subroutine RadTrans_init
  end interface

  interface
     subroutine RadTrans_planckInt(x, p)
       implicit none
       real, intent(in) :: x
       real, intent(out) :: p
     end subroutine RadTrans_planckInt
  end interface

  interface 
     subroutine RadTrans_mgdGetBound(g, b)
       implicit none
       integer, intent(in) :: g
       real, intent(out) :: b
     end subroutine RadTrans_mgdGetBound
  end interface

  interface 
     subroutine RadTrans_mgdSetBound(g, b)
       implicit none
       integer, intent(in) :: g
       real, intent(in) :: b
     end subroutine RadTrans_mgdSetBound
  end interface

  interface
     subroutine RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
       implicit none
       integer, intent(in) :: blockId
       integer, intent(in) :: axis(MDIM)
       real,    intent(in) :: trad
       real,    intent(out), optional :: tradActual
     end subroutine RadTrans_mgdEFromT
  end interface

  interface
     subroutine RadTrans_sumEnergy(ivar, nblk, blklst)
       implicit none
       integer, intent(in) :: ivar
       integer, intent(in) :: nblk
       integer, intent(in) :: blklst(nblk)
     end subroutine RadTrans_sumEnergy
  end interface

  interface
     subroutine RadTrans_mgdSetEnergy(blockId, axis, grpNum, eg)
       implicit none
       integer, intent(in) :: blockId
       integer, intent(in) :: axis(MDIM)
       integer, intent(in) :: grpNum
       real,    intent(in) :: eg
     end subroutine RadTrans_mgdSetEnergy
  end interface

  interface
     subroutine RadTrans_mgdSetBc(ig, bcTypes, bcValues, f, bcType, bcValue)
       implicit none

       integer, intent(in) :: ig
       
       integer, optional, intent(in) :: bcTypes(6)
       real, optional, intent(in) :: bcValues(6)
       
       integer, optional, intent(in) :: f
       integer, optional, intent(in) :: bcType
       real, optional, intent(in) :: bcValue

     end subroutine RadTrans_mgdSetBc
  end interface

  interface 
     subroutine RadTrans_finalize ()
       implicit none
     end subroutine RadTrans_finalize       
  end interface

end module RadTrans_interface
