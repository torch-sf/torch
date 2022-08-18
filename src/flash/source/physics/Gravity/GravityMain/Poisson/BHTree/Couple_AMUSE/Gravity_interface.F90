!!****h* source/physics/Gravity/Gravity_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module Gravity_interface

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  interface
     subroutine Gravity_accelAtCoords (numPoints, iCoords,jCoords,kCoords, accelDir,&
          accel, blockID, &
          potentialIndex)
       integer, intent(IN) :: accelDir, numPoints
       real, dimension(:),INTENT(in) :: iCoords,jCoords,kCoords
       real, dimension(numPoints),INTENT(OUT) :: accel
       integer, intent(IN),optional :: blockID
       integer, intent(IN),optional :: potentialIndex
     end subroutine Gravity_accelAtCoords
  end interface

  interface Gravity_accelListOfBlocks
     subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
          accelIndex, potentialIndex)
       integer,intent(IN)                      :: blockCount
       integer,dimension(blockCount), intent(IN)     :: blockList
       integer, INTENT(in) ::  component
       integer, intent(in) :: accelIndex
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_accelListOfBlocks
  end interface

  interface
     subroutine Gravity_accelOneRow (pos,sweepDir,blockID, numCells, grav, &
           varIndex, extraAccelVars)
       implicit none
       integer, intent(IN) :: sweepDir,blockID,numCells
       integer, dimension(2),INTENT(in) ::pos
       real, dimension(numCells),INTENT(inout) :: grav
       integer, intent(IN), optional :: varIndex 
       integer, intent(IN),OPTIONAL      :: extraAccelVars(MDIM)
     end subroutine Gravity_accelOneRow
  end interface

  interface Gravity_accelOneBlock
     subroutine Gravity_accelOneBlock ( blockID, ngcellcomp, gvec, potentialIndex)
        integer, intent(in)             :: blockID,  ngcellcomp
        real, dimension(:,:,:,:),intent(out)  :: gvec
        integer, intent(in),optional    :: potentialIndex
     end subroutine Gravity_accelOneBlock
  end interface


  interface Gravity_computeDt
     subroutine Gravity_computeDt (blockID, dt_grav, dt_minloc)
       real,intent(OUT)       ::  dt_grav
       integer, intent(IN)    ::  blockID
       integer, intent(INOUT) :: dt_minloc(5)
     end subroutine Gravity_computeDt
  end interface

  interface Gravity_finalize
     subroutine Gravity_finalize()
     end subroutine Gravity_finalize
  end interface

  interface Gravity_init
     subroutine Gravity_init()
       
     end subroutine Gravity_init
  end interface

  interface Gravity_potentialListOfBlocks
     subroutine Gravity_potentialListOfBlocks(blockCount,blockList, &
          potentialIndex)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       integer, intent(IN), optional :: potentialIndex
     end subroutine Gravity_potentialListOfBlocks
  end interface

  interface Gravity_unitTest
     subroutine Gravity_unitTest( fileUnit, perfect)
       implicit none
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine Gravity_unitTest
  end interface

! The following subroutines and functions are needed by the tree solver
  interface Gravity_bhAccBotNode
    subroutine Gravity_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, dimension(:,:,:) :: locCoords
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(IN) :: botnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine Gravity_bhAccBotNode
  end interface

  interface Gravity_bhAccNode
    subroutine Gravity_bhAccNode(subnode, accnode)
      real, dimension(:), intent(IN)  :: subnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine Gravity_bhAccNode
  end interface


  interface Gravity_bhBTInit
    subroutine Gravity_bhBTInit()
    end subroutine Gravity_bhBTInit
  end interface

  interface
    logical function Gravity_bhCheckAccuracy(blockno, point, blkLimits, solnData)
    implicit none
    integer, intent(IN) :: blockno
    integer, dimension(MDIM), intent(IN) :: point
    integer, dimension(2,MDIM), intent(IN)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end function Gravity_bhCheckAccuracy
  end interface

  interface Gravity_bhFillBotNode
    subroutine Gravity_bhFillBotNode(blockno, point, blkLimits, solnData, botnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(INOUT) :: botnode
    end subroutine Gravity_bhFillBotNode
  end interface

  interface Gravity_bhFinalizeIter
    subroutine Gravity_bhFinalizeIter()
    end subroutine Gravity_bhFinalizeIter
  end interface

  interface Gravity_bhNodeContrib
    subroutine Gravity_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      integer, intent(IN) :: trLevel, refLevel, blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhNodeContrib
  end interface

  interface Gravity_bhNormalizeNode
    subroutine Gravity_bhNormalizeNode(smr, node)
      real, dimension(MDIM), intent(IN) :: smr
      real, dimension(:), intent(INOUT) :: node
    end subroutine Gravity_bhNormalizeNode
  end interface

  interface Gravity_bhPostprocNode
    subroutine Gravity_bhPostprocNode(ndSize, node)
      real, intent(IN) :: ndSize
      real, dimension(:), intent(INOUT) :: node
    end subroutine Gravity_bhPostprocNode
  end interface

  interface Gravity_bhStartBlock
    subroutine Gravity_bhStartBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhStartBlock
  end interface

  interface Gravity_bhTreeWalkEnd
    subroutine Gravity_bhTreeWalkEnd()
    end subroutine Gravity_bhTreeWalkEnd
  end interface

  interface Gravity_bhFinalizeBlock
    subroutine Gravity_bhFinalizeBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
    end subroutine Gravity_bhFinalizeBlock
  end interface


  interface Gravity_bhInitFieldVar
    subroutine Gravity_bhInitFieldVar(gpotVar)
      implicit none
      integer, intent(IN) :: gpotVar
    end subroutine Gravity_bhInitFieldVar
  end interface

  interface Gravity_bhMAC
    logical function Gravity_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize2
    integer, dimension(MDIM), intent(IN) :: point
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    end function Gravity_bhMAC
  end interface
  
  interface Gravity_bhPartErr
    subroutine Gravity_bhPartErr(node, ndSize, dr, perr)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize
    real, dimension(MDIM+2), intent(IN) :: dr
    real :: perr
    end subroutine Gravity_bhPartErr
  end interface
  
  interface Gravity_bhSelfContrib
    subroutine Gravity_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
      implicit none
      real, dimension(:), intent(IN) :: node
      real, dimension(MDIM), intent(IN) :: cellsize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine Gravity_bhSelfContrib
  end interface

  interface Gravity_bhEwaldAccV42
    function Gravity_bhEwaldAccV42(x, y, z)
      implicit none
      real,intent(in)    :: x, y, z
      real :: Gravity_bhEwaldAccV42(IAXIS:KAXIS)
    end function Gravity_bhEwaldAccV42
  end interface

  interface Gravity_bhEwaldPotV42
    function Gravity_bhEwaldPotV42(x, y, z)
      implicit none
      real,intent(in)    :: x, y, z
      real :: Gravity_bhEwaldPotV42
    end function Gravity_bhEwaldPotV42
  end interface

  interface Gravity_bhGetNodeStruct
    subroutine Gravity_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
      implicit none
      integer, intent(IN) :: im, ix, iy, iz
      integer, intent(INOUT) :: bnsize, nsize
    end subroutine Gravity_bhGetNodeStruct
  end interface

  interface Gravity_bhEwaldAcc
    function Gravity_bhEwaldAcc(x, y, z, drAbsInv)
      implicit none
      real,intent(in)    :: x, y, z, drAbsInv
      real :: Gravity_bhEwaldAcc(IAXIS:KAXIS)
    end function Gravity_bhEwaldAcc
  end interface

  interface Gravity_bhEwaldPot
    function Gravity_bhEwaldPot(x, y, z, drAbsInv)
      implicit none
      real,intent(in)    :: x, y, z, drAbsInv
      real :: Gravity_bhEwaldPot
    end function Gravity_bhEwaldPot
  end interface


! end of tree solver subroutines/functions

!!! Additions for AMUSE interface. - Josh Wall

  interface Gravity_getAccelAtPoint
    subroutine Gravity_getAccelAtPoint(x, y, z, gaccel_x, gaccel_y, gaccel_z)
      implicit none
      real, intent(IN) :: x, y, z
      real, intent(INOUT) :: gaccel_x, gaccel_y, gaccel_z
    end subroutine Gravity_getAccelAtPoint
  end interface
  
  interface Gravity_getPotentialAtPoint
    subroutine Gravity_getPotentialAtPoint(x, y, z, grav_pot)
      implicit none
      real, intent(IN) :: x, y, z
      real, intent(INOUT) :: grav_pot
    end subroutine Gravity_getPotentialAtPoint
  end interface
  
  interface Gravity_potentialListOfBlocksPDENonly
     subroutine Gravity_potentialListOfBlocksPDENonly(blockCount,blockList, &
          potentialIndex, which_particles)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
       integer, intent(IN), optional :: potentialIndex
       integer, intent(IN), optional :: which_particles
     end subroutine Gravity_potentialListOfBlocksPDENonly
  end interface

!!! End AMUSE

end Module Gravity_interface
