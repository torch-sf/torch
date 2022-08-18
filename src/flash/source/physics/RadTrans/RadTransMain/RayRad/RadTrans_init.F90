!!****if* source/physics/RadTrans/RadTransMain/RadTrans_init
!!
!!  NAME 
!!
!!  RadTrans_init
!!
!!  SYNOPSIS
!!
!!  call RadTrans_init()
!!
!!  DESCRIPTION 
!!    Initialize radiative transfer unit
!!
!! ARGUMENTS
!!
!!
!!***

#include "constants.h"

subroutine RadTrans_init()
  use RadTrans_data
  use rt_interface, ONLY: rt_init
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm
  use Grid_interface, ONLY : Grid_releaseBlkPtr, Grid_getBlkPtr, Grid_getListOfBlocks

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: blockID, thisBlock
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  real, pointer, dimension(:,:,:,:) :: solnData


  call Driver_getMype(MESH_COMM,rt_meshMe)

  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)

  ! Store physical constants:
  call PhysicalConstants_get("speed of light",rt_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",rt_radconst)
  rt_radconst = 4.0 * rt_radconst / rt_speedlt

  call PhysicalConstants_get("Boltzmann", rt_boltz)
  ! not used

  call RuntimeParameters_get("rt_dtFactor", rt_dtFactor)
  call RuntimeParameters_get("meshCopyCount", rt_meshCopyCount)

  call Driver_getMype(MESH_ACROSS_COMM, rt_acrossMe)
  call Driver_getComm(MESH_ACROSS_COMM, rt_acrossComm)
  call Driver_getComm(GLOBAL_COMM,rt_globalComm)

  call rt_init

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do thisBlock = 1, blockCount
    blockID = blockList(thisBlock)
	! Get a pointer to solution data 
    call Grid_getBlkPtr(blockID,solnData)
    solnData(PHIO_VAR,:,:,:) =  0d0
    call Grid_releaseBlkPtr(blockID,solnData)
  enddo ! block loop

  return
end subroutine RadTrans_init
