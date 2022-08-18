!!****h* source/physics/sourceTerms/Heat/Heat_interface
!!
!! This is the header file for the heat module that defines its
!! public interfaces.
!!***
Module Heat_interface
#include "constants.h"
#include "Flash.h"

  interface 
     subroutine Heat (blockCount,blockList,dt,time)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt,time
     end subroutine Heat
  end interface

! This is not called in sourceTerms as Heat, but instead is 
! subcycled over in RadTrans simultaineously with
! ray tracing and ionization calculation for both
! speed and accuracy. -Josh Wall
  interface 
     subroutine RadHeat (blockCount,blockList,dt,time)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt,time
     end subroutine RadHeat
  end interface

  interface Heat_computeDt
     subroutine Heat_computeDt (block_no, &
          x, dx, uxgrid, &
          y, dy, uygrid, &
          z, dz, uzgrid, &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_check, dt_minloc )

       integer, intent(IN) :: block_no
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: dx, uxgrid
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: dy, uygrid
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: dz, uzgrid
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: dx, uxgrid
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: dy, uygrid
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: dz, uzgrid
#endif
       real,INTENT(OUT)    :: dt_check
       integer,INTENT(OUT)    :: dt_minloc(5)
       real, pointer, dimension(:,:,:,:) :: solnData
     end subroutine Heat_computeDt

  end interface

  interface Heat_init
     subroutine Heat_init()
         
     end subroutine Heat_init
  end interface


  interface Heat_finalize
     subroutine Heat_finalize ()
     end subroutine Heat_finalize
  end interface




end Module Heat_interface
