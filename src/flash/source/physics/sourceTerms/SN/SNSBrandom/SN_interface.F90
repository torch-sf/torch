!!****h* source/physics/sourceTerms/Heat/Heat_interface
!!
!! This is the header file for the heat module that defines its
!! public interfaces.
!!***
Module SN_interface
#include "constants.h"
#include "Flash.h"

  interface SN
     subroutine SN (blockCount,blockList,dt,time)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt,time
     end subroutine SN
  end interface

  interface SN_computeDt
     subroutine SN_computeDt (block_no, &
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
     end subroutine SN_computeDt

  end interface

  interface SN_init
     subroutine SN_init()
         
     end subroutine SN_init
  end interface


  interface SN_finalize
     subroutine SN_finalize ()
     end subroutine SN_finalize
  end interface

end Module SN_interface
