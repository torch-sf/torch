!*******************************************************************************

! Module:       RadTrans_TreeCommData

! Description:  Tree communication data types and structures.


module RadTrans_TreeCommData

!===============================================================================

  save

  integer, parameter :: MAX_NDIM = 3
  integer, parameter :: MAX_LOCBLOCKS = MAXBLOCKS, MAX_PROCS = 1024
  integer, parameter :: MAX_CHILD = 2**MAX_NDIM, MAX_FACE = 2*MAX_NDIM
!  integer, parameter :: MAX_MSGS = 3

!  integer, dimension(2,MAX_CHILD,MAX_LOCBLOCKS,0:MAX_PROCS-1) :: g_child
!  integer, dimension(2,MAX_LOCBLOCKS,0:MAX_PROCS-1)           :: g_parent
!  integer, dimension(MAX_LOCBLOCKS,0:MAX_PROCS-1)             :: g_node_type, &
!                                                                 g_lrefine
  integer, dimension(MAX_LOCBLOCKS,0:MAX_PROCS-1)             :: g_lrefine
!  integer, dimension(2,MAX_FACE,MAX_LOCBLOCKS,0:MAX_PROCS-1)  :: g_nbr
!  real, dimension(MAX_NDIM,MAX_LOCBLOCKS,0:MAX_PROCS-1)       :: g_size, g_coord
!  integer, dimension(0:MAX_PROCS-1)                           :: g_lnblocks
  real, dimension(2,MAX_NDIM,MAX_LOCBLOCKS,0:MAX_PROCS-1)     :: g_bndbox

  integer, dimension(0:MAX_PROCS-1)               :: my_n_to_send, partcounts
  integer, dimension(0:MAX_PROCS-1,0:MAX_PROCS-1) :: n_to_send
  
!===============================================================================

end module RadTrans_TreeCommData
