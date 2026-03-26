!! Photoionization coupled to VETTAM for modelling EUV radiation. This is just a stub implementation; see subdirs for actual implementation.
!!
!!***

SUBROUTINE rt_ionise(blockCount_,blockList_,dt,time)
  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt, time

  !Stub implementation
  return
END SUBROUTINE rt_ionise