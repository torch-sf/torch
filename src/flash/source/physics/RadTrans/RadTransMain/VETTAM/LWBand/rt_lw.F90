!! Photodissociation due to LW photons of H2 molecules
!! AUTHOR
!!  Shyam Harimohan Menon (2023)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!  time : the current time
!! 
!!
!!***

SUBROUTINE rt_lw(blockCount_,blockList_,dt)

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt
  !Stub implementation
  return
  
END SUBROUTINE rt_lw