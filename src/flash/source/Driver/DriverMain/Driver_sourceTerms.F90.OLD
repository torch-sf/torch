!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)

  use Polytrope_interface, ONLY : Polytrope
  use Driver_data, ONLY: dr_simTime
  use Flame_interface, ONLY : Flame_step
  use Stir_interface, ONLY : Stir
  use Heat_interface, ONLY : Heat
  use Heatexchange_interface, ONLY : Heatexchange
  use Burn_interface, ONLY : Burn
  use Cool_interface, ONLY : Cool
  use Ionize_interface, ONLY : Ionize
  use Flame_interface, ONLY : Flame_step
  use Turb_interface, ONLY : Turb_calc
  use EnergyDeposition_interface, ONLY : EnergyDeposition
  use Deleptonize_interface, ONLY : Deleptonize
  use SN_interface, ONLY : SN

  implicit none

  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList
  integer, OPTIONAL, intent(IN):: pass

#ifdef DRIVER_EDEARLY
  call EnergyDeposition(blockCount, blockList, dt, dr_simTime, pass)
#endif
  call Polytrope(blockCount, blockList, dt)
  call Stir(blockCount, blockList, dt) 
  call Turb_calc(blockCount, blockList)
  call Flame_step(blockCount, blockList, dt)
  call Burn(blockCount, blockList, dt) 
! SN before heating as it adds a heating rate. 
! The order should be ionising radiation heating, SN heating
! then actual heating/cooling routine.
  call SN(blockCount, blockList, dt, dr_simTime)
  call Heat(blockCount, blockList, dt, dr_simTime) 
  call Heatexchange(blockCount, blockList, dt)
  call Cool(blockCount, blockList, dt, dr_simTime)
  call Ionize(blockCount, blockList, dt, dr_simTime)
#ifndef DRIVER_EDEARLY
  call EnergyDeposition(blockCount, blockList, dt, dr_simTime, pass)
#endif
  call Deleptonize(blockCount, blockList, dt, dr_simTime)

  return
end subroutine Driver_sourceTerms
