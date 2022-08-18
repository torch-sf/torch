!!****if* source/Driver/DriverMain/Driver_initSourceTerms
!!
!! NAME
!!   
!!  Driver_initSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_initSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes all source terms Units by 
!!   calling their respective initialization routines
!!   viz. Stir_init, Burn_init, Heat_init, Cool_init, etc.
!!  
!! ARGUMENTS
!!   myPE - current processor number
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_initSourceTerms( restart)

  use Polytrope_interface, ONLY:  Polytrope_init
  use Burn_interface, ONLY:  Burn_init
  use Stir_interface, ONLY : Stir_init
  use Heat_interface, ONLY : Heat_init
  use Heatexchange_interface, ONLY : Heatexchange_init
  use Cool_interface, ONLY : Cool_init
  use Diffuse_interface, ONLY : Diffuse_init
  use Ionize_interface, ONLY : Ionize_init
  use Flame_interface, ONLY : Flame_init
  use Turb_interface, ONLY : Turb_init
  use RadTrans_interface, ONLY : RadTrans_init
  use EnergyDeposition_interface, ONLY : EnergyDeposition_init
  use Deleptonize_interface, ONLY : Deleptonize_init
  use GridInject_interface, ONLY : GridInject_init
  use SN_interface, ONLY : SN_init

  implicit none
  
  logical, intent(in) :: restart

  call Polytrope_init()
  call Stir_init( restart)
  call Cool_init()
  call Diffuse_init()
  call Heat_init()
  call GridInject_init()
  call SN_init()
  call Heatexchange_init( restart)
  call Ionize_init()
  call Burn_init()
  call Turb_init()
  call Flame_init()
  call RadTrans_init()
  call EnergyDeposition_init()
  call Deleptonize_init()

end subroutine Driver_initSourceTerms
