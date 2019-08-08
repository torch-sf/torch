!!****if* source/Driver/DriverMain/Driver_finalizeSourceTerms
!!
!! NAME
!!   
!!  Driver_finalizeSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_finalizeSourceTerms(logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Finalizes all source terms Units by 
!!   calling their respective termination routines
!!   viz. Stir_finalize, Burn_finalize, Heat_finalize, Cool_finalize, etc.
!!  
!! ARGUMENTS
!!
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_finalizeSourceTerms( restart)

  use Polytrope_interface, ONLY:  Polytrope_finalize
  use Burn_interface, ONLY:  Burn_finalize
  use Stir_interface, ONLY : Stir_finalize
  use Heat_interface, ONLY : Heat_finalize
  use Heatexchange_interface, ONLY: Heatexchange_finalize
  use Diffuse_interface, ONLY : Diffuse_finalize
  use Cool_interface, ONLY : Cool_finalize
  use Flame_interface, ONLY : Flame_finalize
  use Turb_interface, ONLY : Turb_finalize
  use EnergyDeposition_interface, ONLY : EnergyDeposition_finalize
  use Deleptonize_interface, ONLY : Deleptonize_finalize
  use RadTrans_interface, ONLY : RadTrans_finalize
  use SN_interface, ONLY : SN_finalize

  implicit none
  
  logical, intent(in) :: restart

  call Polytrope_finalize()
  call Stir_finalize()
  call Cool_finalize()
  call Diffuse_finalize()
  call SN_finalize()
  call Heat_finalize()
  call Heatexchange_finalize()
  call EnergyDeposition_finalize()
  call Deleptonize_finalize()
  call RadTrans_finalize()

!!$  call Ioniz_finalize()
  call Flame_finalize()
  call Burn_finalize()
  call Turb_finalize()

end subroutine Driver_finalizeSourceTerms
