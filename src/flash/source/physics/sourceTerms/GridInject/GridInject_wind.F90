!!****f* source/physics/sourceTerms/GridInject/GridInject_wind
!!
!! NAME
!!
!!  GridInject_wind
!!
!! SYNOPSIS
!!
!!  call GridInject_wind(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: injectMassIn
!!      real(IN)          :: injectVelocityIn
!!      real(IN)          :: twind
!!      real(IN)          :: dt
!!      real(INOUT)       :: bgDens
!!      optional,logical(IN) :: snap_to_grid
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit a wind onto grid with both kinetic and thermal energy;
!!  winds can be modified (mass loaded) via runtime parameters.
!!
!! ARGUMENTS
!!
!!  xloc              : where to inject
!!  yloc              : where to inject
!!  zloc              : where to inject
!!  injectMassIn      : amount of mass to inject
!!  injectVelocityIn  : velocity of injected mass
!!  twind             : lifetime of injected wind
!!  dt                : timestep for wind injection, to get dm/dt
!!  bgDens            : initial background density of wind.
!!                      if bgDens==0 and runtime param wind_var_radius is set,
!!                      calculate and return to caller, else use input value
!!  snap_to_grid      : for testing/debugging
!!
!!***
subroutine GridInject_wind (xloc, yloc, zloc, injectMassIn, injectVelocityIn, &
                            twind, dt, bgDens, snap_to_grid)
  implicit none

  real, intent(IN)    :: xloc, yloc, zloc
  real, intent(IN)    :: injectMassIn, injectVelocityIn, twind, dt
  real, intent(INOUT) :: bgDens
  logical, optional, intent(IN) :: snap_to_grid

  return
end subroutine GridInject_wind
