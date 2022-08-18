!!****if* source/physics/sourceTerms/cool/molecules+dust/Cool_data
!!
!! NAME
!!  Cool_data
!!
!! SYNOPSIS
!!  Cool_data()
!!
!! DESCRIPTION
!!  Stores the local data for Source Term: cool/molecules+dust
!!
!! PARAMETERS
!!
!!
!! WRITTEN BY
!!  RB + DS 2012
!!
!!***

Module Cool_data

  real, save    :: T_cool_min, nd_cool_min, nd_cool_max, tstep_cool_factor
  real, save    :: T_max, gasConstant, pi, newton, sigma
  real, save    :: kB, gamma, gammam1, mp
  real, save    :: T_max_core, T_max_core_radius

  integer, PARAMETER :: TEMP_PTS = 320
  integer, PARAMETER :: DENS_PTS = 190
  real, save    :: cool_dat(DENS_PTS, TEMP_PTS, 3)
  real, save    :: T_min, nd_min

  logical, save :: useDustCool
  character(len=24), save :: he_int_method

  integer, external  :: get_cooling_data !, find
  real, external     :: get_dust_temperature
!   external      find, get_cooling_data, get_dust_temperature

interface
  integer function find(x, x0, N)
    integer, intent(in) :: N
    real, intent(in)    :: x(N), x0
  end function find
end interface


end Module Cool_data
