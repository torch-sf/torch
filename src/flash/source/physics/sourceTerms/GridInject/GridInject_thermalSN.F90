!!****f* source/physics/sourceTerms/GridInject/GridInject_thermalSN
!!
!! NAME
!!
!!  GridInject_thermalSN
!!
!! SYNOPSIS
!!
!!  call GridInject_thermalSN(
!!      real(IN)          :: xloc
!!      real(IN)          :: yloc
!!      real(IN)          :: zloc
!!      real(IN)          :: energy
!!      real(IN)          :: mass
!!      real(IN)          :: r_min
!!      real(IN)          :: r_max
!!      real(IN)          :: nms
!!      real(OUT)         :: r_exp
!!      real(OUT)         :: m_exp
!!      real(OUT)         :: rho_avg
!!  )
!!
!! DESCRIPTION
!!
!!  Deposit SN onto grid as thermal energy
!!
!! ARGUMENTS
!!
!!  xloc    : where to inject
!!  yloc    : where to inject
!!  zloc    : where to inject
!!  energy  : amount of energy to inject
!!  mass    : how much mass should lie in the explosion sphere?
!!  r_min   : minimum SN explosion radius
!!  r_max   : minimum SN explosion radius
!!  nms     : number of mass shells to search between r_min and r_max
!!
!! RETURNS
!!
!!  r_exp : explosion radius
!!  m_exp : mass within explosion sphere
!!  rho_avg : average density within explosion sphere
!!
!!***
subroutine GridInject_thermalSN (xloc, yloc, zloc, energy, mass, &
    r_min, r_max, nms, r_exp, m_exp, rho_avg)

  implicit none

  real, intent(IN)    :: xloc, yloc, zloc, energy, mass, r_min, r_max
  integer, intent(IN) :: nms
  real, intent(OUT)   :: r_exp, m_exp, rho_avg

  return
end subroutine GridInject_thermalSN
