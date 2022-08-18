subroutine calc_ionization(dt,temp0,eldens,ndens,xh,xhold, phih)
     
! Calculates time dependent ionization state for hydrogen and the electron
! density (cf. Schmidt-Voigt & Koeppen 1987).

! magic numbers ahoy
  use rt_data, only : albpow, bh00, temph0, colh0
  implicit none

  real,    intent(in)    :: dt
  real,    intent(in)    :: temp0
  real,    intent(inout) :: eldens
  real,    intent(in)    :: ndens
  real,    intent(out)   :: xh(0:1)
  real,    intent(in)    :: phih
  real,    intent(in)    :: xhold(0:1)

  real    :: brech0, sqrtt0, acolh0
  real    :: aih0, delth, eqxh0, eqxh1
  real    :: deltht, ee

!-------------------------------------------------------------------------------

! Hydrogen recombination rate at the local temperature.

  brech0 = bh00*(temp0/1.0e4)**albpow
!  brech0 = bh00!*(temp0/1.0e4)**albpow

! Hydrogen collisional ionization rate at the local temperature

  sqrtt0 = sqrt(temp0)
  acolh0 = colh0 * sqrtt0 * exp(-temph0/temp0)

! Calculate the ionization fractions for a constant electron density.

  aih0   = phih + eldens * acolh0  ! kph + n_e*C_H
  delth  = aih0 + eldens * brech0  ! kph + n_e*C_H + n_e*a_b
  eqxh1  = aih0 / delth  ! kph + n_e*C_H / (kph + n_e*C_H + n_e*a_b )
  eqxh0  = eldens * brech0 / delth ! n_e*a_b / ( kph + n_e*C_H + n_e*a_b )
  deltht = delth * dt ! 
  ee     = exp(-deltht)
  xh(1)  = (xhold(1) - eqxh1) * ee + eqxh1
  xh(0)  = (xhold(0) - eqxh0) * ee + eqxh0

  call electrondens(eldens,ndens,xh)

! Determine neutral densities (take care of precision fluctuations)

  if (xh(0).lt.1e-40.and. abs(xh(0)) .lt.1.0e-10) then
#ifdef DEBUG_DORIC
     write(*,*) 'doric_calc_ionization: precision problem:'
     write(*,*) xhold,xh
     write(*,*) phih/(xhold(0)*ndens),eldens * acolh0,eldens * brech0
     write(*,*) temp0
#endif
     xh(0)=1e-40
  endif

  return
end subroutine calc_ionization

!=========================================
! helper functions
!=========================================
subroutine electrondens(eldens,ndens,xh)
!
! Find electron density.
  use rt_data, only : abu_c
  implicit none
!
  real, intent(out) :: eldens
  real, intent(in)  :: ndens
  real, intent(in)  :: xh(0:1)
!
  eldens = ndens*(xh(1)+abu_c)
!
end subroutine electrondens
