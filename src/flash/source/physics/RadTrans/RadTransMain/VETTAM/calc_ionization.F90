module calc_ion

! A re-writing of Christian's original ionization calculation. Here
! we solve the equation implicitly, since it's 1) guaranteed to be 
! stable and 2) allows much larger timesteps. Also, we estimate the
! timestep using the second derivative explicitly, then subcycle over
! this until we hit the hydro timestep. This will make solving for the
! new ionized fractions much, much faster than before.
! - Josh Wall and Stephen McMillan, Drexel, 05-2016

! Note: Christian's original comments lack the -JW at the end.

! Calculates time dependent ionization state for hydrogen and the electron
! density (cf. Schmidt-Voigt & Koeppen 1987).

contains

subroutine calc_ionization(dt, del_t_new, temp0, ndens, xh, xhold, k_ion)

#ifdef ONE_CELL_TESTING
#define debug
#endif
#define debug2

  use rt_data, only : albpow, bh00, temph0, colh0
  implicit none

  real,    intent(in)    :: dt
  real,    intent(in)    :: temp0
  real,    intent(out)   :: del_t_new
!  real,    intent(inout) :: eldens
  real,    intent(in)    :: ndens
  real,    intent(out)   :: xh(0:1)
  real,    intent(in)    :: k_ion ! the artist formerly known as phih
  real,    intent(in)    :: xhold(0:1)

  real    :: alpha_b, sqrtt0, acolh0
  real    :: x0, x1
  real    :: A_coeff, B_coeff, C_coeff
  real    :: del_t, sqroot
  real, parameter :: subcycle_factor=0.8
! min and max valid recombination coefficients. Its valid from 30 K to 30,000 K.
  real, parameter :: max_alpha_b = 2.54e-13*(3d0/1d4)**(-0.8163-0.0208*log(3d0/1d4))
  real, parameter :: min_alpha_b = 2.54e-13*(3d4/1e4)**(-0.8163-0.0208*log(3d4/1e4))
!  real    :: aih0, delth, eqxh0, eqxh1
!  real    :: deltht, ee

!-------------------------------------------------------------------------------

! Current ionization fraction.
  x0 = xhold(1)
  x1 = 0.0
  del_t = 0.0
  del_t_new = 0.0

! Hydrogen recombination rate at the local temperature.

!  alpha_b = bh00*(temp0/1.0e4)**albpow
!  brech0 = bh00*(temp0/1.0e4)**albpow

! From Draine 2008, eqn 14.6
! This should be more accurate (although Draine notes its still
! not totally correct much past 1e4). Difference b/t this and
! Christian's implementation is this has less recombinations at
! high temperature, allowing for more ionization when the gas is
! close to totally ionized, but the recombination amount still increases
! with decreasing temperature. Note also that Draine makes a
! point out of the fact that no single power law can capture
! recombination properly. -JW

  alpha_b = 2.54e-13*(temp0/1e4)**(-0.8163-0.0208*log(temp0/1e4))
  alpha_b = max(min_alpha_b,min(max_alpha_b,alpha_b))


if (alpha_b .gt. 1e10) then
   print*, "a_b too big!"
   call flush(6)
   stop
end if

! Hydrogen collisional ionization rate at the local temperature

  sqrtt0 = dsqrt(temp0)

  acolh0 = colh0 * sqrtt0 * max(0.0d0,dexp(-temph0/temp0))


! Calculate the ionization fractions for a constant electron density.

! Here I reformulate the equation as quadratic in x1, where x1 is the
! new ionization fraction we calculate implicitly at the hydro timestep.
! We then check if this timestep was too large by estimating the timestep
! using the second derivative of x1. If it was, we subcycle at the
! estimated timestep until we hit the hydro dt. - JW

! The solution of the quadratic form of the implict equation, where the
! derivative is evaulated at x1 (not x0!). So the form:
! A*(x1^2)+B*x+C=0 A=(acolh0-alpha_B)*n_H*dt B=k_ion*dt-1 C=x0
! gives the solution (we only consider the positive one): - JW

! These coefficients based on Christian's paper.
!  A_coeff = (acolh0-alpha_B)*ndens*dt
!  B_coeff =  k_ion*dt-1d0
!  C_coeff =  x0

! These coefficients based on my own calculations starting
! with the evolution of neutral hydrogen. - JW

!  A_coeff = -(acolh0+alpha_B)*ndens*dt
!  B_coeff =  (acolh0*ndens-k_ion)*dt-1d0
!  C_coeff =  k_ion*dt + x0
  
  A_coeff =  (acolh0+alpha_B)*ndens*dt
  B_coeff =  (-acolh0*ndens+k_ion)*dt+1d0
  C_coeff =  -k_ion*dt - x0  
  
  sqroot = B_coeff**2.0 - 4.0*A_coeff*C_coeff
  x1 = (-B_coeff + sqrt(sqroot))/(2.0*A_coeff)
  
  if (x1 .lt. 0.0d0) then
    print*, "Heads up, x1 < 0."
    x1 = min(1d-8, x1)
  end if
  
  if (sqroot .lt. 0.0) then
  print*, "About to make imaginary stuff. sqroot =", sqroot
  print*, "x1 =", x1
  call flush(6)
  stop
  end if
  
! Now check the timestep and see if we did something crazy.
! This form is based on f=dx/dt, and del_t = f/f' - JW

! NOTE: k_ion = dN / n_h(1-x)*dx^3*dt and in dx/dt = k_ion*(1-x),
! so in d^2 x/ dt^2 k_ion is a constant and drops out. - JW

!  write(*, '(A, 4ES13.3)') 'acolh0, alpha_b, x1, ndens', acolh0, alpha_b, x1, ndens

  del_t_new = subcycle_factor*1.0/max(abs(acolh0*ndens-2.0*(acolh0+alpha_b)*ndens*x1), 1d-50)
  del_t = del_t_new

!    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
!                            del_t type 1 =", del_t_new


#ifdef debug

!if (temp0 .gt. 1e5) then
!  write(*,'(A,F12.3)') "[calc_ionization]: Inside subcycle. &
!                            New ionization fraction =", x1
!#endif
!#ifdef debug2
  write(*,'(A,ES22.10E3)') "[calc_ionization]: &
                          Old ionization fraction =", x0
  write(*,'(A,ES22.10E3)') "[calc_ionization]: &
                          New ionization fraction =", x1
  write(*,'(A,ES22.10E3)') "[calc_ionization]: &
                          Old neutral fraction =", 1d0 - x0
  write(*,'(A,ES22.10E3)') "[calc_ionization]: &
                          New neutral fraction =", 1d0 - x1
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            A =", A_coeff
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            B =", B_coeff
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            C =", C_coeff
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            k_ion =", k_ion
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            alpha_b =", alpha_b
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            C_col =", acolh0
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                          ndens =", ndens
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                          temp0 =", temp0
  write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            sub_dt =", del_t_new

  call flush(6)
  !stop
!endif
#endif

if (del_t .lt. 0.0) then
  print*, "del_t is negative, stopping. del_t=", del_t
!  stop
end if


if (x1 .lt. 0.0) then
  print*, "WARNING! x1 is less than zero! x1=", x1
  call flush(6)
  stop
  x1 = 1e-40
end if

  xh(1) = x1
  xh(0) = (1d0 - x1)

#ifdef debug2
  if (isnan(xh(0)) .or. isnan(xh(1))) then

    print*, "X1 or X0 is nan!."
    write(*,'(A,F12.3)') "[calc_ionization]: &
                            Old ionization fraction =", x0
    write(*,'(A,F12.3)') "[calc_ionization]: &
                            New ionization fraction =", x1
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            del_t =", del_t
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            A =", A_coeff
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            B =", B_coeff
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            C =", C_coeff
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            k_ion =", k_ion
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            alpha_b =", alpha_b
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            C_col =", acolh0
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            ndens =", ndens
    write(*,'(A,ES12.3E3)') "[calc_ionization]: &
                            temp0 =", temp0

  end if
#endif 
!  aih0   = phih + eldens * acolh0  ! kph + n_e*C_H

!  delth  = aih0 - eldens * brech0  ! kph + n_e*C_H - n_e*a_b
!  eqxh1  = aih0 / delth  ! kph + n_e*C_H / (kph + n_e*C_H + n_e*a_b )
!  eqxh0  = eldens * brech0 / delth ! n_e*a_b / ( kph + n_e*C_H + n_e*a_b )
!  deltht = delth * dt ! 
!  ee     = exp(-deltht)
!  xh(1)  = (xhold(1) - eqxh1) * ee + eqxh1
!  xh(0)  = (xhold(0) - eqxh0) * ee + eqxh0

!  call electrondens(eldens,ndens,xh)

!Determine neutral densities (take care of precision fluctuations)

!  if (xh(0).lt.1e-40.and. abs(xh(0)) .lt.1.0e-10) then
!#ifdef DEBUG_DORIC
!     write(*,*) 'doric_calc_ionization: precision problem:'
!     write(*,*) xhold,xh
!     write(*,*) phih/(xhold(0)*ndens),eldens * acolh0,eldens * brech0
!     write(*,*) temp0
!#endif
!     xh(0)=1e-40
!  endif

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
end module calc_ion
