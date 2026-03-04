!******************************************************************************
!*                       Code generated with sympy 1.0                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                       This file is part of 'project'                       *
!******************************************************************************

ELEMENTAL REAL function Psiu_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2

Psiu_0 = (a**2*g0*tau1**2 - a**2*tau1**2 - 2*a*g0*tau1**2 - a*g0*tau1* &
      tau2 - 2*a*g0*tau1 + 2*a*tau1**2 + a*tau1*tau2 + 2*a*tau1 + b**2* &
      tau1**2 - 2*b*tau1**2 - b*tau1*tau2 - 2*b*tau1 + g0*tau1**2 + g0* &
      tau1*tau2 + 2*g0*tau1 + g0*tau2 + 2*g0)/(tau1*(tau1 + tau2))

end function

ELEMENTAL REAL function Psi0_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2

Psi0_0 = -a**2*g0*tau1/tau2 + a**2*tau1/tau2 + a*g0*tau1/tau2 + a*g0 + 2 &
      *a*g0/tau2 - a*tau1/tau2 - a - 2*a/tau2 - b**2*tau1/tau2 + b*tau1 &
      /tau2 + b + 2*b/tau2 - g0/tau2 - g0/tau1 - 2*g0/(tau1*tau2)

end function

ELEMENTAL REAL function Psid_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2

Psid_0 = (a**2*g0*tau1**2 - a**2*tau1**2 - a*g0*tau1**2 - 2*a*g0*tau1 + &
      a*tau1**2 + 2*a*tau1 + b**2*tau1**2 - b*tau1**2 - 2*b*tau1 + g0* &
      tau1 + 2*g0)/(tau2*(tau1 + tau2))

end function

ELEMENTAL REAL function Psiu_1(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1

Psiu_1 = (2*g0 + tau1**2*(a**2*g0 - a**2 - 2*a*g0 + 2*a + b**2 - 2*b + &
      g0) + 2*tau1*(-a*g0 + a - b + g0))/tau1**2

end function

ELEMENTAL REAL function Psi0_1(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1

Psi0_1 = (-2*g0 + tau1**2*(-a**2*g0 + a**2 + 2*a*g0 - 2*a - b**2 + 2*b) &
      + 2*tau1*(a*g0 - a + b - g0))/tau1**2

end function

ELEMENTAL REAL function Psiu_2(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1

Psiu_2 = (-2*g0 + tau1**2*(-a**2*g0 + a**2 - b**2 + g0) + 2*tau1*(a*g0 - &
      a + b))/tau1**2

end function

ELEMENTAL REAL function Psi0_2(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1

Psi0_2 = (2*g0 + tau1**2*(a**2*g0 - a**2 + b**2) + 2*tau1*(-a*g0 + a - b &
      ))/tau1**2

end function

ELEMENTAL REAL function get_g0(a, b, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: tau1

get_g0 = -exp(tau1*(a - b)) + 1

end function

ELEMENTAL REAL function get_Sc(S0, Sd, Su, tau1, tau2)
implicit none
REAL, intent(in) :: S0
REAL, intent(in) :: Sd
REAL, intent(in) :: Su
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2

get_Sc = S0 - 1.0d0/2.0d0*tau1*(tau1*(-S0 + Sd)/(tau2*(tau1 + tau2)) + &
      tau2*(S0 - Su)/(tau1*(tau1 + tau2)))

end function
