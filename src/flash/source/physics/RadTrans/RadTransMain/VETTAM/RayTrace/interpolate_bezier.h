!******************************************************************************
!*                       Code generated with sympy 1.0                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                       This file is part of 'project'                       *
!******************************************************************************


interface
ELEMENTAL REAL function Psiu_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2
end function
end interface
interface
ELEMENTAL REAL function Psi0_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2
end function
end interface
interface
ELEMENTAL REAL function Psid_0(a, b, g0, tau1, tau2)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2
end function
end interface
interface
ELEMENTAL REAL function Psiu_1(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
end function
end interface
interface
ELEMENTAL REAL function Psi0_1(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
end function
end interface
interface
ELEMENTAL REAL function Psiu_2(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
end function
end interface
interface
ELEMENTAL REAL function Psi0_2(a, b, g0, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: g0
REAL, intent(in) :: tau1
end function
end interface
interface
ELEMENTAL REAL function get_g0(a, b, tau1)
implicit none
REAL, intent(in) :: a
REAL, intent(in) :: b
REAL, intent(in) :: tau1
end function
end interface
interface
ELEMENTAL REAL function get_Sc(S0, Sd, Su, tau1, tau2)
implicit none
REAL, intent(in) :: S0
REAL, intent(in) :: Sd
REAL, intent(in) :: Su
REAL, intent(in) :: tau1
REAL, intent(in) :: tau2
end function
end interface

