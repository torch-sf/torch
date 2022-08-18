module OdeData

! Set double precision type.
integer, parameter :: dp = kind(1.0d0)
! Max number of steps the integrator performs.
integer, save :: MAX_STEPS = 50000

real, save    :: atol, rtol, hmin, hmax

end module OdeData
