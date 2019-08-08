! ODE integrator based on the methods in Numerical Reciepes.

! Written by Joshua Wall
! Drexel University.
!#define debug_int
!#define debug_driver

module OdeInt

use OdeData

contains

! Sets up initial values for the OdeDriver that would be set
! regardless of the stepper (method) we use.
subroutine OdeInit(atol_in, rtol_in, hmin_in, hmax_in, MAX_STEPS_IN)

implicit none

real(dp), intent(in)           :: atol_in, rtol_in, hmin_in
real(dp), optional, intent(in) :: hmax_in
integer, optional, intent(in)  :: MAX_STEPS_IN

! Maximum number of steps the integrator will take.
if (present(MAX_STEPS_IN)) then
    MAX_STEPS = MAX_STEPS_IN
else
    MAX_STEPS = 10000000
end if

if (present(hmax_in)) then
    hmax = hmax_in
else
    hmax = hmin_in*1d5
end if

atol = atol_in ! Absolute tolerance for error.
rtol = rtol_in ! Relative error.
hmin = hmin_in ! Absolute smallest step size (could be zero).
hmax = hmax_in

end subroutine OdeInit

subroutine OdeDriver(f, yinit, x1, x2, hinit, yend)
! Manages the integration of an ODE from x1 to x2.
! Calls the underlying stepper routine to do the integration.

implicit none

real(dp), intent(in) :: yinit, x1, x2, hinit
real(dp), intent(out):: yend

real(dp)             :: x, y, h, x_new, y_new
integer              :: step

interface
    function f(x_par, y_par)
    implicit none
    integer, parameter    :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface

y = yinit
x = x1
h = hinit
step = 0

do while (x < x2)

#ifdef debug_driver
    write(*, '(A,ES13.3E3,X,A,ES13.3E3)') "[OdeDriver]: At x =", x, "y =", y
    write(*, '(A,ES13.3E3,X,A,ES13.3E3)') "[OdeDriver]: x2 =", x2, "hmax =", hmax
    if (hmax < 1e4) stop
#endif
    call OdeStepperDoPr5(f, y, x, x2, h, x_new, y_new)
    
    x = x_new
    y = y_new
    
    if (y < 0.0) then
        write(*,'(A)') "[OdeDriver]: y < 0."
	write(*,'(A,ES13.3E3)') "[OdeDriver]: y=", y
	write(*,'(A,ES13.3E3)') "[OdeDriver]: x=", x
	write(*,'(A,ES13.3E3)') "[OdeDriver]: h=", h
	stop
    end if
    
!    if (x > x2+hmin) then
!        write(*,'(A)') "[OdeDriver]: x now larger than x_end point."
!	write(*,'(A,ES13.3E3)') "[OdeDriver]: y=", y
!	write(*,'(A,ES13.3E3)') "[OdeDriver]: x=", x
!	write(*,'(A,ES13.3E3)') "[OdeDriver]: h=", h
!	stop
!    end if

    if (step > MAX_STEPS) then
	x = x_new
	y = y_new
        write(*,'(A)') "[OdeDriver]: Max steps exceeded."
	write(*,'(A,ES13.3E3)') "[OdeDriver]: y=", y
	write(*,'(A,ES13.3E3)') "[OdeDriver]: x=", x
	write(*,'(A,ES13.3E3)') "[OdeDriver]: h=", h
	stop
	exit
    end if

    step = step + 1
end do

#ifdef debug_driver
write(*,'(A,ES13.3E3)') "[OdeDriver]: y=", y
write(*,'(A,ES13.3E3)') "[OdeDriver]: x=", x
write(*,'(A,ES13.3E3)') "[OdeDriver]: h=", h
#endif
yend = y


end subroutine OdeDriver
  
! Stepper that uses the RK5 method from Numerical Recipes with the coeff
! from Dormand and Price 85 to integrate the equations.
subroutine OdeStepperDoPr5(f, ystart, xx1, xmax, h1, xx2, yend)

implicit none

real(dp), intent(in)    :: ystart, xx1, xmax ! y_initial, x_initial, absolute max x allowed (end of integration).
real(dp), intent(inout) :: h1                ! initial step
real(dp), intent(out)   :: yend, xx2         ! final y answer, final x at which this y occurs.


integer, parameter  :: max_stepper_count = 100
real(dp), parameter :: order = 5.0_dp, beta = 0.4_dp / order
real(dp), parameter :: alpha = 1.0_dp / order - 0.75_dp*beta
real(dp)            :: error, olderr
real(dp)            :: h, y, yn1, x
integer             :: nok, nbad, nvar
! Smallest stepsize.
real(dp) :: s1=0.9, del, oldx, oldy, oldh


interface
    function f(x_par, y_par)
    implicit none
    integer, parameter    :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface

    nbad  = 0
    h     = h1
    y     = ystart
    x     = xx1
    error = 1.0_dp
    
    do while (nbad < max_stepper_count)

! Note it is important that we capture the old timestep before any
! correction to prevent overstepping.
	oldy   = y
	oldx   = x
	oldh   = h
	olderr = error
! Make sure we don't step past the end time or take a step larger
! than the maximum set by the user. Note this assumes that x
! is an increasing variable for this to work properly.
	if ((x+h) > xmax) &
	    h = xmax-x
	    
	call DoPr5(f,x,y,h,yn1,del)

#ifdef debug_int
	write(*,'(A,ES13.3E3)') "[OdeStepper]: y = ", y
	write(*,'(A,ES13.3E3)') "[OdeStepper]: yn1 = ", yn1
	write(*,'(A,ES13.3E3)') "[OdeStepper]: x = ", x
	write(*,'(A,ES13.3E3)') "[OdeStepper]: h = ", h
	write(*,'(A,ES13.3E3)') "[OdeStepper]: del = ", del
	if (yn1 /=yn1) stop
#endif
	!error = abs(del)/(atol+max(abs(y),abs(yn1))*rtol)
	error = abs(del)/(atol+abs(yn1)*rtol)
#ifdef debug_int
	write(*,'(A,ES13.3E3)') "[OdeStepper]: error = ", error
#endif
	if (error == 0.0_dp) then
#ifdef debug_int
	    print*, "[OdeStepper]: Error = 0, h, hmax =", h, hmax
#endif
	    h = hmax
	else
	
	    h = max(s1*oldh*error**(-alpha)*olderr**(beta),hmin)
	
	end if
	!write(*,'(A,ES13.3E3)') "[OdeStepper]: dt = ", h
	
	!stop
	
	if (error <= 1.0_dp) then
	
	    xx2  = x+oldh
	    yend = yn1
	    h1   = h
	    
	        
	    if (yn1 < 0.0) then
		write(*,'(A,ES13.3E3)') "yn1=", yn1
		write(*,'(A,ES13.3E3)') "xx2=", xx2
		write(*,'(A,ES13.3E3)') "xmax=", xmax
		write(*,'(A,ES13.3E3)') "h1=", oldh
		stop
	    end if
#ifdef debug_int
	    write(*,'(A,ES13.3E3)') "[OdeStepper]: Returning. yend=", yend
#endif
	    return
	    
	else
	    
	    y = oldy
	    x = oldx
	    ! Don't allow steps larger than x_end-x.
	    h = min(h,hmax)
	
	end if
	
	nbad = nbad + 1
	   !write(*, '(A,I8)') "[OdeStepper]: nbad =", nbad 
	if (nbad .ge. max_stepper_count) then
	   write(*, '(A,I8)') "[OdeStepper]: nbad so bad! nbad =", nbad
	   write(*,'(A,ES13.3E3)') "[OdeStepper]: y = ", y
	   write(*,'(A,ES13.3E3)') "[OdeStepper]: yn1 = ", yn1
	   write(*,'(A,ES13.3E3)') "[OdeStepper]: x = ", x
	   write(*,'(A,ES13.3E3)') "[OdeStepper]: h = ", h
	   write(*,'(A,ES13.3E3)') "[OdeStepper]: del = ", del
	   stop
	end if
	
    end do

end subroutine OdeStepperDoPr5

subroutine DoPr5(f,x,y,h,yn1,del)

! Input independent and dependent variables.
real(dp),intent(in)  :: x,y,h

! 6th order result of RK y_n+1.
real(dp),intent(out) :: yn1

! error from xn1 - xsn1
real(dp),intent(out) :: del

! c_i coefficients for the actual RK method.
real(dp), parameter :: c2=1.0_dp/5.0_dp, c3=3.0_dp/10.0_dp, c4=4.0_dp/5.0_dp, &
                       c5=8.0_dp/9.0_dp, c6=1.0_dp, c7=1.0_dp
! a_ij coefficients for the actual RK method.
real(dp), parameter :: a21=1.0_dp/5.0_dp, &
		       a31=3.0_dp/40.0_dp, a32=9.0_dp/40.0_dp, &
		       a41=44.0_dp/45.0_dp, a42=-56.0_dp/15.0_dp, a43=32.0_dp/9.0_dp, &
		       a51=19372.0_dp/6561.0_dp, a52=-25360.0_dp/2187.0_dp, &
		       a53=64448.0_dp/6561.0_dp, a54=-212.0_dp/729.0_dp, &
		       a61=9017.0_dp/3168.0_dp, a62=-355.0_dp/33.0_dp, &
		       a63=46732.0_dp/5247.0_dp, a64=49.0_dp/176.0_dp, &
		       a65=-5103.0_dp/18656.0_dp, &
		       a71=35.0_dp/384.0, a72=0.0_dp, a73=500.0_dp/1113.0_dp, &
		       a74=125.0_dp/192.0_dp, a75=-2187.0_dp/6784.0_dp, &
		       a76=11.0_dp/84.0_dp
		       
! b_i coefficients for the final 6th order sum of all the RK terms. Note that since we use FSAL a7i = bi.
real(dp), parameter  :: b1=35.0_dp/384.0, b2=0.0_dp, b3=500.0_dp/1113.0_dp, &
		        b4=125.0_dp/192.0_dp, b5=-2187.0_dp/6784.0_dp, &
		        b6=11.0_dp/84.0_dp, b7=0.0_dp

! bstar_i coefficients for the embedded formula estimation at 5th order.
real(dp), parameter  :: bs1=5179.0_dp/57600.0_dp, bs2=0.0_dp, bs3=7571.0_dp/16695.0_dp, &
		        bs4=393.0_dp/640.0_dp, bs5=-92097.0_dp/339200.0_dp, &
		        bs6=187.0_dp/2100.0_dp, bs7=1.0_dp/40.0_dp
		       
! ki returns from each evaluation
real(dp)             :: k1, k2, k3, k4, k5, k6
real(dp), save       :: k7 = 0.0_dp

! 5th order result of embedded formula y*_n+1.
real(dp)             :: ysn1

logical, save        :: first_call = .false.

interface
    function f(x_par, y_par)
    implicit none
    integer, parameter    :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface

if (first_call) then
    k1 = h*f(x,y)
    first_call = .false.
else
    k1=k7
end if

    k2 = h*f(x+c2*h,y+a21*k1)
    k3 = h*f(x+c3*h,y+a31*k1+a32*k2)
    k4 = h*f(x+c4*h,y+a41*k1+a42*k2+a43*k3)
    k5 = h*f(x+c5*h,y+a51*k1+a52*k2+a53*k3+a54*k4)
    k6 = h*f(x+c6*h,y+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    k7 = h*f(x+c7*h,y+a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    
    ! 6th order result. Note we return the higher order result, even
    ! though the method is only truely 5th order (the so called
    ! 'local extrapolation').
    yn1  = y + b1*k1  + b2*k2  + b3*k3  + b4*k4  + b5*k5  + b6*k6 + b7*k7
    !ysn1  = y + b1*k1  + b2*k2  + b3*k3  + b4*k4  + b5*k5  + b6*k6 + b7*k7
    ! 5th order result derived from same k_i's (embedded formula).
    ysn1 = y + bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4 + bs5*k5 + bs6*k6 + bs7*k7
    !yn1 = y + bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4 + bs5*k5 + bs6*k6 + bs7*k7
    ! Error between the methods.
    del = yn1 - ysn1

    return

end subroutine DoPr5


end module OdeInt
