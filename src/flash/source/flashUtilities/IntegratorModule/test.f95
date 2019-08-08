program test

!use num_methods
!use runge_kutta4
use OdeInt

implicit none

!integer, parameter    :: dp=KIND(1d0)

integer       :: i, num_steps=100
real(dp)      :: x=0.0, t=0.0, lam=0.0, dt = 0.1
real(dp)      :: x_new, t_new, error
character(len=*), parameter    :: frmt='(es15.3,t20,es15.3)'
logical       :: adpt = .true.

call OdeInit(1d-3, 1d-4, 1d-10)

open(unit=10, file='test.dat')

write(10,frmt) x, t

if (adpt) then

do while (t < 10.0)    

!    call ad_rk4(dx_dt, x, t, dt, lam, x_new, t_new)
    t_new = t + 0.1
    call OdeDriver(line, x, t, t_new, dt, x_new)
    x = x_new
    t = t_new
!    write(*,*) x_new, t_new
    write(*,'(A,ES13.3E3)') "[test]: Current dt =", dt
    write(10,frmt) x, t
    !stop
   
    
end do

else

do while (t < 10.0)    

!    call ad_rk4(dx_dt, x, t, dt, lam, x_new, t_new)
    t_new = t + dt
    call DoPr5(dy_dx, t, x, dt, x_new, error)
    x = x_new
    t = t_new
!    write(*,*) x_new, t_new
!    write(*,'(A,ES13.3E3)') "[test]: Current dt =", dt
    write(10,frmt) x, t
    
end do

end if

stop

contains

function dy_dx(x_par, y_par)

    real(dp), intent(in)     :: x_par, y_par
    real(dp)                 :: dy_dx
    
    dy_dx = -y_par
    
end function dy_dx

function line(x_par, y_par)

    real(dp), intent(in)     :: x_par, y_par
    real(dp)                 :: line
    
    line = x_par

end function line

!function rk4(f, x, t, dt, param)

!    implicit none
    
!    real*8      :: rk4, x, t, dt, param
!    real*8      :: F1, F2, F3, F4    
!    interface
!        function f(x_par, t_par, param)
!        implicit none
!        real*8, intent(in) :: x_par, t_par, param
!        real*8              :: f
!        end function f
!    end interface

!    F1 = f(x, t, param)
!    F2 = f(x + 0.5*dt*F1, t + 0.5*dt, param)
!    F3 = f(x + 0.5*dt*F2, t + 0.5*dt, param)
!    F4 = f(x + dt*F3, t + dt, param)

!    rk4 = x + (1.0/6.0)*dt*(F1 + 2.0*F2 + 2.0*F3 + F4)
    
!end function rk4

end program test
