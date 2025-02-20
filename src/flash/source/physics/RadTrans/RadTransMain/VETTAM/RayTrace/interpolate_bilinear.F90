
!-------------------------------------------------------------------------
!             SUBROUTINE FOR BILINEAR INTERPOLATION 
!-------------------------------------------------------------------------
#ifndef CHECKS_RT
elemental &
#endif
function interpolate_bilinear(q11,q21,q12,q22,           &
                              x1,x2,y1,y2,               &
                              x,y)
  !
  !  Linear interpolation between four grid points at (x,y)
  !    
  !                        q12-R2---q22--y2
  !                        |   .      |
  !                        |..........|..y
  !                        |   .      |
  !                        |   .      |
  !                        |   .      |
  !                        q11-R1---q21--y1
  !                        |   .      |
  !                        x1  x      x2
  !
  implicit none
  real :: interpolate_bilinear
  real, intent(in) :: q11,q21,q12,q22
  real, intent(in) :: x1,x2,y1,y2
  real, intent(in) :: x,y
  real :: R1,R2,f,a,b
  real, parameter :: precision=1.e+5*EPSILON(q11)
  !
  ! Check for precision errors
  !
!  if(ABS(x2-x1).LT.precision) then
!     R1 = q11
!     R2 = q12
!  else
     a = MAX(x2-x,0.)
     b = MAX(x-x1,0.)
     f = a / (a+b)
     R1 = f * q11 + (1.-f) * q21
     R2 = f * q12 + (1.-f) * q22

IF(MIN(q11,q21).le.precision*MAX(q11,q21)) R1 = MAX(q11,q21)
IF(MIN(q12,q22).le.precision*MAX(q12,q22)) R2 = MAX(q12,q22)
     !R1 = abs(x2-x)/abs(x2-x1) * q11 + abs(x-x1)/abs(x2-x1) * q21
     !R2 = abs(x2-x)/abs(x2-x1) * q12 + abs(x-x1)/abs(x2-x1) * q22
!  endif
  !
!  if(ABS(y2-y1).lt.precision) then
!     interpolate_bilinear = R1
!  else
     a = MAX(y2-y,0.)
     b = MAX(y-y1,0.)
     f = a / (a+b)
     interpolate_bilinear = f * R1 + (1.-f)* R2
IF(MIN(R1,R2).le.precision*MAX(R1,R2)) interpolate_bilinear=MAX(R1,R2)
     !interpolate_bilinear = abs(y2-y)/abs(y2-y1) * R1 + abs(y-y1)/abs(y2-y1) * R2
!  endif

#ifdef CHECKS_RT
  if(isNaN(interpolate_bilinear).or.(interpolate_bilinear.gt.1.e100.or.interpolate_bilinear < 0.0)) then
    write(*,*) 'interpolate_bilinear',interpolate_bilinear
    write(*,*) 'x1,x2:',x1,x2
    write(*,*) 'y1,y2:',y1,y2
    write(*,*) 'x ,y :',x ,y 
    write(*,*) 'q11,q21 :',q11,q21
    write(*,*) 'q12,q22 :',q12,q22
    call Driver_abortFlash('ERROR: error in function "interpolate_bilinear"') 
  end if
#endif

  return
  !
end function interpolate_bilinear

