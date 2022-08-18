subroutine sphere_and_cell_frac(vol4,R4,x4,y4,z4,cellsize4)
! calculates the volume of intersection between a sphere of radius R,
! and a cube with side length cellsize with a difference of x, y, z
! between their centers.
implicit none
real,intent(in) :: R4,x4,y4,z4,cellsize4
real(kind=8) :: R,x,y,z,cellsize,vol
real,intent(out) :: vol4

! coordinates will be swapped internally so that we integrate over a shallow
! function. x and y in these variables refer to a new coordinate system
real(kind=8) :: xbegin, ybegin, xstep, ystep, xpos, ypos
real(kind=8) :: xnear, xfar, ynear, yfar
real(kind=8) :: taylorvol
integer :: xcount, ycount
integer, parameter :: numxsteps = 60
integer, parameter :: numysteps = 60

R = real(R4,8)/real(cellsize4,8)
x = real(x4,8)/real(cellsize4,8)
y = real(y4,8)/real(cellsize4,8)
z = real(z4,8)/real(cellsize4,8)

! Check if sphere encompasses whole cell
if (R**2 .gt. (abs(x)+5d-1)**2 + (abs(y)+5d-1)**2 + &
        (abs(z)+5d-1)**2) then
    vol = cellsize**3
! Check if sphere does not reach cell
else if (R**2 .lt. (abs(x)-5d-1)**2 + (abs(y)-5d-1)**2 + &
        (abs(z)-5d-1)**2) then
    vol = 0d0
else
    ! Determine coordinate with largest difference. The other two will be the
    ! limits of integration
    if ( abs(x) .gt. abs(y) .and. abs(x) .gt. abs(z) ) then
        xbegin = abs(y)-5d-1
        ybegin = abs(z)-5d-1
    else if (abs(y) .gt. abs(x) .and. abs(y) .gt. abs(z) ) then
        xbegin = abs(x)-5d-1
        ybegin = abs(z)-5d-1
    else
        xbegin = abs(x)-5d-1
        ybegin = abs(y)-5d-1
    end if

end if

xstep = 1d0/numxsteps
ystep = 1d0/numysteps

vol = 0d0
do xcount=0,numxsteps-1
    xpos = xbegin+xcount*xstep 
    do ycount=0,numysteps-1
        ypos = ybegin+ycount*ystep 
        vol = vol + taylorvol(R,xpos,xpos+xstep,ypos,ypos+ystep,-abs(z))
    end do
end do

vol4 = real(vol,4)
return

end subroutine sphere_and_cell_frac


real(kind=8) function sphere(R,x,y,offset)
real(kind=8) :: R,x,y,offset
sphere = 0d0
sphere = sqrt(R**2-x**2-y**2)+offset
return
end function sphere

! integrate the function "sphere" over the limits xn (xmin) to xm (xmax) and
! yn (ymin) to ym (ymax)
real(kind=8) recursive function taylorvol(R,xn,xm,yn,ym,offset) result(res)
real(kind=8) :: a,b,R,xn,xm,yn,ym,dx,dy,offset,sphere
real(kind=8),dimension(4) :: bounds
real(kind=8),dimension(4) :: subvols
real :: minref = 1e-5
a=(xm+xn)/2d0
b=(ym+yn)/2d0
dx = (xm-xn)/2d0
dy = (ym-yn)/2d0
res = 0d0

!primitive stop condition
if (dx .lt. minref .or. dy .lt. minref) then
    res = min(sphere(R,a,b,offset),1d0)*4*dx*dy
    return
end if

bounds(1) = sphere(R,xn,yn,offset)
bounds(2) = sphere(R,xn,ym,offset)
bounds(3) = sphere(R,xm,yn,offset)
bounds(4) = sphere(R,xm,ym,offset)

if (bounds(1) .gt. 0 .and. bounds(2) .gt. 0 .and. bounds(3) .gt. 0 .and. bounds(4) .gt. 0) then

    !f= sphere(R,a,b,offset)
    !fxx= (b**2-R**2)/(R**2-a**2-b**2)**1.5
    !fyy= (a**2-R**2)/(R**2-a**2-b**2)**1.5
    !res = f*(xm-xn)*(ym-yn) &
    !        + fxx/6d0*((xm-a)**3-(xn-a)**3)*(ym-yn) &
    !        + fyy/6d0*((ym-b)**3-(yn-b)**3)*(xm-xn)

    res = 4*dx*dy*sphere(R,a,b,offset)+2d0/3d0*(dx**4*(b**2-R**2)+dy**4*(a**2-R**2))/sphere(R,a,b,0d0)**3
    res = min(res,(xm-xn)*(ym-yn))
    return

else if (bounds(1) .lt. 0 .and. bounds(2) .lt. 0 .and. bounds(3) .lt. 0 .and. bounds(4) .lt. 0) then

    res = 0d0
    return

else 

    subvols(1) = taylorvol(R,xn,a,yn,b,offset)
    subvols(2) = taylorvol(R,xn,a,b,ym,offset)
    subvols(3) = taylorvol(R,a,xm,yn,b,offset)
    subvols(4) = taylorvol(R,a,xm,b,ym,offset)

    res = sum(subvols)
    res = min(res,(xm-xn)*(ym-yn))

    return

end if

end function taylorvol

