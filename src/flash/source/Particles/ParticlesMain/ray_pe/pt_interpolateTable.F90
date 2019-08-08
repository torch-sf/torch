!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   interpolate logarithmically or linearly on some value in some table
!!   interpolation based on InterpolateTable.F90 in WindDrivingMain by Andrea Gatto

subroutine pt_interpolateTable(valTable,intTable,nval,offset,indx,method,outval)

  implicit none

#include "constants.h"
#include "Flash.h"
  
  real, intent (OUT)   :: outval
! log or linear, offset in table, size of passed array
  integer, intent (IN) :: method, indx, nval
  real,    intent (IN) :: offset
! valTable is data to interpolate on
! intTable is spacing on which to interpolate
  real, dimension(nval), intent (IN) :: valTable, intTable

! log
  if(method .eq. 1 ) then
! h is bin size, T elapsed time, t is current time
! this is  [x(t)*x(t+h)/x(t)]^(T - t(t)) 
    if(intTable(indx+1) .eq. intTable(indx) ) then
      outval = valTable(indx+1)
    else
      outval = valtable(indx)*(valtable(indx+1)/valtable(indx))**( (offset-intTable(indx) )/( intTable(indx+1)-intTable(indx)) )
    endif
  endif

! linear 
  if(method .eq. 2 ) then
! is there anything to interpolate?
    if(intTable(indx+1) .eq. intTable(indx) ) then
      outval = valTable(indx+1)
    else
      outval = valTable(indx)+(valTable(indx+1)-valTable(indx))*(offset-intTable(indx))/(intTable(indx+1)-intTable(indx))
    endif    
  endif

  return
end subroutine pt_interpolateTable
