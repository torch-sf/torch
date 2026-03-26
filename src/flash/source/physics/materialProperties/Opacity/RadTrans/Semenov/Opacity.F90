!!****f* source/physics/materialProperties/Opacity/RadTrans/Semenov
!!
!! NAME
!!
!!  Opacity
!!
!! SYNOPSIS
!!
!!  call Opacity (real    (in)  :: soln(:),
!!                real    (out) :: opac)
!!
!! DESCRIPTION
!!
!!  get the gray opacity according to Semenov 2003 
!!
!! ARGUMENTS
!!
!!   soln              : The solution vector for the cell
!!   opac              : the gray (Rosseland and Planck mean) opacity
!!
!!***

subroutine Opacity (soln, opac)

#include "Flash.h"

  use SemenovOpacities, ONLY: getOpacity

  implicit none
  
  real, intent(out) :: opac
  real, intent(in), dimension (:) :: soln
  
  real :: rho, temp, kappa

  rho  = soln(DENS_VAR)
#ifdef TEMP_VAR  
  temp = soln(TEMP_VAR)
#endif
  call getOpacity(temp, rho, kappa)

  if (kappa .lt. 1e-40) then
     write(*,*) 'kappa', kappa
     write(*,*) 'rho, temp = ', rho, temp
     call Driver_abortFlash('Opacity: kappa < 1.e-40. ' &
      // 'Probably rho or temp out of bounds (rho \in [2.e-18,2.e-7], temp \in [5., 1.E+5]).')
  endif

  opac = kappa * soln(DENS_VAR)

  return
end subroutine Opacity
