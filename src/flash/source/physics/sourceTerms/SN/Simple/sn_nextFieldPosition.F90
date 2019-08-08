!!****f* source/physics/sourceTerms/SN/Simple/sn_nextFieldPosition
!!
!! NAME
!!
!!  sn_nextFieldPosition
!!
!!
!! SYNOPSIS
!!
!!  sn_nextFieldPosition (  int(IN) : sn_type
!!                        real(OUT) : x
!!                        real(OUT) : y
!!                        real(OUT) : z )
!!
!! DESCRIPTION
!!
!!  Where will the next field SN blow up?
!!  This subroutine is state-ful; repeat calls may return different locations,
!!  and behavior depends upon value of sn_nSN.
!!
!! ARGUMENTS
!!
!!  sn_type     Supernova type
!!  x, y, z     Three directional coordinates of next supernova
!!
!!***

subroutine sn_nextFieldPosition (sn_type, x, y, z)

#include "Flash.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_abortFlash

  use SN_data, only: sn_imin, sn_imax, sn_jmin, sn_jmax, sn_kmin, sn_kmax, &
                     sn_nSN, &
                     sn_fieldMode, sn_stratifySN, sn_hstar1, sn_hstar2, &
                     sn_single_x, sn_single_y, sn_single_z

  use mtmodSN, ONLY : grndSN

  implicit none

  !character(len=6), intent(IN)  :: mode
  integer, intent(IN)           :: sn_type
  real, intent(OUT)             :: x, y, z

  real :: hscale

  select case (sn_fieldMode)
  case("random")

    x = grndSN()*(sn_imax-sn_imin) + sn_imin
    y = grndSN()*(sn_jmax-sn_jmin) + sn_jmin

    if (sn_stratifySN) then

      if (sn_type == 1) then
        hscale = sn_hstar1
      else if (sn_type == 2) then
        hscale = sn_hstar2
      else
        call Driver_abortFlash("Bad sn_type for stratified random SNe.")
      end if

      ! draw until the SN lands in simulation z domain
      z = sn_kmax
      do while (z .ge. sn_kmax .or. z .le. sn_kmin)
        z = grndSN()*2.0 - 1.0   !-1 to 1
        z = abs(z)/z * log(abs(z)) * hscale ! inverse transform, sgn(z)*ln(|z|)*hscale
      end do

    else

      z = grndSN()*(sn_kmax-sn_kmin) + sn_kmin

    end if

  case("single")

    x = sn_single_x
    y = sn_single_y
    z = sn_single_z

    if (sn_nSN .ge. 2) then
      call Driver_abortFlash("2nd call to SN_nextFieldPosition(mode=single).")
    end if

  case("double")

    if (sn_nSN .eq. 1) then
      x = 0d0
      y = 0d0
      z = 0d0
    else if (sn_nSN .eq. 2) then
      x = 0.6*(sn_imax-sn_imin) + sn_imin
      y = 0.6*(sn_jmax-sn_jmin) + sn_jmin
      z = 0.6*(sn_kmax-sn_kmin) + sn_kmin
    else
      call Driver_abortFlash("3rd call to SN_nextFieldPosition(mode=double).")
    end if

  case default

    call Driver_abortFlash("Bad fieldMode for SN_nextFieldPosition.")

  end select

  return
end subroutine sn_nextFieldPosition
