!!****ih* source/physics/RadTrans/RadTransMain/RayRad/SinksRad/get_ionization.F90
!!
!! NAME
!!
!!   get_ion
!!
!! SYNOPSIS
!!
!!   use get_ion, ONLY get_ionization
!!
!! Migrated from heatCool.F90 on 2019 June 07.
!! get_ionization(...) is only used by RadTrans_combine3.F90
!! and so is not in production at all; only provided so that code can compile.
!! Maybe useful for carrying out Josh's vision of a RadTrans cleanup and
!! merger in the future.  - AT 2019 June 07
!!
!!***

module get_ion

contains

function get_ionization(xHp_in, ndens_in, rho_in, temp_in, PhotoIon_in, dt_in)

use heatCool, ONLY : approx_column_dens
use Heat_data
use calc_ion

real :: get_ionization
real :: xh(2), phih, ion_dt
real, intent(in) :: xHp_in, ndens_in, rho_in, temp_in, PhotoIon_in, dt_in
logical :: fully_ionized

    ! You cannot photoionize more hydrogen than neutral hydrogen exists in the cell!
    fully_ionized = ((xHp_in .ge. 1.0d0) .or. ((1d0-xHp_in) .le. 0.0d0))
    ! Add background cosmic ray ionization rate to the photoionization rate.
    if (he_use_cr_heating .and. .not. fully_ionized) then
        if (he_crIonNH == 1.0 .and. he_crIonExp == 1.0) then
            phih = PhotoIon_in + he_crIonRate ! Uniform background ionization from CRs
        else ! Use eqn 27 from Padovani et. al. 2009
            phih = PhotoIon_in + he_crIonRate* &
                   (max(he_crIonNH,approx_column_dens(rho_in, ndens_in, temp_in)) &
                       / he_crIonNH)**(-he_crIonExp)
        end if
    end if

    !print*, phih, he_crIonRate, min(he_crIonNH,approx_column_dens(rho, ndens, temp))

    xh = (/1.0d0-xHp_in, xHp_in/)
    call calc_ionization(dt_in, ion_dt, temp_in, ndens_in, xh, xh, phih)
    get_ionization = xh(2)
return
end function get_ionization

end module get_ion
