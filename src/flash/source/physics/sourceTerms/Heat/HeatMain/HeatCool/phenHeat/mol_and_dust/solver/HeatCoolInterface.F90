module heatInterface

! Note interfaces don't import on their own, must use import command.
use OdeData, ONLY : dp

implicit none

interface
    subroutine heating(temp_in, ndens_in, rho_in, ephen_in, zpos_in, eplus)
	!use Heat_data ! Use Christian's heating data file.
	import :: dp

	real(dp), intent(in)    :: temp_in, ndens_in, rho_in, ephen_in, zpos_in
	real(dp), intent(inout) :: eplus
    end subroutine heating
end interface

interface
    subroutine pe_heating(Gflux_in, xHp_in, ndens_in, temp_in, rho_in, PeHeat, dust_heating)

	import :: dp

	real(dp), intent(in)  :: Gflux_in, xHp_in, ndens_in, temp_in, rho_in
	real(dp), intent(out) :: PeHeat, dust_heating
    end subroutine pe_heating
end interface

interface
    subroutine cooling(temp_in, ndens_in, conf_in, emin)

	import :: dp

	!use Heat_data ! Use Christian's heating data file.
	!use Cool_data ! Use Robi's cooling data file.

	real(dp), intent(in)  :: temp_in, ndens_in, conf_in
	real(dp), intent(inout) :: emin
    end subroutine cooling
end interface

interface
    subroutine molecular_cooling(ndens_in, temp_in, mu_mol_in, cool_pw)

	import :: dp

	real(dp), intent(in)  :: ndens_in, temp_in, mu_mol_in
	real(dp), intent(out) :: cool_pw
    end subroutine molecular_cooling
end interface

interface
    subroutine dust_cooling(ndens_in, temp_in, tdust_in, rho_in, cool_dg)

	import :: dp

	real, intent(in)  :: ndens_in, temp_in, tdust_in, rho_in
	real, intent(out) :: cool_dg
    end subroutine dust_cooling
end interface

interface
    subroutine Radloss(T, radia)
	
	import  :: dp

	real(dp), intent(IN) :: T
	real(dp), intent(OUT) :: radia
    end subroutine Radloss
end interface

end module heatInterface
