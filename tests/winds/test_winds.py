"""
test_winds
==========

Compares the evolution of an expanding wind bubble to the analytical solution from
Weaver (1977). The tests asserts that the simulation reproduces the analytical
solution within XX%. The simulation should initialize a uniform background with
density at 10 cm**-3 and temperature of 100 K and evolve wind_driving_star.amuse
for 4 kyr. The star has an initial mass of 60 Msol and is placed at the center of
the box.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy.ndimage import gaussian_filter1d

from amuse.io import read_set_from_file
from amuse.units import units
from torch_se import compute_dmdt_vterm

yt.set_log_level(50)


class Weaver1977:
    """
    Analytical solutions for the structure of a wind-driven bubble.

    This class implements the self-similar solutions of Weaver et al. (1977)
    for the characteristic radii of a stellar wind bubble expanding into a
    uniform ambient medium.

    Structure of a wind-driven bubble
                  R1                      RC                          R2
    --------------|-----------------------|---------------------------|------------>
    * Star        |                       |                           |           r
                  |                       |                           |
    Stellar wind  | Shocked stellar wind  | Shocked interstellar gas  | Ambient interstellar gas
    <--------- --><----------------------><--------------------------><------------------------->
    R1 - Wind termination shock
    RC - Contact discontinuity
    R2 - Forward shock

    Reference: Weaver et al., Astrophysical Journal, Vol. 218, p. 377-395 (1977)

    Attributes:
        alpha: dimensionless expansion constant.
        xiC: ratio between the contact discontinuity radius to the outer shock radius. 
        f_therm: thermal fraction of inject wind energy.
    """
    alpha = (250.0/(308.0*np.pi))**(1.0/5.0)
    xiC = 0.86
    f_therm=5.0/11.0

    def wind_termination_shock_radius(self, time, dmdt, vel, rho0):
        """
        Return the radius of the wind termination shock.

        This is the radius where the freely expanding (hypersonic) stellar wind
        is shocked by the hot interior of the wind bubble, following the
        self-similar solution of Weaver et al. (1977).

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        Returns:
            Wind termination shock radius.
        """
        R1 = 0.9 * self.alpha**(3/2) * (dmdt/rho0)**(3/10) * (vel)**(1/10) * (time)**(2/5)
        return R1

    def forward_shock_radius(self, time, dmdt, vel, rho0):
        """
        Return the radius of the outer shock.

        This is the forward shock propagating into the ambient medium, following
        the self-similar solution of Weaver et al. (1977).

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        Returns:
            Outer shock radius.
        """
        Lw = self.wind_luminosity(dmdt, vel)
        R2 = self.alpha * (Lw * time**3 / rho0)**(1/5)
        return R2

    def contact_discontinuity_radius(self, time, dmdt, vel, rho0):
        """
        Return the radius of the contact discontinuity.

        The contact discontinuity separates the shocked stellar wind from the
        swept-up shocked interstellar medium.

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        Returns:
            Contact discontinuity radius.
        """
        RC = self.xiC * self.forward_shock_radius(time, dmdt, vel, rho0)
        return RC

    def wind_luminosity(self, dmdt, vel):
        """
        Mechanical wind luminosity
            Lw = 0.5 * dM/dt * v_w^2

        Arguments:
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
        """
        return 0.5 * dmdt * vel**2

    def wind_momentum_injection_rate(self, dmdt, vel):
        """
        Direct wind momentum injection rate
            p_w = dM/dt * v_w

        Arguments:
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
        """
        return dmdt * vel

    def shell_velocity(self, time, dmdt, vel, rho0):
        """
        Shock/shell velocity
            dR/dt = (3/5) * R_2/t

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        """
        R2 = self.forward_shock_radius(time, dmdt, vel, rho0)
        return (3.0 / 5.0) * R2 / time

    def swept_mass(self, time, dmdt, vel, rho0):
        """
        Swept-up ambient shell mass
            M_sh = (4pi/3) rho_0 * R_2**3

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        """
        R2 = self.forward_shock_radius(time, dmdt, vel, rho0)
        return (4.0 * np.pi / 3.0) * rho0 * R2**3


    def shell_momentum(self, time, dmdt, vel, rho0):
        """
        Momentum of swept-up shell momentum in the energy-conserving limit
            p_sh = M_sh * v_sh

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density. 
        """
        Msh = self.swept_mass(time, dmdt, vel, rho0)
        vsh = self.shell_velocity(time, dmdt, vel, rho0)
        return Msh * vsh


    def direct_wind_momentum(self, time, dmdt, vel):
        """
        Momentum injected directly
            p = dM/dt * v_w * t

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
        """
        return self.wind_momentum_injection_rate(dmdt, vel) * time


    def injected_wind_energy(self, time, dmdt, vel):
        """
        Total mechanical energy
            E_w = L_w * t

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
        """
        return self.wind_luminosity(dmdt, vel) * time


    def injected_thermal_energy(self, time, dmdt, vel):
        """
        Approximate hot-bubble thermal energy in the Weaver solution, commonly
            E_th ~ 5/11 * Lw * t

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
        """
        return self.f_therm * self.injected_wind_energy(time, dmdt, vel)


    def shell_kinetic_energy(self, time, dmdt, vel, rho0):
        """
        Kinetic energy of the swept-up shell
            Ek_w = 1/2 * M_sh * v_sh**2

        Arguments:
            time: times at which to calculate termination shock.
            dmdt: star mass-loss rate.
            vel: wind terminal velocity.
            rho0: initial ambient density.
        """
        Msh = self.swept_mass(time, dmdt, vel, rho0)
        vsh = self.shell_velocity(time, dmdt, vel, rho0)
        return 0.5 * Msh * vsh**2


class TestWinds:
    """
    Compare Torch wind-bubble simulations against Weaver et al. (1977).

    Attributes:
        directory: simulation data directory
        tolerance: asserts agreement within this tolerance
        weaver: analytical solution from Weaver et al. (1977)
        flash_files: plot files from flash
        star_files: amuse particle files
        plot: create plots
        saveto: directory where plot files are save if created
        initial_state: initial ambient medium and stellar wind properties
        wind_evolution: quanties that will be compared to Weaver et al. (1977) solution
    """

    def __init__(self, directory, number_checks = 5, tolerance=0.1, massloss_method='seba_puls', max_gamma=0.8, plot=False):
        """
        Initialize class by loading data and calculating the properties that will be tested.

        Arguments:
            directory: simulation data directory
            number_checks: number of outputs to check
            tolerance: asserts agreement within this tolerance
            massloss_method: mass loss method assumed by Torch
            max_gamma: maximum Eddington factor assumed by Torch
            plot: create plots, all=verbose, tests=show only tests
        """
        self.directory = Path(directory)
        self.number_checks = number_checks
        self.tolerance = tolerance
        self.weaver = Weaver1977()
        self.star_files = sorted(self.directory.glob("stars*.amuse"))[1:]
        self.flash_files = sorted(self.directory.glob("*_hdf5_plt_cnt_*"))
        self.plot = plot

        if self.plot:
            self.saveto = Path("./plots/")
            self.saveto.mkdir(exist_ok=True)

        if not self.flash_files:
            raise FileNotFoundError("No FLASH outputs found.")

        if len(self.star_files) < 2:
            raise FileNotFoundError("At least two stellar snapshots are required.")

        self.rho0 = self.get_rho0()
        self.center, self.dmdt, self.vterm = self.get_star_properties(massloss_method, max_gamma)
        self.wind_evolution = self.get_wind_evolution()

    def get_rho0(self):
        """
        Load initial density from first flash output.

        Returns:
            Median density in first output
        """
        ds = yt.load(str(self.flash_files[0]))
        ad = ds.all_data()

        rho0 = np.median(ad[("gas", "density")].to("g/cm**3").value) | (units.g / units.cm**3)
        return rho0

    def get_star_properties(self, massloss_method, max_gamma):
        """
        Get mean stellar wind properties.

        Arguments:
            massloss_method: massloss_method used in torch_user.py
            max_gamma: max_gamma used in torch_user.py
        Returns
            center, dmdt, vterm
        """
        dmdt_list = []
        vterm_list = []
        age_list = []

        star_prev = read_set_from_file(str(self.star_files[0]), "amuse")[0]
        center = star_prev.position
        age_list.append(star_prev.age.value_in(units.Myr))

        for starfile in self.star_files[1:]:
            star_now = read_set_from_file(starfile)[0]
            dt = star_now.age - star_prev.age

            if dt <= (0 | units.s):
                continue

            dmdt, vterm = compute_dmdt_vterm(
                star_prev.mass,
                star_now.temperature,
                star_now.radius,
                star_now.mass,
                star_now.luminosity,
                dt,
                massloss_method=massloss_method,
                max_gamma=max_gamma)

            dmdt_list.append(dmdt.value_in(units.MSun / units.yr))
            vterm_list.append(vterm.value_in(units.km / units.s))
            age_list.append(star_now.age.value_in(units.Myr))
            star_prev = star_now

        if len(dmdt_list) == 0:
            raise ValueError("No valid stellar evolution intervals.")

        dmdt = np.average(dmdt_list, weights=np.diff(age_list)) | units.MSun / units.yr
        vterm = np.average(vterm_list, weights=np.diff(age_list)) | units.km / units.s

        if self.plot:
            fig=plt.figure()
            plt.ylabel(r"Stellar mass loss [${\rm M}_{\odot}\,{\rm yr}^{-1}$]")
            plt.xlabel(r"Stellar age [Myr]")
            plt.plot(age_list[1:], dmdt_list, ".k")
            plt.axhline(dmdt.value_in(units.MSun / units.yr), color='r')
            plt.savefig(self.saveto / 'dmdt.png', dpi=120, bbox_inches="tight")
            plt.close(fig)

            fig=plt.figure()
            plt.ylabel(r"Wind velocity [${\rm km\, s}^{-1}$]")
            plt.xlabel(r"Stellar age [Myr]")
            plt.plot(age_list[1:], vterm_list, '.k')
            plt.axhline(vterm.value_in(units.km / units.s), color='r')
            plt.savefig(self.saveto / 'vterm.png', dpi=120, bbox_inches="tight")
            plt.close(fig)

        return center, dmdt, vterm

    def get_wind_evolution(self):
        """
        Process all FLASH snapshots and calculate test quantities.

        Returns:
            Dictionary including time, momentum, R1, and R2
        """

        def _get_profile(r, quantity, weights=None, nbins=256, log = True, smooth=2):
            """
            Get radial profile.

            Arguments:
                r: radius
                quantity: quantity evalutated at r
                wegiths: weights applied to quantity
                nbins: number of radial bins to use
                log: return logarithm of profile
                smooth: number of bins used to smooth profile

            Returns:
                center, profile
            """
            edges = np.linspace(0, r.max(), nbins + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])

            if weights is None:
                weights = np.ones(quantity.size)
                w_sum = np.hisogram(r, bins=edges)
            else:
                w_sum, _ = np.histogram(r, bins=edges, weights=weights)

            hist, _ = np.histogram(r, bins=edges, weights=quantity*weights)
            profile = np.divide(hist, w_sum, out=np.full(nbins, np.nan), where=w_sum > 0)

            valid = w_sum > 0
            profile = np.interp(centers, centers[valid], profile[valid])

            if log:
                profile = np.log10(profile)

            if smooth > 0:
                profile = gaussian_filter1d(profile, smooth)

            return centers, profile

        def _bubble_radii_index(rbins, profile, output_nr=None):
            """
            Measure radius of wind termination shock (R1) and forward shock (R2).

            Arguments:
                rbins: radial bins
                profile: gas profile
                output_nr: output number for error messages

            Returns:
                iR1, iR2
            """

            dp = np.gradient(profile, rbins)

            # Termination shock: strongest positive pressure gradient
            try:
                iR1 = np.nanargmax(dp)
            except ValueError: 
                raise ValueError(f"All-NaN slice encountered in output {output_nr:05d}.")

            # Forward shock: strongest negative pressure gradient outside R1
            try:
                iR2 = iR1 + np.nanargmin(dp[iR1 + 1:]) + 1
            except ValueError: 
                raise ValueError(f"All-NaN slice outside R1 encountered in output {output_nr:05d}.")

            return iR1, iR2


        time_arr = np.zeros(self.number_checks)
        momentum_arr = np.zeros(self.number_checks)
        R1_arr = np.zeros(self.number_checks)
        R2_arr = np.zeros(self.number_checks)

        shift = int((len(self.flash_files))/self.number_checks)
        for i in range(self.number_checks):
            filename = self.flash_files[(1+i)*shift-1]
            ds = yt.load(filename)
            ad = ds.all_data()

            time_arr[i] = ds.current_time.to("s").value

            # Hydro quantites
            position = np.column_stack([ad[("gas", f"{axis}")].to("cm").value for axis in "xyz"])
            velocity = np.column_stack([ad[("gas", f"velocity_{axis}")].to("cm/s").value for axis in "xyz"])
            density = ad[("gas", "density")].to("g/cm**3").value
            pressure = ad["gas", "pressure"].to("dyne/cm**2").value
            temperature = ad["gas", "temperature"].to("K").value
            volume = ad[("gas", "cell_volume")].to("cm**3").value

            # Momentum
            radius = np.sqrt(np.sum((position - self.center.value_in(units.cm))**2, axis=1))
            rhat = (position - self.center.value_in(units.cm)) / radius.reshape((-1, 1))
            vr = np.sum(velocity * rhat, axis=1)
            mass = density * volume
            momentum_arr[i] = np.sum(mass * vr)

            # Characteristic radii from pressure profile
            rbins, pres = _get_profile(radius, pressure, weights = volume)
            iR1, iR2 = _bubble_radii_index(rbins, pres, output_nr=i+1)
            R1_arr[i] = rbins[iR1]
            R2_arr[i] = rbins[iR2]

            if self.plot:
                fig, axs = plt.subplots(ncols=1, nrows=5, sharex=True, figsize=(4,6))
                plt.subplots_adjust(hspace=0.02, wspace=0.02)
                for ax in axs:
                    ax.tick_params(axis='both', labelsize=7)
                axs[-1].set_xlabel("Radius [pc]", fontsize=8)
                axs[0].text(0.95, 0.3, f"Time = {ds.current_time.to("kyr").value:.0f} kyr", transform=axs[0].transAxes, ha='right', va='center', fontsize=8)
                axs[0].text(0.95, 0.5, r"Density [${g\ \rm cm}^{-3}$]", transform=axs[0].transAxes, ha='right', va='center', fontsize=8)
                axs[1].text(0.95, 0.5, r"Temperature [K]", transform=axs[1].transAxes, ha='right', va='center', fontsize=8)
                axs[2].text(0.95, 0.5, r"Pressure [${\rm dyne}\, {\rm cm}^{-2}$]", transform=axs[2].transAxes, ha='right', va='center', fontsize=8)
                axs[3].text(0.95, 0.5, r"Velocity [km/s]", transform=axs[3].transAxes, ha='right', va='center', fontsize=8)
                axs[4].text(0.95, 0.5, "Gradients", transform=axs[4].transAxes, ha='right', va='center', fontsize=8)

                _, rho = _get_profile(radius, density, weights = volume)
                _, temp = _get_profile(radius, temperature, weights = volume)
                _, vel = _get_profile(radius, vr, weights = volume)
                axs[0].plot(rbins, rho, '-r')
                axs[1].plot(rbins, temp, '-g')
                axs[2].plot(rbins, pres, '-b')
                axs[3].plot(rbins, vel, '-k')

                axs[4].plot(rbins, np.gradient(temp, rbins), '-r', lw=1)
                axs[4].plot(rbins, np.gradient(rho, rbins), '-g', lw=1)
                axs[4].plot(rbins, np.gradient(pres, rbins), '-b', lw=1)
                axs[4].plot(rbins, np.gradient(vel, rbins), '-k', lw=1)
                for ax in axs:
                    [ax.axvline(rad, linestyle=ls, color='k', alpha=0.6) for rad,ls in zip([rbins[iR1], rbins[iR2]],[':', '--'])]
                plt.savefig(self.saveto / f'profile_{i+shift+1:05d}.png', dpi=120, bbox_inches="tight")
                plt.close(fig)

        return {"time": time_arr | units.s,
                "momentum": momentum_arr | units.g * units.cm / units.s,
                "R1": R1_arr | units.cm,
                "R2": R2_arr | units.cm}

    def _expected(self):
        """
        Calculate analytical predictions at all snapshot times.

        Returns:
            Dictionary containing R1, RC, R2, and momentum of Weaver et al. (1977) solution
        """
        t = self.wind_evolution["time"]

        return {"R1": self.weaver.wind_termination_shock_radius(t, self.dmdt, self.vterm, self.rho0),
                "RC": self.weaver.contact_discontinuity_radius(t, self.dmdt, self.vterm, self.rho0),
                "R2": self.weaver.forward_shock_radius(t, self.dmdt, self.vterm, self.rho0),
                "momentum": self.weaver.shell_momentum(t, self.dmdt, self.vterm, self.rho0)}

    def _assert_agreement(self, measured, expected, name):
        """
        Assert agreement at all times within relative tolerance.

        Arguments:
            measured: array of measured quantities to check
            expected: array of expected value of quantity
            name: name of quantities
        """
        measured_values = measured.value_in(expected.unit)
        expected_values = expected.value_in(expected.unit)
        relative_error = np.abs(measured_values / expected_values - 1)

        print(f"{name}: maximum relative error "
              f"{np.max(relative_error):.1%} "
              f"with a tolerance of {self.tolerance:.1%}")

    def test_time_evolution(self):
        """
        Compare characteristic radii against Weaver (1977).
        """
        expected = self._expected()

        if self.plot:
            fig = plt.figure()
            plt.gca().tick_params(axis='both', labelsize=7)
            plt.xlabel('Time [kyr]', fontsize=10)
            plt.ylabel('Radius [pc]', fontsize=10)
            plt.yscale('log')
            plt.xscale('log')
            plt.plot([],[], '-k', label = 'Simulation')
            plt.plot([],[], '--k', label = 'Weaver et al. (1977)')
            for quantity, c in zip(["R1", "R2"],["r", "b"]):
                plt.plot(self.wind_evolution["time"].value_in(units.kyr), self.wind_evolution[quantity].value_in(units.pc), '-', label=quantity, color=c)
                plt.plot(self.wind_evolution["time"].value_in(units.kyr), expected[quantity].value_in(units.pc), '--', color=c)

            plt.legend(loc='upper left', fontsize=8, frameon=False)
            plt.savefig(self.saveto / "bubble_radii.png", dpi=120, bbox_inches="tight")
            plt.close(fig)

        for quantity in ("R1", "R2"):
            self._assert_agreement(self.wind_evolution[quantity], expected[quantity], quantity)

    def test_momentum_evolution(self):
        """
        Compare radial shell momentum against Weaver (1977).
        """
        expected = self._expected()

        if self.plot:
            fig = plt.figure()
            plt.gca().tick_params(axis='both', labelsize=7)
            plt.xlabel('Time [kyr]', fontsize=10)
            plt.ylabel(r'Momentum [${\rm M}_{\odot}\,{\rm km\, s}^{-1}$]', fontsize=10)
            plt.yscale('log')
            plt.plot(self.wind_evolution["time"].value_in(units.kyr), self.wind_evolution["momentum"].value_in(units.MSun * units.km / units.s), '-k', label="Simulation")
            plt.plot(self.wind_evolution["time"].value_in(units.kyr), expected["momentum"].value_in(units.MSun * units.km / units.s), '--k', label="Weaver et al. (1977)")

            plt.legend(loc='upper left', fontsize=8, frameon=False)
            plt.savefig(self.saveto / "momentum.png", dpi=120, bbox_inches="tight")
            plt.close(fig)

        self._assert_agreement(self.wind_evolution["momentum"], expected["momentum"], "Shell momentum")

def main():

    parser = argparse.ArgumentParser(description="Test Torch wind bubble against Weaver (1977).")

    parser.add_argument("directory", 
        type=Path, 
        help="Simulation directory")

    parser.add_argument("--number_checks",
        type=float,
        default=10,
        help="Number of outputs to check",)
    
    parser.add_argument("--tolerance",
        type=float,
        default=0.1,
        help="Relative tolerance (default: 0.1)",)

    parser.add_argument("--plot",
        action="store_true",
        help="Plot tests (tests) or plot everything (all)")

    args = parser.parse_args()
    test = TestWinds(args.directory, args.number_checks, args.tolerance, plot=args.plot)
    test.test_momentum_evolution()
    test.test_time_evolution()

if __name__ == "__main__":
    main()

