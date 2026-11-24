"""
test_sne
========

Compares the evolution of a blastwave to the analytical Sedov-Taylor solution and
fitting functions from Kim & Ostriker (2015). The tests asserts that the simulation 
reproduces the analytical solution within XX%. The simulation should initialize a 
uniform background with density at 10 cm**-3 and temperature of 100 K and evolve the 
star in sne_test.amuse for 50 kyr. The star has an initial mass of XX Msol, exploded
within a few timestaps and is placed at the center of the box.
"""
import argparse
from pathlib import Path
import re

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import quad
import yt

from amuse.units import constants, units

yt.set_log_level(50)

class KimOstriker2015:
    """Spherical, energy-conserving blastwave in a uniform ambient medium.

    Reference: Kim & Ostriker, ApJ, Vol. 802, Issue 2, article id. 99, 25 pp. (2015).

    Attributes:
        xi0: dimensionless expansion constant.
        gamma: adiabatic index
        mean_particle_mass: mean mass per particle in fully ionized gas.
        f_therm: thermal fraction of the explosion energy.
    """

    xi0 = 1.15167
    gamma = 5.0 / 3.0
    mean_particle_mass = (1.4 / 2.3) * 1.6735575e-24 | units.g
    f_therm = 0.717

    def free_expansion_time(self, energy, ejecta_mass, rho0):
        """
        End of ballistic expansion: 
            t_free = r_free / v_ej                                               [Eq. 1]

        Arguments:
            energy: explosion kinetic energy
            ejecta_mass: mass initially ejected by the explosion
            rho0: uniform ambient mass density

        Returns:
            Time when swept-up mass equals ejecta mass (end of free expansion)
        """
        velocity = (2.0 * energy / ejecta_mass)**0.5
        return self.free_expansion_radius(ejecta_mass, rho0) / velocity

    def free_expansion_radius(self, ejecta_mass, rho0):
        """
        Radius where swept-up mass equals ejecta mass:
            r_free = (3 M_ej / (4 pi rho0))**(1/3)                               [Eq. 2]

        Arguments:
            ejecta_mass: mass initially ejected by the explosion
            rho0: uniform ambient mass density

        Returns:
            Radius at end of free expansion
        """
        return (3.0 * ejecta_mass / (4.0 * np.pi * rho0))**(1.0 / 3.0)

    def free_expansion_momentum(self, energy, ejecta_mass):
        """
        Ballistic ejecta radial momentum
            p_free = M_ej v_ej                                                  [Eq. 15]

        Arguments:
            energy: explosion kinetic energy.
            ejecta_mass: mass initially ejected by the explosion.

        Returns:
            Radial momentum at end of free expansio.
        """
        return (2.0 * ejecta_mass * energy)**0.5

    def sedovtaylor_radius(self, time, energy, rho0):
        """
        Shock radius of Sedov-Taylor blastwave 
            r_ST = alpha (E_SN t**2 / rho0)**(1/5)                          [c.f. Eq. 3]

        Arguments:
            time: elapsed time since the explosion
            energy: explosion energy
            rho0: uniform ambient mass density

        Returns:
            Sedov-Taylor radius, valid between free expansion and before cooling.
        """
        return self.xi0 * (energy * time**2 / rho0)**(1.0 / 5.0)

    def sedovtaylor_velocity(self, time, energy, rho0):
        """ 
        Shock velocity of Sedov-Taylor blastwave
            v_ST = (2/5) r_ST / t                                                [Eq. 4]

        Arguments:
            time: elapsed time since the explosion
            energy: explosion energy
            rho0: uniform ambient mass density

        Returns:
            Sedov-Taylor speed, valid between free expansion and before cooling.
        """
        return (2.0 / 5.0) * self.sedovtaylor_radius(time, energy, rho0) / time

    def sedovtaylor_temperature(self, time, energy, rho0):
        """
        Immediate postshock temperature 
            T_ST = 3 mu v_ST**2 / (16 kB)                                        [Eq. 5]

        Arguments:
            time: elapsed time since the explosion
            energy: explosion energy
            rho0: uniform ambient mass density

        Returns:
            Postshock temperature
        """
        velocity = self.sedovtaylor_velocity(time, energy, rho0)
        return (3.0 / 16.0) * self.mean_particle_mass * velocity**2 / constants.kB

    def sedovtaylor_momentum(self, time, energy, rho0):
        """
        Integrated outward momentum of Sedov-Taylor blastwave 
            p_ST = 2.69 rho0 v_ST r_ST**3                                       [Eq. 16]

        Arguments:
            time: elapsed time since the explosion
            energy: explosion energy
            rho0: uniform ambient mass density

        Returns:
            Sedov-Taylor momentum, valid between free expansion and before cooling.
        """
        velocity = self.sedovtaylor_velocity(time, energy, rho0)
        radius = self.sedovtaylor_radius(time, energy, rho0)
        return 2.69 * rho0 * velocity * radius**3

    def cooling_time(self, time, energy, rho0):
        """Cooling time of freshly shocked gas, not the hot interior 
            t_cool = (2.3/1.2) kB T_ST / ((gamma+1) * n0 * Lambda(T_ST))         [Eq. 6]

        Arguments:
            time: time when gas was shocked
            energy: explosion energy
            rho0: uniform ambient mass density

        Returns:
            Cooling time for shocked gas.
        """

        def _cooling_function(temperature):
            """ Approximate power-law cooling function for the range 1e5 - 1e7.5 K
            """
            _C = 1.1e-22 | units.erg * units.cm**3 / units.s
            _alpha = 0.7
            return _C * (temperature.value_in(units.K) / 1e6)**-_alpha

        temperature = self.sedovtaylor_temperature(time, energy, rho0)
        n0 = rho0 / (1.4 * (1.6735575e-24 | units.g))
        return (2.3 / 1.2) * (constants.kB * temperature) \
            / ((self.gamma + 1) * n0 * _cooling_function(temperature))

    def shell_time(self, energy, rho0):
        """
        Earliest parcel cooling time obtained by finding t_s that minimized t_sf
            t_sf = min[t_s + t_cool(t_s)]                                   [c.f. Eq. 7]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Shell-formation time derived from the adopted cooling law.
        """

        def _shell_time(t):
            t_s = t | units.yr
            return (t_s + self.cooling_time(t_s, energy, rho0)).value_in(units.yr)

        result = minimize_scalar(_shell_time, bounds=(1e2, 1e6), method="bounded")

        if not result.success:
            raise RuntimeError("Failed to determine shell-formation time.")

        return result.fun | units.yr

    def shell_radius(self, energy, rho0):
        """
        Sedov-Taylor radius at shell formation
            r_sh = r_ST(t_sf)                                               [c.f. Eq. 8]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Sedov-Taylor radius evaluated at shell_time(energy, rho0).
        """
        return self.sedovtaylor_radius(self.shell_time(energy, rho0), energy, rho0)

    def shell_velocity(self, energy, rho0):
        """
        Shock speed at shell formation
            v_sh = v_ST(t_sf)                                               [c.f. Eq. 9]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Sedov-Taylor shock speed evaluated at shell formation.
        """
        return self.sedovtaylor_velocity(self.shell_time(energy, rho0), energy, rho0)

    def shell_temperature(self, energy, rho0):
        """
        Postshock temperature at shell formation 
            T_sh = T_ST(t_sf)                                              [c.f. Eq. 10]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Immediate postshock temperature, before that parcel cools.
        """
        return self.sedovtaylor_temperature(self.shell_time(energy, rho0), energy, rho0)

    def shell_mass(self, energy, rho0):
        """
        Swept-up mass at shell formation
            m_sf = (4 pi/3) rho0 r_ST(t_sf)**3                             [c.f. Eq. 11]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Total swept ambient mass, excluding ejecta.
        """
        return (4.0 * np.pi / 3.0) * rho0 * self.shell_radius(energy, rho0)**3

    def shell_momentum(self, energy, rho0):
        """
        Integrated radial momentum at shell formation
            p_sf = p_ST(t_sf)                                              [c.f. Eq. 17]

        Arguments:
            energy: explosion energy.
            rho0: uniform ambient mass density.

        Returns:
            Sedov-Taylor momentum at shell_time.
        """
        return self.sedovtaylor_momentum(self.shell_time(energy, rho0), energy, rho0)

    def pds_radius(self, time, energy, rho0):
        """
        Shock radius during the pressure-driven snowplow phase. Assumes adiabatic 
        expansion of hot interior, see e.g., Cioffi et al. (1988).
        """
        t_sf = self.shell_time(energy, rho0)
        r_sf = self.shell_radius(energy, rho0)

        return r_sf * (time / t_sf)**(2.0 / 7.0)


    def pds_momentum(self, time, energy, rho0):
        """
        Radial momentum during the pressure-driven snowplow phase. Assumes adiabatic 
        expansion of hot interior, implementing Eq. [18] of Kim & Ostriker (2015).
        """
        t_sf = self.shell_time(energy, rho0)
        p_sf = self.shell_momentum(energy, rho0)

        return p_sf * (1.0 + 4.62 * ((time/t_sf)**(1.0/7.0) - 1.0))

    def snr_momentum(self, time, energy, rho0):
        """
        Radial momentum of the supernovae remnant using the effective pressure of the
        interior, implementing Eq. [28] of Kim & Ostriker (2015).
        """
        t_sf = self.shell_time(energy, rho0)
        p_sf = self.shell_momentum(energy, rho0)

        return p_sf * (1.0 + 0.62 * (1.0 - (time/t_sf)**(-6.0/7.0)))


    def _check_implementation(self):
        """ Checking implementation against values in paper.
        """
        energy = 1e51 | units.erg
        rho = 1.4 * 1.6735575e-24 | units.g / units.cm**3
        ejecta_mass = 3 | units.Msun
        t3 = 1000 | units.yr

        print("Free expansion phase")
        print("Time:", self.free_expansion_time(energy=energy, ejecta_mass=ejecta_mass, rho0=rho).value_in(units.yr))
        print("Radius:", self.free_expansion_radius(ejecta_mass=ejecta_mass, rho0=rho).value_in(units.pc))
        print("Momentum:", self.free_expansion_momentum(energy=energy, ejecta_mass=ejecta_mass).value_in(units.Msun * units.km / units.s))
        print("")

        print("Sedov-Taylor expansion")
        print("Radius:", self.sedovtaylor_radius(time=t3, energy=energy, rho0=rho).value_in(units.pc))
        print("Velocity:", self.sedovtaylor_velocity(time=t3, energy=energy, rho0=rho).value_in(units.km / units.s))
        print("Tempertaure:", self.sedovtaylor_temperature(time=t3, energy=energy, rho0=rho).value_in(units.K))
        print("Momentum:", self.sedovtaylor_momentum(time=t3, energy=energy, rho0=rho).value_in(units.Msun * units.km / units.s))
        print("")

        print("Shell formation")
        print("Time:", self.shell_time(energy=energy, rho0=rho).value_in(units.yr))
        print("Radius:", self.shell_radius(energy=energy, rho0=rho).value_in(units.pc))
        print("Velocity:", self.shell_velocity(energy=energy, rho0=rho).value_in(units.km / units.s))
        print("Tempertaure:", self.shell_temperature(energy=energy, rho0=rho).value_in(units.K))
        print("Mass:", self.shell_mass(energy=energy, rho0=rho).value_in(units.Msun))
        print("Momentum:", self.shell_momentum(energy=energy, rho0=rho).value_in(units.Msun * units.km / units.s))
        print("")

class TorchBlastwave:
    """
    Compare Torch blastwave simulations against Kim & Ostriker (2015).

    Attributes:
        directory: simulation data directory
        tolerance: asserts agreement within this tolerance
        number_checks: number of flash outputs to check
        start_checks_time: first output to check
        ko15: analytical solution from Kim & Ostriker (2015)
        flash_files: plot files from flash
        plot: create plots
        saveto: directory where plot files are save if created
        initial_state: initial ambient medium and blastwave properties
        sne_evolution: quanties that will be compared to Kim & Ostriker (2015) solution
    """

    def __init__(self, directory,  tolerance, number_checks = 5, start_checks_time = 1 | units.kyr, sne_energy=1e51 | units.erg, plot=False):
        """
        Initialize class by loading data and calculating the properties that will be tested.

        Arguments:
            directory: simulation data directory
            tolerance: asserts agreement within this tolerance
            number_checks: number of outputs to check
            start_checks_time: first output to check
            sne_energy: energy injection assumed by Torch
            plot: create plots, all=verbose, tests=show only tests
        """
        self.directory = Path(directory)
        self.datadir = self.directory / "data"
        self.tolerance = tolerance
        self.number_checks = number_checks
        self.start_checks_time = start_checks_time
        self.ko15 = KimOstriker2015()
        self.flash_files = sorted(self.datadir.glob("*_hdf5_plt_cnt_*"))
        self.plot = plot

        if self.plot:
            self.saveto = Path("./plots/")
            self.saveto.mkdir(exist_ok=True)

        if not self.flash_files:
            raise FileNotFoundError("No FLASH outputs found.")

        self.rho0 = self.get_rho0()
        self.sne_energy = sne_energy
        self.t0, self.center, self.injected_mass = self.get_sne_properties()

        self.flash_files = [f for f in self.flash_files
            if self.get_time_flash(f) >= self.start_checks_time]

        self.sne_evolution = self.get_sne_evolution()

        if self.plot:
            self.plot_final_output()

    def get_time_flash(self, flash_file):
        """
        Get the timestamp from a flash output.

        Arguments:
            flash_file: flash output file to load
        """

        with h5py.File(flash_file, "r") as file:
            scalars = file["real scalars"][:]

        for args, val in zip(scalars["name"], scalars["value"]):
            if b'time' in args:
                time = val | units.s
                return time

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

    def get_sne_properties(self):
        """ 
        Load properties of supernova injection

        Returns
            center: location of explosion
            injected_mass: mass injected in explosion
        """

        logfile = next(self.directory.glob("slurm-*.out"))
        bridge_pattern = r"Bridge step: it=\d+, t=([\d.eE+-]+), dt=[\d.eE+-]+"
        sne_pattern = r"SN x=([\d.eE+-]+) cm, y=([\d.eE+-]+) cm, z=([\d.eE+-]+) cm, inj_mass=([\d.eE+-]+)"

        with logfile.open() as f:
            for line in f:
                match = re.search(bridge_pattern, line)
                if match:
                    bridge_time = float(match.group(1))

                match = re.search(sne_pattern, line)
                if match:
                    x, y, z, injected_mass = map(float, match.groups())
                    return bridge_time | units.s, (x, y, z) | units.cm, injected_mass | units.Msun

        raise ValueError(f"No SN injection found in {logfile}")

    def plot_final_output(self, 
            xlim = [-12.8, 12.8],
            ylim = [-12.8, 12.8],
            extent = [-12.8, 12.8, -12.8, 12.8],
            vmin = -1.5,
            vmax = 2.5):
        """
        Creates slice plot of final output.
        """

        def _get_frame(ds, quantity, unit, width=(1.0, "unitary"), res=[256, 256]):
            """ Helper function to create slice plot. """
            sl = ds.slice('z', 0.5)
            frb = sl.to_frb(width, res)
            if unit == 'None':
                return np.array(frb[quantity])
            else:
                return np.array(frb[quantity].to(unit))

        fig = plt.figure(figsize=(4.2, 3.6))
        plt.xlabel('x [pc]', fontsize=11)
        plt.ylabel('y [pc]', fontsize=11)
        plt.xlim(*xlim)
        plt.ylim(*ylim)
        plt.gca().set_aspect('equal')
        plt.gca().label_outer()

        ds = yt.load(self.flash_files[-1])
        dens = _get_frame(ds, 'density', 'mp/cm**3')
        im = plt.imshow(np.log10(dens), cmap='viridis', extent=extent, vmin=vmin, vmax=vmax, origin='lower')

        plt.text(0.95, 0.95, f"{ds.current_time.in_units('kyr').value:.2f} kyr", 
                 transform=plt.gca().transAxes, 
                 va='top', ha='right', 
                 fontsize=9)
        plt.text(0.95, 0.05, r"$\Delta x=$" + f"{ds.index.get_smallest_dx().to("pc").value:.2f} pc", 
                 transform=plt.gca().transAxes, 
                 va='bottom', ha='right', 
                 fontsize=9)
        cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.77])
        fig.colorbar(im, cax=cbar_ax, label=r"$\log_(n\ [{\rm cm}^{-3}])$")
        plt.savefig(self.saveto / f'slice.png', dpi=240, bbox_inches="tight")
        plt.close(fig)

    def get_sne_evolution(self):
        """
        Process all FLASH snapshots and calculate test quantities.

        Returns:
            Dictionary including time, momentum, and buuble radius
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

            w_sum, _ = np.histogram(r, bins=edges, weights=weights)

            hist, _ = np.histogram(r, bins=edges, weights=quantity*weights)
            profile = np.divide(hist, w_sum, out=np.full(nbins, np.nan), where=w_sum > 0)

            valid = w_sum > 0
            profile = np.interp(centers, centers[valid], profile[valid])

            if log:
                valid = profile > 0
                profile[valid] = np.log10(profile[valid])
                profile[~valid] = np.nan

            if smooth > 0:
                profile = gaussian_filter1d(profile, smooth)

            return centers, profile

        def _bubble_radii_index(rbins, profile, output_nr=None):
            """
            Measure radial index of forward shock (R).

            Arguments:
                rbins: radial bins
                profile: gas profile
                output_nr: output number for error messages

            Returns:
                iR
            """

            dp = np.gradient(profile, rbins)

            # Forward shock: strongest negative pressure gradient
            try:
                iR = np.nanargmin(dp)
            except ValueError: 
                raise ValueError(f"All-NaN slice encountered in output {output_nr:05d}.")

            return iR

        time_arr = np.zeros(self.number_checks)
        momentum_arr = np.zeros(self.number_checks)
        R_arr = np.zeros(self.number_checks)

        shift = int((len(self.flash_files))/self.number_checks)
        for i in range(self.number_checks):
            ifile = (1+i)*shift-1
            filename = self.flash_files[ifile]
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
            iR = _bubble_radii_index(rbins, pres, output_nr=ifile)
            R_arr[i] = rbins[iR]

            if self.plot:
                fig, axs = plt.subplots(ncols=1, nrows=5, sharex=True, figsize=(4,6))
                plt.subplots_adjust(hspace=0.02, wspace=0.02)
                for ax in axs:
                    ax.tick_params(axis='both', labelsize=7)
                axs[-1].set_xlabel("Radius [pc]", fontsize=8)
                axs[0].text(0.95, 0.3, f"Time = {ds.current_time.to("kyr").value-self.t0.value_in(units.kyr):.1f} kyr", transform=axs[0].transAxes, ha='right', va='center', fontsize=8)
                axs[0].text(0.95, 0.5, r"Density [${g\ \rm cm}^{-3}$]", transform=axs[0].transAxes, ha='right', va='center', fontsize=8)
                axs[1].text(0.95, 0.5, r"Temperature [K]", transform=axs[1].transAxes, ha='right', va='center', fontsize=8)
                axs[2].text(0.95, 0.5, r"Pressure [${\rm dyne}\, {\rm cm}^{-2}$]", transform=axs[2].transAxes, ha='right', va='center', fontsize=8)
                axs[3].text(0.95, 0.5, r"Velocity [km/s]", transform=axs[3].transAxes, ha='right', va='center', fontsize=8)
                axs[4].text(0.95, 0.5, "Gradients", transform=axs[4].transAxes, ha='right', va='center', fontsize=8)

                _, rho = _get_profile(radius, density, weights = volume)
                _, temp = _get_profile(radius, temperature, weights = volume)
                _, vel = _get_profile(radius, vr, weights = volume)
                rbins = rbins | units.cm
                axs[0].plot(rbins.value_in(units.pc), rho, '-r')
                axs[1].plot(rbins.value_in(units.pc), temp, '-g')
                axs[2].plot(rbins.value_in(units.pc), pres, '-b')
                axs[3].plot(rbins.value_in(units.pc), vel, '-k')

                axs[4].plot(rbins.value_in(units.pc), np.gradient(temp, rbins.value_in(units.cm)), '-r', lw=1)
                axs[4].plot(rbins.value_in(units.pc), np.gradient(rho, rbins.value_in(units.cm)), '-g', lw=1)
                axs[4].plot(rbins.value_in(units.pc), np.gradient(pres, rbins.value_in(units.cm)), '-b', lw=1)
                axs[4].plot(rbins.value_in(units.pc), np.gradient(vel, rbins.value_in(units.cm)), '-k', lw=1)
                for ax in axs:
                    ax.axvline(rbins[iR].value_in(units.pc), linestyle='--', color='k', alpha=0.6)
                plt.savefig(self.saveto / f'profile_{ifile:05d}.png', dpi=120, bbox_inches="tight")
                plt.close(fig)

        return {"time": time_arr - self.t0.value_in(units.s) | units.s,
                "momentum": momentum_arr | units.g * units.cm / units.s,
                "R": R_arr | units.cm}

    def _expected(self):
        """
        Calculate analytical predictions at all snapshot times.

        Returns:
            Dictionary containing radius and momentum of Kim & Ostriker (2015) solution
        """
        t = self.sne_evolution["time"]
        R = np.zeros(len(t)) | units.pc
        momentum = np.zeros(len(t)) | (units.Msun * units.km / units.s)

        free_expansion = (t < self.ko15.free_expansion_time(
            energy=self.sne_energy, 
            ejecta_mass=self.injected_mass,
            rho0=self.rho0))
        snowplow = (t > self.ko15.shell_time(
            energy=self.sne_energy,
            rho0=self.rho0))
        sedov_taylor = (~free_expansion) & (~snowplow)

        R[free_expansion] = self.ko15.free_expansion_radius(ejecta_mass=self.injected_mass, rho0=self.rho0)
        R[sedov_taylor] = self.ko15.sedovtaylor_radius(time=t[sedov_taylor], energy=self.sne_energy, rho0=self.rho0)
        R[snowplow] = self.ko15.pds_radius(time=t[snowplow], energy=self.sne_energy, rho0=self.rho0)

        momentum[free_expansion] = self.ko15.free_expansion_momentum(energy=self.sne_energy, ejecta_mass=self.injected_mass)
        momentum[sedov_taylor] = self.ko15.sedovtaylor_momentum(time=t[sedov_taylor], energy=self.sne_energy, rho0=self.rho0)
        momentum[snowplow] = self.ko15.snr_momentum(time=t[snowplow], energy=self.sne_energy, rho0=self.rho0)

        return {"R": R,
                "momentum": momentum,
                "sedov_taylor": sedov_taylor,
                "snowplow": snowplow}

    def _assert_agreement(self, measured, expected, name):
        """
        Assert agreement based on mean relative error and warn about outliers.

        Arguments:
            measured: measured quantity
            expected: expected quantity
            name: name of quantity
        """
        measured_values = np.atleast_1d(measured.value_in(expected.unit))
        expected_values = np.atleast_1d(expected.value_in(expected.unit))

        relative_error = np.abs(measured_values / expected_values - 1)

        mean_error = np.mean(relative_error)
        max_error = np.max(relative_error)

        print(f"{name}: mean relative error {mean_error:.1%}, "
            f"maximum relative error {max_error:.1%}")

        assert mean_error < self.tolerance, (
            f"{name}: mean relative error {mean_error:.1%} "
            f"exceeds tolerance {self.tolerance:.1%}")

    def time_evolution(self):
        """
        Compare characteristic radii against Weaver (1977).
        """
        expected = self._expected()

        if self.plot:
            fig = plt.figure()
            plt.gca().tick_params(axis='both', labelsize=7)
            plt.xlabel('Time [kyr]', fontsize=10)
            plt.ylabel('Radius [pc]', fontsize=10)
            plt.xlim(0, np.max(self.sne_evolution["time"].value_in(units.kyr)))
            plt.yscale('log')
            plt.plot([],[], '-k', label = 'Simulation')
            plt.plot([],[], '--k', label = 'Kim & Ostriker (2015)')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr), self.sne_evolution["R"].value_in(units.pc), '-', color='k')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr)[expected['sedov_taylor']], expected["R"].value_in(units.pc)[expected['sedov_taylor']], '--', color='b')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr)[expected['snowplow']], expected["R"].value_in(units.pc)[expected['snowplow']], '--', color='r')
            plt.plot(self.ko15.free_expansion_time(energy=self.sne_energy, ejecta_mass=self.injected_mass, rho0=self.rho0).value_in(units.kyr),
                     self.ko15.free_expansion_radius(ejecta_mass=self.injected_mass, rho0=self.rho0).value_in(units.pc), 'x', color='g')
            plt.plot(self.ko15.shell_time(energy=self.sne_energy, rho0=self.rho0).value_in(units.kyr),
                     self.ko15.shell_radius(energy=self.sne_energy, rho0=self.rho0).value_in(units.pc), 'x', color='g')
            plt.legend(loc='lower right', fontsize=8, frameon=False)
            plt.savefig(self.saveto / "bubble_radii.png", dpi=120, bbox_inches="tight")
            plt.close(fig)

        for quantity in ("R",):
            self._assert_agreement(self.sne_evolution[quantity], expected[quantity], quantity)

    def momentum_evolution(self):
        """
        Compare radial shell momentum against Weaver (1977).
        """
        expected = self._expected()

        if self.plot:
            fig = plt.figure()
            plt.gca().tick_params(axis='both', labelsize=7)
            plt.xlabel('Time [kyr]', fontsize=10)
            plt.ylabel(r'Momentum [${\rm M}_{\odot}\,{\rm km\, s}^{-1}$]', fontsize=10)
            plt.xlim(0, np.max(self.sne_evolution["time"].value_in(units.kyr)))
            plt.yscale('log')
            plt.plot([],[], '-k', label = 'Simulation')
            plt.plot([],[], '--k', label = 'Kim & Ostriker (2015)')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr), self.sne_evolution["momentum"].value_in(units.MSun * units.km / units.s), '-k')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr)[expected['sedov_taylor']], expected["momentum"].value_in(units.MSun * units.km / units.s)[expected['sedov_taylor']], '--', color='b')
            plt.plot(self.sne_evolution["time"].value_in(units.kyr)[expected['snowplow']], expected["momentum"].value_in(units.MSun * units.km / units.s)[expected['snowplow']], '--', color='r')
            plt.plot(self.ko15.free_expansion_time(energy=self.sne_energy, ejecta_mass=self.injected_mass, rho0=self.rho0).value_in(units.kyr),
                     self.ko15.free_expansion_momentum(energy=self.sne_energy, ejecta_mass=self.injected_mass).value_in(units.MSun * units.km / units.s), 'x', color='g')
            plt.plot(self.ko15.shell_time(energy=self.sne_energy, rho0=self.rho0).value_in(units.kyr),
                     self.ko15.shell_momentum(energy=self.sne_energy, rho0=self.rho0).value_in(units.MSun * units.km / units.s), 'x', color='g')

            plt.legend(loc='lower right', fontsize=8, frameon=False)
            plt.savefig(self.saveto / "momentum.png", dpi=120, bbox_inches="tight")
            plt.close(fig)

        for quantity in ("momentum",):
            self._assert_agreement(self.sne_evolution[quantity], expected[quantity], quantity)

class TestSNe:

    directory = './sne/'
    tolerance = 0.1

    sne = TorchBlastwave(directory, tolerance)

    def test_momentum_evolution(self):
        self.sne.momentum_evolution()

    def test_time_evolution(self):
        self.sne.time_evolution()

def main():

    parser = argparse.ArgumentParser(description="Test Torch supernovae explosion by comparing to Kim & Ostriker (2015).")

    parser.add_argument("directory", 
        type=Path, 
        help="Simulation directory")

    parser.add_argument("--tolerance",
        type=float,
        default=0.1,
        help="Relative tolerance (default: 0.4)")

    parser.add_argument("--number_checks",
        type=float,
        default=10,
        help="Number of outputs to check")

    parser.add_argument("--start_checks",
        type=float,
        default=7,
        help="Time to start checks in kyr")

    parser.add_argument("--plot",
        action="store_true",
        help="Plot tests (tests) or plot everything (all)")

    args = parser.parse_args()
    test = TorchBlastwave(args.directory, args.tolerance, number_checks=args.number_checks, start_checks_time=args.start_checks | units.kyr, plot=args.plot)

    # Tests
    tests = [
        test.time_evolution,
        test.momentum_evolution]

    for test in tests:
        try:
            test()
            print(f"PASSED: {test.__name__}")
        except AssertionError as e:
            print(f"FAILED: {test.__name__}: {e}")

if __name__ == "__main__":
    main()