"""
ic/sphere
==========

Generates initial conditions for Torch simulations
"""
__all__ = [
    "TurbSphere"
]

import numpy as np

from ..numerics.spectral import SpectralFieldGenerator
from ..physics.profiles import GaussianProfile, SchusterProfile
from ..numerics.grid import CartesianGrid
from .io import write_cube, read_cube
from ..physics.constants import cgs
from ..physics.thermodynamics import CoolingCurve
from ..physics.dynamics import VirialAnalysis
from . import plot

PROFILE_REGISTRY = {
    "gaussian": GaussianProfile,
    "schuster": SchusterProfile,
}

class TurbSphere:
    """
    Object to initialize, display, and generate initial condition for
    turbulent sphere, given some setup parameters.

    This class relies on a density profile from ic/profile, spectral
    field generator from ic/spectral, and vairous routines from ic/io
    and ic/plot. 
    
    Turbsph is created using a number of parameters. These are provided
    as keyword arguments, and include
        NCD : list
            Resolution of domain, e.g., [128,128,128]
        CD  : list of lists
            Computational domain [pc], e.g., [[-10,10],[-10,10],[-10,10]] 
        profile : str
            Density profile shape, e.g., "gaussian" 
        klim : list
            Wave numbers for turbulent velocity, e.g., [1,32]
        spectral_slope : float
            Slope of turbulent power spectra, e.g., -5/3 for Kolmogorov
        mu : float
            Mean molecular weight, e.g., 1.3 for cold netural medium
        virial_ratio : float 
            Virial ratio (Ekin/|Epot|), e.g., 0.5 
        amb_rho : float
            Density of ambient medium [cm^-3], e.g., 1e-1
        scale_radius : float
            Scale radius [pc] for density profile, e.g., 7.0
        rho_ratio : float 
            Density ratio for gaussian profile, e.g., 0.3
        mass : float
            Total mass [Msun] of cloud, e.g., 1e4
        seed : int
            Random seed for sampling

    Note
    ----
    All units are in CGS units. Plots are displayed in relevant units.
    """

    def __init__(self, **setup):
        """ 
        Initialized the class.
        
        Arguments
        ---------
            **setup : dict
                Dictonary of intitial condition parameters.
        """
        self.setup = setup
        for k, v in setup.items():
            setattr(self, k, v)

        # Units
        self.mass*=cgs.Msun
        self.scale_radius*=cgs.pc
        self.amb_rho*=cgs.mH*self.mu
        self.CD = np.array(self.CD)*cgs.pc
        
        self.grid = CartesianGrid(self.NCD, self.CD)
        self.profile = self._create_profile()
        self.spectral = SpectralFieldGenerator(
            self.NCD, *self.klim, self.spectral_slope,
            seed=self.seed,
        )
        self.cooling = CoolingCurve(mu=self.mu)

        self.fields = {}

    def _create_profile(self):
        """
        Function generatating density profile from available models.

        Return
        ------
            profile : SphericalDensityProfile
                Class defining the spherical density profile
        
        Notes
        -----
        Profiles available for initial conditions must be included 
        in PROFILE_REGISTRY.
        """
        try:
            cls = PROFILE_REGISTRY[self.profile.lower()]
        except KeyError:
            available = ", ".join(PROFILE_REGISTRY.keys())
            raise ValueError(
                f"Unknown profile '{self.profile}'. Available profiles: {available}"
            )
        return cls(**self._profile_kwargs())
    
    def _profile_kwargs(self):
        """
        Function defining which keyword arguments to be used for different
        denisty profiles.

        Return
        ------
            args : dict
                Dictionary of arguments for specified density profile.
        """
        if self.profile.lower() == "gaussian":
            return dict(mass=self.mass, Rsph=self.scale_radius, rho_ratio=self.rho_ratio)

        elif self.profile.lower() == "schuster":
            return dict(mass=self.mass, Rcore=self.scale_radius, beta=self.beta)

    def write(self, filename):
        """
        Writes initial conditions to file readable by Torch simulation.

        Arguments
        ---------
            filename : str
                Name of initial conditions file.
        """

        data = {"rho":self.fields["rho"],
                "P":self.fields["P"],
                "velx":self.fields["vel"][0],
                "vely":self.fields["vel"][1],
                "velz":self.fields["vel"][2],
                "phi":self.fields['phi']}
        
        write_cube(filename, data=data, NCD=self.NCD)

    def load(self, filename):
        """
        Load Torch initial conditions file and sets the necessary fields.

        Arguments
        ---------
            filename : str
                Name of initial conditions file.
        """
        data, self.NCD = read_cube(filename)
        self.fields["rho"] = data["col0"]
        self.fields["P"] = data["col1"]
        self.fields["vel"] = np.stack([data["col2"], data["col3"], data["col4"]], axis=0)
        self.fields["phi"] = data["col5"]
        
        # Derived fields
        self.fields["T"] = self.cooling.temperature(rho=self.fields["rho"])

    def build_all(self):
        """
        Generates all fields needed to generate initial conditions for Torch
        """
        self.build_density()
        self.build_thermodynamics()
        self.build_velocity()
        self.normalize_virial()

    def build_density(self, Nr=1000):
        """
        Generates density field and derives the gravitational potential
        for that field.

        Arguments
        ---------
            Nr : int
                Number of points used to define 1D radial profile.
        """
        r1d = np.linspace(0, np.max(self.grid.r), Nr)

        rho_1d = self.profile.density(r1d)
        if hasattr(self, "amb_rho"):
            rho_1d = self.profile.match_to_ambient(r1d, rho_1d, self.amb_rho, self.scale_radius)
        phi_1d = self.profile.potential(r1d, rho=rho_1d)

        self.fields["rho"] = self.grid.interp_radial(
            rho_1d, r1d
        )
        self.fields["phi"] = self.grid.interp_radial(
            phi_1d, r1d
        )

    def build_thermodynamics(self):
        """
        Generates pressure and density fields from density.

        Notes
        -----
        If density has not been derived, the function calls self.build_density
        """
        if "rho" not in self.fields:
            self.build_density()

        rho = self.fields["rho"]
        self.fields["P"] = self.cooling.pressure(rho=rho)
        self.fields["T"] = self.cooling.temperature(rho=rho)

    def build_velocity(self):
        """
        Generates initial velocity fields (vx,vy,vz) using the spectral generator.
        """
        vx, vy, vz = self.spectral.generate_vector_field(components=3)
        self.fields["vel"] = np.stack([vx, vy, vz], axis=0)

    def normalize_virial(self):
        """
        Normalize velocities to set specific virial parameter for cloud.

        Notes
        -----
        - Must first initialize density, potential and velocity fields.
        - Assumes Q = E_kin / |E_pot|
        - Should probably include magnetic field in future
        """
        mass = self.fields["rho"] * self.grid.dV
        vel = self.fields["vel"]
        phi = self.fields["phi"]

        vir = VirialAnalysis(mass, vel, phi)

        Q = vir.virial_parameter()
        factor = np.sqrt(self.virial_ratio / Q)

        self.fields["vel"] *= factor
    
    def plot_density_profile(self):
        """
        Shows a density profile.
        """
        plot.profile(self.grid.r/cgs.pc,
                     self.fields['rho']/cgs.mH/self.mu, 
                     ylabel=r'$n\ [{\rm cm}^{-3}$]',
                     xlabel="Radius [pc]",
                     ylog=True,
                     )

    def plot_density(self):
        """
        Shows a slice plot of the density field.
        """
        plot.field_slice(self.fields['rho']/cgs.mH/self.mu, 
                    clim=[1e0,1e3],
                    clabel=r'$n\ [{\rm cm}^{-3}$]',
                    xlabel="x [pc]",
                    ylabel="y [pc]",
                    extent=self.CD/cgs.pc
                    )
        
    def plot_pressure(self):
        """
        Shows a slice plot of the pressure field.
        """
        plot.field_slice(self.fields['press']/cgs.kB, 
                    clim=[1e0,1e3],
                    clabel=r'$P/k_{\rm B}$ [K cm$^{-3}$]',
                    cmap='inferno',
                    xlabel="x [pc]",
                    ylabel="y [pc]",
                    extent=self.CD/cgs.pc
                    )
        
    def plot_temperature(self):
        """
        Shows a slice plot of the temperature field.
        """
        plot.field_slice(self.fields['T'], 
                    clim=[1e0,1e7],
                    cmap='inferno',
                    clabel=r'$T$ [K]',
                    xlabel="x [pc]",
                    ylabel="y [pc]",
                    extent=self.CD/cgs.pc
                    )
        
    def plot_velocity(self, axis='z'):
        """
        Shows a slice plot of the velicity field in spedified direction.

        Arguments
        ---------
            axis : str
                Direction of velocity field.
        """
        axis_map = {'x': 0, 'y': 1, 'z': 2}

        plot.field_slice(self.fields['vel'][axis_map[axis]], 
                    clim=[-100,100],
                    clog=False,
                    clabel=f'$v_{axis}$'+' [km s$^{-1}$]',
                    cmap='seismic',
                    xlabel="x [pc]",
                    ylabel="y [pc]",
                    extent=self.CD/cgs.pc
                    )
    
