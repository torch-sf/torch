"""
physics/profiles
===========

Holds classes for generating spherically symmetric density profiles.
"""
__all__ = [
    "GaussianProfile",
    "SchusterProfile",
]

import numpy as np

from .constants import cgs

class SphericalDensityProfile:
    """
    Base class for spherically symmetric density profiles.

    The class defines a radial density profile rho(r) and provides
    numerical methods to compute derived quantities such as enclosed
    mass and gravitational potential.

    Subclasses must implement the density(r) method.

    The gravitational potential is computed assuming spherical symmetry:
        M(r) = \int 4\pi r^2 ρ(r) dr
        Φ(r) = -\int G M(r) / r^2 dr

    Boundary condition:
        Φ(R) = -G M(R) / R
    """

    def __init__(self, **params):
        """
        Initialize the density profile.

        Arguments
        ---------
        **params : dict
            Dictionary containing profile-specific parameters.
        """
        self.params = params
        for key, value in params.items():
            setattr(self, key, value)

    def __repr__(self):
        """
        Return a string representation of the profile.
        
        The representation includes the class name and the model
        parameters in a format resembling the constructor call:
            ClassName(param1=value1, param2=value2, ...)
        """
        params_str = ", ".join(f"{k}={v}" for k, v in self.params.items())
        return f"{self.__class__.__name__}({params_str})"

    def _validate_radius(self, r):
        """
        Validate radial coordinate array.

        Ensures that r is one-dimensional, non-negative,
        and strictly increasing.

        Arguments
        ---------
        r : ndarray
            Radial coordinate.

        Returns
        -------
        r : ndarray
            Validated radius array.
        """
        r = np.asarray(r)

        if r.ndim != 1:
            raise ValueError("Radius array r must be one-dimensional.")

        if np.any(r < 0):
            raise ValueError("Radius values must be non-negative.")

        if np.any(np.diff(r) <= 0):
            raise ValueError("Radius array must be strictly increasing.")

        return r

    def match_to_ambient(self, r, field, background, Rsph, f_trunc=0.05):
        """
        Smoothly match a radial field to a background value.

        Parameters
        ----------
        r : ndarray
            Radial coordinate.

        field : ndarray
            Field defined on r (e.g. density, pressure).

        background : float
            Background (ambient) value.

        Rsph : float
            Characteristic radius where transition occurs.

        f_trunc : float
            Controls width of transition region.

        Returns
        -------
        field : ndarray
            Modified field with smooth transition to background.

        Notes
        -----
        - Matching is performed in log-space.
        - Ensures field >= background everywhere.
        """
        r = self._validate_radius(r)
        field = np.asarray(field)

        # Floor
        field = np.maximum(field, background)

        # Smooth transition kernel
        kernel = 0.5 * (np.tanh((Rsph - r) / (f_trunc * Rsph)) + 1.0)

        # Log-space blending
        field = np.exp(
            (np.log(field) - np.log(background)) * kernel + np.log(background)
        )

        return field

    def density(self, r):
        """
        Evaluate the density profile rho(r).

        Arguments
        ---------
        r : ndarray
            Radial coordinate.

        Returns
        -------
        rho : ndarray
            Density evaluated at r.

        Notes
        -----
        This method relies on self.shape, defined in subclasses.
        """
        r = self._validate_radius(r)

        shape = self.shape(r)

        dr = np.gradient(r)
        dM = 4 * np.pi * r**2 * shape * dr
        M = np.cumsum(dM)
        A = self.mass / M[-1]

        return A * shape

    def enclosed_mass(self, r, rho=None):
        """
        Compute enclosed mass M(r).

        Arguments
        ---------
        r : ndarray
            Radial coordinate.
        rho : ndarray, optional
            Density profile if modified from class function

        Returns
        -------
        M : ndarray
            Enclosed mass profile.
        """
        r = self._validate_radius(r)
        dr = np.gradient(r)
        if rho is None:
            rho = self.density(r)

        dM = 4 * np.pi * r**2 * rho * dr
        return np.cumsum(dM)

    def potential(self, r, rtol=1e-12, rho=None):
        """
        Compute gravitational potential Φ(r).

        The potential is calculated by integrating the gravitational
        acceleration inward from the outer boundary:

            Φ(R) = -G M(R) / R

        Arguments
        ---------
        r : ndarray
            Radial coordinate.
        rtol : float, optional
            Small regularization term to avoid division by zero at r=0.
        rho : ndarray, optional
            Density profile if modified from class function

        Returns
        -------
        Phi : ndarray
            Gravitational potential profile.
        """
        r = self._validate_radius(r)
        dr = np.gradient(r)
        M = self.enclosed_mass(r, rho=rho)

        g = cgs.G * M / (r**2 + rtol)

        Phi = np.zeros_like(r)
        Phi[-1] = -cgs.G * M[-1] / r[-1]

        dPhi = g[:-1] * dr[:-1]

        Phi[:-1] = Phi[-1] - np.cumsum(dPhi[::-1])[::-1]

        return Phi

class GaussianProfile(SphericalDensityProfile):
    """
    Gaussian density profile:

        \rho(r) \propto \exp(-r^2 / \sigma^2)

    where \sigma is set such that \rho(Rsph)/\rho(0) = rho_ratio.
    """

    def __init__(self, mass, Rsph, rho_ratio):
        super().__init__(mass=mass, Rsph=Rsph, rho_ratio=rho_ratio)
        self.sigma = Rsph / np.sqrt(-np.log(rho_ratio))

    def shape(self, r):
        return np.exp(-(r**2) / self.sigma**2)
    
class SchusterProfile(SphericalDensityProfile):
    """
    Schuster density profile:

        \rho(r) \propto (1 + (r/Rcore)^2)^(-\beta)
    """

    def __init__(self, mass, Rcore, beta):
        super().__init__(mass=mass, Rcore=Rcore, beta=beta)

    def shape(self, r):
        return 1.0 / (1.0 + (r / self.Rcore)**2)**self.beta
    