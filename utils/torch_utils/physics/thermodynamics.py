"""
physics/thermodynamics
==========

Holds functions for gas thermodynamics
"""

from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d

from .constants import cgs

class CoolingCurve:
    """
    Tabulated cooling curve.

    This class loads a precomputed equilibrium cooling table and provides
    interpolated thermodynamic quantities as a function of number density.
    """

    filename = Path(__file__).resolve().parent / "data" / "hAc_b_2.0E-17_e_0.021_FUV_1.69.dat"

    def __init__(self, interpolation='linear', fill_value="extrapolate", mu = 1.3):
        """
        Load cooling curve from file and build interpolators.

        Parameters
        ----------
        interpolation : str
            Interpolation method passed to scipy.interpolate.interp1d.

        fill_value : str or float
            Behavior outside interpolation range.

        mu : float
            Mean molecular weight for computing number density

        Notes
        -----
        - Assumes ndens is monotonic.
        """
        self.mu = mu
        data = np.loadtxt(self.filename, unpack=True)

        _, self.ndens, self.temp, _, self.pk, _, _, _ = data

        self.P_of_n = interp1d(
            self.ndens,
            self.pk,
            kind=interpolation,
            fill_value=fill_value
        )

        self.T_of_n = interp1d(
            self.ndens,
            self.temp,
            kind=interpolation,
            fill_value=fill_value
        )

    def number_density(self, rho):
        """
        Return number density as a function of mass density.

        Parameters
        ----------
        rho : float or ndarray
            Mass density.

        Returns
        -------
        n : float or ndarray
            Number density.
        """
        return rho/self.mu/cgs.mH

    def pressure(self, *, n=None, rho=None):
        """
        Return pressure as a function of number or mass density.

        Parameters
        ----------
        n : float or ndarray, optional
            Number density.
        rho : float or ndarray, optional
            Mass density.

        Returns
        -------
        P : float or ndarray
            Pressure.

        Notes
        -----
        Exactly one of (n, rho) must be provided.
        """
        if (n is None) == (rho is None):
            raise ValueError("Provide exactly one of 'n' or 'rho'.")

        if rho is not None:
            n = self.number_density(rho)

        return self.P_of_n(n)

    def temperature(self, *, n=None, rho=None):
        """
        Return temperature as a function of number or mass density.

        Parameters
        ----------
        n : float or ndarray, optional
            Number density.
        rho : float or ndarray, optional
            Mass density.

        Returns
        -------
        T : float or ndarray
            Interpolated temperature.
        
        Notes
        -----
        Exactly one of (n, rho) must be provided.
        """
        if (n is None) == (rho is None):
            raise ValueError("Provide exactly one of 'n' or 'rho'.")

        if rho is not None:
            n = self.number_density(rho)

        return self.T_of_n(n)
