"""
physics/dynamics
==========

Holds functions for energy related quantities
"""
__all__ = [
    "VirialAnalysis"
]

import numpy as np

class VirialAnalysis:
    """
    Compute energy components and virial parameter for a system.

    The class operates on fields defined on a grid, where
    each cell carries a mass and velocity.

    Notes
    -----
    - The input mass is assumed to be the mass per cell (i.e. rho * dV).
    - The velocity field must have shape (3, ...), where the first axis
      corresponds to (vx, vy, vz).
    - Magnetic energy is not included in this class.
    """

    def __init__(self, mass, vel, pot, mask=None):
        """
        Initialize the virial analysis.

        Arguments
        ---------
        mass : ndarray
            Mass per cell.
        vel : ndarray
            Velocity field with shape (3, ...).
        pot : ndarray
            Gravitational potential field.
        mask : ndarray, optional
            Mask selecting region of interest (e.g. sphere).
        """
        self.mass = mass
        self.vel = vel
        self.pot = pot
        self.mask = mask if mask is not None else np.ones_like(mass)

    def kinetic_energy(self):
        """
        Compute kinetic energy.

        Returns
        -------
        Ekin : float
            Total kinetic energy.
        """
        v2 = np.sum(self.vel**2, axis=0)
        return 0.5 * np.sum(self.mask * self.mass * v2)

    def potential_energy(self):
        """
        Compute gravitational potential energy.

        Returns
        -------
        Epot : float
            Total potential energy.
        """
        return 0.5 * np.sum(self.mask * self.mass * self.pot)

    def virial_parameter(self):
        """
        Compute virial parameter:

            Q = E_kin / |E_pot|

        Returns
        -------
        Q : float
            Virial parameter.
        """
        return self.kinetic_energy() / np.abs(self.potential_energy())
