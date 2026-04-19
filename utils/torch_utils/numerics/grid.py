"""
numerics/grid
=======

Holds functions for creating grids and mapping data onto the grid.
"""
__all__ = [
    "CartesianGrid"
]

import numpy as np

class CartesianGrid:
    """
    Uniform Cartesian grid with cell-centered coordinates.

    The grid spans a rectangular domain and provides coordinate
    arrays and geometric quantities such as cell size and volume.
    """

    def __init__(self, NCD, CD):
        """
        Initialize the grid.

        Parameters
        ----------
        NCD : tuple of int
            Number of grid cells (Nx, Ny, Nz).

        CD : tuple of tuples
            Domain bounds:
                ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        """
        self.NCD = NCD
        self.CD = CD

        self.nx, self.ny, self.nz = NCD

        # Cell sizes
        self.dx = (CD[0][1] - CD[0][0]) / self.nx
        self.dy = (CD[1][1] - CD[1][0]) / self.ny
        self.dz = (CD[2][1] - CD[2][0]) / self.nz

        self.dV = self.dx * self.dy * self.dz

        # Build coordinates
        self._build_coordinates()

    def _build_coordinates(self):
        """
        Construct cell-centered coordinate arrays.
        """
        ax = np.linspace(self.CD[0][0] + 0.5*self.dx, self.CD[0][1], self.nx)
        ay = np.linspace(self.CD[1][0] + 0.5*self.dy, self.CD[1][1], self.ny)
        az = np.linspace(self.CD[2][0] + 0.5*self.dz, self.CD[2][1], self.nz)

        self.x, self.y, self.z = np.meshgrid(ax, ay, az, indexing="ij")
        self.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def interp_radial(self, field, r, fill_value=0.0):
        """
        Interpolate a radial profile onto the grid.

        Parameters
        ----------
        field : ndarray
            Values defined as function of radius.
        r : ndarray
            Radius.
        fill_value : float, optional
            Value outside interpolation range.

        Returns
        -------
        field3d : ndarray
            Interpolated field on the grid.
        """
        if np.any(np.diff(r) <= 0):
            raise ValueError("Radius array must be strictly increasing.")

        return np.interp(self.r, r, field, left=fill_value, right=fill_value)