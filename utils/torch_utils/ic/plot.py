"""
ic/plot
==========

Holds plotting routines for initial conditions.
"""
__all__ = [
    "field_slice",
    "profile"
]

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import numpy as np
from scipy.stats import binned_statistic

def field_slice(field,
                axis='z', slice_at=0.0,
                clog=True, clim=None, clabel=None, cmap='viridis',
                xlabel=None, ylabel=None,
                xlim=None, ylim=None, extent=None):
    """
    Plot a field along a slice along a given axis.

    Arguments
    ---------
        field : ndarray
            Field 3D data
        axis : str
            Axis perpendicular to slice plane
        slice_at : float
            Slice plane position along the axis in terms of extent
        clog : bool
            Set logscale for field colorbar
        clim : list or tupple
            Limits on colorbar, if not None
        clabel : str
            Label on colorbar, if not None
        cmap : str
            Colormap for colorbar
        xlabel : str
            Label on x-axis
        ylabel : str
            Label on y-axis
        xlim : list or tupple
            Limis on x-axis
        ylim : list or tupple
            Limits on y-axis
        extent : list of lists
            Extent of field [[xmax,xmax],[ymax,ymax],[zmax,zmax]]. 
            If None, assumes boxlen centered on 0.
    """
    if extent is None:
        extent = [[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]]

    clabel = clabel or "Field"

    if clim is None:
        clim = (np.min(field), np.max(field))

    if clog:
        clim = (max(clim[0], 1e-30), clim[1])
        norm = LogNorm(*clim)
    else:
        norm = Normalize(*clim)

    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_idx = axis_map[axis]
    amin, amax = extent[axis_idx]
    n = field.shape[axis_idx]

    slice_idx = int((slice_at - amin) / (amax - amin) * (n - 1))
    slice_idx = np.clip(slice_idx, 0, n - 1)

    if axis_idx == 0:
        img = field[slice_idx, :, :]
        extent = [*extent[1], *extent[2]]
        xlabel = xlabel or "y"
        ylabel = ylabel or "z"
    elif axis_idx == 1:
        img = field[:, slice_idx, :]
        extent = [*extent[0], *extent[2]]
        xlabel = xlabel or "x"
        ylabel = ylabel or "z"
    else:
        img = field[:, :, slice_idx]
        extent = [*extent[0], *extent[1]]
        xlabel = xlabel or "x"
        ylabel = ylabel or "y"

    # Plot
    plt.figure()
    plt.gca().set_aspect('equal')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    plt.imshow(img.T,
               norm=norm,
               extent=extent,
               origin="lower",
               cmap=cmap)

    plt.colorbar(label=clabel, pad=0.01)
    plt.show()

def profile(r, field,
            bins = 128,
            statistic="mean",
            ylabel = "Field",
            xlabel = "Radius",
            xlim = None,
            ylim = None,
            xlog = False,
            ylog = False):

    plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    if isinstance(bins, int):
        if xlog:
            rbins = np.logspace(np.min(r), np.max(r), bins)
        else:
            rbins = np.linspace(0, np.max(r), bins)

    profile, rbins, _  = binned_statistic(r.flatten(), field.flatten(), statistic=statistic, bins=rbins)
    plt.step(rbins[1:], profile, color='k', lw=2)
    plt.show()