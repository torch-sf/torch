""" Generates passive scalar fields and adds them to existing initial conditions file.

    Usage
        python add_scalars.py -f cube128 -mean 0.5 -std 0.05 -s -cp -ps -ds 1 2 4 8
    Arguments
        -f (--filename) str                 - Initial condition filename.
        -ow (--overwrite): bool = False     - Overwrite existing file.
        -mean (--mean_value): float = 0.5   - Mean value of scalar field.
        -std (--standard_deviation): float = 0.05 - Standard deviation of scalar field.
        -ds (--down_sample): list = [1]     - Downsample the resolution of the pattern by a 
                                              factor in base 2
        -b (--box_side): float = 10         - Distance from the center of the cube to the
                                              edge of the cube (half the length of a side).
                                              Only relevant for plotting.
        -s (--scatter): bool = False.       - Generate pattern at uniform resolution where 
                                              the value of each pattern is drawn from a 
                                              normal distribution
    
        -cb (--checkerboard): bool = False  - Generate checkerboard pattern at uniform resolution. 
                                              Pattern is rescaled to mean +- std.
        -ps (--powerspectrum): bool = False - Generate pattern from power spectrum. Pattern is 
                                              normalized to the mean and standard deviation
        -kmin (--kmin): int = 1             - Smallest wavenumber of the turbulence
        -kmax (--kmax): int = 64            - Largest wavenumber of the turbulence
        -e (--turb_exp): float = 5/3        - Exponent of the energy spectrum for the turbulence
        -np (--no_plots): bool = False      - Don't show plots

    Functions

        load_icfile(str: filename, float: boxlen=1.0)
            Load initial conditions from file (e.g., cube128)
        save_icfile(str: filename, array: data, str: header="128 128 128")
            Save initial conditions to file (e.g., cube128)
"""


import numpy as np

import os

def load_icfile(filename, boxlen=1.0):
    """ Load initial condition from file and organizes them for plotting.

        Arguments
            str: filename - Filename of Torch initial conditions file to load.
            float: boxlen - Physical size of initial condition (used to determine positions)

        Return
            array (3, ncell): data - Data in original structure (i,j,k)
            array (3, 2): CD - Computational domain in size given by boxlen
            list (3): NCD - Size of grid (nx, ny, nz)
            array (nx, ny, nz): ii - Index in x
            array (nx, ny, nz): jj - Index in y
            array (nx, ny, nz): kk - Index in z
            array (nx, ny, nz): mx - Position in x [boxlen]
            array (nx, ny, nz): my - Position in y [boxlen]
            array (nx, ny, nz): mz - Position in z [boxlen]
            array (nx, ny, nz): rho - Density [cgs]
            array (nx, ny, nz): pres - Pressure [cgs]
            array (nx, ny, nz): velx - Velocity in x [cgs]
            array (nx, ny, nz): vely - Velocity in y [cgs]
            array (nx, ny, nz): velz - Velocity in z [cgs]
            array (nx, ny, nz): gpot - Gravitational potential [cgs]

    """

    with open(filename, 'r') as f:
        NCD = [int(res) for res in f.readline().split()[1:]]

    CD = np.array(((-boxlen, boxlen), (-boxlen, boxlen), (-boxlen, boxlen)))
    dx = (CD[0][1] - CD[0][0]) / NCD[0]
    dy = (CD[1][1] - CD[1][0]) / NCD[1]
    dz = (CD[2][1] - CD[2][0]) / NCD[2]
    ax = np.arange(CD[0][0]+0.5*dx, CD[0][1], dx)
    ay = np.arange(CD[1][0]+0.5*dy, CD[1][1], dy)
    az = np.arange(CD[2][0]+0.5*dz, CD[2][1], dz)
    (mx, my, mz) = np.meshgrid(ax, ay, az)

    data = np.loadtxt(filename)

    ii = data.T[0].reshape(NCD)
    jj = data.T[1].reshape(NCD)
    kk = data.T[2].reshape(NCD)
    rho = data.T[3].reshape(NCD)
    pres = data.T[4].reshape(NCD)
    velx = data.T[5].reshape(NCD)
    vely = data.T[6].reshape(NCD)
    velz = data.T[7].reshape(NCD)
    gpot = data.T[8].reshape(NCD)

    return data, CD, NCD, ii, jj, kk, mx, my, mz, rho, pres, velx, vely, velz, gpot


def save_icfile(filename, data, header="128 128 128"):
    """ Load initial condition from file and organizes them for plotting.

        Arguments
            str: filename - Filename to save data.
            array (3, ncell): data - Data to save.
            str: header - Header for initial conditions file.

    """

    if os.path.isfile(filename):
        raise FileExistsError(filename+" already exists.")

    fmt = '%3d %3d %3d '
    nvar = data.shape[1]
    for i in range(nvar-3):
        fmt += '%15.7e '

    np.savetxt(filename, data,
               fmt=(fmt),
               header=header)


def scatter(mean, sigma, CD, NCD, downsample=1, lim=[0, 1]):
    """ Generate pattern at uniform resolution where the value of each pattern 
        is drawn from a normal distribution.

        Arguments
            float: mean - Mean value of normal distribution
            float: sigma - Standard devation of normal distribution
            list (3, 2): CD - Computational domain
            list (3): NCD - Resolution of domain
            int: downsample - Number of times resolution is downgraded for pattern
            list (2): lim - Upper and lower limit for values in pattern

        Return
            array (nx, ny, nz): met - Pattern
    """

    # Partial down sampling will not work.
    if NCD[0] % downsample != 0:
        raise ValueError("NCD must be divisible with downsample")

    # Target resolution of pattern
    res = (int(NCD[0]/downsample), int(NCD[1] /
           downsample), int(NCD[2]/downsample))

    # Create 3D structure for scattering values.
    dx = (CD[0][1] - CD[0][0]) / res[0]
    dy = (CD[1][1] - CD[1][0]) / res[1]
    dz = (CD[2][1] - CD[2][0]) / res[2]
    ax = np.arange(CD[0][0]+0.5*dx, CD[0][1], dx)
    ay = np.arange(CD[1][0]+0.5*dy, CD[1][1], dy)
    az = np.arange(CD[2][0]+0.5*dz, CD[2][1], dz)
    (mx, _, _) = np.meshgrid(ax, ay, az)

    # Scatter values around normal distribuion.
    met = np.random.normal(mean, sigma, size=mx.size)
    met[met <= lim[0]] = lim[0]
    met[met >= lim[1]] = lim[1]
    met = met.reshape(res)

    # Perform downsample if warranted.
    if downsample != 1:
        met = met.repeat(downsample, axis=0).repeat(
            downsample, axis=1).repeat(downsample, axis=2)

    return met


def checkerboard(NCD, downsample=1):
    """ Generate checkerboard pattern at uniform resolution. Pattern has values
        0 and 1.

        Arguments
            list (3): NCD - Resolution of domain
            int: downsample - Number of times resolution is downgraded for pattern

        Return
            array (nx, ny, nz): plist - Pattern

    """

    # Partial down sampling will not work.
    if NCD[0] % downsample != 0:
        raise ValueError("NCD must be divisible with downsample")

    # Target resolution of pattern.
    res = (int(NCD[0]/downsample), int(NCD[1] /
           downsample), int(NCD[2]/downsample))

    # Create pattern.
    pattern = np.tile(np.array([[0, 1], [1, 0]]),
                      (int(res[0]/2), int(res[1]/2))).astype(bool)
    plist = [pattern]
    for i in range(res[2]-1):
        if i % 2:
            plist.append(pattern)
        else:
            plist.append(~pattern)

    plist = np.stack(plist)

    # Perform downsample if warranted.
    if downsample != 1:
        plist = plist.repeat(downsample, axis=0).repeat(
            downsample, axis=1).repeat(downsample, axis=2)

    return plist


def distribute_with_powerlaw_spectrum(NCD, slope, freq_min, freq_max):
    """ Generate pattern from power spectrum. Pattern is normalized
        to its mean and standard deviation.

        Arguments
            list (3): NCD - Resolution of domain
            float: slope - Slope of power law
            float: freq_min - Minimim in frequence space
            float: freq_max - Maximum in frequence space

        Return
            array (nx, ny, nz): space - Pattern

    """

    def _kspace(NCD, alpha=-3, kmin=1, kmax=32):

        # Uniform random spectrum covering entire domain [-sqrt(2)/2, sqrt(2)/2],
        # resulting in [-1,1] outside Fourier space.
        space = np.zeros(NCD, dtype=complex)
        space.real = np.sqrt(2)*(np.random.rand(*NCD) - 0.5)
        space.imag = np.sqrt(2)*(np.random.rand(*NCD) - 0.5)

        # Create space with all possible modes.
        [ax, ay, az] = [NCD[i]/2-np.abs(np.arange(NCD[i])-NCD[i]/2)
                        for i in range(len(NCD))]
        (mx, my, mz) = np.meshgrid(ax, ay, az)

        modes = np.sqrt(mx**2 + my**2 + mz**2)
        modes[(modes < kmin)] = 0
        modes[(modes > kmax)] = 0
        # Scale all positive modes by powerlaw to shift power between different scales.
        mask = (modes > 0)
        modes[mask] *= modes[mask]**(alpha)

        return space*modes

    space = np.fft.ifftn(
        _kspace(NCD, alpha=slope, kmin=freq_min, kmax=freq_max)).real
    # rescale to represent distribution
    space = (space - np.mean(space))/np.std(space)
    return space


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser(description=
    """ Script to add passive scalars with different distributions
        to an existing initial conditions file (e.g., cube128)
    """)
    parser.add_argument("-f", "--filename", required=True,
                    help="Initial condition filename.")
    parser.add_argument("-ow", "--overwrite", action='store_true', 
                        required=False, help="Overwrite exisiting file, otherwise create new file")
    
    parser.add_argument("-mean", "--mean_value", default=0.5, required=False, type=float,
                    help="Mean value of scalar field")
    parser.add_argument("-std", "--standard_deviation", default=0.05, required=False, type=float,
                    help="Standard deviation of scalar field")
    parser.add_argument("-ds", "--down_sample", nargs="+", type=float, default=[1],
                        help="""Downsample the resolution of the pattern by a factor in base 2.
                        """)
    parser.add_argument("-b", "--box_side", default=10, required=False, type=float,
                    help="""Distance from the center of the cube to the
                            edge of the cube (half the length of a side).
                            Only relevant for plotting.""")
    
    # Scatter
    parser.add_argument("-s", "--scatter", action='store_true',
                    help="""Generate pattern at uniform resolution where 
                            the value of each pattern is drawn from a 
                            normal distribution.""")
    
    # Checkerboard
    parser.add_argument("-cb", "--checkerboard", action='store_true',
                    help="""Generate checkerboard pattern at uniform resolution. 
                            Pattern is rescaled to mean +- std.""")
    
    # Power-spectrum
    parser.add_argument("-ps", "--powerspectrum", action='store_true',
                    help="""Generate pattern from power spectrum. Pattern is 
                    normalized to the mean and standard deviation.""")
    parser.add_argument("-kmin", "--kmin", default=1, required=False, type=int,
                    help="Smallest wavenumber of the turbulence")
    parser.add_argument("-kmax", "--kmax", default=64, required=False, type=int,
                    help="Largest wavenumber of the turbulence")
    parser.add_argument("-e", "--turb_exp", default=5./3., required=False, type=float,
                    help="Exponent of the energy spectrum for the turbulence")

    parser.add_argument('-np','--no_plots', action='store_true', default=False,
                        help="Don't produce plots")

    args = parser.parse_args()

    data, CD, NCD, i, j, k, x, y, z, dens, pres, vx, vy, vz, pot = load_icfile(
        args.filename, boxlen=args.box_side)

    value = args.mean_value
    sigma = args.standard_deviation
    downsample = args.down_sample
    
    met = []
    # Normal distribution
    if args.scatter:
        for ds in downsample:
            sample = scatter(value, sigma, CD, NCD, downsample=ds)
            met.append(sample)

    # Checkerboard
    if args.checkerboard:
        for ds in downsample:
            pattern = checkerboard(NCD, downsample=ds)
            sample = np.ones(NCD)
            sample[pattern] = value-sigma
            sample[~pattern] = value+sigma
            met.append(sample)

    # Powerspecturm (Kolmogorov scaling)
    if args.powerspectrum:
        for ds in downsample:
            freq = int(args.kmax/ds)
            space = distribute_with_powerlaw_spectrum(
                NCD, slope=-args.turb_exp, freq_min=args.kmin, freq_max=freq)
            sample = (1.0+space*2*sigma)*value
            met.append(sample)

    if not args.no_plots:
        for i in range(len(met)):
            print("Mean, sigma = ", np.mean(met[i].flatten()), np.std(
                met[i].flatten()))
            figsize = plt.figure(figsize=[5, 2.6])
            plt.hist(met[i].flatten(), bins=40, histtype='step', lw=2)
            plt.xlim(0, 1)
            plt.xlabel('Passive scalar')
            plt.ylabel('Number of cells')
            plt.yscale('log')
            plt.show()

            figsize = plt.figure()
            plt.imshow(met[i][:, :, int(NCD[2]/2)], 
                    extent=[-args.box_side, args.box_side, -args.box_side, args.box_side], cmap='seismic', vmin=0.4, vmax=0.6)
            plt.colorbar(label=r'Passive scalar', pad=0.01)
            plt.xlabel('x [pc]')
            plt.ylabel('y [pc]')
            plt.show()

    # Add to exisiting initial condition
    for i in range(len(met)):
        data = np.c_[data, met[i].flatten()]

    if args.overwrite:
        save_icfile(args.filename, data)
    else:
        save_icfile(args.filename+f'_{len(met)}scalars', data)
