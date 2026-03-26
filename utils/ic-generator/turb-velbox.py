#!/usr/bin/env python

from __future__ import division, print_function

from argparse import ArgumentParser
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import ifftn
from os import path

# manually extract values from yt and hardcode
# CCC 16/07/2024 - following AT, 2019 August 13, turb-sphere.py
MSun = 1.98841586e+33
cmpc = 3.08567758e+18
G    = 6.67384e-08
kB   = 1.3806488e-16
mH   = 1.67373522e-24

def create_spspace_box(NCD, kmin, kmax, Eslp):
    """
    NCD = (nx, ny, nz)
    kmin, kmax = wavenumber min/max
    Eslp = slope of energy spectrum
    """
    ssbox = np.zeros(NCD, dtype=complex)
    # first, create uniform complex random spectrum in whole array
    # range of coeffs is (-sqrt(2)/2, sqrt(2)/2), it results in
    # range (-1,1) for coefficient magnitudes
    ssbox.real = np.sqrt(2.)*(np.random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)
    ssbox.imag = np.sqrt(2.)*(np.random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)

    # create mask that filters selected modes and gives weights according to the
    # given spectrum (Eslp)
    # each array runs from 0, 1, 2, ..., ndim/2-1, ndim/2, ndim/2-1, ..., 2, 1
    # (i.e., standard fft frequency layout)
    ax = NCD[0]/2 - abs(np.arange(NCD[0],dtype=np.float64) - NCD[0]/2)
    ay = NCD[1]/2 - abs(np.arange(NCD[1],dtype=np.float64) - NCD[1]/2)
    az = NCD[2]/2 - abs(np.arange(NCD[2],dtype=np.float64) - NCD[2]/2)

    (mx, my, mz) = np.meshgrid(ax,ay,az)

    mask = mx*mx + my*my + mz*mz  # k^2
    mask[np.where(mask < kmin*kmin)] = 0
    mask[np.where(mask > kmax*kmax)] = 0
    # E(k) ~ k^Eslp; number of modes grows with k^2 => Eslp-2
    # v ~ sqrt(E) => 0.5*Eslp-1, k^2 is in mask => 0.25*Eslp-0.5
    mask[np.where(mask>0)] = mask[np.where(mask>0)]**(0.25*Eslp-0.5)

    ssbox *= mask

    return ssbox


def kolmogorov_vel(NCD, kmin, kmax, Eslp):
    """
    Create velocity field in spectral space (sp_vel[xyz])
    and FFT it to configuration space (vel[xyz])
    """
    sp_velx = create_spspace_box(NCD, kmin, kmax, Eslp)
    sp_vely = create_spspace_box(NCD, kmin, kmax, Eslp)
    sp_velz = create_spspace_box(NCD, kmin, kmax, Eslp)

    velx = ifftn(sp_velx).real  # np.fft.ifftshift is not needed
    vely = ifftn(sp_vely).real
    velz = ifftn(sp_velz).real

    return velx, vely, velz


def scale_vel(velx, vely, velz, target_v):
    """
    Scale the velocity field so that the velocity dispersion
    is equal to target_v in cm/s
    """

    vel_disp = np.sqrt((velx*velx+vely*vely+velz*velz).sum())
    print('vel disp', vel_disp)
    velx *= target_v / vel_disp
    vely *= target_v / vel_disp
    velz *= target_v / vel_disp
    
    return velx, vely, velz


def pressure(NCD, rho_arr, T, mu):
    """
    Create uniform pressure distribution from
    volume density and temperature
    """
    p_arr = (kB/mu/mH)*rho_arr*T
    
    return p_arr


def density(NCD, L, sigma):
    """
    Create volume density distribution from
    surface density in solar masses per pc^2
    """
    densbox = np.ones(NCD)
    
    Lcm = L*cmpc     # Box side in cm
    Sgg = sigma*MSun*(cmpc)**(-2) # Mass in g, box side in cm
    
    density = Sgg / Lcm
    
    densbox *= density
    
    return densbox
    


if __name__ == '__main__':

    parser = ArgumentParser(description="Build a cube data set with \
                                         turbulent velocity spectrum")
    parser.add_argument("-n", "--nblock", default=128, required=False, type=int,
                        help="elements per cube side. Default is 128.")
    parser.add_argument("-L", "--box_side", default=4, required=False, type=float,
                        help="Box side in parsecs.")
    parser.add_argument("-kmin", "--kmin", default=1, required=False, type=int,
                        help="Smallest wavenumber of the turbulence. \
                              Default is 1.")
    parser.add_argument("-kmax", "--kmax", default=32, required=False, type=int,
                        help="Largest wavenumber of the turbulence. \
                              Default is 32.")
    parser.add_argument("-e", "--turb_exp", default=5./3., required=False, type=float,
                        help="Exponent of the energy spectrum for the turbulence. \
                              Default is Kolmogorov (-5/3). Note that you pass the positive \
                              value and the code tacks on the negative sign.")
    parser.add_argument("-S", "--sigma", default=100, required=False, type=int,
                        help="Surface density in MSun/pc^2. Default is 100.")
    parser.add_argument("-T", "--temp", default=100, required=False, type=int,
                        help="Gas temperature in K. Default is 100.")
    parser.add_argument("-mu", "--mu", default=1.3, required=False, type=int,
                        help="Mean particle mass in mH. Default is 1.3.")
    parser.add_argument("-v", "--target_v", default=1e5, required=False, type=int,
                        help="Target velocity dispersion in cm/s. Default is 1 km/s.")
    parser.add_argument("-s", "--seed", default=0, required=False, type=int,
                        help="Random seed.")
    parser.add_argument("-f", "--filename", default=None, required=False,
                        help="Output filename.")
    parser.add_argument("-np","--no_plots", action='store_true', default=False,
                        help="Don't make slice plots of velx, pres, temp, dens. \
                        Default is true (write plot files).")
    parser.add_argument("-c", "--clobber", action='store_true',
                        help="Overwrite existing filename?")

    args = parser.parse_args()
    
    if args.seed != -1:
        np.random.seed(args.seed)

    NCD = (args.nblock, args.nblock, args.nblock)

    print("NCD=", NCD)
    print("kmin=", args.kmin)
    print("kmax=", args.kmax)
    print("Eslp=", -1*args.turb_exp)
    print("filename=", args.filename)
    velx, vely, velz = kolmogorov_vel(NCD, args.kmin, args.kmax,
                                        -1*args.turb_exp)
    velx, vely, velz = scale_vel(velx, vely, velz, args.target_v)
    
    rho_arr = density(NCD, args.box_side, args.sigma)
    print("density=", rho_arr[0, 0, 0])
    
    p_arr = pressure(NCD, rho_arr, args.temp, args.mu)
    print("pressure=", p_arr[0, 0, 0])
    print("speed of sound=", np.sqrt((5./3)*p_arr[0, 0, 0]/rho_arr[0, 0, 0]))
    
    # Potential is zero everywhere, since the density is uniform
    pot_arr = np.zeros(NCD)

    if args.filename is not None:
        # Write out the initial conditions file.
        print("writing to", args.filename)

        if path.exists(args.filename):
            if args.clobber:
                print("File exists, clobbering...")
            else:
                raise Exception("File exists, not clobbering")

        started = datetime.now()
        with open(args.filename, 'w') as f:
            f.write('# {} {} {} \n'.format(NCD[0], NCD[1], NCD[2]) )
            # Low-priority TODO: rewrite this by flattening arrays
            for i in range(NCD[0]):
                for j in range(NCD[1]):
                    for k in range(NCD[2]):
                        arr = np.array((i, j, k, rho_arr[i,j,k], p_arr[i,j,k], velx[i,j,k], vely[i,j,k], velz[i,j,k], pot_arr[i,j,k]))
                        arr = arr.reshape(1,9)
                        np.savetxt(f, arr, fmt=('%3d %3d %3d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e'))
                    #f.write('\n')  # not sure what hte purpose of this newline is, or how it works with fortran read - AT 20190322
        print("File write elapsed", datetime.now() - started)
        
        if not args.no_plots: # Added plots following turb-sphere.py, CCC 16/07/2024

            # plot a slice of the velocity field
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            from matplotlib import cm
            from matplotlib import colors

            def symshow(ax, x, vbnd=None, **kwargs):
                """mpl.plt.imshow wrapper, force vmin/vmax symmetric about zero"""
                assert 'vmin' not in kwargs
                assert 'vmax' not in kwargs
                if vbnd is None:
                    filt = np.isfinite(x)
                    vbnd = np.max(np.abs(x[filt]))
                if 'cmap' not in kwargs:
                    kwargs['cmap'] = 'RdBu_r'
                return ax.imshow(x, vmin=-vbnd, vmax=+vbnd, **kwargs)

            fig, ax = plt.subplots()
            im = symshow(ax, velx[:,:,64], aspect='auto')
            fig.colorbar(im)
            plt.savefig(args.filename+'velx.png')
            plt.clf()

            #fig, ax = plt.subplots()
            #im = ax.imshow(p_arr[:,:,64]/kB, aspect='auto', cmap='viridis',
            #               norm=mpl.colors.LogNorm())
            #fig.colorbar(im)
            #plt.savefig(args.filename+'pres.png')
            #plt.clf()

            #fig, ax = plt.subplots()
            #im = ax.imshow(p_arr[:,:,64]/(((kB/musph/mH)*mask[:,:,64] + (kB/mu_amb/mH)*(1.0-mask[:,:,64]))*rho_arr[:,:,64]),
            #               aspect='auto', cmap='viridis', norm=mpl.colors.LogNorm())
            #fig.colorbar(im)
            #plt.savefig(args.filename+'temp.png')
            #plt.clf()

            fig, ax = plt.subplots()
            im = ax.imshow(rho_arr[:,:,64], aspect='auto', cmap='viridis',
                           norm=mpl.colors.LogNorm())
            fig.colorbar(im)
            plt.savefig(args.filename+'dens.png')
            plt.clf()

    print("Done!")
