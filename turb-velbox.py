#!/usr/bin/env python

from __future__ import division, print_function

from argparse import ArgumentParser
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import ifftn
from os import path

def create_spspace_box(NCD, kmin, kmax, Eslp):
    """
    NCD = (nx, ny, nz)
    kmin, kmax = wavenumber min/max
    Eslp = slope of energy spectrum
    """
    ssbox = np.zeros(NCD, dtype=np.complex)
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

    return (velx, vely, velz)


if __name__ == '__main__':

    parser = ArgumentParser(description="Build a cube data set with \
                                         turbulent velocity spectrum")
    parser.add_argument("-n", "--nblock", default=128, required=False, type=int,
                        help="elements per cube side. Default is 128.")
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
    parser.add_argument("-f", "--filename", default=None, required=False,
                        help="Output filename.")
    parser.add_argument("-c", "--clobber", action='store_true',
                        help="Overwrite existing filename?")

    args = parser.parse_args()

    NCD = (args.nblock, args.nblock, args.nblock)

    print("NCD=", NCD)
    print("kmin=", args.kmin)
    print("kmax=", args.kmax)
    print("Eslp=", -1*args.turb_exp)
    print("filename=", args.filename)
    (velx, vely, velz) = kolmogorov_vel(NCD, args.kmin, args.kmax,
                                        -1*args.turb_exp)

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
                        arr = np.array((i, j, k, velx[i,j,k], vely[i,j,k], velz[i,j,k]))
                        arr = arr.reshape(1,6)
                        np.savetxt(f, arr, fmt=('%3d %3d %3d %15.7e %15.7e %15.7e'))
                    #f.write('\n')  # not sure what hte purpose of this newline is, or how it works with fortran read - AT 20190322
        print("File write elapsed", datetime.now() - started)

    print("Done!")
