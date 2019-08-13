#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

from numpy import *
#import yt
import sys
import argparse
import numpy as np
import os
from scipy.interpolate import interp1d as ip1d

########################
### Units
########################

#MSun = yt.units.Msun.in_cgs().v
#cmpc = yt.units.pc.in_cgs().v
#G    = yt.physical_constants.G.v
#kB   = yt.physical_constants.kb.v
#mH   = yt.physical_constants.mass_hydrogen_cgs.v

# manually extract values from yt and hardcode - AT, 2019 August 13
# yt import is slow.
MSun = 1.98841586e+33
cmpc = 3.08567758e+18
G    = 6.67384e-08
kB   = 1.3806488e-16
mH   = 1.67373522e-24

#######################
### Command line parser
#######################

parser = argparse.ArgumentParser(description='Builds a data set in a cube \
                                 with a Gaussian density profile and a \
                                 turbulent velocity spectrum. Note the default \
                                 background density and temperature assumes the \
                                 background is cool neutral medium (CNM).')

parser.add_argument("-m", "--mass", default=None, required=False, type=float,
                   help="Total mass of the sphere in MSun.")
parser.add_argument("-B", "--magnetic_field", default=0.0, required=False, type=float,
                   help="Background magnetic field magnitude in Gauss.")
parser.add_argument("-n", "--number_density", default=1.25, required=False, type=float,
                   help="Background number density of H in cm^-3.")
parser.add_argument("-r", "--radius", default=None, required=False, type=float,
                   help="Radius of the sphere in parsecs.")
parser.add_argument("-v", "--virial_ratio", default=0.4, required=False, type=float,
                   help="Virial ratio of the sphere. CUrrently only \
                         calculates virial ratio using kinetic and gravitational \
                         energy. NOTE: equilibrium is v=0.5!")
parser.add_argument("-b", "--box_side", default=None, required=False, type=float,
                   help="Distance from the center of the cube to the \
                         edge of the cube (half the length of a side). \
                         If not specified, is Rsph*1.25 by default.")

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

parser.add_argument("-mus", "--musph", default=1.3, required=False, type=float,
                   help="Mean particle mass of sphere gas in mH.")
parser.add_argument("-mua", "--muamb", default=1.3, required=False, type=float,
                   help="Mean particle mass of ambient gas in mH.")
parser.add_argument("-Ts", "--sphere_temperature", default=30., required=False, type=float,
                   help="Sphere gas temperature in Kelvin. \
                         Default is 30 Kelvin.")
parser.add_argument("-Tb", "--background_temperature", default=8e3, required=False, type=float,
                   help="Background gas temperature in Kelvin. \
                         Default is 8e3 Kelvin.")
parser.add_argument("--Ts_from_cool_curve", action='store_true',
                   help="Compute Ts for given density from equilibrium \
                         cooling curve. Overrides input sphere temperature.")
parser.add_argument("--rho_match", action='store_true',
                   help="Match density at sphere edge to ambient medium dens \
                         by tanh smoothing or flooring. \
                         This alters the sphere mass.")

parser.add_argument("-f", "--filename", default=None, required=False,
                   help="Output filename.")
parser.add_argument('-nd','--no_data', action='store_true', default=False,
                    help="Don't write a data output file, just do the calculations. \
                          Default is false (write data).")
parser.add_argument('-wp','--with_plots', action='store_true', default=True,
                    help="Make slice plots of velx, pres, temp, dens")
parser.add_argument('-rf','--read_file', default=None,
                    help="Read all input data from file, formatted in the following order: \
                          Mass	Radius	box	virial_ratio	NumDensSph	Tsph	Tamb	kmin	kmax	Eslp	Bmag	filename \
                          Pass in the input filename.")

#############################
### Definitions
#############################

def create_spspace_box(kmin, kmax, NCD, Eslp):
    ssbox = zeros(NCD, dtype=complex)

    # first, create uniform complex random spectrum in whole array
    # range of coeffs is (-sqrt(2)/2, sqrt(2)/2), it results in
    # range (-1,1) for coeffiecient magnitudes
    ssbox.real = sqrt(2.)*(random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)
    ssbox.imag = sqrt(2.)*(random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)

    # create mask that filters selected modes and gives weights according to the
    # given spectrum (Eslp)
    ax = NCD[0]/2-abs(arange(NCD[0], dtype=float64)-NCD[0]/2)
    ay = NCD[1]/2-abs(arange(NCD[1], dtype=float64)-NCD[1]/2)
    az = NCD[2]/2-abs(arange(NCD[2], dtype=float64)-NCD[2]/2)

    (mx, my, mz) = meshgrid(ax,ay,az)

    mask = mx*mx + my*my + mz*mz
    mask[where(mask < kmin*kmin)] = 0
    mask[where(mask > kmax*kmax)] = 0
    # E(k) ~ k^Eslp; number of modes grows with k^2 => Eslp-2
    # v ~ sqrt(E) => 0.5*Eslp-1, k^2 is in mask => 0.25*Eslp-0.5
    mask[where(mask>0)] = mask[where(mask>0)]**(0.25*Eslp-0.5)

    ssbox *= mask

    return ssbox


def kolmogorov_vel(NCD, kmin, kmax, Eslp):
    # create velocity field in spectral space (sp_vel[xyz])
    # and FFT it to configuration space (vel[xyz])
    sp_velx = create_spspace_box(kmin, kmax, NCD, Eslp)
    sp_vely = create_spspace_box(kmin, kmax, NCD, Eslp)
    sp_velz = create_spspace_box(kmin, kmax, NCD, Eslp)

    velx = fft.ifftn(sp_velx).real
    vely = fft.ifftn(sp_vely).real
    velz = fft.ifftn(sp_velz).real

    return (velx, vely, velz)


# numerically calculates grav. potential of a spherically symmetric density
# distribution given by a radial profile
def calc_sphsym_pot(r_arr, rho_arr):
    C4pio3 = 4*pi/3.0
    dr_arr       = correlate(r_arr, [-1,1])     # thicknesses of shells between r-points
    rmid_arr     = convolve(r_arr, [0.5,0.5])   # mid-points, rmid[0] = 0
    rmid_arr[-1] = r_arr[-1]                    # last rmid segment end at the surface
    drmid_arr    = correlate(rmid_arr, [-1,1])  # thicknesses of shells between mid-points
    r3mid_arr    = rmid_arr**3                  # r^3 at mid-points
    dr3mid_arr   = correlate(r3mid_arr, [-1,1]) # volumes of shells between mid-points
    dMrmid_arr   = rho_arr*dr3mid_arr*C4pio3    # masses of shells -"-
    Mrmid_arr    = dMrmid_arr.cumsum()          # enclosed mass below mid-points
    Mrmid_arr    = insert(Mrmid_arr, 0, 0.0)    # central point
    Frmid_arr    = -G*Mrmid_arr/(rmid_arr**2+1e-99) # grav. force at mid-points
    dPhi_arr     = Frmid_arr[1:-1] * dr_arr     # prepare for integration of Fr*dr
    pot_arr      = zeros(r_arr.shape[0], dtype=float64)
    pot_arr[-1]  = -G*Mrmid_arr[-1]/r_arr[-1]   # -G*M/R at the sphere surface
    pot_arr[:-1] = dPhi_arr                     # copy Fr*dr into the pot_array
    pot_arr      = cumsum(pot_arr[::-1])[::-1]  # integrate from surface to the centre
    #print 'G = ', G
    #print 'dr_arr = ', dr_arr, dr_arr.shape
    #print 'rmid_arr = ', rmid_arr, rmid_arr.shape
    #print 'drmid_arr = ', drmid_arr, drmid_arr.shape
    #print 'r3mid_arr = ', r3mid_arr, r3mid_arr.shape
    #print 'dr3mid_arr = ', dr3mid_arr, dr3mid_arr.shape
    #print 'dMrmid_arr = ', dMrmid_arr, dMrmid_arr.shape
    #print 'Mrmid_arr = ', Mrmid_arr, Mrmid_arr.shape
    #print 'Frmid_arr = ', Frmid_arr, Frmid_arr.shape
    #print 'pot_arr = ', pot_arr, pot_arr.shape
    return (pot_arr, Mrmid_arr, Frmid_arr)

def gauss_dens_prof(Rsph, Msph, rho_rat, Nr):
    rarr   = np.linspace(0.0, Rsph, Nr+1)
    sig_R  = Rsph / sqrt(-log(rho_rat)) # characteristic radius
    rho_rarr = exp(-rarr*rarr/(sig_R*sig_R))
    # calculate mass of the sphere with temporary density profile
    (pot_rarr, Mrmid_rarr, Frmid_rarr) = calc_sphsym_pot(rarr, rho_rarr)
    # scale the density profile to get correct sphere mass
    rho_rarr *= Msph / Mrmid_rarr[-1]
    return (rarr, rho_rarr)


def schuster_dens_prof(Rsph, Msph, Rcore, beta, Nr):
    rarr   = arange(0.0, Rsph+Rsph/Nr, Rsph/Nr)
    rho_rarr = 1.0 / (1.0 + (rarr*rarr)/(Rcore*Rcore))**beta
    # calculate mass of the sphere with temporary density profile
    (pot_rarr, Mrmid_rarr, Frmid_rarr) = calc_sphsym_pot(rarr, rho_rarr)
    # scale the density profile to get correct sphere mass
    rho_rarr *= Msph / Mrmid_rarr[-1]
    return (rarr, rho_rarr)


def dens_pot_3darr(Rsph, Msph, rarr, rho_rarr, rho_amb, Nr, NCD, CD):
    # calculate gravitational potential
    (pot_rarr, Mrmid_rarr, Frmid_rarr) = calc_sphsym_pot(rarr, rho_rarr)

    # create computational domain
    dx = (CD[0][1] - CD[0][0]) / NCD[0]
    dy = (CD[1][1] - CD[1][0]) / NCD[1]
    dz = (CD[2][1] - CD[2][0]) / NCD[2]
    ax = arange(CD[0][0]+0.5*dx, CD[0][1], dx)
    ay = arange(CD[1][0]+0.5*dy, CD[1][1], dy)
    az = arange(CD[2][0]+0.5*dz, CD[2][1], dz)
    (mx, my, mz) = meshgrid(ax,ay,az)
    rarr = sqrt(mx*mx + my*my + mz*mz)
    rind_arr = (rarr * Nr / Rsph).astype(int)
    CDSize = sqrt((CD[0][1]-CD[0][0])**2 + (CD[1][1]-CD[1][0])**2 \
    +             (CD[2][1]-CD[2][0])**2)

    # extend radius, density and potential 1D arrays to cover whole CD
    rarr_ext = arange(0.0, CDSize, Rsph/Nr)
    rho_rarr_ext = ones((rarr_ext.shape[0]), dtype=float64)*rho_amb
    rho_rarr_ext[0:Nr+1] = rho_rarr
    pot_rarr_ext = -G*Msph/(rarr_ext+1e-99)
    pot_rarr_ext[0:Nr+1] = pot_rarr

    # fill 3D arrays of density and potential by linear interpolation
    r3d = rarr_ext[rind_arr]
    weight_arr = (rarr - r3d) * Nr / Rsph
    rho_arr = (1-weight_arr) * rho_rarr_ext[rind_arr] \
    +            weight_arr  * rho_rarr_ext[rind_arr+1]
    pot_arr = (1-weight_arr) * pot_rarr_ext[rind_arr] \
    +            weight_arr  * pot_rarr_ext[rind_arr+1]

#    for i in range(NCD[0]):
#        for j in range(NCD[1]):
#            for k in range(NCD[2]):
#                print '%3d %3d %3d %15.7e %15.7e %15.7e %15.7e %15.7e' \
#                % (i,j,k, rarr[i,j,k], r3d[i,j,k] \
#                , weight_arr[i,j,k], rho_arr[i,j,k], pot_arr[i,j,k])
#            print
#        print

    return (rarr, rho_arr, pot_arr)


def dens_pot_3darr_noextrap(r_rarr, rho_rarr, NCD, CD):
    """No extrapolation w/rho_amb, just apply density profile to cube domain"""
    # calculate gravitational potential
    (pot_rarr, Mrmid_rarr, Frmid_rarr) = calc_sphsym_pot(r_rarr, rho_rarr)

    # create computational domain
    dx = (CD[0][1] - CD[0][0]) / NCD[0]
    dy = (CD[1][1] - CD[1][0]) / NCD[1]
    dz = (CD[2][1] - CD[2][0]) / NCD[2]
    ax = arange(CD[0][0]+0.5*dx, CD[0][1], dx)
    ay = arange(CD[1][0]+0.5*dy, CD[1][1], dy)
    az = arange(CD[2][0]+0.5*dz, CD[2][1], dz)

    (mx, my, mz) = meshgrid(ax,ay,az)
    r_arr = sqrt(mx*mx + my*my + mz*mz)

    assert np.amax(r_arr) <= r_rarr[-1]

    rho_arr = np.interp(r_arr, r_rarr, rho_rarr)
    pot_arr = np.interp(r_arr, r_rarr, pot_rarr)

    return (r_arr, rho_arr, pot_arr)


# The actual code to make the data file.
def make_data_cube(Msph, Rsph, box, n0, Tsph, T_amb, musph, mu_amb, vir_rat,
                   kmin, kmax, Eslp, Bmag, filename, write_data,
                   Ts_from_cool_curve=False,
                   cool_curve='hAc_b_2.0E-17_e_0.021_FUV_1.69.dat',
                   rho_match=False,
                   with_plots=True):

    rho_rat = 1.0/3.0  # density ratio between border and centre
    Nr      = 10000    # Number of points in 1-D
    f_trunc = 0.05      # rho-matching, hardcoded constant for tanh function, reasonable for NCD=(128,128,128)

    rho_amb = n0*mH*mu_amb # Ambient density

    # computational domain
    CD   = array(((-box, box), (-box, box), (-box, box)), dtype=float64)
    NCD  = (128,128,128)

    if rho_match:  # new method - AT, 2019 Aug 13

        r_rarr = linspace(0, 1.74*box, Nr+1)  # 1.74 is just above sqrt(3)
        rho_rarr = exp(-r_rarr**2/Rsph**2 * -1*log(rho_rat))
        indRsph = np.searchsorted(r_rarr, Rsph)  # index of smallest r satisfying r >= Rsph

        # calculate mass of the sphere with temporary density profile
        (pot_rarr, Mrmid_rarr, Frmid_rarr) = calc_sphsym_pot(r_rarr[:indRsph], rho_rarr[:indRsph])
        # scale the density profile to get correct sphere mass
        rho_rarr *= Msph / Mrmid_rarr[-1]

        # smoothly match radial density profile to ambient medium.
        # this alters the enclosed mass.  gaussian dens profile concentrates
        # mass at center, so tweaking edge dens should not alter mass too much.
        if rho_rarr[indRsph-1] > rho_amb:

            rho_rarr[rho_rarr < rho_amb] = rho_amb  # BEFORE applying kernel, enforce rho >= rho_amb everywhere

            # kernel(r=0) -> 1, kernel(r=Rsph) = 0.5, kernel(r->infty) -> 0
            kernel = 0.5*(np.tanh((Rsph-r_rarr)/f_trunc/Rsph)+1)
            rho_rarr = np.exp( (np.log(rho_rarr) - np.log(rho_amb))*kernel + np.log(rho_amb) )

            # pressure match based on cooling curve.
            # I'm too lazy to write "general" solution. - AT, 2019 Aug 13
            assert Ts_from_cool_curve  # not tested without --Ts_from_cool_curve, so just require it
            assert cool_curve == 'hAc_b_2.0E-17_e_0.021_FUV_1.69.dat'
            # density regime with dP/dn < 0 at thermal equilibrium
            # is unstable and will evolve towards warm or cold ISM.
            # enforce pressure continuity and skip the evolve.
            unstable_ndens_min = 0.5888  # P/kB = 4.711E+03  # also skip the stable range n=0.5888 to 1.871
            #unstable_ndens_min = 1.871  # P/kB = 1.156E+04
            unstable_ndens_max = 23.50  # P/kB = 4.712E+03

            assert rho_rarr[0]/mH/musph > unstable_ndens_max
            assert n0 < unstable_ndens_min
            rho_rarr[:indRsph] = np.maximum(rho_rarr[:indRsph], unstable_ndens_max*mH*musph)
            rho_rarr[indRsph:] = np.minimum(rho_rarr[indRsph:], unstable_ndens_min*mH*mu_amb)

        elif rho_rarr[indRsph-1] < rho_amb:
            # apply a density floor
            print("Warning: sphere edge density is below ambient density, flooring...")
            rho_rarr[rho_rarr < rho_amb] = rho_amb

        if (with_plots):
            import matplotlib.pyplot as plt
            plt.plot(r_rarr/Rsph, rho_rarr, '-k')
            plt.yscale('log')
            plt.xlabel('r/Rsph')
            plt.ylabel('rho (g/cm3)')
            plt.savefig(filename+'profile.png')
            plt.clf()

        (r_arr, rho_arr, pot_arr) = dens_pot_3darr_noextrap(r_rarr, rho_rarr, NCD, CD)

    else:

        # calculate density, pressure and potential fields
        (r_rarr, rho_rarr) = gauss_dens_prof(Rsph, Msph, rho_rat, Nr)
        (r_arr, rho_arr, pot_arr) = dens_pot_3darr(Rsph, Msph, r_rarr, rho_rarr, rho_amb, Nr, NCD, CD)

    mask = (r_arr <= Rsph).astype(float)

    if Ts_from_cool_curve:
        # Lets just set the temperature initially from the density in the sphere
        # (which is likely in the unstable regime) directly from the equilibrium
        # cooling curve calculated previously.
        numdens_arr   = rho_arr/musph/mH*mask + rho_arr/mu_amb/mH*(1.0-mask)
        data_P, T_arr = get_P_and_T_from_Eq_Cooling_Curve(numdens_arr, data_file=cool_curve)
        p_arr = (kB/musph/mH)*mask*rho_arr*T_arr + (kB/mu_amb/mH)*(1.0-mask)*rho_arr*T_arr
    else:

        # Wunsch's old method
        p_arr = (kB*Tsph/musph/mH)*mask*rho_arr + (kB*T_amb/mu_amb/mH)*(1.0-mask)*rho_arr

        # My method where I pressure match at the boundary with the sphere. - JW
        # NOTE: This just doesn't matter. Anything more massive than a few 100 solar masses
        #       is just Jeans unstable and collapsing, its not in "pressure equ" with the
        #       surrounding medium. And adding any amount of density doesn't make it in
        #       equ with the surrounding, because its collapsing. If anything, its likely
        #       *accreting from the surroundings* (see Vazquez-Semadeni 2009).
        #       Therefore we're going back to separately setting the ambient and
        #       core temps and pressures independently.  - JW 8/30/18
        #p_arr = (kB*Tsph/musph/mH)*mask*rho_arr + (kB*Tsph/musph/mH)*(1.0-mask)*rho_rarr[-1] # Psph edge = Pamb
        # Invert to get the ambient density.

    # this OVERRIDES the rho_amb passed to dens_pot_3darr(...)
    if not rho_match:
        rho_arr = rho_arr*mask + p_arr/(kB*T_amb/mu_amb/mH)*(1.0-mask)

    # calculate turbulent velocity field
    (velx, vely, velz) = kolmogorov_vel(NCD, kmin, kmax, Eslp)

    # calculate total kinetic and potential energy
    dx = (CD[0][1] - CD[0][0]) / NCD[0]
    dy = (CD[1][1] - CD[1][0]) / NCD[1]
    dz = (CD[2][1] - CD[2][0]) / NCD[2]
    Ekin = 0.5*(mask*rho_arr*(velx*velx+vely*vely+velz*velx)).sum()*dx*dy*dz
    Epot = 0.5*(mask*rho_arr*pot_arr).sum()*dx*dy*dz
    Emag = 0.5*(4.0/3.0*np.pi*Rsph**3.0*(Bmag**2.0/4.0/np.pi))
    Qvir = (Ekin) / np.abs(Epot)
    #Qvir = (Ekin+Emag*0.5) / np.abs(Epot)
    print 'Ekin, Epot, Emag/2, Qvir = ', Ekin, Epot, 0.5*Emag, Qvir
    # rescale velocity field to get pre-set virial ratio
    velx *= sqrt(vir_rat / Qvir)*mask
    vely *= sqrt(vir_rat / Qvir)*mask
    velz *= sqrt(vir_rat / Qvir)*mask

    # recalculate kinetic energy
    Ekin = 0.5*(mask*rho_arr*(velx*velx+vely*vely+velz*velx)).sum()*dx*dy*dz
    Qvir = (Ekin) / np.abs(Epot)
    #Qvir = (Ekin+Emag*0.5) / np.abs(Epot)
    print 'Ekin, Epot, Emag/2, Qvir = ', Ekin, Epot, 0.5*Emag, Qvir
    print 'Msph, Rsph = ', Msph/MSun, Rsph/cmpc

    if (with_plots):

        # plot a slice of the velocity field
        import matplotlib.pyplot as plt
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
        plt.savefig(filename+'velx.png')
        plt.clf()

        fig, ax = plt.subplots()
        im = ax.imshow(p_arr[:,:,64]/kB, aspect='auto', cmap='viridis',
                       norm=mpl.colors.LogNorm())
        fig.colorbar(im)
        plt.savefig(filename+'pres.png')
        plt.clf()

        fig, ax = plt.subplots()
        im = ax.imshow(p_arr[:,:,64]/(((kB/musph/mH)*mask[:,:,64] + (kB/mu_amb/mH)*(1.0-mask[:,:,64]))*rho_arr[:,:,64]),
                       aspect='auto', cmap='viridis', norm=mpl.colors.LogNorm())
        fig.colorbar(im)
        plt.savefig(filename+'temp.png')
        plt.clf()

        fig, ax = plt.subplots()
        im = ax.imshow(rho_arr[:,:,64], aspect='auto', cmap='viridis',
                       norm=mpl.colors.LogNorm())
        fig.colorbar(im)
        plt.savefig(filename+'dens.png')
        plt.clf()

    if (write_data):

        print "writing metadata to "+filename+".dat"
        with open(filename+'.dat', 'w') as f:

            f.write("# {:<9} {:<35}\n".format('filename:', filename))
            if Ts_from_cool_curve:
                f.write("# cooling curve for Tsph: {}\n".format(cool_curve))
            f.write("\n")

            def write_pars(key, unit, val):
                assert len(key) == len(unit)
                assert len(key) == len(val)
                # construct format strings
                # hdr: "# {:>8} {:>10} {:>10} . . ."
                # val: "{:>10.2E} {:>10.2E} {:>10.2E} . . ."
                hdrfmt = "# {:>8}" + " {:>10}"*(len(key)-1) + "\n"
                valfmt = " ".join(["{:>10.2E}"]*len(key)) + "\n"
                f.write(hdrfmt.format(*key))
                f.write(hdrfmt.format(*unit))
                f.write(valfmt.format(*val))

            write_pars(
                ['Msph', 'Rsph', 'x/y/z max', 'Ekin', 'Epot', 'Emag', 'Qvir', 'B0', 'n0'],
                ['(Msun)', '(pc)', '(pc)', '(erg)', '(erg)', '(erg)', '(-)', '(Gauss)', '(cm^-3)'],
                [Msph/MSun, Rsph/cmpc, box/cmpc, Ekin, Epot, Emag, Qvir, Bmag, n0],
            )

            f.write('\n')

            write_pars(
                ['B0', 'n0', 'Tsph', 'Tamb', 'musph', 'muamb', 'kmin', 'kmax', 'turb_exp'],
                ['(Gauss)', '(cm^-3)', '(K)', '(K)', '(mH)', '(mH)', '(-)', '(-)', '(-)'],
                [Bmag, n0, Tsph, T_amb, musph, mu_amb, kmin, kmax, Eslp],
            )

        # Write out the initial conditions file.
        print "writing data file."
        f = open(filename, 'w')
        f.write('# {} {} {} \n'.format(NCD[0], NCD[1], NCD[2]) )
        for i in range(NCD[0]):
            for j in range(NCD[1]):
                for k in range(NCD[2]):
                   savetxt(f, array((i,j,k, rho_arr[i,j,k], p_arr[i,j,k]\
                    , velx[i,j,k], vely[i,j,k], velz[i,j,k], pot_arr[i,j,k])).reshape(1,9), \
                    fmt=('%3d %3d %3d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e') )
                f.write('\n')

    f.close()


def get_P_and_T_from_Eq_Cooling_Curve(numdens, data_file='hAc_b_2.0E-17_e_0.021_FUV_1.69.dat'):

    [time, ndens, temp, ei, pk, xHp, mu_mol, tdust] = np.loadtxt(
                                                   data_file, unpack=True)

    p_from_n = ip1d(ndens, pk, kind='linear')
    T_from_n = ip1d(ndens, temp, kind='linear')

    P = p_from_n(numdens)
    T = T_from_n(numdens)

    return P, T

################################
### Main code.
################################

args = parser.parse_args()

infile     = args.read_file
write_data = not args.no_data

if (infile is None):


    if (args.radius is None and args.box_side is None):

        if   (args.mass == 1e3):
            args.radius     =  5.0
            args.box_side   =  7.0
        elif (args.mass == 1e4):
            args.radius     = 10.0
            args.box_side   = 12.5
        elif (args.mass == 1e5):
            args.radius     = 50.0
            args.box_side   = 55.0

        print "Using the default factors for Rsph and box size:"
        print "--radius", args.radius, "--box_size", args.box_side

    elif (args.radius is not None and args.box_side is None):

        args.box_side = args.radius * 1.25

        print "Using default box size +/-(1.25x Rsph):"
        print "--box_size", args.box_side


    if (args.radius is None or args.mass is None or args.virial_ratio is None or args.filename is None):
        raise Exception("Error: You must either pass sphere radius, mass,"
                        + " virial ratio and output filename or you must pass"
                        + " in an input file with the proper values specified!")


    # sphere params
    Rsph       = [args.radius*cmpc]
    Msph       = [args.mass*MSun]
    box        = [args.box_side*cmpc]
    vir_rat    = [args.virial_ratio]
    kmin       = [args.kmin]
    kmax       = [args.kmax]
    Eslp       = [-args.turb_exp]
    Tsph       = [args.sphere_temperature]
    T_amb      = [args.background_temperature]
    musph      = [args.musph]
    mu_amb     = [args.muamb]
    filename   = [args.filename]
    n0         = [args.number_density]
    Bmag       = [args.magnetic_field]
    num_runs   = 1


else: # load from the file

    #print np.shape(np.genfromtxt(infile, dtype=None, skip_header=2, usecols=range(8), unpack=True))

    # Not tested/working currently. 2019 Mar 22, AT

    Msph, Rsph, box, vir_rat, n0, Tsph, T_amb, kmin, kmax, Eslp, Bmag = np.loadtxt(
            infile, dtype=np.float, skiprows=2, usecols=range(11), unpack=True)
    filename = np.loadtxt(infile, dtype=np.str, skiprows=2, usecols=11, unpack=True)

    num_runs = len(Msph)
    Rsph *= cmpc
    Msph *= MSun


for i in range(num_runs):

    make_data_cube(Msph[i], Rsph[i], box[i], n0[i], Tsph[i], T_amb[i],
                   musph[i], mu_amb[i],
                   vir_rat[i], kmin[i], kmax[i], Eslp[i], Bmag[i],
                   filename[i], write_data,
                   Ts_from_cool_curve=args.Ts_from_cool_curve,
                   rho_match=args.rho_match,
                   with_plots=args.with_plots)

