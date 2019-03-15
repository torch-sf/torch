#!/usr/bin/env python

from numpy import *
import yt
import fftw3
import sys
import argparse
import numpy as np
import os
from scipy.interpolate import interp1d as ip1d


#filename = 'clouds/cubeM3V04E'

########################
### Units
########################

MSun = yt.units.Msun.in_cgs().v
cmpc = yt.units.pc.in_cgs().v
G    = yt.physical_constants.G.v
kB   = yt.physical_constants.kb.v
mH   = yt.physical_constants.mass_hydrogen_cgs.v

#######################
### Command line parser
#######################

parser = argparse.ArgumentParser(description='Builds a data set in a cube \
                                 with a Gaussian density profile and a \
                                 turbulent velocity spectrum. Note the default \
                                 background density and temperature assumes the \
                                 background is cool neutral medium (CNM).')

parser.add_argument("-m", "--mass", default=None, required=False, type=float,
                   help="Input the total mass of the sphere in MSun.")
#parser.add_argument("-mus", "--musph", default=2.3, required=False, type=float,
#                   help="Atomic mass of gas in the sphere.")
parser.add_argument("-B", "--magnetic_field", default=0.0, required=False, type=float,
                   help="Input the background magnetic field magnitude in Gauss.")
parser.add_argument("-n", "--number_density", default=1.25, required=False, type=float,
                   help="Input the background number density of H in cm^-3.")
#parser.add_argument("-mua", "--muamb", default=1.3, required=False, type=float,
#                   help="Atomic mass of gas in the sphere.")
parser.add_argument("-r", "--radius", default=None, required=False, type=float,
                   help="Input the radius of the sphere in parsecs.")
parser.add_argument("-v", "--virial_ratio", default=0.4, required=False, type=float,
                   help="Input the virial ratio of the sphere. CUrrently only \
                         calculates virial ratio using kinetic and gravitational \
                         energy. NOTE: equilibrium is v=0.5!")
parser.add_argument("-b", "--box_side", default=None, required=False, type=float,
                   help="Input the distance from the center of the cube to the \
                         edge of the cube (half the length of a side). \
                         If not specified, is Rsph*1.25 by default.")
parser.add_argument("-f", "--filename", default=None, required=False,
                   help="Output filename.")
parser.add_argument("-kmin", "--kmin", default=1, required=False, type=int,
                   help="Input the smallest wavenumber of the turbulence. \
                         Default is 1.")
parser.add_argument("-kmax", "--kmax", default=32, required=False, type=int,
                   help="Input the largest wavenumber of the turbulence. \
                         Default is 32.")
parser.add_argument("-e", "--turb_exp", default=5./3., required=False, type=float,
                   help="Input the exponent of the energy spectrum for the turbulence. \
                         Default is Kolmogorov (-5/3). Note that you pass the positive \
                         value and the code tacks on the negative sign.")
parser.add_argument("-Ts", "--sphere_temperature", default=30., required=False, type=float,
                   help="Input the sphere temperature of the gas in Kelvin. \
                         Default is 30 Kelvin.")
parser.add_argument("-Tb", "--background_temperature", default=8e3, required=False, type=float,
                   help="Input the background temperature of the gas in Kelvin. \
                         Default is 8e3 Kelvin.")
parser.add_argument('-nd','--no_data', action='store_true', default=False,
                    help="Don't write a data output file, just do the calculations. \
                          Default is false (write data).")
parser.add_argument('-rf','--read_file', default=None,
                    help="Read all input data from file, formatted in the following order: \
                          Mass	Radius	box	virial ratio	NumDensSph	Tsph	Tamb	filename \
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

    # multiply random coefficients with mask
    ssbox *= mask

    #ssbox = zeros(NCD, dtype=complex)
    #ssbox[63,63,63] = 1.0 + 0.0j
    #return ssbox
    return ssbox


def kolmogorov_vel(NCD, kmin, kmax, Eslp):
    # prepare fftw plan
    inarr  = zeros(NCD, dtype=complex)
    outarr = zeros(NCD, dtype=complex)
    fft = fftw3.Plan(inarr, outarr, direction='backward', flags=['measure'])

    # arrays for velocity field
    velx  = zeros(NCD, dtype=float64)
    vely  = zeros(NCD, dtype=float64)
    velz  = zeros(NCD, dtype=float64)

    # create velocity field in spectral space (sp_vel[xyz])
    # and FFT it to configuration space (vel[xyz])
    sp_velx = create_spspace_box(kmin, kmax, NCD, Eslp)
    inarr[:,:,:] = sp_velx
    fft.execute()
    velx[:,:,:] = outarr.real
    sp_vely = create_spspace_box(kmin, kmax, NCD, Eslp)
    inarr[:,:,:] = sp_vely
    fft.execute()
    vely[:,:,:] = outarr.real
    sp_velz = create_spspace_box(kmin, kmax, NCD, Eslp)
    inarr[:,:,:] = sp_velz
    fft.execute()
    velz[:,:,:] = outarr.real


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
    rarr   = arange(0.0, Rsph+Rsph/Nr, Rsph/Nr)
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

def window_avg(x,n):

    from scipy.signal import convolve


    """
    A function that averages over the elements of an array x.

    Attributes
    ----------
    x       : numpy array (float)
              Array that you want to average over.
    n       : iteger
              Number of elements that you want to average over (window radius).
    avg_all : Whether you want to average over all elements within the window (True)
              or only those at either end of the window boundary (False). Default is True.
    """

    # Numpy version does 1 d arrays only.
    #return np.convolve(x, np.ones(n)*1.0/float(n), 'full') #[(n-1):-(n-1)]
    # Scipy version can do any d array. Nice.
    return convolve(x, np.ones((n,n,n))*1.0/float(n**3.0), mode='same')

# The actual code to make the data file.
def make_data_cube(Msph, Rsph, box, n0, Tsph, T_amb, vir_rat, kmin, kmax, Eslp, Bmag, filename, write_data):

    #Msph    = 1e4*MSun # Sphere mass
    #Tsph    = 30.      # Sphere temperature
    musph   = 1.3 #2.3      # molecular mass inside sphere
    rho_rat = 1.0/3.0  # density ratio between border and centre
    Nr      = 1000     # Number of points in 1-D

    # ambient gas
    mu_amb  = 1.3      # molecular mass of ambient gas
    rho_amb = n0*mH*mu_amb # Ambient density
    #T_amb   = 100.     # Ambient temperature

    #box     = 7.0     # dist from center to edge of box in pc
    #box    = 12.5
    #box     = 55.0     # dist from center to edge of box in pc

    # energy spectrum
    #Eslp = -5./3.      # spectrum exponent (maybe also try Burger spectrum?)
    #vir_rat = 0.4      # virial ratio
    # per M-MML 99 make kmax=2, kmin=1
    #kmin = 1           # longest wave number of turbulence
    #kmax = 32          # shorest wave number of turbulence

    # computational domain
    CD   = array(((-box, box), (-box, box), (-box, box)), dtype=float64)
    NCD  = (128,128,128)

    # calculate density, pressure and potential fields
    (r_rarr, rho_rarr) = gauss_dens_prof(Rsph, Msph, rho_rat, Nr)
    (rarr, rho_arr, pot_arr) = dens_pot_3darr(Rsph, Msph, r_rarr, rho_rarr, rho_amb, Nr, NCD, CD)
    mask = (rarr <= Rsph).astype(float)

    # Lets just set the temperature initially from the density in the sphere
    # (which is likely in the unstable regime) directly from the equilibrium
    # cooling curve calculated previously.
    #numdens_arr   = rho_arr/musph/mH*mask + rho_arr/mu_amb/mH*(1.0-mask)
    #data_P, T_arr = get_P_and_T_from_Eq_Cooling_Curve(numdens_arr)
    #p_arr = (kB/musph/mH)*mask*rho_arr*T_arr + (kB/mu_amb/mH)*(1.0-mask)*rho_arr*T_arr

    # Wunsch's old method
    p_arr = (kB*Tsph/musph/mH)*mask*rho_arr + (kB*T_amb/mu_amb/mH)*(1.0-mask)*rho_arr

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

    #velx *= sqrt((vir_rat*np.abs(Epot)-Emag*0.5)/Ekin)*mask
    #vely *= sqrt((vir_rat*np.abs(Epot)-Emag*0.5)/Ekin)*mask
    #velz *= sqrt((vir_rat*np.abs(Epot)-Emag*0.5)/Ekin)*mask

    # recalculate kinetic energy
    Ekin = 0.5*(mask*rho_arr*(velx*velx+vely*vely+velz*velx)).sum()*dx*dy*dz
    Qvir = (Ekin) / np.abs(Epot)
    #Qvir = (Ekin+Emag*0.5) / np.abs(Epot)
    print 'Ekin, Epot, Emag/2, Qvir = ', Ekin, Epot, 0.5*Emag, Qvir
    print 'Msph, Rsph = ', Msph/MSun, Rsph/cmpc

    # plot a slice of the velocity field
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colors

    fig, ax = plt.subplots()
    cmap = plt.get_cmap('viridis')
    #bounds = [velx.min(), velx.max()]
    #print "bounds", bounds
    #cNorm = colors.BoundaryNorm(velx[:,:,64], cmap.N)
    #print cNorm
    im = ax.imshow(velx[:,:,64], aspect='auto', cmap=cmap) #, norm=cNorm)
    fig.colorbar(im)
    plt.savefig(filename+'velx.png')
    plt.clf()

    fig, ax = plt.subplots()
    #bounds = [rho_arr[np.where(rho_arr > 0.0)[0]].min(), rho_arr.max()]
    #print "bounds", bounds
    #cNorm = colors.BoundaryNorm(rho_arr[:,:,64], cmap.N)
    im = ax.imshow(np.log10(p_arr[:,:,64]/kB), aspect='auto', cmap=cmap) #, norm=cNorm)
    fig.colorbar(im)
    plt.savefig(filename+'pres.png')
    plt.clf()

    fig, ax = plt.subplots()
    #bounds = [rho_arr[np.where(rho_arr > 0.0)[0]].min(), rho_arr.max()]
    #print "bounds", bounds
    #cNorm = colors.BoundaryNorm(rho_arr[:,:,64], cmap.N)
    im = ax.imshow(np.log10(p_arr[:,:,64]/(((kB/musph/mH)*mask[:,:,64] + (kB/mu_amb/mH)*(1.0-mask[:,:,64]))*rho_arr[:,:,64])), aspect='auto', cmap=cmap) #, norm=cNorm)
    fig.colorbar(im)
    plt.savefig(filename+'temp.png')
    plt.clf()

    fig, ax = plt.subplots()
    #bounds = [rho_arr[np.where(rho_arr > 0.0)[0]].min(), rho_arr.max()]
    #print "bounds", bounds
    #cNorm = colors.BoundaryNorm(rho_arr[:,:,64], cmap.N)
    im = ax.imshow(np.log10(rho_arr[:,:,64]), aspect='auto', cmap=cmap) #, norm=cNorm)
    fig.colorbar(im)
    plt.savefig(filename+'dens.png')
    plt.clf()

    #print '# ', NCD[0], NCD[1], NCD[2]
    #for i in range(NCD[0]):
        #for j in range(NCD[1]):
            #for k in range(NCD[2]):
                #print '%3d %3d %3d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e' \
                #% (i,j,k, rho_arr[i,j,k], p_arr[i,j,k]\
                #, velx[i,j,k], vely[i,j,k], velz[i,j,k], pot_arr[i,j,k])
            #print
        #print

    # Write info about what was set when making this data cube.
    f = open(filename+'.dat', 'w')
    f.write("{:<9} {:<35}  \n \n".format('filename:', filename))
    f.write("{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n".format('Msph', 'Rsph', 'x/y/z max', 'Ekin', 'Epot', 'Emag', 'Qvir', 'B0', 'n0'))
    f.write("{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} \n".format('(M)', '(pc)', '(pc)', '(ergs)', '(ergs)', '(ergs)', ' ', '(Gauss)', '(cm^-3)'))
    f.write("{:>10.2E} {:>10.2E} {:>10.2E} {:>10.2E} {:>10.2E} {:>10.2E} {:>10.2E} {:>10.2E} \n".format(Msph/MSun, Rsph/cmpc, box, Ekin, Epot, Emag, Qvir, Bmag))
    f.close()

    if (write_data):
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
    filename   = [args.filename]
    n0         = [args.number_density]
    Bmag       = [args.magnetic_field]
    num_runs   = 1


else: # load from the file

    #print np.shape(np.genfromtxt(infile, dtype=None, skip_header=2, usecols=range(8), unpack=True))

    Msph, Rsph, box, vir_rat, n0, Tsph, T_amb, kmin, kmax, Eslp, Bmag = np.loadtxt(
            infile, dtype=np.float, skiprows=2, usecols=range(11), unpack=True)
    filename = np.loadtxt(infile, dtype=np.str, skiprows=2, usecols=11, unpack=True)

    num_runs = len(Msph)
    Rsph *= cmpc
    Msph *= MSun


for i in range(num_runs):

    make_data_cube(Msph[i], Rsph[i], box[i], n0[i], Tsph[i], T_amb[i],
                   vir_rat[i], kmin[i], kmax[i], Eslp[i], Bmag[i],
                   filename[i], write_data)

