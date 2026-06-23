#!/usr/bin/env python

######################
#
# Modification of turb-sphere.py in order to build a cube with uniform density 
# of the specified density and temperature. Developed by Eric Andersson, based 
# on the turb-sphere.py script.
# 
# Eventually this should also be set up to generate a uniform cube with turbulent 
# initial velocities, but this doesn't work yet.
#  -SA  20241209
#
######################


import numpy as np
from scipy.interpolate import interp1d

import argparse
import os

# Constants
mH   = 1.67373522e-24
kB   = 1.3806488e-16

parser = argparse.ArgumentParser(description='Creates cube file with uniform density.\
                                              Can include turbulent velocity field')

parser.add_argument("-rho", "--density", type=float,
                   help="Particle number density in cube [cm^-3].")
parser.add_argument("-temp", "--temperature", default=-1, type=float,
                   help="Temperature in cube [K]. If temp<0, derive pressure from cooling curve.")
parser.add_argument("-res", "--resolution", default=128, type=int,
                   help="Resolution of data.")
parser.add_argument("-turb", "--use_turb", default=False, type=bool,
                   help="Add turbulent velocity field")
parser.add_argument("-rms", "--turb_rms_velocity", default=1.0, type=float,
                   help="Turbulent velocity field [km/s]")
parser.add_argument("-kmi", "--kmin", default=1, type=int,
                   help="Smallest wavenumber of the turbulence.")
parser.add_argument("-kma", "--kmax", default=32, type=int,
                   help="Largest wavenumber of the turbulence.")
parser.add_argument("-exp", "--turb_exp", default=5./3., type=float,
                   help="Exponent of the energy spectrum for the turbulence. \
                         Default is Kolmogorov (-5/3). Note that you pass the positive \
                         value and the code tacks on the negative sign.")
parser.add_argument("-vel", "--vel_disp", default=20, required=False, type=float,
                   help="Velocity dispersion of the turbulent velocities [cm/s].")
parser.add_argument("-s", "--seed", default=-1, type=int,
                   help="Random seed.")
parser.add_argument("-o", "--filename", default="cube", type=str,
                   help="Filename for output.")

def write(filename, cube):
    
    # Write out the initial conditions file.
    
    with open(filename, 'w') as file:
        file.write('# {} {} {} \n'.format(*cube.shape[1:]) )
        for i in range(cube.shape[1]):
            for j in range(cube.shape[2]):
                for k in range(cube.shape[3]):
                    np.savetxt(file, np.asarray([i,j,k, cube[0,i,j,k], cube[1,i,j,k],
                        cube[2,i,j,k], cube[3,i,j,k], cube[4,i,j,k], cube[5,i,j,k]]).reshape(1,9),
                        fmt=('%3d %3d %3d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e') )


def create_spspace_box(NCD, kmin, kmax, exp):
    ssbox = np.zeros(NCD, dtype=complex)

    # first, create uniform complex random spectrum in whole array
    # range of coeffs is (-sqrt(2)/2, sqrt(2)/2), it results in
    # range (-1,1) for coeffiecient magnitudes
    ssbox.real = np.sqrt(2.)*(np.random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)
    ssbox.imag = np.sqrt(2.)*(np.random.rand(NCD[0],NCD[1],NCD[2]) - 0.5)

    # create mask that filters selected modes and gives weights according to the
    # given spectrum
    ax = NCD[0]/2-np.abs(np.arange(NCD[0], dtype='float64')-NCD[0]/2)
    ay = NCD[1]/2-np.abs(np.arange(NCD[1], dtype='float64')-NCD[1]/2)
    az = NCD[2]/2-np.abs(np.arange(NCD[2], dtype='float64')-NCD[2]/2)

    (mx, my, mz) = np.meshgrid(ax,ay,az)

    mask = mx*mx + my*my + mz*mz
    mask[np.where(mask < kmin*kmin)] = 0
    mask[np.where(mask > kmax*kmax)] = 0
    # E(k) ~ k^exp; number of modes grows with k^2 => exp-2
    # v ~ sqrt(E) => 0.5*exp-1, k^2 is in mask => 0.25*exp-0.5
    mask[np.where(mask>0)] = mask[np.where(mask>0)]**(0.25*exp-0.5)

    ssbox *= mask

    return ssbox


def kolmogorov_vel(NCD, kmin, kmax, exp):
    # create velocity field in spectral space (sp_vel[xyz])
    # and FFT it to configuration space (vel[xyz])
    sp_velx = create_spspace_box(NCD, kmin, kmax, exp)
    sp_vely = create_spspace_box(NCD, kmin, kmax, exp)
    sp_velz = create_spspace_box(NCD, kmin, kmax, exp)

    velx = np.fft.ifftn(sp_velx).real
    vely = np.fft.ifftn(sp_vely).real
    velz = np.fft.ifftn(sp_velz).real

    return velx, vely, velz

def get_temp_from_cooling_curve(dens, data_file='hAc_b_2.0E-17_e_0.021_FUV_1.69.dat'):

    [_, dens_cooling, temp_cooling, _, _, _, _, _] = np.loadtxt(
                                                   data_file, unpack=True)

    interp = interp1d(dens_cooling, temp_cooling, kind='linear')

    return interp(dens)

if __name__ == "__main__":
    args = parser.parse_args()
    
    # Initialize cube
    cube = np.ones([6]+[args.resolution]*3) # [Density, Pressure, 3*Velocity, Potential]
    
    # Density
    cube[0] *= args.density*1.3*mH # factor 1.3 is mean particle mass
    
    # Pressure
    if args.temperature<0:
        args.temperature = get_temp_from_cooling_curve(args.density)
    cube[1] *= args.density*kB*args.temperature # Ideal gas law

    # Velocity
    if args.use_turb:        
        NCD = (args.resolution,args.resolution,args.resolution)
        # calculate turbulent velocity field
        (velx, vely, velz) = kolmogorov_vel(NCD, args.kmin, args.kmax, args.turb_exp)
        v_mag = np.sqrt(velx**2 + vely**2 + velz**2)
        print('Raw velocity dispersion', np.std(v_mag))
        vel_ratio = args.vel_disp/np.std(v_mag)

        # scale velocities to achieve desired dispersion
        velx *= vel_ratio
        vely *= vel_ratio
        velz *= vel_ratio

        v_mag = np.sqrt(velx**2 + vely**2 + velz**2)
        print('Scaled velocity dispersion', np.std(v_mag))

        cube[2], cube[3], cube[4] = velx, vely, velz
    else:
        cube[2:5] *= 0.0

    # Potential
    cube[5] *= 0.0
    

    filename = args.filename 
    if os.path.isfile(filename):
        raise FileExistsError(f"File with name {filename} already exisits")
    
    write(filename, cube)
