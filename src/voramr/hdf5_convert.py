# voramr_convert.py
#
# Contains utility functions to extract data from
# a Voronoi data structure output to manipulate
# and resave as a file of data recognizable by FLASH.

import h5py
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units

def extract_data(file_name, apply_consts=True):
    pctocm, kmtocm, msuntog, scale0, scale1  = 1, 1, 1, 1, 1
    if (apply_consts):
        # AREPO uses different units than FLASH, these are the conversions.
        # https://www.illustris-project.org/data/docs/specifications/
        length, mass, velocity, hubble = 3.08567759e+21, 1.989e43, 1.0e5, 0.7
    f = h5py.File(file_name, 'r')
    ds = f['PartType0']
    c = ds['Coordinates'][:]*length*hubble
    d = ds['Density'][:]*mass*(1./hubble**2)*(1./length**3)
    m = ds['Masses'][:]*mass*hubble
    ie = ds['InternalEnergy'][:]*velocity**2
    v = ds['Velocities'][:]*velocity
    gpot = ds['Potential'][:]*velocity**2
    
    coords = np.array([c[:,0], c[:,1], c[:,2]]).T
    vels = np.array([v[:,0], v[:,1], v[:,2]]).T
    f.close()
    return coords, vels, d, m, ie, gpot

def rescale_coords_vels(coords, vels, masses, apply_consts=True, use_com_coords=False):
    pctocm, kmtocm, msuntog, scale0, scale1  = 1, 1, 1, 1, 1
    u_coord = units.pc
    u_vels = units.km/units.s
    if (apply_consts):
        pctocm, kmtocm, msuntog, scale0, scale1 = 3.08567759e+18, 1.0e5, 1.989e33, 0.7e3, 0.7e10
        u_coord = units.cm
        u_vels = units.cm/units.s
    x_cor = (coords[:,0].max()+coords[:,0].min())/2
    y_cor = (coords[:,1].max()+coords[:,1].min())/2
    z_cor = (coords[:,2].max()+coords[:,2].min())/2
    
    parts = Particles(len(vels[:,0]))
    parts.x, parts.y, parts.z = coords[:,0] | u_coord, coords[:,1] | u_coord, coords[:,2] | u_coord
    parts.vx, parts.vy, parts.vz = vels[:,0] | u_vels, vels[:,1] | u_vels, vels[:,2] | u_vels
    parts.mass = masses
    
    if (use_com_coords):
        com_coords = parts.center_of_mass().value_in(u_coord)
        x_cor, y_cor, z_coor = com_coords[0], com_coords[1], com_coords[2]

    com_vels = parts.center_of_mass_velocity().value_in(u_vels)
    vx_cor = com_vels[0]
    vy_cor = com_vels[1]
    vz_cor = com_vels[2]
    
    coords_cor = coords - np.array([x_cor, y_cor, z_cor]).reshape(1,3)
    vels_cor = vels - np.array([vx_cor, vy_cor, vz_cor]).reshape(1,3)
    
    return coords_cor, vels_cor

def write_corrected_file(output_filename, coords, vels, dens, masses, ie, gpot):
    f = h5py.File(output_filename, 'w')
    group = f.create_group('PartType0')
    dset = group.create_dataset('Coordinates', data=coords, dtype='d')
    dset = group.create_dataset('Velocities', data=vels, dtype='d')
    dset = group.create_dataset('Density', data=dens, dtype='d')
    dset = group.create_dataset('Masses', data=masses, dtype='d')
    dset = group.create_dataset('InternalEnergy', data=ie, dtype='d')
    dset = group.create_dataset('Potential', data=gpot, dtype='d')
    f.close()
