# voramr_convert.py
#
# Contains utility functions to extract data from
# a Voronoi data structure output to manipulate
# and resave as a file of data recognizable by FLASH.

import h5py
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
from voramr.voramr_stdout import vprint

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

    #Extract star dataset
    sds = f['PartType4']
    c = sds['Coordinates'][:]*length*hubble
    sm = sds['Masses'][:]*mass*hubble
    v = sds['Velocities'][:]*velocity
    im = sds['GFM_InitialMass'][:]*mass*hubble
    a = sds['GFM_StellarFormationTime'][:]
    smet = sds['GFM_Metallicity'][:]

    scoords = np.array([c[:,0], c[:,1], c[:,2]]).T
    svels = np.array([v[:,0], v[:,1], v[:,2]]).T
    f.close()
    return coords, vels, d, m, ie, gpot, scoords, svels, sm, im, a, smet

def rescale_coords_vels(coords, vels, masses, scoords, svels, apply_consts=True, use_com_coords=False):
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

    scoords_cor = scoords - np.array([x_cor, y_cor, z_cor]).reshape(1,3)
    svels_cor = svels - np.array([vx_cor, vy_cor, vz_cor]).reshape(1,3)
    return coords_cor, vels_cor, scoords_cor, svels_cor

def write_corrected_file(output_filename, coords, vels, dens, masses, ie, gpot,
                         scoords, svels, smass, sinitmass, sfmtime, smetal, local_ref=None):
    # Write all gas data to file to be included in interpolation kdtree regardless if we
    # are refining on a region of interest.
    f = h5py.File("kdtree-"+output_filename, 'w')
    group = f.create_group('PartType0')
    dset = group.create_dataset('Coordinates', data=coords, dtype='d')
    dset = group.create_dataset('Velocities', data=vels, dtype='d')
    dset = group.create_dataset('Density', data=dens, dtype='d')
    dset = group.create_dataset('Masses', data=masses, dtype='d')
    dset = group.create_dataset('InternalEnergy', data=ie, dtype='d')
    dset = group.create_dataset('Potential', data=gpot, dtype='d')
    f.close()
    vprint("Wrote all gas field values to", "kdtree-"+output_filename)
    
    #f = h5py.File(output_filename, 'w')
    # Recreate gas dataset
    #group = f.create_group('PartType0')
    
    if(local_ref):
        vprint("DOING LOCALIZED REFINEMENT. Limiting gas particles written. Reading ",output_filename)
        # open file to fill with region-of-interest gas only --> FLASH refinement
        f = h5py.File(output_filename, 'w')
        group = f.create_group('PartType0')
        vprint("locx = ", local_ref[0])
        vprint("locy = ", local_ref[1])
        vprint("locz = ", local_ref[2])
        vprint("locr = ", local_ref[3])
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[3]
        diffr = np.sqrt((coords[:,0]-locx)**2 + (coords[:,1]-locy)**2 + (coords[:,2]-locz)**2)
        ind = np.where(diffr < locr)
        vprint("INDICIES < locr:", ind)
        vprint("coords shape: ", coords[ind].shape)
        #vprint("masses shape: ", masses[ind].shape)
        dset = group.create_dataset('Coordinates', data=coords[ind], dtype='d')
        #dset = group.create_dataset('Velocities', data=vels[ind], dtype='d')
        #dset = group.create_dataset('Density', data=dens[ind], dtype='d')
        #dset = group.create_dataset('Masses', data=masses[ind], dtype='d')
        #dset = group.create_dataset('InternalEnergy', data=ie[ind], dtype='d')
        #dset = group.create_dataset('Potential', data=gpot[ind], dtype='d')   
    else:
        vprint("USING ALL GAS PARTICLES, NO LOCAL REFINEMENT.")
        # open file to fill with ALL gas data --> FLASH refinement
        f = h5py.File(output_filename, 'w')
        group = f.create_group('PartType0')
        vprint("coords shape: ", coords.shape)
        vprint("masses shape: ", masses.shape)
        dset = group.create_dataset('Coordinates', data=coords, dtype='d')
        #dset = group.create_dataset('Velocities', data=vels, dtype='d')
        #dset = group.create_dataset('Density', data=dens, dtype='d')
        #dset = group.create_dataset('Masses', data=masses, dtype='d')
        #dset = group.create_dataset('InternalEnergy', data=ie, dtype='d')
        #dset = group.create_dataset('Potential', data=gpot, dtype='d')

    # Recreate stars dataset
    vprint("Including all stars")
    group_s = f.create_group('PartType4')
    dset = group_s.create_dataset('Coordinates', data=scoords, dtype='d')
    dset = group_s.create_dataset('Velocities', data=scoords, dtype='d')
    dset = group_s.create_dataset('Masses', data=smass, dtype='d')
    dset = group_s.create_dataset('GFM_InitialMass', data=sinitmass, dtype='d')
    dset = group_s.create_dataset('GFM_StellarFormationTime', data=sfmtime, dtype='d')
    dset = group_s.create_dataset('GFM_Metallicity', data=smetal, dtype='d')

    vprint("Wrote refinement gas and stars to", output_filename)
    f.close()
def write_voramr_data_to_txt_file(voramr_txt_filename, coords, local_ref=None):
    vprint("~~ Writing input gas coordinate data to text file ~~")
    f = open(voramr_txt_filename, 'w')
    if(local_ref):
        vprint("Doing local refinement")
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[3]
        diffr = np.sqrt((coords[:,0]-locx)**2 + (coords[:,1]-locy)**2 + (coords[:,2]-locz)**2)
        ind = np.where(diffr < locr)
        vprint("INDICIES < locr:", ind)
        vprint("coords shape: ", coords[ind].shape)
        np.savetxt(f, coords[ind], fmt=('%15.7e'))
    else:
        vprint("Not doing local refinement")
        vprint("coords shape: ", coords.shape)
        np.savetxt(f, coords, fmt=('%15.7e'))
    f.close()

    vprint("~~ Done writing to text file ~~")
