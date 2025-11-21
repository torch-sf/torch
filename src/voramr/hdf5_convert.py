# voramr_convert.py
#
# Contains utility functions to extract data from
# a Voronoi data structure output to manipulate
# and resave as a file of data recognizable by FLASH.
#
# Written by Sean C. Lewis - Drexel University

import h5py
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
from voramr.voramr_stdout import vprint
from torch_param import FlashPar

def extract_data(file_name, apply_consts=True):
    #pctocm, kmtocm, msuntog, scale0, scale1  = 1, 1, 1, 1, 1
    if (apply_consts):
        # AREPO uses different units than FLASH, these are the conversions.
        # https://www.illustris-project.org/data/docs/specifications/
        length, mass, velocity, hubble = 3.08567759e+21, 1.989e43, 1.0e5, 0.7
    f = h5py.File(file_name, 'r')
    ds = f['PartType0']
    c = ds['Coordinates'][:]*length*(1./hubble)
    d = ds['Density'][:]*mass*(hubble**2)*(1./length**3)
    m = ds['Masses'][:]*mass*(1./hubble)
    ie = ds['InternalEnergy'][:]*velocity**2
    v = ds['Velocities'][:]*velocity
    gpot = ds['Potential'][:]*velocity**2
    
    coords = np.array([c[:,0], c[:,1], c[:,2]]).T
    vels = np.array([v[:,0], v[:,1], v[:,2]]).T

    sim_time = f["/Header"].attrs['Time']
    #Extract star dataset
    sds = f['PartType4']
    c = sds['Coordinates'][:]*length*(1./hubble)
    sm = sds['Masses'][:]*mass*(1./hubble)
    smassive = sds['MassiveStarMass'][:] # units in solar masses
    im = sds['GFM_InitialMass'][:]*mass*(1./hubble)
    v = sds['Velocities'][:]*velocity
    a = (sim_time - sds['GFM_StellarFormationTime'][:])*length/velocity/hubble
    smet = sds['GFM_Metallicity'][:]

    scoords = np.array([c[:,0], c[:,1], c[:,2]]).T
    svels = np.array([v[:,0], v[:,1], v[:,2]]).T
    f.close()
    return coords, vels, d, m, ie, gpot, scoords, svels, sm, im, a, smet, smassive

def rescale_coords_vels(coords, vels, masses, scoords, svels, use_com_coords=False):
    #pctocm, kmtocm, msuntog, scale0, scale1  = 1, 1, 1, 1, 1
    #u_coord = units.pc
    #u_vels = units.km/units.s
    #if (apply_consts):
        #pctocm, kmtocm, msuntog, scale0, scale1 = 3.08567759e+18, 1.0e5, 1.989e33, 0.7e3, 0.7e10
    u_coord = units.cm
    u_vels = units.cm/units.s
    # This is a pretty crude method for centering the domain at (0,0,0)cm but it works for now.
    x_cor = (coords[:,0].max()+coords[:,0].min())/2
    y_cor = (coords[:,1].max()+coords[:,1].min())/2
    z_cor = (coords[:,2].max()+coords[:,2].min())/2
    
    parts = Particles(len(vels[:,0]))
    parts.x, parts.y, parts.z = coords[:,0] | u_coord, coords[:,1] | u_coord, coords[:,2] | u_coord
    parts.vx, parts.vy, parts.vz = vels[:,0] | u_vels, vels[:,1] | u_vels, vels[:,2] | u_vels
    parts.mass = masses
    
    if (use_com_coords):
        # Set the center of mass as (0,0,0)cm
        com_coords = parts.center_of_mass().value_in(u_coord)
        x_cor, y_cor, z_coor = com_coords[0], com_coords[1], com_coords[2]

    # Scale velocities such that center of mass of system is (0,0,0)cm/s
    com_vels = parts.center_of_mass_velocity().value_in(u_vels)
    vx_cor = com_vels[0]
    vy_cor = com_vels[1]
    vz_cor = com_vels[2]

    coords_cor = coords - np.array([x_cor, y_cor, z_cor]).reshape(1,3)
    vels_cor = vels - np.array([vx_cor, vy_cor, vz_cor]).reshape(1,3)

    # Same for stars if we have them
    scoords_cor = scoords - np.array([x_cor, y_cor, z_cor]).reshape(1,3)
    svels_cor = svels - np.array([vx_cor, vy_cor, vz_cor]).reshape(1,3)
    return coords_cor, vels_cor, scoords_cor, svels_cor

def write_corrected_file(output_filename, coords, vels, dens, masses, ie, gpot,
                         scoords, svels, smass, sinitmass, sage, smetal,
                         use_localRef=False, local_ref=[0.0,0.0,0.0], recenter_coords=False):
    # Write all gas data to file to be included in interpolation kdtree regardless if we
    # are refining on a region of interest. Include all field values.
    #f = h5py.File("kdtree-"+output_filename, 'w')
    f = h5py.File("interp-data.hdf5", 'w') 
    group = f.create_group('PartType0')
    coords_i = coords.copy()
    if (local_ref and recenter_coords):
        # Only want to limit what we write to interpolation file if we are rescaling coords,
        # since then we are presumably zooming in on ROI and do not need to interp
        # to full computational domain.
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[3]
        diffr = np.sqrt((coords[:,0]-locx)**2 + (coords[:,1]-locy)**2 + (coords[:,2]-locz)**2)
        ind = np.where(diffr < locr)
        # Set particle coord array to be only those within ROI.
        #print(coords*3.24078e-19,coords.shape)
        coords_i = coords.copy()[ind]
        #print(coords*3.24078e-19,coords.shape)
        vels = vels[ind]
        masses = masses[ind]
        vprint("Shifting coordinates and velocities of interpolation file for recentered local refinement")
        # RESCALING SCOORDS AND SVELS FOR ROI NOT IMPLEMENTED YET - SCL 04/15/23
        coords_cor, vels_cor, scoords_cor, svels_cor = rescale_coords_vels(coords_i, vels, masses, scoords, svels, use_com_coords=False)
        coords_i = coords_cor
        #print(coords*3.24078e-19, coords.shape)
        #for c in coords:
        #    print(c)
        #print(a)
        vels = vels_cor
        dens = dens[ind]
        ie = ie[ind]
        gpot = gpot[ind]
    vprint("{} particles saved to interp-data.hdf5".format(len(coords)))
    dset = group.create_dataset('Coordinates', data=coords_i, dtype='d')
    dset = group.create_dataset('Velocities', data=vels, dtype='d')
    dset = group.create_dataset('Density', data=dens, dtype='d')
    dset = group.create_dataset('Masses', data=masses, dtype='d')
    dset = group.create_dataset('InternalEnergy', data=ie, dtype='d')
    dset = group.create_dataset('Potential', data=gpot, dtype='d')
    f.close()
    #vprint("Wrote all gas field values to", "kdtree-"+output_filename)
    vprint("Wrote all gas field values to", "interp-data.hdf5")
    
    #f = h5py.File(output_filename, 'w')
    # Recreate gas dataset
    #group = f.create_group('PartType0')
    
    if(use_localRef):
        vprint("DOING LOCALIZED REFINEMENT. Limiting gas particles written. Opening",output_filename)
        # open file to fill with region-of-interest gas only --> FLASH refinement
        # therefore only need coordinate data, commented out all other field values
        # to reduce file size.
        f = h5py.File(output_filename, 'w')
        group = f.create_group('PartType0')
        
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[3]
        vprint("locx = ", local_ref[0])
        vprint("locy = ", local_ref[1])
        vprint("locz = ", local_ref[2])
        vprint("locr = ", local_ref[3])

        #diffr = np.sqrt((coords[:,0]-locx)**2 + (coords[:,1]-locy)**2 + (coords[:,2]-locz)**2)
        #ind = np.where(diffr < locr)
        # Set particle coord array to be only those within ROI.
        #coords = coords[ind] 
        #vprint("INDICIES < locr:", ind)
        #vprint("coords shape: ", coords.shape)

        # New addition 04/18/23 extracting particles inside cube of side 2/sqrt(2)*locr centered at loc{x,y,z}
        # From https://stackoverflow.com/questions/42352622/finding-points-within-a-bounding-box-with-numpy
        vprint("Creating new mask for particles in bounding cube of side 2*locr/sqrt(2)")
        sqrt2 = np.sqrt(2)
        bound_x = np.logical_and(coords[:, 0] > locx-locr/sqrt2, coords[:, 0] < locx+locr/sqrt2)
        bound_y = np.logical_and(coords[:, 1] > locy-locr/sqrt2, coords[:, 1] < locy+locr/sqrt2)
        bound_z = np.logical_and(coords[:, 2] > locz-locr/sqrt2, coords[:, 2] < locz+locr/sqrt2)
        bb_filter = np.logical_and(np.logical_and(bound_x, bound_y), bound_z)
        coords = coords[bb_filter]
        vprint("len boundx:", len(np.where(bound_x==True)[0]))
        vprint("len boundy:", len(np.where(bound_y==True)[0]))
        vprint("len boundz:", len(np.where(bound_z==True)[0]))
        vprint("len boundx and boundy:", len(np.where(np.logical_and(bound_x, bound_y) == True)[0]))
        vprint("INDICIES < locr:", len(np.where(bb_filter==True)[0]))
        vprint("coords shape: ", coords.shape)

        if (recenter_coords):
            x_cor = (coords[:,0].max()+coords[:,0].min())/2
            y_cor = (coords[:,1].max()+coords[:,1].min())/2
            z_cor = (coords[:,2].max()+coords[:,2].min())/2
            coords = coords - np.array([x_cor, y_cor, z_cor]).reshape(1,3)
        dset = group.create_dataset('Coordinates', data=coords, dtype='d')
        #dset = group.create_dataset('Velocities', data=vels[ind], dtype='d')
        #dset = group.create_dataset('Density', data=dens[ind], dtype='d')
        #dset = group.create_dataset('Masses', data=masses[ind], dtype='d')
        #dset = group.create_dataset('InternalEnergy', data=ie[ind], dtype='d')
        #dset = group.create_dataset('Potential', data=gpot[ind], dtype='d')   
    else:
        vprint("USING ALL GAS PARTICLES, NO LOCAL REFINEMENT.")
        # open file to fill with ALL gas data --> FLASH refinement
        # also would only need coordinate data.
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

    # Removed for now, as we do not port AREPO stars into FLASH,
    # but this in principle could be done. - SCL
    # Recreate stars dataset
    #vprint("Including all stars")
    #group_s = f.create_group('PartType4')
    #dset = group_s.create_dataset('Coordinates', data=scoords, dtype='d')
    #dset = group_s.create_dataset('Velocities', data=scoords, dtype='d')
    #dset = group_s.create_dataset('Masses', data=smass, dtype='d')
    #dset = group_s.create_dataset('GFM_InitialMass', data=sinitmass, dtype='d')
    #dset = group_s.create_dataset('GFM_StellarFormationTime', data=sfmtime, dtype='d')
    #dset = group_s.create_dataset('GFM_Metallicity', data=smetal, dtype='d')

    vprint("Wrote FLASH refinement gas and stars to", output_filename)
    f.close()

def make_background_sinks(scoords, svels, smass, sinitmass, sage, smet, smassive, age_cut=5.0|units.Myr, apply_roi=True):
    """
    Make sink particles in Torch from the already formed young star particles in Arepo sims.
    We only form the star particles that are young and have a massive star particle inside.
    Initializes sinks in the hydro worker using amuse functions. 

    Inputs: star coordinates and velocities in the Torch simulation reference frame (i.e. after
    running rescale_coords_vels), star masses, star initial masses, stellar formation time,
    stellar metallicity, and the amount of mass in massive stars in the stellar particle, age 
    limit of young stars to include as sinks, boolean that determines whether to only include stars
    inside the region of interest defined by the derefinement region.

    Outputs: none. 

    """

    x = scoords[:,0] | units.cm
    y = scoords[:,1] | units.cm
    z = scoords[:,2] | units.cm
    vx = svels[:,0] | units.cm/units.s
    vy = svels[:,1] | units.cm/units.s
    vz = svels[:,2] | units.cm/units.s
    m = smass | units.g
    massive = smassive | units.MSun
    age = sage | units.s 

    sink_filter = [True]*len(x)

    # step one: check that these stars are in the derefinement region as described in the flash.par.
    if apply_roi:
        flashp = FlashPar("flash.par")
        deref_xl = flashp['deref_xl'] | units.cm
        deref_xr = flashp['deref_xr'] | units.cm
        deref_yl = flashp['deref_yl'] | units.cm
        deref_yr = flashp['deref_yr'] | units.cm
        deref_zl = flashp['deref_zl'] | units.cm
        deref_zr = flashp['deref_zr'] | units.cm

        sink_filter = (
            (x > deref_xl) & (x < deref_xr) &
            (y > deref_yl) & (y < deref_yr) &
            (z > deref_zl) & (z < deref_zr)
        )

    # step two: apply age cut and massive star mass cut
    sink_filter = sink_filter & (smassive > 8.0 | units.MSun) & (age<age_cut)

    # step three: make the sinks in flash. 
    hydro.set_particle_pointers('sink')
    sink_tags = hydro.make_sink(sink_tag, x[sink_filter], y[sink_filter], z[sink_filter])
    hydro.set_particle_velocity(sink_tag, vx[sink_filter], vy[sink_filter], vz[sink_filter])
    hydro.set_particle_mass(sink_tag, m[sink_filter])
    hydro.set_particle_extr(sink_tag, [1]*len(m))
    hydro.set_particle_creation_time(sink_tag, -age[sink_filter])
    hydro.set_particle_pointers('mass')

    return None

def write_voramr_data_to_txt_file(voramr_txt_filename, coords, use_localRef=False, local_ref=[0.0,0.0,0.0]):
    # Meant to provide a way to get away from serial hdf5 read in FLASH.
    # This routine writes the file properly, but the FLASH side is still
    # in dev as of Apr 1st, 2023 - SCL
    vprint("~~ Writing input gas coordinate data to text file ~~")
    f = open(voramr_txt_filename, 'w')
    if(use_localRef):
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
