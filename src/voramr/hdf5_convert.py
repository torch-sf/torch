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
from imf_sample import sample_stellar_mass
from torch_sf import random_three_vector

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
    xion = 1.0-ds['NeutralHydrogenAbundance'][:]
    
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
    formtime = sds['GFM_StellarFormationTime'][:]
    vprint('sim time and max formtime raw', sim_time, max(formtime))
    a = (sim_time - sds['GFM_StellarFormationTime'][:])*length/velocity/hubble
    vprint('max converted age in seconds', max(a))
    smet = sds['GFM_Metallicity'][:]

    scoords = np.array([c[:,0], c[:,1], c[:,2]]).T
    svels = np.array([v[:,0], v[:,1], v[:,2]]).T
    f.close()
    return coords, vels, d, m, ie, gpot, xion, scoords, svels, sm, im, a, smet, smassive

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

def write_corrected_file(output_filename, coords, vels, dens, masses, ie, gpot, xion,
                         scoords, svels, smass, sinitmass, sage, smetal,
                         use_localRef=False, local_ref=[0.0,0.0,0.0], recenter_coords=False):
    # Write all gas data to file to be included in interpolation kdtree regardless if we
    # are refining on a region of interest. Include all field values.
    #f = h5py.File("kdtree-"+output_filename, 'w')
    f = h5py.File("interp-data.hdf5", 'w') 
    group = f.create_group('PartType0')
    coords_i = coords.copy()
    scoords_i = scoords.copy()
    if (local_ref and recenter_coords):
        # Only want to limit what we write to interpolation file if we are rescaling coords,
        # since then we are presumably zooming in on ROI and do not need to interp
        # to full computational domain.
        locx, locy, locz, locvx, locvy, locvz, locr = local_ref
        loc_pos = [locx, locy, locz]
        loc_vel = [locvx, locvy, locvz]
        diffr = np.sqrt((coords[:,0]-locx)**2 + (coords[:,1]-locy)**2 + (coords[:,2]-locz)**2)
        ind = np.where(diffr < locr)
        # Set particle coord array to be only those within ROI.
        #print(coords*3.24078e-19,coords.shape)
        coords_i = coords.copy()[ind] - loc_pos
        #print(coords*3.24078e-19,coords.shape)
        vels = vels[ind] - loc_vel
        masses = masses[ind]
        scoords = scoords - loc_pos
        svels = svels - loc_vel
        dens = dens[ind]
        ie = ie[ind]
        gpot = gpot[ind]
        xion = xion[ind]
    vprint("{} particles saved to interp-data.hdf5".format(len(coords)))
    dset = group.create_dataset('Coordinates', data=coords_i, dtype='d')
    dset = group.create_dataset('Velocities', data=vels, dtype='d')
    dset = group.create_dataset('Density', data=dens, dtype='d')
    dset = group.create_dataset('Masses', data=masses, dtype='d')
    dset = group.create_dataset('InternalEnergy', data=ie, dtype='d')
    dset = group.create_dataset('Potential', data=gpot, dtype='d')
    dset = group.create_dataset('IonizationFraction', data=xion, dtype='d')
    f.close()
    #vprint("Wrote all gas field values to", "kdtree-"+output_filename)
    vprint("Wrote all gas field values to", "interp-data.hdf5")

    if(use_localRef):
        vprint("DOING LOCALIZED REFINEMENT. Limiting gas particles written. Opening",output_filename)
        # open file to fill with region-of-interest gas only --> FLASH refinement
        # therefore only need coordinate data, commented out all other field values
        # to reduce file size.
        f = h5py.File(output_filename, 'w')
        group = f.create_group('PartType0')
        
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[-1]
        vprint("locx = ", local_ref[0])
        vprint("locy = ", local_ref[1])
        vprint("locz = ", local_ref[2])
        vprint("locr = ", local_ref[3])

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
    else:
        vprint("USING ALL GAS PARTICLES, NO LOCAL REFINEMENT.")
        # open file to fill with ALL gas data --> FLASH refinement
        # also would only need coordinate data.
        f = h5py.File(output_filename, 'w')
        group = f.create_group('PartType0')
        vprint("coords shape: ", coords.shape)
        vprint("masses shape: ", masses.shape)
        dset = group.create_dataset('Coordinates', data=coords, dtype='d')

    vprint("Wrote FLASH refinement gas and stars to", output_filename)
    f.close()

    return scoords, svels

def make_background_sinks(hydro, scoords, svels, smass, sinitmass, sage, smet, smassive, age_cut=5.0|units.Myr, sink_rad=None, apply_roi=True,
                          num_bins=100, min_samp_mass=0.08|units.MSun, max_samp_mass=100.0|units.MSun, sum_small=False, m_small=1.0|units.MSun):
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

    x = scoords[:,0]
    y = scoords[:,1]
    z = scoords[:,2]
    vx = svels[:,0]
    vy = svels[:,1]
    vz = svels[:,2]
    m = smass | units.g
    massive = smassive | units.MSun
    age = sage | units.s 

    roi_filter = x==x

    # step one: check that these stars are in the derefinement region as described in the flash.par.
    if apply_roi:
        flashp = FlashPar("flash.par")
        deref_xl = flashp['deref_xl']
        deref_xr = flashp['deref_xr']
        deref_yl = flashp['deref_yl']
        deref_yr = flashp['deref_yr']
        deref_zl = flashp['deref_zl']
        deref_zr = flashp['deref_zr']
        roi_filter = (
            (x > deref_xl) & (x < deref_xr) &
            (y > deref_yl) & (y < deref_yr) &
            (z > deref_zl) & (z < deref_zr)
        )

    # step two: apply age cut and massive star mass cut
    sink_filter = roi_filter & (massive > 8.0 | units.MSun) & (age<age_cut)

    x = x[sink_filter]
    y = y[sink_filter]
    z = z[sink_filter]
    vx = vx[sink_filter]
    vy = vy[sink_filter]
    vz = vz[sink_filter]
    m = m[sink_filter]
    age = age[sink_filter]

    # step three: make the sinks in flash. 
    hydro.set_particle_pointers('mass')

    for i in range(len(x)):
        # get imf list with arepo star mass
        spawn_masses = sample_stellar_mass(
                            m[i].value_in(units.MSun),
                            num_bins,
                            min_samp_mass=min_samp_mass.value_in(units.MSun),
                            max_samp_mass=max_samp_mass.value_in(units.MSun),
                            sum_small=sum_small,
                            m_small=m_small.value_in(units.MSun))
        nnew = len(spawn_masses)
        if nnew == 0:
            continue
        sink_pos = [x[i], y[i], z[i]] | units.cm
        sink_vel = [vx[i], vy[i], vz[i]] | units.cm/units.s
        spawn_vel = 1.0e5

        star          = Particles(nnew)
        star.mass = spawn_masses | units.MSun
        star.position = sink_pos + (sink_rad.value_in(units.cm)*np.random.rand(nnew,1)*random_three_vector(nnew) | units.cm)
        star.velocity = sink_vel + (np.random.normal(scale=spawn_vel, size=(nnew,3)) | units.cm/units.s)
        star.age = -age[i]
        star_tag = hydro.add_particles(star.x, star.y, star.z)
        hydro.set_particle_mass(star_tag, star.mass)
        hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
        hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
        hydro.set_particle_creation_time(star_tag, star.age)
        vprint('spawning ',nnew, 'stars from arepo')

    return None

def write_voramr_data_to_txt_file(voramr_txt_filename, coords, use_localRef=False, local_ref=[0.0,0.0,0.0]):
    # Meant to provide a way to get away from serial hdf5 read in FLASH.
    # This routine writes the file properly, but the FLASH side is still
    # in dev as of Apr 1st, 2023 - SCL
    vprint("~~ Writing input gas coordinate data to text file ~~")
    f = open(voramr_txt_filename, 'w')
    if(use_localRef):
        vprint("Doing local refinement")
        locx, locy, locz, locr = local_ref[0], local_ref[1], local_ref[2], local_ref[-1]
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
