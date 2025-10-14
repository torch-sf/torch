"""
Torch code to do star formation

Currently, just implements scheme for creating stars from sinks.
Could be made more general in the future.

Joshua Wall, Drexel University

Modified to account for primordial binaries (CCC, 05/2020, 02/2021, 11/2021)
"""



import numpy as np

from amuse.datamodel import Particles
from amuse.units import units

from torch_stdout import tprint
from imf_sample import sample_stars, sample_binaries


def add_particles_to_grav(state, hydro, grav, mult, se):
    """
    Send prtl from hydro to grav + AMUSE

    This gets called in two cases
    1. restarting with prtl,
    2. immediately after making new stars from sinks

    Separating hydro->grav update from sink->amuse->star->hydro update
    allows for possibility of hydro creating its own stars.

    Treat hydro as main fount of knowledge; copy stars from hydro to
    AMUSE and grav.

    se (SeBa or other stellar evolution worker) is only used to get
    correct stellar type for restarts; newborn stars are assumed to be on ZAMS

    postcondition:
        stars updated
        grav updated
    """
    add_parts_restart = False
    num_new_parts = hydro.get_number_of_new_tags()

    if num_new_parts > 0:

        newtags = hydro.get_new_tags(range(1,num_new_parts+1))

    else:
        
        tprint("add_particles_to_grav: assuming restart because Flash reports no new particles!")
        tprint("add_particles_to_grav: sync all stars from Flash to grav.")
        add_parts_restart = True
        num_new_parts = hydro.get_number_of_particles()
        newtags = hydro.get_particle_tags(range(1,num_new_parts+1))

    newtags.sort()

    position = hydro.get_particle_position(newtags)
    velocity = hydro.get_particle_velocity(newtags)
    mass     = hydro.get_particle_mass(newtags)
    initMass = hydro.get_particle_oldmass(newtags)

    # Get SeBa properties from checkpoint
    relMass  = hydro.get_particle_rel_mass(newtags)
    relAge   = hydro.get_particle_rel_age(newtags)
    COcoreM  = hydro.get_particle_co_corem(newtags)
    coreM    = hydro.get_particle_corem(newtags)
    sType    = hydro.get_particle_stype(newtags)
    radius   = hydro.get_particle_radius(newtags)
    
    # Make AMUSE particles for grav code.
    add_star = Particles(num_new_parts)
    add_star.mass = mass
    add_star.x    = position[:,0]
    add_star.y    = position[:,1]
    add_star.z    = position[:,2]
    add_star.vx   = velocity[:,0]
    add_star.vy   = velocity[:,1]
    add_star.vz   = velocity[:,2]

    # Add saved SeBa properties to AMUSE particles - CCC 25/04/2024, 06/11/2024
    add_star.relative_mass = relMass
    add_star.relative_age  = relAge
    add_star.COcore_mass   = COcoreM
    add_star.core_mass     = coreM
    add_star.age           = relAge

    add_star.tag  = newtags  # AMUSE stars know their FLASH tags
    add_star.initial_mass = initMass # for SE/SN uses
    # Set stellar type and radius - CCC 07/08/2024
    # If restart or user ICs, take values from FLASH, otherwise use sensible guess
    add_star.stellar_type = sType | units.stellar_type
    add_star.radius       = radius
    # For new stars
    _new_stars = np.where(sType == 0)[0]
    add_star[_new_stars].stellar_type = 1 | units.stellar_type # ZAMS star
    # Initial guess for the radius if running with user ICs - CCC 12/05/2023
    # It must be somewhat realistic in case there is a contact system
    # Empirical relation from https://articles.adsabs.harvard.edu/pdf/1991Ap%26SS.181..313D
    # Use linear MRR for upper mass range
    # Note that radius now denotes a physical radius and not a collisional radius
    add_star[_new_stars].radius = (1.01 * (add_star[_new_stars].mass / (1 | units.MSun)) ** 0.57) | units.RSun
        
    # only used by ph4... without this, ph4 complains about reused user IDs
    add_star.id = state.stars_next_id + np.arange(num_new_parts)
    state.stars_next_id += num_new_parts

    state.stars.add_particles(add_star)
    state.stars = state.stars.sorted_by_attribute('tag')

    grav.particles.add_particles(add_star)
    
    #Add particles to stellar evolution, CCC 10/05/2024
    se.particles.add_particles(add_star)

    if mult is not None:
        mult._inmemory_particles.add_particles(add_star)
        # Multiples module needs an "id" attribute for internal book-keeping.
        # AMUSE example scripts set "id" directly; we use "index_in_code".
        mult.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    if add_parts_restart:
        hydro.set_starting_local_tag_numbers()

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    return

def remove_particles_outside_bndbox(overwrite, state, hydro, grav, mult, se):
    """
    Remove any particles that have left the simulation.
    WARNING: assumes a box-shaped domain specified by xmin, xmax, etc. in
    FLASH runtime parameters.

    Note: if any star in multiple is outside bndbox, remove the entire multiple
    system, even if rest of system is inside bndbox.

    Arguments:
        state = TorchState(...)
        hydro = FLASH worker code instance
        grav = N-body gravity code instance
        mult = Multiples worker code instance, OR None
    """
    p = state.stars
    if len(p) == 0:
        return

    xmin = hydro.get_runtime_parameter('xmin') | units.cm
    xmax = hydro.get_runtime_parameter('xmax') | units.cm
    ymin = hydro.get_runtime_parameter('ymin') | units.cm
    ymax = hydro.get_runtime_parameter('ymax') | units.cm
    zmin = hydro.get_runtime_parameter('zmin') | units.cm
    zmax = hydro.get_runtime_parameter('zmax') | units.cm

    outside = np.logical_or.reduce([
        p.x >= xmax, p.x <= xmin,
        p.y >= ymax, p.y <= ymin,
        p.z >= zmax, p.z <= zmin,
    ])

    stars_rem = p[outside]
    grav_rem  = stars_rem.copy()
    se_rem    = stars_rem.copy()

    if len(stars_rem) > 0:

        tprint("Removing", len(stars_rem), "star(s) outside bndbox")

        if mult is not None:

            root_rem = Particles(0)

            # remove the entire tree if any leaf outside bndbox
            for root, tree in mult.root_to_tree.items():
                leaves = tree.get_leafs_subset()
                leaves_outside = stars_rem.get_intersecting_subset_in(leaves)
                leaves_inside = leaves - leaves_outside
                if leaves_outside:
                    for leaf in leaves_inside:
                        stars_rem.add_particle(leaf.as_particle_in_set(state.stars))
                    grav_rem.remove_particles(leaves_outside)
                    grav_rem.add_particle(root)
                    root_rem.add_particle(root)

            # stars_rem contains all single stars outside bndbox, and all
            # leaves for multiple-star systems straddling/outside bndbox.
            # grav_rem contains all single stars outside bndbox, and all root
            # particles for multiple-star systems straddling/outside bndbox.
            tprint("Removing", len(root_rem), "multiple system(s) on/outside bndbox")
            tprint("Removing", len(stars_rem), "star(s) on/outside bndbox")

            for root in root_rem:
                del mult.root_to_tree[root]

        # hydro requires sorted tags for removal
        # only the stars particle set has a tag attribute.
        t = stars_rem.tag
        t = np.sort(np.array(t).flatten())
        
        stars_rem.escape_time = hydro.get_time()
        state.out_escaped_stars(stars_rem, overwrite)

        hydro.remove_particles(t)
        state.stars.remove_particles(stars_rem)
        grav.particles.remove_particles(grav_rem)
        se.particles.remove_particles(se_rem)
        if mult is None:
            grav.particles.synchronize_to(state.stars)
        else:
            mult._inmemory_particles.remove_particles(grav_rem)
            grav.particles.synchronize_to(mult._inmemory_particles)
        
    return


def queue_stars(state, hydro, min_imf_mass=None, max_imf_mass=None,
                sample_imf_mass=10000|units.MSun, 
                sum_small=False, m_small=1.0|units.MSun,
                binaries=False, sample_imf_bins=100, mult_frac='field',
                pdist='field', qdist='field', edist='field'):

    """Check hydro for new sinks, queue stars for spawning"""

    new_sink_ = False # CCC 27/04/2023, to save original sink list
    
    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()
    if num_sinks == 0:
        hydro.set_particle_pointers('mass')
        return

    sink_tags = hydro.get_particle_tags(list(range(1,num_sinks+1)))  # does not work with empty list

    # Josh wrote efficient code to update new sinks using cached, sorted list
    # of old sinks, which I (AT) removed for brevity.
    # Simple for-loop should work fine for up to few thousand sinks...
    
    # Only set system_masses, all_positions and all_velocities if binaries=True - CCC 18/06/2025
    if binaries:
        
        for sink_tag in sink_tags:

            if sink_tag not in state.all_masses:
                state.all_masses[sink_tag]     = np.array([])
                state.system_masses[sink_tag]  = np.array([])    # Added by CCC for binaries, 04/04/2020
                state.all_positions[sink_tag]  = np.empty([0,3]) # Added by CCC for binaries, 04/04/2020
                state.all_velocities[sink_tag] = np.empty([0,3]) # Added by CCC for binaries, 04/04/2020
                tprint("... new sink tag {}".format(sink_tag))
                new_sink_ = True # save original sink list

            while np.sum(state.all_masses[sink_tag]) | units.MSun <= hydro.get_particle_mass(sink_tag):
                new_masses, new_system_masses, new_positions, new_velocities = sample_binaries(sample_imf_mass.value_in(units.MSun),
                                                                                               num_bins=sample_imf_bins,
                                                                                               min_samp_mass=min_imf_mass.value_in(units.MSun),
                                                                                               max_samp_mass=max_imf_mass.value_in(units.MSun),
                                                                                               mult_frac=mult_frac, pdist=pdist,
                                                                                               qdist=qdist, edist=edist
                                                                                               )
            
                tprint("... sink tag {}".format(sink_tag), end='')
                print(" queued {} stars,".format(len(new_masses)), end='')
                print(" mass {},".format(np.sum(new_masses)), end='')
                print(" max mass {}".format(np.amax(new_masses)))

                state.all_masses[sink_tag]     = np.concatenate((state.all_masses[sink_tag], new_masses))
                state.system_masses[sink_tag]  = np.concatenate((state.system_masses[sink_tag], new_system_masses))
                state.all_positions[sink_tag]  = np.concatenate((state.all_positions[sink_tag], new_positions))
                state.all_velocities[sink_tag] = np.concatenate((state.all_velocities[sink_tag], new_velocities))
                
    else:
        
        for sink_tag in sink_tags:

            if sink_tag not in state.all_masses:
                state.all_masses[sink_tag]     = np.array([])
                tprint("... new sink tag {}".format(sink_tag))
                new_sink_ = True # save original sink list

            while np.sum(state.all_masses[sink_tag]) | units.MSun <= hydro.get_particle_mass(sink_tag):
                new_masses = sample_stars(sample_imf_mass.value_in(units.MSun),
                                          num_bins=sample_imf_bins,
                                          min_samp_mass=min_imf_mass.value_in(units.MSun),
                                          max_samp_mass=max_imf_mass.value_in(units.MSun),
                                          sum_small=sum_small
                                          )
            
                tprint("... sink tag {}".format(sink_tag), end='')
                print(" queued {} stars,".format(len(new_masses)), end='')
                print(" mass {},".format(np.sum(new_masses)), end='')
                print(" max mass {}".format(np.amax(new_masses)))

                # Only set masses
                state.all_masses[sink_tag]     = np.concatenate((state.all_masses[sink_tag], new_masses))

    hydro.set_particle_pointers('mass')
    
    return new_sink_ # save original sink list


def make_stars_from_sinks(state, hydro, sink_rad=None, binaries=False):
    """
    Given an initial sampling of the IMF, distribute the stars randomly
    as sinks accrete the required mass to form them.

    Post-condition:
    * hydro: new star prtl(s) created, sink mass decremented
    * AMUSE: sink queue updated.  particle set NOT updated.

    Return: True if formed stars, False otherwise
    """
    assert sink_rad is not None  # required kwarg

    formed_stars = False

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()
    if num_sinks == 0:
        # can't get sink tags w/ empty list so need to exit early
        hydro.set_particle_pointers('mass')
        return formed_stars

    sink_tags = hydro.get_particle_tags(list(range(1,num_sinks+1)))  # does not work with empty list
    sink_tags.sort()  # is this necessary?

    for sink_tag in sink_tags:

        hydro.set_particle_pointers('sink')
        sink_mass = hydro.get_particle_mass(sink_tag)
        sink_pos = hydro.get_particle_position(sink_tag)
        sink_vel = hydro.get_particle_velocity(sink_tag)
        sink_cs  = hydro.get_sink_mean_cs(sink_tag)
        
        if binaries:
            # get all the stars that we can form now
            csum = np.cumsum(state.system_masses[sink_tag]) #Changed from all_masses to accomodate binaries
            i = np.searchsorted(csum, sink_mass.value_in(units.MSun), side='left')
            assert i < len(csum)  # ensure csum[-1] = sum(queue) > sink_mass
        else:
            # get all the stars that we can form now
            csum = np.cumsum(state.all_masses[sink_tag])
            i = np.searchsorted(csum, sink_mass.value_in(units.MSun), side='left')
            assert i < len(csum)  # ensure csum[-1] = sum(queue) > sink_mass

        spawn_masses     = state.all_masses[sink_tag][:i]
        if binaries:
            spawn_systems    = state.system_masses[sink_tag][:i]
            spawn_positions  = state.all_positions[sink_tag][:i]
            spawn_velocities = state.all_velocities[sink_tag][:i]
        
        nnew = len(spawn_masses)
        if binaries:
            nbin = nnew - np.count_nonzero(spawn_systems)
            nsin = nnew - 2 * nbin

        if nnew == 0:

            tprint("... sink tag {} did not spawn stars".format(sink_tag))


        else:

            tprint("... sink tag {} spawned".format(sink_tag), end='')
            print(" {} stars".format(nnew), end='')
            if binaries:
                print(" ({} single stars".format(nsin), end='')  #Added by CCC, May 9, 2020 to account for binaries
                print(" and {} binaries),".format(nbin), end='')
            print(" total mass {},".format(np.sum(spawn_masses)), end='')
            print(" max mass {}".format(np.amax(spawn_masses)))
            formed_stars = True

            # Remove newly-created stars from sink's queue
            state.all_masses[sink_tag]     = state.all_masses[sink_tag][nnew:]
            if binaries:
                state.system_masses[sink_tag]  = state.system_masses[sink_tag][nnew:]
                state.all_positions[sink_tag]  = state.all_positions[sink_tag][nnew:]
                state.all_velocities[sink_tag] = state.all_velocities[sink_tag][nnew:]

            # Remove the mass from the sink.
            sink_mass = sink_mass - (np.sum(spawn_masses)|units.MSun)
            hydro.set_particle_mass(sink_tag, sink_mass)

            star          = Particles(nnew)
            star.mass     = spawn_masses | units.MSun

            # For-loop to use the same random position/velocity for stars in a binary
            # COM positions come from an isothermal spherical distribution
            # COM velocities come from a Gaussian distribution satisfying <vx**2> = sink_cs**2
            # so that stars' specific energy 1/2 <v**2> = (3/2)*sink_cs**2                                           
            # matches gas specific energy P/rho/(gamma-1) for gamma=5/3                                              
            # with cs = sqrt(P/rho) from Particles_sinkCreateAccrete.F90

            if binaries:
                for j in range(len(star)):
                    spawn_position = spawn_positions[j]
                    spawn_velocity = spawn_velocities[j]
                    if spawn_systems[j] == 0:
                        star[j].position = sink_pos + random_pos + (spawn_position | units.cm)
                        star[j].velocity = sink_vel + random_vel + (spawn_velocity | units.cm/units.s)
                    else:
                        random_pos = sink_rad*np.random.rand()*random_three_vector()
                        print(random_pos)
                        random_vel = np.random.normal(scale=sink_cs.value_in(units.cm/units.s), size=3) | units.cm/units.s
                        star[j].position = sink_pos + random_pos + (spawn_positions[j] | units.cm)
                        star[j].velocity = sink_vel + random_vel + (spawn_velocity | units.cm/units.s)
            else:
                star.position = sink_pos + sink_rad*np.random.rand(nnew,1)*random_three_vector(nnew)
                if np.isnan(sink_cs.value_in(units.cm/units.s)):
                    # with the maximum value corresponding to gas at T=100K
                    star.velocity = sink_vel + (np.random.normal(117200.0, size=(nnew,3)) | units.cm/units.s)
                else:
                    star.velocity = sink_vel + (np.random.normal(scale=sink_cs.value_in(units.cm/units.s), size=(nnew,3)) | units.cm/units.s)

            # Create new stars in FLASH
            hydro.set_particle_pointers('mass')
            star_tag = hydro.add_particles(star.x, star.y, star.z)
            hydro.set_particle_mass(star_tag, star.mass)
            hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
            hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.

    # if we made no stars, need to reset pointers
    hydro.set_particle_pointers('mass')

    return formed_stars

def random_three_vector(n=1):
    """
    Generates a unit vector by randomly sampling points uniformly
    distributed on the surface of a sphere with radius 1.
    """
    if n <=0:
        raise ValueError('n must be larger than 0')
    if n == 1:
        vec = np.random.normal(size=3)
        return vec/np.sqrt(np.sum(vec**2))
    else:
        vec = np.random.normal(size=[n, 3])
        return vec/np.sqrt(np.sum(vec**2, axis=1))[:,None]

    
if __name__ == '__main__':
    pass
