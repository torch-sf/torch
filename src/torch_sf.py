"""
Torch code to do star formation

Currently, just implements scheme for creating stars from sinks.
Could be made more general in the future.

Joshua Wall, Drexel University
"""

from __future__ import division, print_function

import numpy as np
import pprint

from amuse.datamodel import Particles
from amuse.units import units

from torch_stdout import tprint
from imf_sample import sample_stellar_mass


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
    angMom   = hydro.get_particle_ang_mom(newtags) # add angular momentum code -SA 20230301

    # Get SeBa properties from checkpoint - CCC 25/04/2024, 06/11/2024
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
    # Set stellar type and radius - CCC 06/11/2024
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
    add_star.initial_mass = initMass # for SE/SN uses

    #add_star.ang_mom = angMom # add angular mometum code -SA 20230301
    # Removing the line to add ang_mom to the AMUSE particle set 
    # as the amuse particle set doesn't need it and breaks if it's added. -SA 20250226
    
    # only used by ph4... without this, ph4 complains about reused user IDs
    add_star.id = state.stars_next_id + np.arange(num_new_parts)
    state.stars_next_id += num_new_parts

    # Many print statements to check issue with add_particles call.
    # Commented out 20250304 but could be removed entirely soon. -SA
    #tprint("Check add_star before adding particles: id, tag, mass, stellar_type")
    #print(add_star.id, add_star.tag, add_star.mass, add_star.stellar_type)
    #tprint("Print x, y, z, vx, vy, vz")
    #print(add_star.x, add_star.y, add_star.z, add_star.vx, add_star.vy, add_star.vz)
    #tprint("Print radius, initMass")
    #print(add_star.radius, add_star.initial_mass) #, add_star.ang_mom)
    #tprint("Number of new particles to be added: ")
    #print(num_new_parts)
    #print(len(add_star.id))
    #tprint("Done printing properties. ")

    #tprint("Test pprint option.")
    #pprint.pp(add_star)
    #tprint("Finished pprint.")

    state.stars.add_particles(add_star)
    tprint("Finished stars.add_particles()")
    state.stars = state.stars.sorted_by_attribute('tag')
    #tprint("Finished sorted_by_attribute")

    grav.particles.add_particles(add_star)
    tprint("Finished grav.particles.add_particles()")
   
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

    tprint("Completed add_particles_to_grav")
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
    tprint("Entering remove_particles_outside_bndbox")
    p = state.stars
    if len(p) == 0:
        return

    xmin = hydro.get_runtime_parameter('xmin') | units.cm
    xmax = hydro.get_runtime_parameter('xmax') | units.cm
    ymin = hydro.get_runtime_parameter('ymin') | units.cm
    ymax = hydro.get_runtime_parameter('ymax') | units.cm
    zmin = hydro.get_runtime_parameter('zmin') | units.cm
    zmax = hydro.get_runtime_parameter('zmax') | units.cm

    tprint("Acquired runtime_parameters. Now looking for stars outside bndbox.")
    outside = np.logical_or.reduce([
        p.x >= xmax, p.x <= xmin,
        p.y >= ymax, p.y <= ymin,
        p.z >= zmax, p.z <= zmin,
    ])

    tprint("Now removing stars.")
    stars_rem = p[outside]
    grav_rem  = stars_rem.copy()
    se_rem    = stars_rem.copy()

    if len(stars_rem) > 0:

        if mult is None:

            tprint("Removing", len(stars_rem), "star(s) outside bndbox")

        else:

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
    tprint("Completed remove_particles_outside_bndbox.")
        
    return


def queue_stars(state, hydro, min_imf_mass=None, max_imf_mass=None,
                sample_imf_mass=10000|units.MSun, sum_small=False, m_small=1.0|units.MSun,
                sample_imf_bins=10, jet_fraction=0.0,  #Add default value of jet_fraction -SA 20220819
                minimum_jet_mass=100|units.MSun, maximum_jet_mass=0.01|units.MSun): #Add jet mass range -SA 20230728
    """Check hydro for new sinks, queue stars for spawning"""

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()
    if num_sinks == 0:
        hydro.set_particle_pointers('mass')
        return

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))  # does not work with empty list

    # Josh wrote efficient code to update new sinks using cached, sorted list
    # of old sinks, which I (AT) removed for brevity.
    # Simple for-loop should work fine for up to few thousand sinks...
    for sink_tag in sink_tags:

        if sink_tag not in state.all_masses:
            state.all_masses[sink_tag] = np.array([])
            state.starjet_masses[sink_tag] = np.array([]) # Added starjet masses -SA 20220819
            tprint("... new sink tag {}".format(sink_tag))

        while np.sum(state.all_masses[sink_tag]) | units.MSun <= hydro.get_particle_mass(sink_tag):
            new_masses, new_starjet_masses = sample_stellar_mass(
                            sample_imf_mass.value_in(units.MSun),
                            num_bins=sample_imf_bins,
                            min_samp_mass=min_imf_mass.value_in(units.MSun),
                            max_samp_mass=max_imf_mass.value_in(units.MSun),
                            sum_small=sum_small, 
                            m_small=m_small.value_in(units.MSun),
                            jet_fraction=jet_fraction,  # Added jet_fraction -SA 20220819
                            minimum_jet_mass=minimum_jet_mass, 
                            maximum_jet_mass=maximum_jet_mass,  #Added jet mass range -SA 20230728
            )

            tprint("... sink tag {}".format(sink_tag), end='')
            print(" queued {} stars,".format(len(new_masses)), end='')
            print(" mass {},".format(np.sum(new_masses)), end='')
            print(" max mass {}".format(np.amax(new_masses)))

            state.all_masses[sink_tag] = np.concatenate((state.all_masses[sink_tag], new_masses))
            state.starjet_masses[sink_tag] = np.concatenate((state.starjet_masses[sink_tag], new_starjet_masses))
            #  Added starjet state - SA 20220819

    hydro.set_particle_pointers('mass')

    return


def make_stars_from_sinks(state, hydro, sink_rad=None):
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

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))  # does not work with empty list
    sink_tags.sort()  # is this necessary?

    for sink_tag in sink_tags:

        hydro.set_particle_pointers('sink')
        sink_mass = hydro.get_particle_mass(sink_tag)
        sink_pos = hydro.get_particle_position(sink_tag)
        sink_vel = hydro.get_particle_velocity(sink_tag)
        sink_cs  = hydro.get_sink_mean_cs(sink_tag)
        sink_angMom = hydro.get_particle_ang_mom(sink_tag)  # Added -SA 20230405

        # Add quick print statement to double check sink accretion method for jets -SA 20221012
        tprint("Sink mass before forming stars for sink tag {}".format(sink_tag), "mass is {}".format(sink_mass))  
        tprint("Sink angular momentum before forming stars for sink tag {}".format(sink_tag), "ang momentum is {}".format(sink_angMom))

        # get all the stars that we can form now
        csum = np.cumsum(state.all_masses[sink_tag])
        i = np.searchsorted(csum, sink_mass.value_in(units.MSun), side='left')
        assert i < len(csum)  # ensure csum[-1] = sum(queue) > sink_mass

        spawn_masses = state.all_masses[sink_tag][:i]
        spawn_starjet = state.starjet_masses[sink_tag][:i] # Added starjet masses - SA 20220819
        nnew = len(spawn_masses)

        if nnew == 0:

            tprint("... sink tag {} did not spawn stars".format(sink_tag))

        elif np.isnan(sink_cs.value_in(units.cm/units.s)):

            tprint("... sink tag {} blocked from spawning".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total starjet masses {},".format(np.sum(spawn_starjet)), end='') # Added print statement -SA 20220819
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" due to absence of nearby cold gas")

        else:

            tprint("... sink tag {} spawned".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total starjet masses {},".format(np.sum(spawn_starjet)), end='') # Added print statement -SA 20220819
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" max mass {:.2f}".format(np.amax(spawn_masses)))

            formed_stars = True

            # Remove newly-created stars from sink's queue
            state.all_masses[sink_tag] = state.all_masses[sink_tag][nnew:]
            state.starjet_masses[sink_tag] = state.starjet_masses[sink_tag][nnew:]  # Added -SA 20220819

            # Remove the mass from the sink and jet, as well as the corresponding angular momentum fraction. -SA updated 20230712
            remaining_mass_frac = (sink_mass - (np.sum(spawn_starjet)|units.MSun)) /sink_mass #fraction of mass that remains (for ang_mom decrease) -SA 20230405
            sink_mass = sink_mass - (np.sum(spawn_starjet)|units.MSun) # Update the sink mass -SA 20230405 (updated to remove spawn_starjet instead of span_masses - SA 20230712)
            #print("sink angular momentum before ang_mom removal: ", sink_angMom) # SA 20230407
            
            hydro.set_particle_mass(sink_tag, sink_mass)
            #Remove ang_mom from sink, based on mass of stars formed, for each direction -SA 20230405
            sink_angMom_x = (sink_angMom[0]*remaining_mass_frac).as_quantity_in(units.cm**2.0 * units.g / units.s)
            sink_angMom_y = (sink_angMom[1]*remaining_mass_frac).as_quantity_in(units.cm**2.0 * units.g / units.s) 
            sink_angMom_z = (sink_angMom[2]*remaining_mass_frac).as_quantity_in(units.cm**2.0 * units.g / units.s)
            tprint("Sink ang momentum x component before updating: ", sink_angMom_x)
            tprint("Sink ang momentum all before updating: ", sink_angMom)
            hydro.set_particle_ang_mom(sink_tag, sink_angMom_x, sink_angMom_y, sink_angMom_z, 1 ) #Assign new sink ang mom


            star          = Particles(nnew)
            star.mass     = spawn_masses | units.MSun # Leave this as spawn_masses to form correct stellar mass - SA 20220819
            # Isothermal spherical distribution.
            star.position = sink_pos + sink_rad*np.random.rand(nnew,1)*random_three_vector(nnew)
            # Gaussian distribution satisfying <vx**2> = sink_cs**2
            # so that stars' specific energy 1/2 <v**2> = (3/2)*sink_cs**2
            # matches gas specific energy P/rho/(gamma-1) for gamma=5/3
            # with cs = sqrt(P/rho) from Particles_sinkCreateAccrete.F90
            star.velocity = sink_vel + (np.random.normal(scale=sink_cs.value_in(units.cm/units.s), size=(nnew,3)) | units.cm/units.s)

            # Calculate star angular momentum -SA 20230103
            # First, we need to check whether the sink angular momentum is 0, which is rare but possible. -SA 20241121
            if all(j==(0.0 | units.cm**2.0*units.g/units.s) for j in sink_angMom): #counts to see if all directions are zero
                tprint("Sink angular momentum is zero: ", sink_angMom, "So stars ang_mom will be random")
                star_angMom = random_three_vector(nnew) | units.cm**2.0*units.g/units.s #use uniform distribution defined below
                tprint("Star ang momentum given sink_angMom=0, before norm: ", star_angMom)
                # Normalization can proceed the same since it's based on the individual star's total ang. mom.
            else: # Assume sink ang_mom is nonzero in at least one direction
                ## Add a small, random variation to the direction of the ang_mom for each star. -SA 20230405
                size_dir_vary = 0.1 # standard deviation of random fluctuations to star ang_mom orientation, currently chosen arbitrarily
                star_angMom_vary = np.random.normal(loc = 1.0, scale = size_dir_vary, size = (nnew,3))
                star_angMom = (sink_angMom * star_angMom_vary) #3d array of momentum vector for each star - does it have units?
                tprint("Star ang momentum after vary, before norm: ", star_angMom)
            
            ## We only need the unit vector of the sink's ang_mom since we want to set the direction of the star's ang. mom. without
            ## handling the full ang_mom conservation. -SA 20230405
            star_angMom_mag = [(np.sqrt(star_angMom[i,0]**2 + star_angMom[i,1]**2 + star_angMom[i,2]**2)) for i in range(nnew)]
            tprint("Star ang momentum magnitude for normalizing: ", star_angMom_mag)
            ## Now set the star angular momentum
            # NOTE the star angular momentum is normalized to a magnitude of 1 and therefore a dimensionless quantity
            # however, the set_particle_ang_mom expects units so the star.ang_mom currently keeps the units. -SA
            star.ang_mom = [ (star_angMom[i]/star_angMom_mag[i].value_in(units.cm**2.0 * units.g / units.s)) for i in range(nnew)]
            tprint("Star ang momentum after norm, with units? : ", star.ang_mom)

            # Create new stars in FLASH
            hydro.set_particle_pointers('mass')
            star_tag = hydro.add_particles(star.x, star.y, star.z)
            hydro.set_particle_mass(star_tag, star.mass)
            hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
            hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
            #print("Setting star angular momentum now - e.g., ", star.ang_mom[0])
            hydro.set_particle_ang_mom(star_tag, star.ang_mom[:,0], star.ang_mom[:,1], star.ang_mom[:,2], nnew )  # Add angular momentum -SA 20230425

        # Add quick print statement (2 of 2) to double check sink accretion method for jets -SA 20221012
        tprint("Sink mass after checking whether to form stars for sink tag {}".format(sink_tag), "mass is {}".format(sink_mass))

    # Indicate we're done spawning stars
    tprint("Finished checking whether to spawn stars from sink particles.")

    # if we made no stars, need to reset pointers
    hydro.set_particle_pointers('mass')

    return formed_stars


def random_three_vector(n=1):
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """
    three_vector = np.zeros((n,3))

    phi = np.random.uniform(0,np.pi*2,n)
    costheta = np.random.uniform(-1,1,n)

    theta = np.arccos( costheta )
    three_vector[:,0] = np.sin( theta) * np.cos( phi )
    three_vector[:,1] = np.sin( theta) * np.sin( phi )
    three_vector[:,2] = np.cos( theta )
    return three_vector


if __name__ == '__main__':
    pass
