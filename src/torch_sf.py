"""
Torch code to do star formation

Currently, just implements scheme for creating stars from sinks.
Could be made more general in the future.

Joshua Wall, Drexel University
"""

from __future__ import division, print_function

import numpy as np

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

    # Get SeBa properties from checkpoint - CCC 25/04/2024
    relMass  = hydro.get_particle_rel_mass(newtags)
    relAge   = hydro.get_particle_rel_age(newtags)
    COcoreM  = hydro.get_particle_co_corem(newtags)
    coreM    = hydro.get_particle_corem(newtags)
    
    # Make AMUSE particles for grav code.
    add_star = Particles(num_new_parts)
    add_star.mass = mass
    add_star.x    = position[:,0]
    add_star.y    = position[:,1]
    add_star.z    = position[:,2]
    add_star.vx   = velocity[:,0]
    add_star.vy   = velocity[:,1]
    add_star.vz   = velocity[:,2]

    # Add saved SeBa properties to AMUSE particles - CCC 25/04/2024
    add_star.relative_mass = relMass
    add_star.relative_age  = relAge
    add_star.COcore_mass   = COcoreM
    add_star.core_mass     = coreM

    add_star.tag  = newtags  # AMUSE stars know their FLASH tags
    add_star.stellar_type = 1 | units.stellar_type # ZAMS star
    add_star.radius = 100 | units.AU # initial collision radius
    add_star.initial_mass = initMass # for SE/SN uses
# don't need to carry this around because we don't need history
# just update directly in hydro
    #if with_lyc:
#    add_star.nion = 0.0 | units.s**-1 # ionizing flux
#    add_star.eion = 0.0 | units.erg # ionizing energy *OVER* 13.6 eV
#    add_star.sigh = 0.0 | units.cm**2 # ionizing cross section.
    #if with_pe_heat:
#    add_star.npe   = 0.0 | units.s**-1 # PE photon flux
#    add_star.epe   = 0.0 | units.erg # PE photon energy (should be around 8 eV)
#    add_star.sigpe = 0.0 | units.cm**2 # dust cross section per hydrogen atom
    #if with_wind:
#    add_star.dm_dt = 0.0 | units.g/units.s
#    add_star.vterm = 0.0 | units.cm/units.s

    # fast-forward stellar evolution to get current stellar type, because
    # torch_sf looks for change in stellar type to decide when to deposit SN
    # Commented out, should not be necesasry if we're restarting with SeBa properties, CCC 10/05/2024
    #if add_parts_restart:
    #    t_evol = hydro.get_time() - hydro.get_particle_creation_time(newtags)
    #    # TODO hardcoded solar metallicity Z=0.02 should be chosen by user.  -AT, 2019oct14
    #    _tmp = se.evolve_star(add_star.initial_mass, t_evol, 0.02)
    #    se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp
    #    add_star.stellar_type = se_type
    #    # Temporary - Set new properties here - CCC, 25/04/2024
    #    add_star.relative_mass = relMass
    #    add_star.relative_age  = relAge
    #    add_star.COcore_mass   = COcoreM
    #    add_star.core_mass     = coreM
        
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


def remove_particles_outside_bndbox(state, hydro, grav, mult):
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
    grav_rem = stars_rem.copy()

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

        hydro.remove_particles(t)
        state.stars.remove_particles(stars_rem)
        grav.particles.remove_particles(grav_rem)
        if mult is None:
            grav.particles.synchronize_to(state.stars)
        else:
            mult._inmemory_particles.remove_particles(grav_rem)
            grav.particles.synchronize_to(mult._inmemory_particles)

    return


def queue_stars(state, hydro, min_imf_mass=None, max_imf_mass=None,
                sample_imf_mass=10000|units.MSun, sum_small=False,
                sample_imf_bins=10):
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
            tprint("... new sink tag {}".format(sink_tag))

        while np.sum(state.all_masses[sink_tag]) | units.MSun <= hydro.get_particle_mass(sink_tag):
            new_masses = sample_stellar_mass(
                            sample_imf_mass.value_in(units.MSun),
                            num_bins=sample_imf_bins,
                            min_samp_mass=min_imf_mass.value_in(units.MSun),
                            max_samp_mass=max_imf_mass.value_in(units.MSun),
                            sum_small=sum_small,
            )

            tprint("... sink tag {}".format(sink_tag), end='')
            print(" queued {} stars,".format(len(new_masses)), end='')
            print(" mass {},".format(np.sum(new_masses)), end='')
            print(" max mass {}".format(np.amax(new_masses)))

            state.all_masses[sink_tag] = np.concatenate((state.all_masses[sink_tag], new_masses))

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

        # get all the stars that we can form now
        csum = np.cumsum(state.all_masses[sink_tag])
        i = np.searchsorted(csum, sink_mass.value_in(units.MSun), side='left')
        assert i < len(csum)  # ensure csum[-1] = sum(queue) > sink_mass

        spawn_masses = state.all_masses[sink_tag][:i]
        nnew = len(spawn_masses)

        if nnew == 0:

            tprint("... sink tag {} did not spawn stars".format(sink_tag))

        elif np.isnan(sink_cs.value_in(units.cm/units.s)):

            tprint("... sink tag {} blocked from spawning".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" due to absence of nearby cold gas")

        else:

            tprint("... sink tag {} spawned".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" max mass {:.2f}".format(np.amax(spawn_masses)))

            formed_stars = True

            # Remove newly-created stars from sink's queue
            state.all_masses[sink_tag] = state.all_masses[sink_tag][nnew:]

            # Remove the mass from the sink.
            sink_mass = sink_mass - (np.sum(spawn_masses)|units.MSun)
            hydro.set_particle_mass(sink_tag, sink_mass)

            star          = Particles(nnew)
            star.mass     = spawn_masses | units.MSun
            # Isothermal spherical distribution.
            star.position = sink_pos + sink_rad*np.random.rand(nnew,1)*random_three_vector(nnew)
            # Gaussian distribution satisfying <vx**2> = sink_cs**2
            # so that stars' specific energy 1/2 <v**2> = (3/2)*sink_cs**2
            # matches gas specific energy P/rho/(gamma-1) for gamma=5/3
            # with cs = sqrt(P/rho) from Particles_sinkCreateAccrete.F90
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
