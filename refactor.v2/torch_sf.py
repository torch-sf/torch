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


def add_particles_to_grav(state, hydro, grav):
    """
    Send prtl from hydro to grav + AMUSE

    This gets called in two cases
    1. restarting with prtl,
    2. immediately after making new stars from sinks

    Separating hydro->grav update from sink->amuse->star->hydro update
    allows for possibility of hydro creating its own stars.

    Treat hydro as main fount of knowledge; copy stars from hydro to
    AMUSE and grav.

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
    #age      = hydro.get_time() - hydro.get_particle_creation_time(newtags)

    # Make AMUSE particles for grav code.
    add_star = Particles(num_new_parts)
    add_star.mass = mass
    #add_star.age  = age
    add_star.x    = position[:,0]
    add_star.y    = position[:,1]
    add_star.z    = position[:,2]
    add_star.vx   = velocity[:,0]
    add_star.vy   = velocity[:,1]
    add_star.vz   = velocity[:,2]
    add_star.tag  = newtags  # AMUSE stars know their FLASH tags
    add_star.stellar_type = 1 | units.stellar_type # ZAMS star
    add_star.radius = 100 | units.AU # initial collision radius
    add_star.initial_mass = initMass # for SE/SN uses
    #if with_lyc:
    add_star.nion = 0.0 | units.s**-1 # ionizing flux
    add_star.eion = 0.0 | units.erg # ionizing energy *OVER* 13.6 eV
    add_star.sigh = 0.0 | units.cm**2 # ionizing cross section.
    #if with_pe_heat:
    add_star.npe   = 0.0 | units.s**-1 # PE photon flux
    add_star.epe   = 0.0 | units.erg # PE photon energy (should be around 8 eV)
    add_star.sigpe = 0.0 | units.cm**2 # dust cross section per hydrogen atom
    #if with_wind:
    add_star.dm_dt = 0.0 | units.g/units.s
    add_star.vterm = 0.0 | units.cm/units.s

    state.stars.add_particles(add_star)
    state.stars = state.stars.sorted_by_attribute('tag')

    grav.particles.add_particles(add_star)

    if add_parts_restart:
        hydro.set_starting_local_tag_numbers()

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    return


def remove_particles_outside_bndbox(state, hydro, grav):
    """
    Remove any particles that have left the simulation.
    WARNING: assumes a box-shaped domain specified by {x,y,z}{min,max}
    """
    # AMUSE stars have tags to FLASH, but grav stars don't have those tags.
    # so more convenient to operate on AMUSE star particles.
    p = state.stars
    if len(p) == 0:
        return False

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

    prem = p[outside]

    if len(prem) > 0:

        tprint("Removing", len(prem), "particles outside bndbox")

        t = prem.tag
        t = np.sort(np.array(t).flatten())  # hydro requires sorted tags for removal

        hydro.remove_particles(t)
        grav.particles.remove_particles(prem)
        state.stars.remove_particles(prem)

        grav.particles.synchronize_to(state.stars)

        return True

    return False


def queue_stars(state, hydro, min_imf_mass=None, max_imf_mass=None,
                sample_imf_mass=10000|units.MSun, sum_small=False,
                sample_imf_bins=10):
    """Check hydro for new sinks, queue stars for spawning"""

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()
    if num_sinks == 0:
        return

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))  # does not work with empty list

    # Josh wrote efficient code to update new sinks using cached, sorted list
    # of old sinks, which I (AT) removed for brevity.
    # Simple for-loop should work fine for up to few thousand sinks...
    for sink_tag in sink_tags:
        if sink_tag not in state.all_masses:
            new_masses = sample_stellar_mass(
                            sample_imf_mass.value_in(units.MSun),
                            num_bins=sample_imf_bins,  # TODO fugly naming and unit handling
                            min_samp_mass=min_imf_mass.value_in(units.MSun),
                            max_samp_mass=max_imf_mass.value_in(units.MSun),
                            sum_small=sum_small,
            )
            state.all_masses[sink_tag] = new_masses

            tprint("... sink tag {} queued".format(sink_tag), end='')
            print(" {} stars,".format(len(new_masses)), end='')  # note mixing tprint(...) and print(...)
            print(" total mass {},".format(np.sum(new_masses)), end='')
            print(" max mass {}".format(np.amax(new_masses)))

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
        spawn_masses = []
        for i in range(len(state.all_masses[sink_tag])):
            m = state.all_masses[sink_tag][i] | units.MSun
            if sink_mass > m:
                sink_mass -= m
                spawn_masses.append(m)

        if spawn_masses:
            tprint("... sink tag {} created {} new stars".format(sink_tag, len(spawn_masses)))
            formed_stars = True
            hydro.set_particle_mass(sink_tag, sink_mass)

            # Remove newly-created stars from sink's queue
            n = len(spawn_masses)
            state.all_masses[sink_tag] = state.all_masses[sink_tag][n:]

            star          = Particles(n)
            star.mass     = spawn_masses | units.MSun
            # Isothermal spherical distribution.
            star.position = sink_pos
            star.position = star.position + ( sink_rad * np.random.rand() * random_three_vector(n) )
            star.velocity = sink_vel
            star.velocity = star.velocity + sink_cs*np.random.uniform(-1,+1,size=(n,3))

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
