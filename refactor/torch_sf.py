"""
Torch code to do star formation

Currently, just implements scheme for creating stars from sinks.
Could be made more general in the future.

Joshua Wall, Drexel University
"""

from __future__ import division

import numpy as np

from imf_sample import sample_stellar_mass


def make_stars_from_sinks2(state, hydro, min_imf_mass, max_imf_mass,
                           sample_imf_mass=10000 | units.MSun,
                           local_sfe=1.0, sum_small=False):
    """
    Given an initial sampling of the IMF, distribute the stars randomly
    as sinks accrete the required mass to form them.

    Post-condition:
    * INSERT new sink, sampled stars INTO state.all_masses
    * new star prtl(s) created in hydro
    * sink mass updated in hydro

    Return: True if formed stars, False otherwise
    """

    formed_stars = False

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()

    if num_sinks < 1:
        hydro.set_particle_pointers('mass')
        return formed_stars

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))
    sink_tags.sort()

    # Each new sink needs a list of star masses.
    #
    # Josh wrote efficient code to update new sinks using cached, sorted list
    # of old sinks, which I (AT) removed for brevity.
    # Simple for-loop should work fine for up to few thousand sinks...
    for sink_tag in sink_tags:
        if sink_tag not in state.all_masses:
            new_masses = sample_stellar_mass(
                            sample_imf_mass,
                            num_bins=10,
                            min_samp_mass=min_imf_mass.value_in(units.MSun),
                            max_samp_mass=max_imf_mass.value_in(units.MSun),
                            eff=local_sfe,
                            sum_small=sum_small
            )
            state.all_masses[sink_tag] = new_masses

    sink_masses = hydro.get_particle_mass(sink_tags).value_in(units.Msun)

    for sink_tag, sink_mass in zip(sink_tags, sink_masses):
        # Does this sink have enough mass to make the next star in queue?
        while sink_mass > state.all_masses[sink_tag][0]:

            sink_pos = hydro.get_particle_position(sink_tag)
            sink_vel = hydro.get_particle_velocity(sink_tag).value_in(units.cm/units.s)
            sink_cs  = hydro.get_sink_mean_cs(sink_tag).value_in(units.cm/units.s)

            # Create new star in AMUSE
            star          = Particles(1)
            star.mass     = state.all_masses[sink_tag][0]  | units.MSun
            star.velocity = np.random.uniform(sink_vel, np.ones(3)*sink_cs) | units.cm/units.s
            # Singular isothermal spherical distribution.
            rvec = (random_three_vector(1)[:,:]*(np.random.rand(1))[:,None]*SINK_RAD)
            star.position = np.add(rvec.value_in(units.cm), sink_pos.value_in(units.cm)) | units.cm

            # Make the new star particle in FLASH.
            hydro.set_particle_pointers('mass')
            star_tag = hydro.add_particles(star.x, star.y, star.z)
            hydro.set_particle_mass(star_tag, star.mass)
            hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
            hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
            # Switch back to sinks to continue the loop.
            hydro.set_particle_pointers('sink')

            # Remove mass from the sink.
            sink_mass = sink_mass - new_star_mass
            hydro.set_particle_mass(sink_tag, sink_mass|units.MSun)
            # Remove newly-created star from sink's queue
            state.all_masses[sink_tag] = np.delete(state.all_masses[sink_tag], 0)

            # Tell the main code we made a star.
            formed_stars = True

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
