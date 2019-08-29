"""
Torch initial conditions module

Joshua Wall, Drexel University
"""
#from amuse.units import units
from amuse.lab import *
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

def make_single_star_in_hydro(hydro, x, y, z, mass, initMass=0.0|units.MSun,
                              age=0.0|units.Myr, vx=0.0|units.cm/units.s,
                              vy=0.0|units.cm/units.s,
                              vz=0.0|units.cm/units.s):
    """Place a single star into the hydro simulation"""

    if (initMass.value_in(units.MSun) == 0.0):
        initMass = mass
    creation_time = hydro.get_time() - age
    tag = hydro.add_particles(x,y,z)
    hydro.set_particle_velocity(tag, vx, vy, vz)
    hydro.set_particle_mass(tag, mass)
    hydro.set_particle_oldmass(tag, initMass)
    hydro.set_particle_creation_time(tag, creation_time)

    return tag


def make_cluster_in_hydro(hydro, cluster, xinit=0.0|units.cm,
                          yinit=0.0|units.cm, zinit=0.0|units.cm):
    """
    Place an entire star cluster into the hydro simulation (after it was made
    in AMUSE).
    """

    x = cluster.x + initial_x
    y = cluster.y + initial_y
    z = cluster.z + initial_z

    tag = hydro.add_particles(x,y,z)
    hydro.set_particle_velocity(tag, cluster.vx, cluster.vy, cluster.vz)
    hydro.set_particle_mass(tag, cluster.mass)

    return tag



def make_cluster(convert, nm_part, bndbox, fractal=False, equal_mass=False, eq_mass=1.0 | units.MSun):
    """
    Make a fractal cluster that is contained within a bounding box
    (generally the hydro box size) with an AMUSE particle set.
    """

    stars_out=True

    while (stars_out):

        if (fractal):
            cluster = new_kroupa_mass_distribution(nm_part, mass_max=100.0|units.MSun)
            cluster = new_fractal_cluster_model(masses=cluster, convert_nbody=convert, do_scale=False, virial_ratio=5.0)
        else:
            cluster = new_plummer_sphere(nm_part, convert_nbody=convert, do_scale=False)

        if (equal_mass):
            cluster.mass = eq_mass
        else:
            #cluster.mass = new_kroupa_mass_distribution(nm_part, mass_max = (100.0 |units.MSun))
            cluster.mass = new_salpeter_mass_distribution(nm_part,
                                mass_min=0.1|units.MSun,
                                mass_max=100.0|units.MSun
            )

        remove_stars = cluster.select(lambda r: bndbox < max(abs(r)), ["position"])
        stars_out = len(remove_stars) > 0

    return cluster


if __name__ == '__main__':
    pass
