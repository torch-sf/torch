from amuse.units import units
from amuse.datamodel import Particles, Particle
from amuse.community.seba import Seba
import numpy as np
from matplotlib import pyplot as plt
from amuse.io import write_set_to_file

def evolve_star_from_sf(
    mass = 20 | units.Msun,
    t_end = 5 | units.Myr,
    t_step = 0.1 | units.kyr,
):

    star = Particles(mass = mass)
    stellar = Seba()
    stellar.particles.add_particle(star)

    snapshots = Particles()

    time = 0 | units.Myr
    while time < t_end:
        time += t_step
        stellar.evolve_model(time)
        snapshots.add_particles(stellar.particles)

    stellar.stop()

    return snapshots

evolved_star = evolve_star_from_sf(mass = 20 | units.Msun, t_end = 10 | units.Myr)
write_set_to_file(evolved_star, 'data/seba_amuse_20Msun.amuse', overwrite_file=True)
evolved_star = evolve_star_from_sf(mass = 40 | units.Msun, t_end = 10 | units.Myr)
write_set_to_file(evolved_star, 'data/seba_amuse_40Msun.amuse', overwrite_file=True)
evolved_star = evolve_star_from_sf(mass = 80 | units.Msun)
write_set_to_file(evolved_star, 'data/seba_amuse_80Msun.amuse', overwrite_file=True)