#!/usr/bin/env python
"""
User file for torch star formation code.
You must define the methods:

    user_initial_conditions(state, hydro)
    user_parameters()

User parameters should have AMUSE units attached, where appropriate.

See the main torch code (torch.py) to understand how this all works.

Design inspired by TRISTAN-MP, Athena++ architecture.
"""

from __future__ import division, print_function

from amuse.datamodel import Particles
from amuse.units import units

from torch_param import FlashPar
from torch_mainloop import run_torch

def get_ntasks_from_run_script(name="run.sh"):
    """formally -n is --ntasks, de facto same as nprocs"""
    n = None
    with open(name) as f:
        for line in f:
            w = line.split()
            if len(w) >= 3 and w[0] == '#SBATCH' and w[1] == '-n':
                assert n is None  # throw error if #SBATCH -n occurs >1x
                n = int(w[2])
    assert n is not None
    return n

def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.
    Usually, this means adding star particles to the hydro code.

    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """

    # ------------------------------------------------------------------------
    # Star removal test: plop a single star into FLASH that exits domain
    # quickly.

#    flashp = FlashPar("flash.par")
#
#    star          = Particles(1)
#    star.mass     = 1 | units.MSun
#    star.position = [flashp['xmax'] - 200, 0, 0] | units.cm
#    star.velocity = [100, 0, 0] | units.cm/units.s
#
#    star_tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(star_tag, star.mass)
#    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.

    # ------------------------------------------------------------------------
    # Multiples test: plop a binary system

#    star        = Particles(2)
#    star.mass   = 1. | units.MSun
#    star.x      = 0.0 | units.cm
#    star.y      = 0.0 | units.cm
#    star.z      = 0.0 | units.cm
#    star.vx     = 0.0 | units.cm/units.s
#    star.vy     = 0.0 | units.cm/units.s
#    star.vz     = 0.0 | units.cm/units.s
#
#    star[0].x = 1.5e16 | units.cm  # 1000 AU away
#    star[1].vy = 1.0e4 | units.cm/units.s  # sqrt(GM/R) = 9.42e4 cm/s ...
#
#    creation_time = hydro.get_time()  # comes with AMUSE units
#
#    tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(tag, star.mass)
#    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(tag, star.mass) # for SE code
#    hydro.set_particle_creation_time(tag, creation_time)

    # ------------------------------------------------------------------------
    # Multiples test: plop a binary system that exits domain QUICKLY

    star        = Particles(2)
    star.mass   = 1. | units.MSun
    star.x      = 0.0 | units.cm
    star.y      = 0.0 | units.cm
    star.z      = 0.0 | units.cm
    star.vx     = 0.0 | units.cm/units.s
    star.vy     = 0.0 | units.cm/units.s
    star.vz     = 0.0 | units.cm/units.s

    # make bound binary
    star[0].x = -1.5e16 | units.cm  # 1000 AU away
    star[1].vy = 1.0e4 | units.cm/units.s  # sqrt(GM/R) = 9.42e4 cm/s

    # place system on exit trajectory
    flashp = FlashPar("flash.par")
    star.x = star.x + ((flashp['xmax'] - 1e17) | units.cm)
    star.vx = (5.0e17/1.0e12) | units.cm/units.s  # 6e4 cm/s

    creation_time = hydro.get_time()  # comes with AMUSE units

    tag = hydro.add_particles(star.x, star.y, star.z)
    hydro.set_particle_mass(tag, star.mass)
    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
    hydro.set_particle_oldmass(tag, star.mass) # for SE code
    hydro.set_particle_creation_time(tag, creation_time)

    # ------------------------------------------------------------------------
    # Multiples test: plop a binary system that exits domain SLOWLY,
    # such that one star lies outside bndbox while COM lies inside bndbox

    return

def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """

    p = {}
    flashp = FlashPar("flash.par")

    # <bridge>

    #p['npy_seed'] = None
    #p['refresh_with_new_rng'] = False
    #p['restart_with_ics'] = False  # meant for testing
    p['npy_seed'] = 103180  # no effect if (restart && refresh_rng=False)
    p['restart_with_new_rng'] = True  # TODO return to False 2019nov01
    p['restart_with_user_ics'] = True  # meant for testing  # TODO return to False 2019nov01
    p['with_bridge'] = True
    p['with_multiples'] = True  # adds two workers: kepler, smalln
    p['with_se'] = True

    # <timestepping>

    p['hy_dt_factor'] = 1.5  # pin bridge timestep to <= hy_dt_factor*(hydro timestep)

    # <star/n-body gravity>

    p['with_ph4'] = False  # use ph4 or Hermite
    p['epsilon'] = 15.0 | units.RSun  # N-body softening = actual radius of a massive star

    # <stellar evolution>

    p['with_lyc'] = True
    p['with_pe_heat'] = True
    p['with_sn'] = True
    p['with_winds'] = True
    p['massloss_method'] = 'puls'
    p['min_feedback_mass'] = 7.0 | units.MSun

    # <star particle creation>

    p['min_imf_mass'] = 0.08 | units.MSun
    p['max_imf_mass'] = 150.0 | units.MSun
    p['sample_imf_mass'] = 10000.0 | units.MSun
    p['sample_imf_bins'] = 10
    p['sink_rad'] = flashp['sink_accretion_radius'] | units.cm
    p['sum_small'] = False

    # <job>

    ntasks = get_ntasks_from_run_script("submit")

    p['num_grav_workers'] = 1
    p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 1  # amuse

    if p['with_se']:
        p['num_hy_workers'] -= 1

    if p['with_multiples']:
        p['num_hy_workers'] -= 2  # SmallN, Kepler

    return p

# ============================================================================

if __name__ == '__main__':
    run_torch(
        user_initial_conditions,
        user_parameters,
    )
