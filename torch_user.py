#!/usr/bin/env python
"""
User file for torch star formation code.
You must define the methods:

    user_initial_conditions(state, hydro)
    user_parameters()

User parameters should have AMUSE units attached, where appropriate.

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
#    star.position = [flashp['xmax'] - 1e10, 0, 0] | units.cm
#    star.velocity = [1e5, 0, 0] | units.cm/units.s
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

#    star        = Particles(2)
#    star.mass   = 1. | units.MSun
#    star.x      = 0.0 | units.cm
#    star.y      = 0.0 | units.cm
#    star.z      = 0.0 | units.cm
#    star.vx     = 0.0 | units.cm/units.s
#    star.vy     = 0.0 | units.cm/units.s
#    star.vz     = 0.0 | units.cm/units.s
#
#    # make bound binary
#    star[0].x = -1.5e16 | units.cm  # 1000 AU away
#    star[1].vy = 1.0e4 | units.cm/units.s  # sqrt(GM/R) = 9.42e4 cm/s
#
#    # place system on exit trajectory
#    flashp = FlashPar("flash.par")
#    star.x = star.x + ((flashp['xmax'] - 1e17) | units.cm)
#    star.vx = (5.0e17/1.0e12) | units.cm/units.s  # 6e4 cm/s
#
#    creation_time = hydro.get_time()  # comes with AMUSE units
#
#    tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(tag, star.mass)
#    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(tag, star.mass) # for SE code
#    hydro.set_particle_creation_time(tag, creation_time)

    # ------------------------------------------------------------------------
    # Multiples test: plop a binary system that exits domain SLOWLY,
    # such that one star lies outside bndbox while COM lies inside bndbox.
    # This is a bit unphysical because the star's position while tracked by
    # Multiples is not well defined.

#    star        = Particles(2)
#    star.mass   = 1. | units.MSun
#    star.x      = 0.0 | units.cm
#    star.y      = 0.0 | units.cm
#    star.z      = 0.0 | units.cm
#    star.vx     = 0.0 | units.cm/units.s
#    star.vy     = 0.0 | units.cm/units.s
#    star.vz     = 0.0 | units.cm/units.s
#
#    # make bound binary with stars initially along x-axis.
#    star[0].x = -1.5e16 | units.cm  # 1000 AU away
#    star[0].vy = -1.0e3 | units.cm/units.s  # sqrt(GM/R) = 9.42e4 cm/s
#    star[1].vy = 1.0e3 | units.cm/units.s  # balance so that COM vy~0
#
#    # first bridge dt = 4.5e12 sec, dx = 4.5e14cm.  Nominally,
#    # system moves 1/3rd of binary sep, so COM in domain, one star past xmax.
#    # next bridge dt, COM exits domain.
#    flashp = FlashPar("flash.par")
#    star.x = star.x + ((flashp['xmax'] - 1e12) | units.cm)
#    star.y = star.y + ((flashp['ymax'] - 1e15) | units.cm)  # be careful about float precision
#    star.z = star.z + ((flashp['zmax'] - 1e15) | units.cm)
#    star.vx = star.vx + (1e3 | units.cm/units.s)  # COM vx=1e2 cm/s too small, gets perturbed by grav
#
#    # The stars fall towards each other and perturb starting positions before
#    # ph4 throws stopping condition, so the setup is inexact.
#    #
#    # Nevertheless, with current (2019 nov 02) example torch settings, after
#    # two bridge steps the binary straddles domain boundary with COM inside and
#    # one star outside.
#    #
#    # This triggers crash in hydro.particles_sort() with old star-removal
#    # algorithm that only looks for ph4 particles outside domain.
#    # New algorithm, which removes entire tree if any leaf outside, works.
#
#    creation_time = hydro.get_time()  # comes with AMUSE units
#
#    tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(tag, star.mass)
#    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(tag, star.mass) # for SE code
#    hydro.set_particle_creation_time(tag, creation_time)

    # ------------------------------------------------------------------------
    # Start with a cluster.  BEWARE: properties are not very carefully chosen,
    # e.g., current initialization scheme probably adds subtle biases etc...

    # The fractal cluster model requires an extra worker, a bit wasteful.
    # It calls "stop()" and should probably release/kill the worker...
    # If we re-initialize grav, delay grav init to after IC setup, or only load
    # cluster from file instead of worker, maybe we can save a process.
    # -AT, 2019 Nov 25
#    from amuse.community.fractalcluster.interface import new_fractal_cluster_model
#    from amuse.ic.brokenimf import new_kroupa_mass_distribution
#    from amuse.ic.plummer import new_plummer_sphere
#    from amuse.ic.salpeter import new_salpeter_mass_distribution
#    from amuse.io import write_set_to_file, read_set_from_file
#    from amuse.units import nbody_system
#    import numpy as np
#
#    def make_cluster(converter, nm_part, bndbox, fractal=False):
#        stars_out = True
#        n = 0
#        while stars_out:
#            if fractal:
#                cluster = new_kroupa_mass_distribution(nm_part, mass_max=(150.0|units.MSun))
#                cluster = new_fractal_cluster_model(masses=cluster, convert_nbody=converter, do_scale=False, virial_ratio=1.0)
#            else:
#                cluster = new_plummer_sphere(nm_part, convert_nbody=converter, do_scale=False)
#                cluster.mass = new_kroupa_mass_distribution(nm_part, mass_min=(0.08|units.MSun), mass_max=(150.0|units.MSun))
#            remove_stars = cluster.select(lambda r: bndbox < max(abs(r)), ["position"])
#            stars_out = len(remove_stars) > 0
#            n += 1
#        print("Made cluster in", n, "attempts.")
#        return cluster
#
#    def make_cluster_in_hydro(cluster, bndbox):
#
#        tag = hydro.add_particles(cluster.x, cluster.y, cluster.z)
#        hydro.set_particle_velocity(tag, cluster.vx, cluster.vy, cluster.vz)
#        hydro.set_particle_mass(tag, cluster.mass)
#        hydro.set_particle_oldmass(tag, cluster.mass)  # for SE code
#        hydro.set_particle_creation_time(tag, hydro.get_time())
#
#        return tag
#
#    flashp = FlashPar("flash.par")
#    xmax = flashp['xmax'] | units.cm
#
#    # create new cluster from scratch...
#    conv_cluster = nbody_system.nbody_to_si(3.0|units.parsec, 300.0|units.MSun)
#    cluster = make_cluster(conv_cluster, 100, xmax, fractal=True)
#    #write_set_to_file(cluster, 'starting_cluster.hdf5', 'hdf5')
#
#    # load cluster from file...
#    #cluster = read_set_from_file('starting_cluster.hdf5', 'hdf5')
#
#    make_cluster_in_hydro(cluster, xmax)

    # ------------------------------------------------------------------------

    return

def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """

    p = {}
    flashp = FlashPar("flash.par")

    # <bridge>

    p['npy_seed'] = None  # no effect if (restart && restart_with_new_rng=False)
    p['restart_with_new_rng'] = False
    p['restart_with_user_ics'] = False  # meant for testing

    p['evolve_async'] = True
    p['with_bridge'] = True
    p['with_multiples'] = True  # adds two workers: kepler, smalln
    p['with_se'] = True

    # <timestepping>

    p['hy_dt_factor'] = 0.99999  # pin bridge timestep to <= hy_dt_factor*(hydro timestep)

    # <star/n-body gravity>

    p['with_ph4'] = True  # use ph4 or Hermite
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
    #p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 2  # if using fractal cluster IC, need extra worker

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
