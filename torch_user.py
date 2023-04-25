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
    nodes = None
    cores = None
    import os
    with open(name) as f:
        for line in f:
            w = line.split()
            if len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('--ntasks-per-node'):
                assert cores is None  # throw error if #SBATCH -n occurs >1x
                cores = int(''.join(char for char in w[1] if char.isdigit()))
                cores = int(os.getenv("SLURM_NTASKS_PER_NODE"))
                #print('Requesting ', cores, 'cores on each of')
            elif len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('--nodes'):
                assert nodes is None  # throw error if #SBATCH -n occurs >1x
                nodes =	int(''.join(char for char in w[1] if char.isdigit()))
                nodes = int(os.getenv("SLURM_JOB_NUM_NODES"))
                #print(nodes, 'nodes')
    assert n is None
    n = nodes*cores
    assert n is not None
    return int(os.getenv("SLURM_NTASKS"))
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
    # SN with SE test: plop a star that goes SN within 5e11 seconds

#    star        = Particles(1)
#    star.mass   = 3.09698e+34 | units.g
#    star.x      = 0.0 | units.cm
#    star.y      = 0.0 | units.cm
#    star.z      = 0.0 | units.cm
#    star.vx     = 0.0 | units.cm/units.s
#    star.vy     = 0.0 | units.cm/units.s
#    star.vz     = 0.0 | units.cm/units.s
#
#    oldmass = 5.10964e+34 | units.g  # about 25.5 MSun
#    creation_time = hydro.get_time() - (2.3861e+14|units.s)  # 7.5611 Myr old
#    # goes SN between 7.5611 and 7.5763 Myr (2.3861e14 to 2.3909e14 s)
#
#    tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(tag, star.mass)
#    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(tag, oldmass) # for SE code
#    hydro.set_particle_creation_time(tag, creation_time)

    # ------------------------------------------------------------------------
#    # Multiples test: plop a binary system

#    star        = Particles(2)
#    star.mass   = 1. | units.MSun
#    star.x      = 0.0 | units.cm
#    star.y      = 0.0 | units.cm
#    star.z      = 0.0 | units.cm
#    star.vx     = 0.0 | units.cm/units.s
#    star.vy     = 0.0 | units.cm/units.s
#    star.vz     = 0.0 | units.cm/units.s

#    # make bound binary                                                                                                                                                                                                                                                        
#    star[0].x = -1.5e14 | units.cm  # 10 AU away
#    star[1].vy = 1.0e5 | units.cm/units.s  # sqrt(GM/R) = 9.42e5 cm/s 
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

    p['npy_seed'] = 0  # random seed for numpy RNG. no effect if (restart && restart_with_new_rng=False)
    p['restart_with_new_rng'] = False  # refresh numpy random seed upon restart?
    p['restart_with_user_ics'] = False  # meant for testing
    p['restart_from_stall'] = False # did PeTar stall and exit? Sets r_out = r_bin for first Torch loop
    p['test_binary'] = False # meant for testing
    
    p['evolve_async'] = True  # evolve hydro (Flash), N-body workers in parallel? (using AMUSE async requests)
    p['with_bridge'] = True  # use bridge leapfrog to evolve posiions and velocities? Warning: "False" is not well tested / supported
    p['with_multiples'] = True  # adds two workers: kepler, smalln
    p['with_se'] = True  # do stellar evolution for individual stars?

    # <timestepping>

    p['hy_dt_factor'] = 0.99999  # pin bridge timestep to <= hy_dt_factor*(hydro timestep)

    # <star/n-body gravity>

    p['with_ph4'] = True  # use ph4 or Hermite
    p['epsilon'] = 15.0 | units.RSun  # N-body softening = actual radius of a massive star

    # <star/n-body gravity & binaries>

    p['with_petar'] = True
    p['r_bin'] = 1.496e15 | units.cm # 100AU
    p['set_timeout'] = 300 | units.s # Set timeout stopping condition to 5 minutes, to allow hydro to finish before timeout, CCC 09/03/2023
    
    # <stellar evolution>

    p['with_lyc'] = True  # ionizing radiation, via ray-tracing from stars
    p['with_pe_heat'] = True  # photoelectric heating from stellar radiation (ray-traced); this is SEPARATE from background diffuse photoelectric heating
    p['with_sn'] = True  # allow stars to deposit SNe at end of life
    p['with_winds'] = True  # allow stars to deposit hot winds. NOTE: if winds are off and the radiation pressure on, timesteps won't be limited enough for velocities from radiation pressure and may cause unphysically high velocities -BP 25Jan23
    p['massloss_method'] = 'puls'
    p['min_feedback_mass'] = 8.0 | units.MSun

    # <star particle creation>

    p['binaries'] = True
    #Not used if binaries is false, can leave to default values                                                                
    p['mult_frac'] = 'field'  #Currently accepted method is 'field'. TO DO: Add fraction. 
    p['pdist'] = 'field' #Currently accepted methods are 'field' and 'inner'. TO DO: Add lognormal. 
    p['qdist'] = 'field' #Currently accepted method is 'field'. TO DO: Add random.
    p['edist'] = 'field' #Currently accepted method is 'field'. TO DO: Add thermal.
    p['min_imf_mass'] = 0.08 | units.MSun
    p['max_imf_mass'] = 100.0 | units.MSun
    p['sample_imf_mass'] = 10000.0 | units.MSun
    p['sample_imf_bins'] = 100 # Number of log-space bins from which we Poisson sample the Kroupa IMF. Value of 10 was used for Wall+19 and Wall+20. Value of 100 used in Cournoyer-Cloutier+21. https://groups.google.com/g/torch-users/c/BB4qsaxJoig
    p['sink_rad'] = flashp['sink_accretion_radius'] | units.cm
    p['sum_small'] = False # agglomerate low-mass stars into particles with mass >= m_small Msun?
    p['m_small'] = 1.0 # agglomerate mass in Msun

    # <amuse file overwrite>

    p['overwrite'] = True # <True> Passes flag to AMUSE write_set_to_file(); allows .amuse files to be overwritten without warning.

    # <job>

    ntasks = get_ntasks_from_run_script()


    p['num_grav_workers'] = 2 # must be power of 2 for PeTar 
    p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 1  # amuse
    #p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 2  # if using fractal cluster IC, need extra worker

    if p['with_petar']:
        p['with_ph4'] = False
        p['with_multiples'] = False
        p['epsilon'] = 0 | units.RSun

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
