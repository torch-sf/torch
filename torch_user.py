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
    # Single star  test: plop a single star at the center of the box
    # ------------------------------------------------------------------------

#    flashp = FlashPar("flash.par")
#
#    star          = Particles(1)
#    star.mass     = 1 | units.MSun
#    star.position = [0, 0, 0] | units.cm
#    star.velocity = [0, 0, 0] | units.cm/units.s
#    star.ang_mom  = [[0, 0, 1.0]] | units.cm**2.0 * units.g / units.s  #To set angular momentum of star -SA 20    230718
#
#    star_tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(star_tag, star.mass)
#    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
#    hydro.set_particle_ang_mom(star_tag, star.ang_mom[:,0], star.ang_mom[:,1], star.ang_mom[:,2])  #To set ang momentum of star - SA 20230718

    # NOTE: This only works with a ZAMS star -- other parameters are needed to start 
    # with an evolved star.
    # ------------------------------------------------------------------------

    return

def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """

    p = {}
    flashp = FlashPar("flash.par")

    # <VorAMR>

    
    try:
        p['with_voramr'] = flashp['use_voramr']
    except KeyError:
        p['with_voramr'] = False
    if p['with_voramr']:    
        p['source_file'] = flashp['voramr_source']
        p['convert_file'] = True
        p['use_localRef'] = flashp['use_localRef']
        p['local_ref'] = [flashp['localRef_x'], flashp['localRef_y'], flashp['localRef_z'], flashp['localRef_r']]
        #None #[3.20621187e+20, 6.24367575e+20, -1.51873194e+20, 1.543e+20] # Restrict particles included in input hdf5 file by defining spherical region. None or [center_x, center_y, center_z, radius] (cm)
        p['center_local_ref'] = flashp['center_localRef']
        p['input_file'] = flashp['voramr_input']
        p['pickle_kdtree'] = False
        p['pickle_file_name'] = "kdtree.pickle"
        p['numBlocks'] = 15000 #345
        p['cellsPerBlock'] = 16
        
    # <bridge>

    p['npy_seed'] = 0  # random seed for numpy RNG. no effect if (restart && restart_with_new_rng=False)
    p['restart_with_new_rng'] = False  # refresh numpy random seed upon restart?
    p['restart_with_user_ics'] = False  # meant for testing
    
    p['evolve_async'] = True  # evolve hydro (Flash), N-body workers in parallel? (using AMUSE async requests)
    p['with_bridge'] = True  # use bridge leapfrog to evolve posiions and velocities? Warning: "False" is not well tested / supported
    p['with_multiples'] = True  # adds two workers: kepler, smalln
    p['with_se'] = True  # do stellar evolution for individual stars?
    p['remove_merged'] = True and p['with_se'] # remove merged stars, only available if running with stellar evolution

    # <timestepping>

    p['hy_dt_factor'] = 0.99999  # pin bridge timestep to <= hy_dt_factor*(hydro timestep)

    # <star/n-body gravity>

    p['with_ph4'] = True  # use ph4 or Hermite
    p['epsilon'] = 15.0 | units.RSun  # N-body softening = actual radius of a massive star

    # <star/n-body gravity & binaries>

    p['with_petar'] = True
    p['petar_rout'] = 0.001 | units.pc # outer radius for tree 

    # <stellar evolution>

    p['with_lyc'] = True  # ionizing radiation, via ray-tracing from stars
    p['with_pe_heat'] = True  # photoelectric heating from stellar radiation (ray-traced); this is SEPARATE from background diffuse photoelectric heating
    p['sigd'] = flashp['sigDust'] # Cross section of dust per hydrogen nulcei
    p['with_sn'] = True  # allow stars to deposit SNe at end of life
    p['with_winds'] = True  # allow stars to deposit hot winds. NOTE: if winds are off and the radiation pressure on, timesteps won't be limited enough for velocities from radiation pressure and may cause unphysically high velocities -BP 25Jan23
    p['massloss_method'] = 'seba' #'puls' 

    p['min_sn_mass'] = 7.0 | units.MSun #Minumum mass for injecting supernovae  -SA 20231007
    p['min_rad_mass'] = 7.0 | units.MSun #Minimum mass for radiation feedback (ionizing and heating) -SA 20231007
    p['minimum_wind_mass'] = flashp['min_wind_mass'] | units.g #Minimum mass for injecting winds -SA 20231007
    
    
    # <set Jets parameters >  - SA 20230808
    # If jets params are in flash.par, we use those values. If not, we assume jets are off and set defaults accordingly.

    if 'min_jet_mass' in flashp:
        p['minimum_jet_mass'] = flashp['min_jet_mass'] | units.g #1.0 | units.MSun  #Minimum mass for producing protostellar jets -SA 20230718
    else: 
        p['minimum_jet_mass'] = 100 | units.MSun   #set defaults to no jets
    if 'max_jet_mass' in flashp:
        p['maximum_jet_mass'] = flashp['max_jet_mass'] | units.g  #7.0 | units.MSun  #Stars at masses equal to or greater than this mass won't produce jets -SA 20230718
    else:
        p['maximum_jet_mass'] = 0.01 | units.MSun #set defaults to no jets (without using 0 in case that causes issues)
    # To ensure a single star only produces either jets OR winds, make sure 'min_wind_mass' and 'max_jet_mass' are equal.
    # However, if a stars mass allows both jets and winds, the jets will be produced at the beginning of the stars life (for the length
    # of the 'jet_lifetime' set below) and then will produce winds.
    # To produce winds but never produce jets, set 'min_jet_mass' to be greater than 'max_jet_mass'. - SA 20230718

    if 'jet_mass_fraction' in flashp:
        p['jet_fraction'] = flashp['jet_mass_fraction']  #0.33  # Set default to 0.0 if no jets
    else:
        p['jet_fraction'] = 0.0  #set default to no jets
    if 'jet_time' in flashp:
        p['jet_lifetime'] = flashp['jet_time'] | units.s  #1e5 | units.yr
    else:
        p['jet_lifetime'] = 0 | units.yr  #set default to no jets
    if 'jet_vel_fraction' in flashp:
        p['jet_vel_frac'] = flashp['jet_vel_fraction'] #1
    else:
         p['jet_vel_frac'] = 1 

    # <star particle creation>

    p['min_imf_mass'] = 0.08 | units.MSun
    p['max_imf_mass'] = 100.0 | units.MSun
    p['sample_imf_mass'] = 10000.0 | units.MSun
    p['sample_imf_bins'] = 100 # Number of log-space bins from which we Poisson sample the Kroupa IMF. Value of 10 was used for Wall+19 and Wall+20. Value of 100 used in Cournoyer-Cloutier+21. https://groups.google.com/g/torch-users/c/BB4qsaxJoig
    p['sink_rad'] = flashp['sink_accretion_radius'] | units.cm
    p['sum_small'] = False # agglomerate low-mass stars into particles with mass >= m_small Msun?
    p['m_small'] = 1.0 | units.MSun # agglomerate mass in Msun

    # <amuse file overwrite>

    p['overwrite'] = True # <True> Passes flag to AMUSE write_set_to_file(); allows .amuse files to be overwritten without warning.

    # <job>

    ntasks = get_ntasks_from_run_script("run.sh")

    p['num_grav_workers'] = 1 # must be power of 2 for PeTar 
    p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 1  # amuse
    #p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 2  # if using fractal cluster IC, need extra worker

    if p['with_petar']:
        p['with_ph4'] = False
        p['with_multiples'] = False

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
