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

from torch_param import WriteOnceDict, FlashPar
from torch_stdout import tprint

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
    return

def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """

    p = {} # WriteOnceDict()
    flashp = FlashPar("flash.par")

    # <bridge>

    p['npy_seed'] = None
    #p['npy_seed'] = 103180  # no effect if (restart && refresh_rng=False)
    p['refresh_rng'] = False
    p['with_bridge'] = True
    p['with_multiples'] = False  # adds three workers: kepler, smalln, multiples
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
