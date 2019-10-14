#!/usr/bin/env python

from __future__ import division, print_function

from amuse.datamodel import Particles
from amuse.units import units

from torch_param import WriteOnceDict, FlashPar
from torch_stdout import tprint

def user_initial_conditions(state, hydro):

    if state.restart:
        return

    # Plop a single star into FLASH - test star remove algorithm.
#
#    star        = Particles(1)
#    star.mass   = 1 | units.MSun
#    star.x      = BNDBOX - (10 | units.cm)
#    star.y      = 0 | units.cm
#    star.z      = 0 | units.cm
#    star.vx     = 100 | units.cm/units.s
#    star.vy     = 0 | units.cm/units.s
#    star.vz     = 0 | units.cm/units.s
#
#    hydro.set_particle_pointers('mass')
#    star_tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(star_tag, star.mass)
#    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
#
#    # Plop another star into FLASH
#
#    star            = Particles(1)
#    star.mass       = 1 | units.MSun
#    star.position   = [0,0,0] | units.cm
#    star.velocity   = [0,0,0] | units.cm/units.s
#
#    hydro.set_particle_pointers('mass')
#    star_tag = hydro.add_particles(star.x, star.y, star.z)
#    hydro.set_particle_mass(star_tag, star.mass)
#    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
#    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
#
    return

def user_parameters():

    flashp = FlashPar("flash.par")

    p = WriteOnceDict()

    # <main controls>
    p['npy_seed'] = None
    #p['npy_seed'] = 103180  # no effect if (restart && refresh_rng=False)
    p['num_grav_workers'] = 4
    p['num_hy_workers'] = 9
    p['refresh_rng'] = False
    p['with_bridge'] = True
    p['with_multiples'] = False  # adds three workers: kepler, smalln, multiples
    p['with_se'] = True

    # <timestepping>
    p['hy_dt_factor'] = 1.5  # pin bridge timestep to <= hy_dt_factor*(hydro timestep)

    # <star/n-body gravity>
    p['epsilon'] = 15.0 | units.RSun  # N-body softening = actual radius of a massive star

    # <stellar evolution>
    p['with_radiation'] = True  # stellar evolution switches
    p['with_pe_heat'] = True
    p['with_sn'] = True
    p['with_winds'] = True
    p['with_massloss'] = True
    p['massloss_method'] = 'puls'
    p['min_feedback_mass'] = 7.0 | units.MSun

    # <star particle creation>
    p['min_imf_mass'] = 0.08 | units.MSun
    p['max_imf_mass'] = 150.0 | units.MSun
    p['sample_imf_mass'] = 10000.0 | units.MSun
    p['sample_imf_bins'] = 10
    p['sink_rad'] = flashp['sink_accretion_radius'] | units.cm
    p['sum_small'] = False

    return p

