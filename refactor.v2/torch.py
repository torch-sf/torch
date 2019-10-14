#!/usr/bin/env python
"""
Fresh start.  No stellar evolution, no multiples, no initial conditions yet.

Coding principles:
* workers and evolution loops should be near top-level, to allow easy hooks / editing
* do not leave debugging hooks lying around.

* "state" object should be as general as possible.
  don't stuff methods into it; try not to use its workers.

* "single source of truth", to the most extent possible.

KNOWN ERRORS:
* if all stars escape domain during bridge loop, second bridge kick fails
  while trying to recompute BGPT_VAR.  -AT, 2019oct11

Aaron Tran
Started 2019 October 09
"""

from __future__ import division, print_function

#import datetime
#import glob
import numpy as np
#import os
#import pickle
#import sys
#import time

#from scipy.integrate import *
np.set_printoptions(precision=3)
#np.random.seed(103180)  # Set initial random seed for testing/debugging.  # TODO DEBUGGING
np.random.seed(203180)  # Set initial random seed for testing/debugging.  # TODO DEBUGGING

from amuse.lab import *
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
from amuse.community.flash import josh_multiples as multiples
#from amuse.rfi.channel import AsyncRequestsPool

from torch_sf import add_particles_to_grav, remove_particles_outside_bndbox, make_stars_from_sinks, queue_stars
from torch_state import TorchState
from torch_stdout import tprint

#import ionizingflux as ion

# ============================================================================
# Multiples boilerplate

SMALLN = None

def init_smalln(converter):
    global SMALLN
    SMALLN = SmallN(convert_nbody=converter)
    SMALLN.initialize_code()

def new_smalln():
    global SMALLN
    SMALLN.reset()
    return SMALLN

def stop_smalln():
    global SMALLN
    SMALLN.stop()

# ============================================================================

def evolve_bridge(state, hydro, grav, mult):

    # FLASH time stepping
    hy_dt           = hydro.get_timestep()
    hy_step         = hydro.get_current_step()
    hy_time         = hydro.get_time()
    hy_max_steps    = hydro.get_max_num_steps()
    hy_max_time     = hydro.get_end_time()

    # grav time stepping
    grav.parameters.begin_time  = hy_time
    grav.parameters.sync_time   = hy_time
    gr_time = grav.get_time()

    # bridge time stepping
    it = 1
    tt = hy_time  # tt = "torch time" or "time in torch"
    dt = min(1.5*hy_dt, hy_max_time-tt)  # TODO MAGIC NUMBERS -AT,2019oct10

    # more bridge loop variables
    made_stars = False
    num_stars = hydro.get_number_of_particles()
    if num_stars > 0:  # restart or user initial conditions
        add_particles_to_grav(state, hydro, grav)
        made_stars = True

    while tt < hy_max_time and hy_step < hy_max_steps:

        tprint("Bridge step: it={}, tt={}, dt={}".format(it, tt, dt))#Current simulation time:", tt, "dt:", dt)
        tprint("... Hydro time:", hy_time)
        tprint("... Grav time:", gr_time)
        tprint("... num_stars=", num_stars)
        tprint("... made_stars=", made_stars)

        if num_stars > 0:

            ### ------------------
            ### Stellar evolution.
            ### ------------------

            # TODO

            ### -----------
            ### First kick.
            ### -----------

            tprint("Evolving hydro with grav for dt:", dt, "to reach t =", tt+dt)

            if WITH_BRIDGE:
                tprint("First kick")
                kick_number = 1
                if made_stars:  # must recompute grav pot (BGPT), accel (BGA{X,Y,Z}) induced by new star prtl
                    kick_number = 2
                hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
                tprint("... grid kicked")
                hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
                tprint("... stars kicked")
                # sync workers
                state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)
                state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])

            ### --------------
            ### Evolve models.
            ### --------------
            tprint("Advance grav")
            grav.evolve_model(tt+dt)
            tprint("Advance hydro")
            hydro.evolve_model(tt+dt)

            # sync workers
            state.grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])
            hydro.set_particle_position(state.stars.tag, grav.particles.x, grav.particles.y, grav.particles.z)
            hydro.set_particle_velocity(state.stars.tag, grav.particles.vx, grav.particles.vy, grav.particles.vz)

            remove_particles_outside_bndbox(state, hydro, grav)
            # sort and also remove stars outside domain, though
            # remove_particles_outside_bndbox(...) should have us covered
            hydro.particles_sort()

            # TODO can we move star creation here?
            # this saves one poisson solve when new stars are being formed.
            # edge case of restarts can get weird
            # -AT,2019Oct10

            ### ------------
            ### Second kick.
            ### ------------
            if WITH_BRIDGE:
                tprint("Second kick")
                kick_number = 2  # update grav pot (BGPT), accel (BGA{X,Y,Z}) for star prtl new positions
                hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
                tprint("... grid kicked")
                hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
                tprint("... stars kicked")
                # sync workers
                state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)
                state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])

        else: # num_stars == 0

            tprint("Evolving hydro without grav for dt:", dt, "to reach t =", tt+dt)
            hydro.evolve_model(tt + dt)

            # two possible cases:
            # 1. no stars yet, so never called grav.evolve_model(...)
            # 2. had stars, but they all escaped
            # not sure if below code works with case 2 of stars -> no stars
            hy_time = hydro.get_time()
            grav.parameters.begin_time  = hy_time
            grav.parameters.sync_time   = hy_time

        # if any new sinks, draw a list of stars from Kroupa IMF
        tprint("Star formation check")
        queue_stars(state, hydro, MIN_SF_MASS, MAX_SF_MASS,
                    sample_imf_mass=SAMPLE_MASS, sum_small=SUM_SMALL,
                    sample_imf_bins=SAMPLE_IMF_BINS)  # TODO magic numbers

        made_stars = make_stars_from_sinks(state, hydro, sink_rad=SINK_RAD)  # in hydro  # TODO fugly parameter handling
        if made_stars:
            add_particles_to_grav(state, hydro, grav)  # push stars hydro->amuse, hydro->grav

        num_stars = hydro.get_number_of_particles()  # loop variable
        assert num_stars == len(state.stars)
        assert num_stars == len(grav.particles)  # only true without multiples

        # write output iff it's time to do so
        tprint("Output check")
        state.output()

        # update bridge loop variables
        tt += dt
        it += 1

        gr_time = grav.get_time()
        hy_step = hydro.get_current_step()
        hy_time = hydro.get_time()

        assert abs((hy_time - gr_time).value_in(units.s)) <= 1e4

        hy_dt = hydro.get_timestep()
        dt = min(1.5*hy_dt, hy_max_time-tt, 2*dt)  # cannot increase dt more than 2x  # TODO MAGIC NUMBERS -AT,2019oct10

    return

# ============================================================================

def user_init(state, hydro):
    """user initialization, just add stuff to hydro"""

    if state.restart:
        return

    # Plop a single star into FLASH
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

# ============================================================================

if __name__ == '__main__':

    # --------------------
    # User configuration

    global EPS
    global NUM_GRAV_WORKERS
    global NUM_HY_WORKERS
    global REFRESH_RNG

    global WITH_BRIDGE
    global WITH_MULTIPLES
    global WITH_RADIATION
    global WITH_PE_HEAT
    global WITH_SE
    global WITH_SN
    global WITH_WINDS
    global WITH_MASSLOSS
    global MASSLOSS_METHOD

    global MIN_MASS
    global SINK_RAD

    global MIN_SF_MASS
    global MAX_SF_MASS
    global SAMPLE_MASS
    global SAMPLE_IMF_BINS
    global SUM_SMALL

    EPS = 15.0 | units.RSun  # N-body softening = actual radius of a massive star
    MULT_DEBUG_LEVEL = 1
    NUM_GRAV_WORKERS = 4
    NUM_HY_WORKERS = 9
    REFRESH_RNG     = False

    # OVERALL code toggles/logic
    WITH_BRIDGE     = True
    WITH_MULTIPLES  = False  # adds three workers: kepler, smalln, multiples
    WITH_SE         = True

    WITH_RADIATION  = True  # stellar evolution switches
    WITH_PE_HEAT    = True
    WITH_SN         = True
    WITH_WINDS      = True
    WITH_MASSLOSS   = True
    MASSLOSS_METHOD = 'puls'

    MIN_MASS    = 7.0 | units.MSun
    SINK_RAD    = 9.16156e+18 | units.cm # TODO hydro.get_runtime_parameter('sink_accretion_radius') | units.cm

    MIN_SF_MASS = 0.08 | units.MSun
    MAX_SF_MASS = 150.0 | units.MSun
    SAMPLE_MASS = 10000 | units.MSun
    SAMPLE_IMF_BINS = 10
    SUM_SMALL   = False

    # --------------------
    # Internal configuration

    # Converter for the N-body code.
    convert = nbody.nbody_to_si(1.0|units.parsec, 1000.0|units.MSun)
    # Converter for the hydro code.
    convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(1.0|units.cm, 1.0|units.g, 1|units.s)

    # --------------------
    # Worker init

    se = SeBa()
    se.initialize_code()

    grav = ph4(convert, number_of_workers=NUM_GRAV_WORKERS, mode='cpu', redirection="none")
    grav.parameters.set_defaults()
    grav.parameters.epsilon_squared = EPS**2.0
    grav.parameters.force_sync = 1  # end exactly at requested time
    grav.parameters.timestep_parameter = 0.14  # timestep accuracy

    mult = None

    if WITH_MULTIPLES:

        grav.parameters.epsilon_squared = 0.0|units.cm**2.0
        grav.stopping_conditions.collision_detection.enable()

        init_smalln(convert)

        kep = Kepler(unit_converter=convert)
        kep.initialize_code()

        mult = multiples.Multiples(grav, new_smalln, kep, constants.G)
        mult.global_debug                = 1
        mult.neighbor_veto               = True
        mult.check_tidal_perturbation    = True
        mult.neighbor_perturbation_limit = 0.05
        mult.wide_perturbation_limit     = 0.08

    hydro = Flash(unit_converter=convert2, number_of_workers=NUM_HY_WORKERS, redirection='none')
    hydro.initialize_code()
    hydro.set_particle_pointers('mass')  # code convention: hydro should point to star prtl by default

    # --------------------
    # AMUSE framework state init (Particles set, channels, sink lists, etc)

    state = TorchState(hydro, grav, mult, refresh=REFRESH_RNG)
    state.initialize() # loads restart files if needed

    # --------------------
    # Apply user initial conditions to hydro

    user_init(state, hydro)

    try:

        evolve_bridge(state, hydro, grav, mult)

    finally:
        pass
        #hydro.timer_summary()
        #hydro.cleanup_code()
        #grav.stop()
        #kep.stop()
        #stop_smalln()
        #del multiples
