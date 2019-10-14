#!/usr/bin/env python
"""
Rewrite of Josh's bridge_multiples.py, to be more readable and extensible.

CODING PRINCIPLES:

* workers and evolution loop should be near top-level, so the user can easily
  edit and add hooks to the evolution loop.

* do not leave debugging hooks lying around, unless you think
  (1) many users / developers will use hooks,
  (2) a single dev will use hooks many times.

* the "state" object is a container for global state, and to help with I/O.
  don't stuff too many methods into it; try not to use its workers.

* "single source of truth", to the most extent possible.
  In general, let hydro hold truth.  It already must do so for e.g. restarts.

* Try not to pass Torch parameter struct into deeper methods.
  Keep abstraction layers well isolated.

"""

from __future__ import division, print_function

import numpy as np
np.set_printoptions(precision=3)

from amuse.lab import *
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
from amuse.community.flash import josh_multiples as multiples
#from amuse.rfi.channel import AsyncRequestsPool

from torch_se import stellar_evolution
from torch_sf import add_particles_to_grav, remove_particles_outside_bndbox, make_stars_from_sinks, queue_stars
from torch_state import TorchState
from torch_stdout import tprint
from torch_user import user_initial_conditions, user_parameters

# ============================================================================
# Multiples boilerplate - required as of Oct 2019, see AMUSE book.

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

def initialize_workers():

    # Converter for the N-body code.
    convert = nbody.nbody_to_si(1.0|units.parsec, 1000.0|units.MSun)
    # Converter for the hydro code.
    convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(1.0|units.cm, 1.0|units.g, 1|units.s)

    grav = ph4(convert, number_of_workers=USER['num_grav_workers'], mode='cpu', redirection="none")
    grav.parameters.set_defaults()
    grav.parameters.epsilon_squared = USER['epsilon']**2.0
    grav.parameters.force_sync = 1  # end exactly at requested time
    grav.parameters.timestep_parameter = 0.14  # timestep accuracy # TODO how was this chosen?! -AT,2019oct13

    mult = None

    if USER['with_multiples']:

        grav.parameters.epsilon_squared = 0.0|units.cm**2.0
        grav.stopping_conditions.collision_detection.enable()

        init_smalln(convert)

        kep = Kepler(unit_converter=convert)
        kep.initialize_code()

        mult = multiples.Multiples(grav, new_smalln, kep, constants.G)
        mult.global_debug                = 1
        mult.neighbor_veto               = True
        mult.check_tidal_perturbation    = True
        mult.neighbor_perturbation_limit = 0.05 # TODO how was this chosen?! -AT,2019oct13
        mult.wide_perturbation_limit     = 0.08

    se = None

    if USER['with_se']:

        se = SeBa()
        se.initialize_code()

    hydro = Flash(unit_converter=convert2, number_of_workers=USER['num_hy_workers'], redirection='none')
    hydro.initialize_code()
    hydro.set_particle_pointers('mass')  # code convention: hydro should point to star prtl by default

    return hydro, grav, mult, se

# ============================================================================

def evolve(state, hydro, grav, mult):

    # FLASH loop control
    hy_dt           = hydro.get_timestep()
    hy_step         = hydro.get_current_step()
    hy_time         = hydro.get_time()
    hy_max_steps    = hydro.get_max_num_steps()
    hy_max_time     = hydro.get_end_time()

    # grav loop control
    grav.parameters.begin_time  = hy_time
    grav.parameters.sync_time   = hy_time
    gr_time = grav.get_time()

    # stellar evolution timestep (hack for SN)
    # TODO this really shuld be handled by HYDRO and not torch -AT, 2019Oct14
    se_dt = 1e99 | units.s

    # bridge loop control
    it = 1
    tt = hy_time  # tt = "torch time" or "time in torch"
    dt = min(1.5*hy_dt, se_dt, hy_max_time-tt)  # TODO MAGIC NUMBERS -AT,2019oct10

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
            if USER['with_se']:
                tprint("Do stellar evolution")
                se_dt = stellar_evolution(
                    tt+dt, dt, state, hydro, se,
                    with_lyc          = USER['with_lyc'],
                    with_pe_heat      = USER['with_pe_heat'],
                    with_winds        = USER['with_winds'],
                    with_sn           = USER['with_sn'],
                    massloss_method   = USER['massloss_method'],
                    min_feedback_mass = USER['min_feedback_mass'],
                )

            ### -----------
            ### First kick.
            ### -----------

            tprint("Evolving hydro with grav for dt:", dt, "to reach t =", tt+dt)

            if USER['with_bridge']:
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
            if USER['with_bridge']:
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

        tprint("Star formation check")
        queue_stars(state, hydro,
            min_imf_mass=USER['min_imf_mass'],
            max_imf_mass=USER['max_imf_mass'],
            sample_imf_mass=USER['sample_imf_mass'],
            sample_imf_bins=USER['sample_imf_bins'],
            sum_small=USER['sum_small'],
        )

        made_stars = make_stars_from_sinks(state, hydro, sink_rad=USER['sink_rad'])  # in hydro
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
        dt = min(USER['hy_dt_factor']*hy_dt, se_dt, hy_max_time-tt)

    return

# ============================================================================

if __name__ == '__main__':

    global USER
    USER = user_parameters()

    if USER['npy_seed'] is not None:
        np.random.seed(USER['npy_seed'])

    hydro, grav, mult, se = initialize_workers()

    state = TorchState(hydro, grav, mult)

    state.initial_io(refresh=USER['refresh_rng'])

    user_initial_conditions(state, hydro)

    try:

        evolve(state, hydro, grav, mult)

    finally:
        pass
        #hydro.timer_summary()
        #hydro.cleanup_code()
        #grav.stop()
        #kep.stop()
        #stop_smalln()
        #del multiples
