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

* In "evolution" methods, update as few workers as possible, to keep modular.
  For example:

    stellar evolution: update {hydro, AMUSE set} only
    bridge kick: update hydro only
    bridge evolve: update {hydro, grav} only

  Sync data between workers via explicit calls in top loop.
  Do not hide sync in "evolution" methods.

* hydro worker points to either "mass" (star) or "sink" particles at any given
  time.  We choose to point at star particles by default.  Anytime you access
  sink particles, don't forget to unpoint when done.

* (1) Write your comments now.  You won't have time to do it later.
  (2) If you change something, update the comments NOW.  Wrong comments are
  worse than no comments.
  -- Adapted from ENZO Developer's Guide
  https://enzo.readthedocs.io/en/latest/developer_guide/ProgrammingGuide.html

* Some general principles: http://google.github.io/styleguide/pyguide.html

"""

from __future__ import division, print_function

import numpy as np
np.set_printoptions(precision=3)

from amuse.lab import *
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
#from amuse.couple import multiples
import multiples_aaron as multiples # TODO -AT,2019oct30, edits to fold into AMUSE repo after testing
#from amuse.rfi.channel import AsyncRequestsPool

from torch_se import stellar_evolution
from torch_sf import (
    add_particles_to_grav,
    remove_particles_outside_bndbox,
    remove_particles_outside_bndbox_mult,
    make_stars_from_sinks,
    queue_stars,
)
from torch_state import TorchState
from torch_stdout import tprint

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

    if USER['with_ph4']:
        grav = ph4(convert, number_of_workers=USER['num_grav_workers'], mode='cpu', redirection="none")
        grav.parameters.set_defaults()
        grav.parameters.epsilon_squared = USER['epsilon']**2.0
        grav.parameters.force_sync = 1  # end exactly at requested time
        grav.parameters.timestep_parameter = 0.14  # timestep accuracy # TODO how was this chosen?! -AT,2019oct13
    else:
        grav = Hermite(convert, number_of_workers=USER['num_grav_workers'])
        grav.parameters.end_time_accuracy_factor = 0.0  # end exactly at requested time
        grav.parameters.dt_param = 0.02  # timestep size control, default 0.03

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

def evolve(state, hydro, grav, mult, se):

    # FLASH loop control
    hy_dt           = hydro.get_timestep()
    hy_step         = hydro.get_current_step()
    hy_time         = hydro.get_time()
    hy_max_steps    = hydro.get_max_num_steps()
    hy_max_time     = hydro.get_end_time()

    # grav loop control
    if USER['with_ph4']:
        grav.parameters.begin_time  = hy_time
        grav.parameters.sync_time   = hy_time
    else:
        grav.parameters.begin_time  = hy_time
        grav.evolve_model(hy_time)
    gr_time = grav.get_time()

    # stellar evolution timestep (hack for SN)
    # TODO this really shuld be handled by HYDRO and not torch -AT, 2019Oct14
    se_dt = 1e99 | units.s

    # bridge loop control
    it = 1
    tt = hy_time  # tt = "torch time" or "time in torch"
    dt = min(USER['hy_dt_factor']*hy_dt, se_dt, hy_max_time-tt)
    made_stars = False
    num_stars = hydro.get_number_of_particles()

    if num_stars > 0:  # restart or user initial conditions
        add_particles_to_grav(state, hydro, grav, mult)
        made_stars = True

    while tt < hy_max_time and hy_step < hy_max_steps:

        tprint("Bridge step: it={}, tt={}, dt={}".format(it, tt, dt))#Current simulation time:", tt, "dt:", dt)
        tprint("... Hydro time:", hy_time)
        tprint("... Grav time:", gr_time)
        tprint("... Num stars in hydro:", num_stars)
        if USER['with_multiples']:
            tprint("... Num in grav:", len(grav.particles))
            tprint("... Num in mult.root_to_tree:", len(mult.root_to_tree))
        tprint("... made_stars:", made_stars)

        if num_stars > 0:

            ### ------------------
            ### Stellar evolution.
            ### ------------------

            if USER['with_se']:
                tprint("Do stellar evolution")
                # update both stars set and hydro properties
                se_dt = stellar_evolution(
                    tt+dt, dt, state, hydro, se,
                    with_lyc          = USER['with_lyc'],
                    with_pe_heat      = USER['with_pe_heat'],
                    with_winds        = USER['with_winds'],
                    with_sn           = USER['with_sn'],
                    massloss_method   = USER['massloss_method'],
                    min_feedback_mass = USER['min_feedback_mass'],
                )
                tprint("... dt from stellar evol:", se_dt)

                # sync mass to gravity code(s) from stars
                state.stars_to_grav.copy_attributes(["mass"])  # AMUSE -> grav singles
                if USER['with_multiples']:
                    mult.channel_from_code_to_memory.copy() # grav  -> multiples
                    state.stars_to_mult_grav_copy("mass")   # AMUSE -> multiples, grav COM

            ### -----------
            ### First kick.
            ### -----------

            tprint("Evolving hydro with grav for dt:", dt, "to reach t =", tt+dt)

            if USER['with_bridge']:
                tprint("First kick")
                kick_number = 1
                if made_stars:  # update grav pot (BGPT), accel (BGA{X,Y,Z}) to account for new star prtl
                    kick_number = 2

                hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
                tprint("... grid kicked")
                hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
                tprint("... stars kicked")

                # sync velocity to stars + gravity code(s) from hydro
                state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)  # hydro -> AMUSE
                state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])              # AMUSE -> grav singles
                if USER['with_multiples']:
                    mult.channel_from_code_to_memory.copy()     # grav  -> multiples
                    state.stars_to_mult_grav_copy("velocity")   # AMUSE -> multiples, grav COM

            ### --------------
            ### Evolve models.
            ### --------------

            tprint("Advance grav")
            if USER['with_multiples']:
                mult.evolve_model(tt+dt)
            else:
                grav.evolve_model(tt+dt)

            tprint("Advance hydro")
            hydro.evolve_model(tt+dt)

            # sync position & velocity to stars + hydro from gravity code(s)
            state.grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])  # grav singles -> AMUSE
            if USER['with_multiples']:
                mult.update_leaves_pos_vel()  # grav COM -> multiples; this synchronizes full tree, i.e. also updates root, tree.particle
                mult.stars.copy_values_of_attributes_to(["x", "y", "z", "vx", "vy", "vz"], state.stars)  # grav singles AND multiples -> AMUSE
            hydro.set_particle_position(state.stars.tag, state.stars.x,  state.stars.y,  state.stars.z)  # AMUSE -> hydro
            hydro.set_particle_velocity(state.stars.tag, state.stars.vx, state.stars.vy, state.stars.vz)

            # this updates all of grav,stars,hydro
            if USER['with_multiples']:
                remove_particles_outside_bndbox_mult(state, hydro, grav, mult)
            else:
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

            num_stars = hydro.get_number_of_particles()

            if USER['with_bridge'] and num_stars > 0:  # in case star exited domain
                tprint("Second kick")
                kick_number = 2  # update grav pot (BGPT), accel (BGA{X,Y,Z}) for star prtl new positions

                hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
                tprint("... grid kicked")
                hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
                tprint("... stars kicked")

                # sync velocity to stars + gravity code(s) from hydro
                state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)  # hydro -> AMUSE
                state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])              # AMUSE -> grav singles
                if USER['with_multiples']:
                    mult.channel_from_code_to_memory.copy()     # grav singles -> multiples
                    state.stars_to_mult_grav_copy("velocity")   # AMUSE -> multiples, grav COM

        else: # num_stars == 0

            ### --------------
            ### Evolve models.
            ### --------------

            tprint("Evolving hydro without grav for dt:", dt, "to reach t =", tt+dt)
            hydro.evolve_model(tt + dt)

            # two possible cases:
            # 1. no stars yet, so never called grav.evolve_model(...)
            # 2. had stars, but they all escaped
            # not sure if below code works with case 2 of stars -> no stars
            hy_time = hydro.get_time()
            if USER['with_ph4']:
                grav.parameters.begin_time  = hy_time
                grav.parameters.sync_time   = hy_time
            else:
                grav.parameters.begin_time  = hy_time
                grav.evolve_model(hy_time)

        ### --------------------------------
        ### Queue and create star particles.
        ### --------------------------------

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
            add_particles_to_grav(state, hydro, grav, mult)  # push stars hydro->amuse, hydro->grav

        ### ---------------------------------------------
        ### Output FLASH and Torch plot,checkpoint files.
        ### ---------------------------------------------

        tprint("Output check")
        state.output()

        # FLASH loop control
        hy_dt = hydro.get_timestep()
        hy_step = hydro.get_current_step()
        hy_time = hydro.get_time()

        # grav loop control
        gr_time = grav.get_time()

        # bridge loop control
        it += 1
        tt += dt
        dt = min(USER['hy_dt_factor']*hy_dt, se_dt, hy_max_time-tt)
        num_stars = hydro.get_number_of_particles()  # loop variable

        assert abs((hy_time - gr_time).value_in(units.s)) <= 1e4
        assert num_stars == len(state.stars)
        if USER['with_multiples']:
            assert num_stars == len(mult.stars)
        else:
            assert num_stars == len(grav.particles)

    return

# ============================================================================

def run_torch(user_initial_conditions, user_parameters):
    """
    Run a Torch simulation.  This is called from a user script, which provides
    initial conditions and parameters for the desired problem set up.

    Arguments: requires two methods as input.

        user_initial_conditions(state, hydro)

            method that alters "state", "hydro" objects
            to set initial conditions for simulation.

        user_parameters()

            method that returns a dict of Torch configuration parameters

    Result: spawn the necessary FLASH, gravity, stellar evolution, etc. workers
    to run a Torch simulation.  Attempt to run the simulation to completion.
    """

    global USER
    USER = user_parameters()

    tprint("Num hydro workers: {:d}".format(USER['num_hy_workers']))
    tprint("Num grav workers: {:d}".format(USER['num_grav_workers']))

    if USER['npy_seed'] is not None:
        np.random.seed(USER['npy_seed'])

    hydro, grav, mult, se = initialize_workers()

    state = TorchState(hydro, grav, mult)

    state.initial_io(refresh=USER['restart_with_new_rng'])

    if not state.restart:
        user_initial_conditions(state, hydro)
    elif state.restart and USER['restart_with_user_ics']:
        # massage the hydro particle structures so that particles from user ICs
        # look like they came from restart checkpoint file.
        hydro.set_starting_local_tag_numbers()
        user_initial_conditions(state, hydro)
        hydro.clear_new_tags()

    try:

        evolve(state, hydro, grav, mult, se)

    finally:
        pass
        #hydro.timer_summary()
        #hydro.cleanup_code()
        #grav.stop()
        #kep.stop()
        #stop_smalln()
        #del multiples
