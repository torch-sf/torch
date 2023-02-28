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

import time
import numpy as np
np.set_printoptions(precision=3)

from amuse.lab import *
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
from amuse.community.petar.interface import Petar
from amuse.couple import multiples

from torch_se import stellar_evolution
from torch_sf import (
    add_particles_to_grav,
    remove_particles_outside_bndbox,
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
    convert = nbody.nbody_to_si(1.0|units.kyr, 1000.0|units.MSun)
    # Converter for the hydro code.
    convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(1.0|units.cm, 1.0|units.g, 1|units.s)

    if USER['with_ph4']:
        grav = ph4(convert, number_of_workers=USER['num_grav_workers'], mode='cpu', redirection='none')
        grav.parameters.set_defaults()
        grav.parameters.epsilon_squared = USER['epsilon']**2.0
        grav.parameters.force_sync = 1  # end exactly at requested time
        grav.parameters.timestep_parameter = 0.14  # timestep accuracy # TODO how was this chosen?! -AT,2019oct13
    elif USER['with_petar']:
        grav = Petar(convert, number_of_workers=USER['num_grav_workers'], mode='cpu', redirection='none')
        grav.parameters.epsilon_squared = USER['epsilon']**2.0
        grav.parameters.r_bin = 1.496e15 | units.cm # 100AU
        #aveStarMass = 1.234e33 | units.g
        #velDisp = 1.7e5 | units.cm/units.s
        G = 6.67428e-8 | units.cm**3 / units.g / units.s**2
        #grav.parameters.r_out = 12.5*grav.parameters.r_bin
        #grav.parameters.dt_soft = (np.pi/8.0)*np.sqrt(((grav.parameters.r_out/2.0)**3)/(2*G*aveStarMass))
        #grav.parameters.r_search_min = grav.parameters.r_out + 3.0*grav.parameters.dt_soft*velDisp
    else:
        grav = Hermite(convert, number_of_workers=USER['num_grav_workers'], redirection='none')
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
        mult.global_debug                = 0
        mult.neighbor_veto               = True
        mult.check_tidal_perturbation    = False # Default: False. True: outputs diagnostics for highest perturbers. - SCL,2021oct5
        mult.neighbor_perturbation_limit = 0.05 # TODO how was this chosen?! -AT,2019oct13
        mult.wide_perturbation_limit     = 0.08

    se = None

    if USER['with_se']:

        se = SeBa()
        se.initialize_code()

    if USER['evolve_async']:
        hydro = Flash(
            unit_converter=convert2,
            number_of_workers=USER['num_hy_workers'],
            redirection='file',
            redirect_stdout_file='flash_worker.out',
            redirect_stderr_file='flash_worker.err',
        )
    else:
        hydro = Flash(
            unit_converter=convert2,
            number_of_workers=USER['num_hy_workers'],
            redirection='none',
        )

    hydro.initialize_code()
    hydro.set_particle_pointers('mass')  # code convention: hydro should point to star prtl by default

    return hydro, grav, mult, se

# ============================================================================

def evolve(state, hydro, grav, mult, se):

    time_file = open("grav_timer.txt",'w')

    # FLASH loop control
    hy_dt           = hydro.get_timestep()
    hy_step         = hydro.get_current_step()
    hy_time         = hydro.get_time()
    hy_max_steps    = hydro.get_max_num_steps()
    hy_max_time     = hydro.get_end_time()

    # stellar evolution timestep (hack for SN)
    # TODO this really shuld be handled by HYDRO and not torch -AT, 2019Oct14
    se_dt = 1e99 | units.s

    # bridge loop control
    it = 1
    dt = min(USER['hy_dt_factor']*hy_dt, se_dt, hy_max_time-hy_time)
    # set initial hydro dt to a power of 2 so PeTar can sync times
    if USER['with_petar']:
        print("nbody time = ",nbody.time)
        dt_nbody = pow(2., np.floor(np.log2(dt.value_in(units.kyr)))) | units.kyr
        dt = dt_nbody

    num_stars = hydro.get_number_of_particles()

    if not USER['with_petar']: # only initialize PeTar if there are stars
        grav.parameters.begin_time  = hy_time
        grav.evolve_model(hy_time)
        gr_time = grav.get_time()

    # worker setup
    if num_stars > 0:  # restart or user initial conditions
        # if this is a restart, FLASH may still have all the
        # particles mis-sorted in the particles array. -JW
        hydro.particles_sort()
        add_particles_to_grav(state, hydro, grav, mult, se)

    if USER['evolve_async']:
        from amuse.rfi.async_request import AsyncRequestsPool
        pool = AsyncRequestsPool()
        pool_table_hydro = []
        pool_table_grav = []
        def handle_result(request, name, i):
            assert request.is_result_available()
            if name == "hydro":
                pool_table_hydro.append(i)
            elif name == "grav":
                pool_table_grav.append(i)



    first_star = 0

    while hy_time < hy_max_time and hy_step < hy_max_steps:

        tprint("Bridge step: it={}, t={:e}, dt={:e}".format(
            it, hy_time.value_in(units.s), dt.value_in(units.s),
        ))
        tprint("... Hydro step:", hy_step)
        if USER['with_multiples']:
            tprint("... Num stars: {:d} (singles {:d}, multiples {:d})".format(
                    num_stars,
                    len(grav.particles) - len(mult.root_to_tree),
                    len(mult.root_to_tree)
            ))
        else:
            tprint("... Num stars:", num_stars)

        if num_stars > 0:

            # initialize PeTar once more than 1!!! star forms
            if num_stars > 1 and first_star == 0:
                first_star = 1
                if USER['with_petar']:
                    tprint("First stars have formed. Initializing PETAR.")
                    grav.parameters.begin_time = hy_time
                    grav.evolve_model(hy_time)
                    print(grav.parameters)

            tprint("Evolving hydro with grav to reach t =", hy_time+dt)

            ### ------------------
            ### First bridge kick.
            ### ------------------
            remove_particles_outside_bndbox(state, hydro, grav, mult)
            hydro.particles_sort()  # also checks for stars outside domain

            if USER['with_bridge']:
                tprint("First bridge kick")
                kick_number = 1  # tell FLASH to NOT recompute grav pot (BGPT), accel (BGA{X,Y,Z}) from stars
                if it == 1:  # but, do calculate BGPT/etc for first time if stars on grid at simulation init
                    tprint("... first bridge step, recompute BGPT_VAR")
                    kick_number = 2
                hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
                tprint("... grid kicked")
                hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
                tprint("... stars kicked")

                # sync velocity to stars + gravity code(s) from hydro
                state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)  # hydro -> AMUSE

                if num_stars > 1: # don't run N-body with only 1 star
                    state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])              # AMUSE -> grav singles
                    if USER['with_multiples']:
                        mult.channel_from_code_to_memory.copy()     # grav  -> multiples
                        state.stars_to_mult_grav_copy("velocity")   # AMUSE -> multiples, grav COM

                remove_particles_outside_bndbox(state, hydro, grav, mult)
                hydro.particles_sort()  # also checks for stars outside domain

            ### ------------------
            ### Stellar evolution.
            ### ------------------

            if USER['with_se']:
                tprint("Do stellar evolution")
                # update both stars set and hydro properties
                se_dt = stellar_evolution(
                    hy_time+dt, dt, state, hydro, se,
                    with_lyc          = USER['with_lyc'],
                    with_pe_heat      = USER['with_pe_heat'],
                    with_winds        = USER['with_winds'],
                    with_sn           = USER['with_sn'],
                    massloss_method   = USER['massloss_method'],
                    min_feedback_mass = USER['min_feedback_mass'],
                )
                tprint("... dt from stellar evol:", se_dt)  # IF we keep this python-level dt management, this probably should enter hydro dt right away... -AT, 2019 nov 26

                # sync mass to gravity code(s) from stars
                if num_stars > 1:
                    state.stars_to_grav.copy_attributes(["mass"])  # AMUSE -> grav singles
                    #print(grav.particles.radius)
                    state.stars_to_grav.copy_attributes(["radius"])
                    #print(grav.particles.radius)
                    if USER['with_multiples']:
                        mult.channel_from_code_to_memory.copy() # grav  -> multiples
                        state.stars_to_mult_grav_copy("mass")   # AMUSE -> multiples, grav COM

            ### --------------
            ### Evolve models.
            ### --------------

            if num_stars > 1:
                if USER['evolve_async']:
                # Example async request code:
                # amuse/src/amuse/test/suite/compile_tests/test_python_implementation.py
                    tprint("Advance grav and hydro asynchronously")

                    if USER['with_multiples']:
                        req_hydro = hydro.evolve_model.asynchronous(hy_time+dt)
                        pool.add_request(req_hydro, handle_result, ["hydro", it])
                        tprint("... hydro submitted")

                        # Multiples is not a worker code, so we can't send it to
                        # the AsyncRequestsPool.
                        mult.evolve_model(hy_time+dt)
                        tprint("... grav advanced")

                        pool.wait()
                        tprint("... both grav and hydro advanced")

                    else:
                        req_hydro = hydro.evolve_model.asynchronous(hy_time+dt)
                        if USER['with_petar']:
                            grav.parameters.dt_soft = dt
                        req_grav = grav.evolve_model.asynchronous(hy_time+dt)
                        pool.add_request(req_hydro, handle_result, ["hydro", it])
                        pool.add_request(req_grav, handle_result, ["grav", it])

                        pool.wait()
                        if pool_table_hydro and pool_table_hydro[-1] == it:
                            tprint("... hydro advanced")
                        elif pool_table_grav and pool_table_grav[-1] == it:
                            tprint("... grav advanced")

                        pool.wait()
                        tprint("... both grav and hydro advanced")

                else:  # evolve models sequentially

                    tprint("Advance grav")
                    if USER['with_multiples']:
                        mult.evolve_model(hy_time+dt)
                    else:
                        if USER['with_petar']:
                            grav.parameters.dt_soft = dt
                        start_t = time.time()
                        grav.evolve_model(hy_time+dt)
                        gr_evolve_time = time.time()-start_t
                        time_file.write(str(gr_evolve_time)+" "+str(num_stars)+" "+str(hy_time+dt)+"\n") 
                        time_file.flush()
                    tprint("Advance hydro")
                    hydro.evolve_model(hy_time+dt)

                if (grav.get_time()-hydro.get_time() >= 1e4|units.s):
                    tprint("Evolving hydro further to sync with PeTar")
                    tprint("grav-hydro time = ",grav.get_time()-hydro.get_time())
                    hydro.evolve_model(grav.get_time())

                # sync position & velocity to stars + hydro from gravity code(s)
                state.grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])  # grav singles -> AMUSE
                if USER['with_multiples']:
                    mult.update_leaves_pos_vel()  # grav COM -> multiples; updates tree.particle and leaves (but not root, weirdly)
                    mult.stars.copy_values_of_attributes_to(["x", "y", "z", "vx", "vy", "vz"], state.stars)  # multiples AND grav singles -> AMUSE
                hydro.set_particle_position(state.stars.tag, state.stars.x,  state.stars.y,  state.stars.z)  # AMUSE -> hydro
                hydro.set_particle_velocity(state.stars.tag, state.stars.vx, state.stars.vy, state.stars.vz)

            else: # num_stars=1

                tprint("Evolving hydro without grav to reach t =", hy_time+dt)

                ### --------------
                ### Evolve models.
                ### --------------

                hydro.evolve_model(hy_time+dt)
                hy_time = hydro.get_time()


        else: # num_stars == 0

            tprint("Evolving hydro without grav to reach t =", hy_time+dt)

            ### --------------
            ### Evolve models.
            ### --------------

            hydro.evolve_model(hy_time+dt)

            # two possible cases:
            # 1. no stars yet, so never called grav.evolve_model(...)
            # 2. had stars, but they all escaped
            # not sure if below code works with case 2 of stars -> no stars
            hy_time = hydro.get_time()
            if not USER['with_petar']: # PeTar cannot be evolved with 0 stars
                grav.parameters.begin_time  = hy_time
                grav.evolve_model(hy_time)

        ### --------------------------------
        ### Queue and create star particles.
        ### --------------------------------

        ### ----------------------------
        ### Remove stars outside domain.
        ### ----------------------------
        # updates all of grav,stars,hydro,mult; can accept mult=None
        remove_particles_outside_bndbox(state, hydro, grav, mult)
        hydro.particles_sort()  # also checks for stars outside domain

        tprint("Star formation check")
        queue_stars(state, hydro,
            min_imf_mass=USER['min_imf_mass'],
            max_imf_mass=USER['max_imf_mass'],
            sample_imf_mass=USER['sample_imf_mass'],
            sample_imf_bins=USER['sample_imf_bins'],
            sum_small=USER['sum_small'],
            binaries=USER['binaries']
        )
        made_stars = make_stars_from_sinks(state, hydro, sink_rad=USER['sink_rad'])  # in hydro
        if made_stars:
            add_particles_to_grav(state, hydro, grav, mult, se)  # push stars hydro->amuse, hydro->grav

        ### ----------------------------
        ### Remove stars outside domain.
        ### ----------------------------
        # updates all of grav,stars,hydro,mult; can accept mult=None
        remove_particles_outside_bndbox(state, hydro, grav, mult)
        hydro.particles_sort()  # also checks for stars outside domain

        ### -------------------
        ### Second bridge kick.
        ### -------------------

        num_stars = hydro.get_number_of_particles()
        if num_stars > 0 and USER['with_bridge']:  # in case all stars exited domain
            tprint("Second bridge kick")
            kick_number = 2  # recompute grav pot (BGPT), accel (BGA{X,Y,Z}) from stars
            hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # star->gas, star->sink kick
            tprint("... grid kicked")
            hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # gas->star, sink->star kick
            tprint("... stars kicked")

            # sync velocity to stars + gravity code(s) from hydro
            state.stars.velocity = hydro.get_particle_velocity(state.stars.tag)  # hydro -> AMUSE

            if num_stars > 1: # Don't run N-body with one star
                state.stars_to_grav.copy_attributes(["vx", "vy", "vz"])              # AMUSE -> grav singles
                if USER['with_multiples']:
                    mult.channel_from_code_to_memory.copy()     # grav  -> multiples
                    state.stars_to_mult_grav_copy("velocity")   # AMUSE -> multiples, grav COM

        ### ---------------------------------------------
        ### Output FLASH and Torch plot,checkpoint files.
        ### ---------------------------------------------

        tprint("Output check")
        state.output(overwrite=USER['overwrite'])

        ### ----------------------
        ### Prepare for next loop.
        ### ----------------------

        # FLASH loop control
        hy_dt = hydro.get_timestep()
        hy_step = hydro.get_current_step() + 1  # need +1 because AMUSE coupling changes FLASH nstep logic
        hy_time = hydro.get_time()

        # grav loop control
        gr_time = grav.get_time()

        # bridge loop control
        it += 1
        dt = min(USER['hy_dt_factor']*hy_dt, se_dt, hy_max_time-hy_time)
        # set initial hydro dt to a power of 2 so PeTar can sync times
        if USER['with_petar']:
            dt_nbody = pow(2., np.floor(np.log2(dt.value_in(units.kyr)))) | units.kyr
            dt = dt_nbody
        num_stars = hydro.get_number_of_particles()  # loop variable

        if USER['with_petar']:
            # only assert time-sync with PeTar if stars have formed
            if first_star==1:
                assert abs(hy_time - gr_time) <= (1e4|units.s)
                print("hydro-grav time = ",hy_time - gr_time)
        else:
            assert abs(hy_time - gr_time) <= (1e4|units.s)
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
    tprint("AMUSE overwrite: {}".format(USER['overwrite']))

    if USER['npy_seed'] is not None:
        np.random.seed(USER['npy_seed'])

    hydro, grav, mult, se = initialize_workers()

    state = TorchState(hydro, grav, mult)

    state.initial_io(overwrite=USER['overwrite'], refresh=USER['restart_with_new_rng'])

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

