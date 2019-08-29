"""
Josh's bridge_multiples.py code, rewritten and stripped down for clarity and
simplicity.  AT, 2019 June 30
"""

def main():
    """Main bridge loop"""

    initialize()

    try:
        bridge_loop()
    finally:
        hydro.timer_summary()

    cleanup()


def initialize():
    """Start AMUSE workers"""

    init_smalln(...)
    kep = Kepler(...)
    hydro = Flash(...)

    load_rnd_state_files(...)

    make_single_star_in_hydro(...)  # this should be user-configurable


def bridge_loop():
    """Do big loop"""

    while t < tmax and curr_hy_step < max_hy_steps:

        dt = min(dtmax, 1.5*hy_dt, (tmax-t), 2*dt_old)

        dt_old = dt
        t_old = t

        make_stars_from_sinks(hydro, msf_min, msf_max)

        if num_particles == 0:

            t = t + dt
            hydro.evolve_model(t)

        else:

            do_stellar_evolution()

            if bridge:
                bridge_kick(hydro, stars, eps, dt, step=1)

            t = t + dt
            grav.evolve_model(t)
            hydro.evolve_model(t)

            if bridge:
                bridge_kick(hydro, stars, eps, dt, step=2)

        # ALTERNATIVELY, write in such a way that
        # passing data structures with stars=[]
        # to gravity, bridge, etc causes a NO-OP


def cleanup():
    """Clean up"""
    whatever(...)


if __name__ == '__main__':
    main()
