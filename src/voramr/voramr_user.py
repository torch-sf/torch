from src.voramr_mainloop import (run_flash, get_ntasks_from_run_script) 

def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.
    Usually, this means adding star particles to the hydro code.

    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """
    # Uncomment the following to pass star placement commands to FLASH durinng
    # hydrodynamical initialization. This code can be modified to include any
    # number of stars hand placed, systematically generated, or imported from
    # other simulations data - Sean C. Lewis
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

    return
def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """
    p = {}

    p['source_file'] = "snapshot_550_9.hdf5"#"voramr_test.hdf5"
    p['convert_file'] = True
    p['input_file'] = "voramr_input.hdf5"
    p['pickle_kdtree'] = True
    p['pickle_file_name'] = "kdtree.pickle"
    p['numBlocks'] = 15000 #345
    p['cellsPerBlock'] = 16
    p['num_hy_workers'] = get_ntasks_from_run_script("submit") - 1
    
    return p

if __name__ == '__main__':
    run_flash(
        user_initial_conditions,
        user_parameters
        )
