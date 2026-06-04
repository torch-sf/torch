from torch_user import *
from functools import wraps

# User ics for test 
def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.

    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """

    # ------------------------------------------------------------------------
    # Stromgren sphere tests

    star          = Particles(1)
    star.mass     = 30 | units.MSun
    star.position = [0, 0, 0] | units.cm
    star.velocity = [0, 0, 0] | units.cm/units.s

    nion    = 1e48 | units.s**-1
    eion    = 4.47337409e-12 | units.erg

    star_tag = hydro.add_particles(star.x, star.y, star.z)
    hydro.set_particle_mass(star_tag, star.mass)
    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
    hydro.set_particle_nion(star_tag, nion)
    hydro.set_particle_eion(star_tag, eion)

    return

# Create a decorator for user_parameters to update only the necessary parameters
def user_parameter_decorator(func):
    """
    Decorates the user_parameter function to update the torch parameters
    needed for the test. This prevents tests from breaking when torch_user.py 
    gets updated.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        p = func()
        
        p['with_bridge'] = False  
        p['with_multiples'] = False 
        p['with_se'] = False 
        p['with_ph4'] = False
        p['with_lyc'] = False 
        p['with_pe_heat'] = False 
        p['with_winds'] = False
        p['remove_merged'] = False

        return p
    return wrapper


# Apply the decorator manually to the imported function
user_parameters = user_parameter_decorator(user_parameters)

# run torch
if __name__ == '__main__':
    run_torch(
        user_initial_conditions,
        user_parameters,
    )
