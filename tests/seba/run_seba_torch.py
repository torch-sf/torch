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
    # Stellar evolution test

    stars          = Particles(3)
    stars.mass     = [20, 40, 80] | units.MSun
    stars.position = [[2, -2, 0], [0, 2, -2], [-2, 0, 2]] | units.pc
    stars.velocity = [[0, 0, 0], [0, 0, 0], [0, 0, 0]] | units.pc/units.Myr

    star_tags = hydro.add_particles(stars.x, stars.y, stars.z)
    hydro.set_particle_mass(star_tags, stars.mass)
    hydro.set_particle_velocity(star_tags, stars.vx, stars.vy, stars.vz)
    hydro.set_particle_oldmass(star_tags, stars.mass)

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

        p['with_se'] = True
        p['with_petar'] = True
        p['with_bridge'] = True
        p['with_multiples'] = False 

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
