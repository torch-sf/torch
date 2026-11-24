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

    # SN with SE test: plop a star that goes SN within 5e11 seconds
    add_stars = read_set_from_file('test_sne.amuse')
    
    creation_time = hydro.get_time() - add_stars.age
    tag = hydro.add_particles(add_stars.x, add_stars.y, add_stars.z)
    hydro.set_particle_mass(tag, add_stars.mass)
    hydro.set_particle_velocity(tag, add_stars.vx, add_stars.vy, add_stars.vz)                                                                                                                 
    hydro.set_particle_oldmass(tag, add_stars.initial_mass)
    hydro.set_particle_creation_time(tag, creation_time)
    hydro.set_particle_rel_age(tag, add_stars.relative_age)
    hydro.set_particle_rel_mass(tag, add_stars.relative_mass)
    hydro.set_particle_corem(tag, add_stars.core_mass)
    hydro.set_particle_co_corem(tag, add_stars.COcore_mass)
    hydro.set_particle_stype(tag, add_stars.stellar_type.value_in(units.stellar_type))
    hydro.set_particle_radius(tag, add_stars.radius)

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
        p['with_sn'] = True

        p['with_winds'] = False
        p['with_bridge'] = False  
        p['with_multiples'] = False 
        p['with_ph4'] = False
        p['with_lyc'] = False 
        p['with_pe_heat'] = False 
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
