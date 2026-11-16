from torch_user import *
from functools import wraps
from amuse.ext.orbital_elements import generate_binaries

# User ics for test 
def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.

    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """

    # ------------------------------------------------------------------------
    # Test for a binary within PeTar's regularization regime.

    primary, secondary = generate_binaries(
                            primary_mass = 0.5 | units.MSun,
                            secondary_mass = 0.5 | units.MSun,
                            semi_major_axis = 55.0 | units.au,
    )

    # The above returns two particlesets; make a particleset from the two particles
    stars = Particles()
    stars.add_particle(primary[0])   
    stars.add_particle(secondary[0])

    star_tag = hydro.add_particles(stars.x, stars.y, stars.z)
    hydro.set_particle_mass(star_tag, stars.mass)
    hydro.set_particle_velocity(star_tag, stars.vx, stars.vy, stars.vz)
    hydro.set_particle_oldmass(star_tag, stars.mass)

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
        
        p['with_bridge'] = True  
        p['with_multiples'] = False   
        p['with_ph4'] = False

        p['binaries'] = False
        p['with_be'] = False

        p['with_se'] = False 
        p['with_lyc'] = False 
        p['with_pe_heat'] = False 
        p['with_winds'] = False
        p['remove_merged'] = False

        p['with_petar'] = True
        p['r_bin'] = 50. | units.au
        p['r_out'] = 0.003 | units.pc #Change with below 
        p['dt_soft'] = 4. | units.yr

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
