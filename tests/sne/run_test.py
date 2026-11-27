from functools import wraps

from amuse.community.seba.interface import SeBa

from torch_user import *
from torch_stdout import tprint

# User ics for test 
def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.

    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """

    # ------------------------------------------------------------------------
    # Supernova blastwave test
    # Single star with initial mass and target age selected to result in
    # a supernova shortly after starting the simulation. Star is initialized at 
    # center of box and evolved by SeBa to reach the appropriate stage.

    # Create star
    initial_mass = 50.0 | units.Msun
    explosion_age = 4.7613 | units.Myr
    target_age = explosion_age - (500 | units.yr)

    star = Particles(1)
    star.mass = initial_mass
    star.initial_mass = initial_mass
    star.position = [0, 0, 0] | units.cm
    star.velocity = [0, 0, 0] | units.cm/units.s
    
    # Evolve star to target age using a temporary SeBa worker.
    seba = SeBa()
    try:
        seba.parameters.metallicity = 0.02
        seba.particles.add_particles(star)
        seba.evolve_model(target_age)

        # Copy evolved stellar properties.
        for _attribute in seba.particles.get_attribute_names_defined_in_store():
            seba.particles.new_channel_to(star).copy_attributes([_attribute])
    finally:
        seba.stop()

    # Store evolved star in hydro
    tag = hydro.add_particles(star.x, star.y, star.z)
    creation_time = hydro.get_time() - star.age
    hydro.set_particle_mass(tag, star.mass)
    hydro.set_particle_velocity(tag, star.vx, star.vy, star.vz)
    hydro.set_particle_oldmass(tag, star.mass)
    hydro.set_particle_creation_time(tag, creation_time)
    hydro.set_particle_rel_age(tag, star.relative_age)
    hydro.set_particle_rel_mass(tag, star.relative_mass)
    hydro.set_particle_corem(tag, star.core_mass)
    hydro.set_particle_co_corem(tag, star.COcore_mass)
    hydro.set_particle_stype(tag, star.stellar_type.value_in(units.stellar_type))
    hydro.set_particle_radius(tag, star.radius)

    tprint("Initialized sne test problem. The following star has been added:")
    for attribute, unit in zip(
        ["age", "mass", "radius", "temperature", "luminosity", "stellar_type"],
        [units.Myr, units.Msun, units.Rsun, units.K, units.Lsun, None]):
        if unit is None:
            tprint(attribute+": ", getattr(star[0], attribute), unit)
        else:
            tprint(attribute+": ", getattr(star[0], attribute).value_in(unit), unit)

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
