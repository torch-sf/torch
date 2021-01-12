"""
Utility functions to handle protoplanetary disks
"""

def  add_particles_to_grav_and_ppds(state, hydro, grav, mult, se, ppds):
    """
    Send prtl from hydro to grav + AMUSE + ppds

    This gets called in two cases
    1. restarting with prtl,
    2. immediately after making new stars from sinks

    Separating hydro->grav update from sink->amuse->star->hydro update
    allows for possibility of hydro creating its own stars.

    Treat hydro as main fount of knowledge; copy stars from hydro to
    AMUSE and grav.

    se (SeBa or other stellar evolution worker) is only used to get
    correct stellar type for restarts; newborn stars are assumed to be on ZAMS

    postcondition:
        stars updated
        grav updated
    """
    add_parts_restart = False
    num_new_parts = hydro.get_number_of_new_tags()

    if num_new_parts > 0:

        newtags = hydro.get_new_tags(range(1,num_new_parts+1))

    else:

        tprint("add_particles_to_grav: assuming restart because Flash reports no new particles!")
        tprint("add_particles_to_grav: sync all stars from Flash to grav.")
        add_parts_restart = True
        num_new_parts = hydro.get_number_of_particles()
        newtags = hydro.get_particle_tags(range(1,num_new_parts+1))

    newtags.sort()

    position = hydro.get_particle_position(newtags)
    velocity = hydro.get_particle_velocity(newtags)
    mass     = hydro.get_particle_mass(newtags)
    initMass = hydro.get_particle_oldmass(newtags)

    # Make AMUSE particles for grav code.
    add_star = Particles(num_new_parts)
    add_star.mass = mass
    add_star.x    = position[:,0]
    add_star.y    = position[:,1]
    add_star.z    = position[:,2]
    add_star.vx   = velocity[:,0]
    add_star.vy   = velocity[:,1]
    add_star.vz   = velocity[:,2]

    add_star.tag  = newtags  # AMUSE stars know their FLASH tags
    add_star.stellar_type = 1 | units.stellar_type # ZAMS star
    add_star.radius = 0.02 | units.pc # initial collision radius, to handle tidal truncations
    add_star.initial_mass = initMass # for SE/SN uses
# don't need to carry this around because we don't need history
# just update directly in hydro
    #if with_lyc:
#    add_star.nion = 0.0 | units.s**-1 # ionizing flux
#    add_star.eion = 0.0 | units.erg # ionizing energy *OVER* 13.6 eV
#    add_star.sigh = 0.0 | units.cm**2 # ionizing cross section.
    #if with_pe_heat:
#    add_star.npe   = 0.0 | units.s**-1 # PE photon flux
#    add_star.epe   = 0.0 | units.erg # PE photon energy (should be around 8 eV)
#    add_star.sigpe = 0.0 | units.cm**2 # dust cross section per hydrogen atom
    #if with_wind:
#    add_star.dm_dt = 0.0 | units.g/units.s
#    add_star.vterm = 0.0 | units.cm/units.s

    # fast-forward stellar evolution to get current stellar type, because
    # torch_sf looks for change in stellar type to decide when to deposit SN
    if add_parts_restart:
        t_evol = hydro.get_time() - hydro.get_particle_creation_time(newtags)
        # TODO hardcoded solar metallicity Z=0.02 should be chosen by user.  -AT, 2019oct14
        _tmp = se.evolve_star(add_star.initial_mass, t_evol, 0.02)
        se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp
        add_star.stellar_type = se_type

    # only used by ph4... without this, ph4 complains about reused user IDs
    add_star.id = state.stars_next_id + np.arange(num_new_parts)
    state.stars_next_id += num_new_parts

    state.stars.add_particles(add_star)
    state.stars = state.stars.sorted_by_attribute('tag')

    grav.particles.add_particles(add_star)

    if mult is not None:
        mult._inmemory_particles.add_particles(add_star)
        # Multiples module needs an "id" attribute for internal book-keeping.
        # AMUSE example scripts set "id" directly; we use "index_in_code".
        mult.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    if ppds is not None:
        # Add stars to ppd code
        ppds.add_star_particles(add_star)

        # Assign total mass
        state.stars.total_mass = ppds.star_particles.total_mass

    if add_parts_restart:
        hydro.set_starting_local_tag_numbers()

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    return
