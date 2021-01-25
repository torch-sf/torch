"""
Utility functions to handle protoplanetary disks
"""

def make_and_add_stars_with_ppds (state, hydro, grav, se, ppds, sink_rad=None):
    '''
    Replacement of make_stars_from_sinks and add_particles_to_grav for runs with
    ppds. Note that this is only for new stars and disks; for old stars and disks,
    use reinitialize_ppds.
    '''

    assert sink_rad is not None  # required kwarg

    formed_stars = False

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()
    if num_sinks == 0:
        # can't get sink tags w/ empty list so need to exit early
        hydro.set_particle_pointers('mass')
        return formed_stars

    add_star = Particles(0)

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))  # does not work with empty list
    sink_tags.sort()  # is this necessary?

    for sink_tag in sink_tags:

        hydro.set_particle_pointers('sink')
        sink_mass = hydro.get_particle_mass(sink_tag)
        sink_pos = hydro.get_particle_position(sink_tag)
        sink_vel = hydro.get_particle_velocity(sink_tag)
        sink_cs  = hydro.get_sink_mean_cs(sink_tag)

        # get all the stars that we can form now
        # disk material must also be present! gas, and 1% dust
        csum = np.cumsum(state.all_masses[sink_tag] + \
            initial_disk_mass(state.all_masses[sink_tag]*1.01)
        i = np.searchsorted(csum, sink_mass.value_in(units.MSun), side='left')
        assert i < len(csum)  # ensure csum[-1] = sum(queue) > sink_mass

        spawn_masses = state.all_masses[sink_tag][:i]
        nnew = len(spawn_masses)

        if nnew == 0:

            tprint("... sink tag {} did not spawn stars".format(sink_tag))

        elif np.isnan(sink_cs.value_in(units.cm/units.s)):

            tprint("... sink tag {} blocked from spawning".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" due to absence of nearby cold gas")

        else:

            tprint("... sink tag {} spawned".format(sink_tag), end='')
            print(" {:d} stars,".format(nnew), end='')
            print(" total mass {:.2f},".format(np.sum(spawn_masses)), end='')
            print(" max mass {:.2f}".format(np.amax(spawn_masses)))

            formed_stars = True

            star          = Particles(nnew)
            star.mass     = spawn_masses | units.MSun
            # Isothermal spherical distribution.
            star.position = sink_pos + sink_rad*np.random.rand(nnew,1)*random_three_vector(nnew)
            # Gaussian distribution satisfying <vx**2> = sink_cs**2
            # so that stars' specific energy 1/2 <v**2> = (3/2)*sink_cs**2
            # matches gas specific energy P/rho/(gamma-1) for gamma=5/3
            # with cs = sqrt(P/rho) from Particles_sinkCreateAccrete.F90
            star.velocity = sink_vel + (np.random.normal(scale=sink_cs.value_in(units.cm/units.s), size=(nnew,3)) | units.cm/units.s)

            start = len(ppds.star_particles)
            ppds.add_star_particles(star)
            total_disk_mass = ppds.star_particles[start:].disk_mass.sum()

            # Remove newly-created stars from sink's queue
            state.all_masses[sink_tag] = state.all_masses[sink_tag][nnew:]

            # Remove the mass from the sink.
            sink_mass = sink_mass - (np.sum(spawn_masses)|units.MSun) - total_disk_mass

            # Prescribed disk mass is not exactly actual disk mass because
            # discretization, which might lead to negative or zero sink masses.
            # The error is small though, so we readd it to the sink and remove it
            # from the last disk -MW
            if sink_mass < 1e-6 | units.MSun:
                ppds.disks[-1].evaporate_mass( (1e-6|units.MSun) - sink_mass )
                sink_mass = 1e-6 | units.MSun
            hydro.set_particle_mass(sink_tag, sink_mass)

            star.total_mass = ppds.star_particles[start:].total_mass

            # Create new stars in FLASH
            hydro.set_particle_pointers('mass')
            star_tag = hydro.add_particles(star.x, star.y, star.z)
            hydro.set_particle_mass(star_tag, star.total_mass)
            hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
            hydro.set_particle_oldmass(star_tag, star.total_mass)

            add_star.add_particles(star)

    # if we made no stars, need to reset pointers
    hydro.set_particle_pointers('mass')


    num_new_parts = hydro.get_number_of_new_tags()
    newtags = hydro.get_new_tags(range(1,num_new_parts+1))

    newtags.sort()

    add_star.tag  = newtags  # AMUSE stars know their FLASH tags
    add_star.stellar_type = 1 | units.stellar_type # ZAMS star
    add_star.radius = 0.02 | units.pc # initial collision radius
    add_star.initial_mass = add_star.mass # for SE/SN uses

    # only used by ph4... without this, ph4 complains about reused user IDs
    add_star.id = state.stars_next_id + np.arange(num_new_parts)
    state.stars_next_id += num_new_parts

    state.stars.add_particles(add_star)
    state.stars = state.stars.sorted_by_attribute('tag')

    start = len(grav.particles)
    grav.particles.add_particles(add_star)
    grav.particles[start:].mass = add_star.total_mass

    hydro.clear_new_tags()

    return


def add_particles_to_grav_and_ppd (state, hydro, grav, se, ppds):

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
    add_star.radius = 0.02 | units.pc # initial collision radius
    add_star.initial_mass = initMass # for SE/SN uses


    if add_parts_restart:
        # fast-forward stellar evolution to get current stellar type, because
        # torch_sf looks for change in stellar type to decide when to deposit SN
        t_evol = hydro.get_time() - hydro.get_particle_creation_time(newtags)
        # TODO hardcoded solar metallicity Z=0.02 should be chosen by user.  -AT, 2019oct14
        _tmp = se.evolve_star(add_star.initial_mass, t_evol, 0.02)
        se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp
        add_star.stellar_type = se_type


    else:

        ppds.add_star_particles(add_star)

    add_star.total_mass = ppds.star_particles.total_mass

    hydro.set_particle_mass(newtags, add_star.total_mass)
    hydro.set_particle_oldmass(newtags, add_star.total_mass)


    # only used by ph4... without this, ph4 complains about reused user IDs
    add_star.id = state.stars_next_id + np.arange(num_new_parts)
    state.stars_next_id += num_new_parts

    state.stars.add_particles(add_star)
    state.stars = state.stars.sorted_by_attribute('tag')

    grav.particles.add_particles(add_star)
    grav.particles.mass = add_star.total_mass

    if add_parts_restart:
        hydro.set_starting_local_tag_numbers()

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    return


def reinitialize_ppds (hydro, ppd_index, rad_field_method, num_viscous_workers):

    if rad_field_method == 'rad_trans':

        ppds = restart_population(hydro.get_output_dir(), ppd_index, 
            number_of_workers=num_viscous_workers, grid_hydro=hydro)

    elif rad_field_method == 'geometric':

        ppds = restart_population(hydro.get_output_dir(), ppd_index, 
            number_of_workers=num_viscous_workers)

    return ppds
