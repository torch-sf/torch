#!/usr/bin/env python

### Gravity bridge implementation for
### the Flash MHD code and a N-body solver.

### Joshua Wall
### Drexel University

# ASSUME:
#  . . . not first step . . .
#  . . . we have particles . . .
# with_multiples = False
# with_se = False
# with_sn = False
# with_winds = False
# with_massloss = False
#
# STRICTLY TRYING TO UNDERSTAND THE GRAVITY/PARTICLE BRIDGE CODE NOW
# AARON TRAN 2019 MARCH 26

max_hy_steps = hydro.get_max_num_steps() # max number of iterations.
curr_hy_step = hydro.get_current_step()

hydro.set_particle_pointers('mass')

while ((t < tmax) and (curr_hy_step < max_hy_steps)):
    i = i + 1

    ### Check for proper bridge timestep based on hydro timestep and crossing time.
    ### Have to write a proper routine to get timestep from hydro.
    hy_dt = hydro.get_timestep()
    dt = min(dtmax, 1.5*hy_dt, (tmax-t), 2.0*dt_old)

    dt_old = dt

    # (AT) if FLASH created a new sink prtl that's "massive enough",
    # allocate a new massive prtl.  Josh uses hydro.set_particle_pointers
    # to toggle between sink and massive particles in FLASH
    made_stars, made_massive_star = make_stars_from_sinks2(hydro, min_sf_mass, max_sf_mass)

    tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
    print "Num particles in grav:", len(grav.particles)
    print "Num particles in tree:", len(tree.particles)
    num_particles = check_particles

    hydro_time = hydro.get_time()
    grav_time  = grav.get_time()

    print "Starting the gravity bridge."

    stars_to_grav.copy()
    hydro.set_particle_mass(tags_keys[:,0], stars.mass)

    # Set the bridge timestep.
    t = t + dt
    print "I'm about to evolve hydro and grav for :" , dt, "to evolve to t =", t

    if (with_bridge):

        print "First kick."
        kick_number = 1
        if (gridChanged):
            kick_number = 2

        hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # this gets gravity AND kicks gas velocities
        print "Grid kicked."
        hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # calculate gravity at star locations
        stars.velocity = hydro.get_particle_velocity(tags_keys[:,0])  # update gravity code with kicked velocity
        print "Stars kicked."

        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
        print "Grav updated."

    print "Evolving models."

    grav.evolve_model(t)

    hydro.evolve_model(t)
    if (with_bridge):

        print "Second kick."
        kick_number = 2

        hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # this gets gravity AND kicks gas velocities
        print "Grid kicked."
        hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # calculate gravity at star locations
        stars.velocity = hydro.get_particle_velocity(tags_keys[:,0])  # update gravity code with kicked velocity
        print "Stars kicked."

        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
        print "Grav updated."

