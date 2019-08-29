"""
Torch module for managing code state.  This includes...

* input: loading simulation from restart files

* output: dumping simulation state

* synchronizing: doing major operations to push stars around between hydro,
  grav, AMUSE particle set, etc...



Joshua Wall, Drexel University.
State load/write logic re-organized by Aaron Tran, aaron.tran@columbia.edu,
2019 July 01.
"""

import os
import pickle

from amuse.io import write_set_to_file


class TorchIOState(object):
    """Holds RNG state and other information for Torch code.
    Handles all file I/O, including restart.
    """
    def __init__(self, hydro, stars, mult_grav, pdir='.', refresh=False):
        """Initialize Torch state,
        if restart, load stuff from file...

        IMPORTANT: if not restart, need to call this AFTER FLASH has already
        dumped first chk/plt files.
        """

        self.hydro = hydro
        self.stars = stars
        self.mult_grav = mult_grav
        self.pdir = pdir
        self.refresh = False

        self.all_masses = {}  # formerly global var in bridge_multiples.py

        # TODO enhancement - read from FLASH's own RuntimeParameter interface,
        # instead of duplicating the flash.par file parsing and default case
        # behavior.  Needs new interface code, see hydro.get_output_dir()
        # (AT) - 2019 July 01
        self.restart = False  # flash defaults
        self.chknum = 0
        self.pltnum = 0
        self.basename = "flash_"

        with open("flash.par") as f:
            for line in f:
                words = line.strip().split()
                if words[0] == '#' or len(words) <= 1:
                    # Comment or unparse-able
                    continue

                if (words[0].lower() == 'restart'):
                    if (words[2].lower() in ['.true.', 't']):
                        self.restart = True
                if (words[0].lower() == 'checkpointfilenumber'):
                    self.chknum = words[2]
                if (words[0].lower() == 'plotfilenumber'):
                    self.pltnum = words[2]
                if (words[0].lower() == 'basenm'):
                    self.basename = words[2].strip("\"")

        self.output_dir = hydro.get_output_dir()

        # flash logic is as follows.
        # restart:
        #   load chk_{chknum}
        #   ...
        #   write chk_{chknum+1}
        #   write plt_{pltnum+1}
        #   ...
        # not restart:
        #   write chk_{chknum}
        #   write plt_{pltnum}
        #   ...
        #   write chk_{chknum+1}
        #   write plt_{pltnum+1}
        #   ...
        #
        # mirror flash logic in our writing/loading of RNG state.
        # and tracking of output mtimes

        # restart: should get chk mtime but not plt mtime
        # non-restart: should get both chk and plt mtime
        self.chk_mtime = -1
        self.plt_mtime = -1
        chk = self.next_chk_file()
        plt = self.next_plt_file()
        if os.path.isfile(chk):
            self.chk_mtime = os.stat(chk).st_mtime
        if os.path.isfile(plt):
            self.plt_mtime = os.stat(plt).st_mtime

        if restart:

            if self.refresh:

                print "WARNING: Refreshing random state with a new seed."

            else:

                rstatefile = os.path.join(output_dir,
                    'rnd_state{:04d}.pickle'.format(self.chknum))
                massesfile = os.path.join(output_dir,
                    'all_masses{:04d}.pickle'.format(self.chknum))

                # if restart and not refresh, demand that files exist;
                # force user to set refresh=T if RNG files are missing.
                assert os.path.isfile(rstatefile)
                assert os.path.isfile(massesfile)

                with open(rstatefile, 'r') as f:
                    rnd_state = pickle.load(f)
                np.random.set_state(rnd_state)
                print "Random state set with file # "+rstatefile

                with open(massesfile, 'r') as f:
                    self.all_masses = pickle.load(f)
                print "Loaded all_masses dictionary from file # "+massesfile

        else:  # not restart

            out_rnd()
            out_mass()
            # no prtl sets yet if not restart, so don't need to write

        # All done with restart loading or non-restart file write,
        # prepare to write new chk/plt files.

        self.chknum = self.chknum + 1
        self.pltnum = self.pltnum + 1


    def next_chk_file():
        return os.path.join(self.output_dir,
            "{:s}hdf5_chk_{:04d}".format(self.basename, self.chknum))


    def next_plt_file():
        return os.path.join(self.output_dir,
            "{:s}hdf5_plt_cnt_{:04d}".format(self.basename, self.pltnum))


    def out(write_psets=False):
        """Write full Torch state to disk"""
        hydro.IO_out('chk')
        hydro.IO_out('pltpart')

        chk = self.next_chk_file()  # just-written files
        plt = self.next_plt_file()

        # Check whether FLASH output anything.  Must check mtime; chk/plt file
        # existence is necessary but not sufficient to prove write success.
        # Old files oft linger in restarted runs.
        # Using short-circuiting 'and' logic here.
        if (os.path.isfile(chk) and os.stat(chk).st_mtime > self.chk_mtime):

            self.out_rnd()
            self.out_mass()

            self.chknum = self.chknum + 1
            self.chk_mtime = os.stat(chk).st_mtime

        if (os.path.isfile(plt) and os.stat(plt).st_mtime > self.plt_mtime):

            if write_psets:
                self.out_psets()

            self.pltnum = self.pltnum + 1
            self.plt_mtime = os.stat(plt).st_mtime


    def out_psets():
        """Write star particle sets to AMUSE-format files"""
        stars_file = os.path.join(self.pdir,
            "stars{:04d}.amuse".format(self.pltnum))
        mult_file = os.path.join(self.pdir,
            "mult{:04d}.amuse".format(self.pltnum))
        multstars = mult_grav.stars.copy_to_new_particles()
        write_set_to_file(stars,    stars_file)
        write_set_to_file(multstars, mult_file)


    def out_rnd():
        """Write current random number state to pickle"""
        fname = os.path.join(self.output_dir,
            "rnd_state{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(np.random.get_state(), f)


    def out_mass():
        """Write dict with all future stars to pickle"""
        fname = os.path.join(self.output_dir,
            "all_masses{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.all_masses, f)


def add_particles_to_grav(state, mult_grav, stars, tree_exists):
    """
    This gets called in two cases
    1. restarting with prtl,
    2. immediately after making new stars from sinks

    Separating hydro->grav update from sink->amuse->star->hydro update
    allows for possibility of hydro creating its own stars.

    there is no point in pretending this method is at all functional in
    paradigm, it has loads of side effects.

    postcondition:
        stars updated
        tree updated
        state.all_masses updated, IFF we reinitialized
        mult_grav channel setup
    """
    # NOTE: state argument is only used for reinitialize_particles...
    # maybe we can get rid of it somehow? I dunno
    # -AT, 2019 July 01

    add_parts_restart = False

    num_new_parts = hydro.get_number_of_new_tags()

    if num_new_parts == 0:
        print "Assuming this is a restart, since Flash reports no new particles!"
        print "Fetching number of particles total as number of new particles for grav."
        num_new_parts = hydro.get_number_of_particles()
        newtags = hydro.get_particle_tags(range(1,num_new_parts+1))
        add_parts_restart = True
    else:
        newtags = hydro.get_new_tags(range(1,num_new_parts+1))

    newtags.sort()

    position = hydro.get_particle_position(newtags)
    velocity = hydro.get_particle_velocity(newtags)
    mass     = hydro.get_particle_mass(newtags)
    initMass = hydro.get_particle_oldmass(newtags)
    age      = hydro.get_time() - hydro.get_particle_creation_time(newtags)

    # Make AMUSE particles for grav code.
    add_star = Particles(num_new_parts)

    for kk in range(num_new_parts):     # Add the properties from hydro to
                                        # the new particles. This will be added to
        add_star[kk].mass= mass[kk]     # stars and the grav code.
        add_star[kk].age = age[kk]
        add_star[kk].x   = position[kk][0]
        add_star[kk].y   = position[kk][1]
        add_star[kk].z   = position[kk][2]
        add_star[kk].vx  = velocity[kk][0]
        add_star[kk].vy  = velocity[kk][1]
        add_star[kk].vz  = velocity[kk][2]
        add_star[kk].tag = newtags[kk]  # AMUSE stars know their FLASH tags
        add_star[kk].stellar_type = 1 |units.stellar_type # ZAMS star.
        add_star[kk].radius = 100 | units.AU # set all the stars initial collision radius.
        if (initMass[kk].value_in(units.MSun) <= 0.001):
            add_star[kk].initial_mass = mass[kk] # record the initial mass of the star for SE/SN uses.
            hydro.set_particle_oldmass(newtags[kk], mass[kk])
        else:
            add_star[kk].initial_mass = initMass[kk] # record the initial mass of the star for SE/SN uses.

        # Initialize radiation parameters to zero if this
        # is a radiation run.
        if use_radiation:
            add_star[kk].nion = 0.0 # ionizing flux
            add_star[kk].eion = 0.0 # ionizing energy *OVER* 13.6 eV
            add_star[kk].sigh = 0.0 # ionizing cross section.

    stars.add_particles(add_star)
    stars = stars.sorted_by_attribute('tag')

    if tree_exists:
        tree.particles = stars.sorted_by_attribute('tag')

    ### I wonder if we can get away without removing all the particles
    ### in the gravity code now that we've switched to ph4?? It would
    ### make multiples integration easier if we don't first remove all particles
    ### then immediately readd them all...
    grav.particles.add_particles(add_star)

    if mult_grav is not None and not add_parts_restart:
        mult_grav._inmemory_particles.add_particles(add_star)
        mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")

    if with_ph4 and add_parts_restart:
        grav.commit_particles()

    if add_parts_restart:
        hydro.set_starting_local_tag_numbers()
        if (clear_particles_on_restart):
            reinitialize_all_particles_from_stars(state, stars, hydro, grav)

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    if add_star.mass.value_in(units.MSun).any() > min_mass.value_in(units.MSun):
        hydro.new_source_flag(True)

    return


def reinitialize_all_particles_from_stars(state, stars, hydro, grav):
    """
    This function first deletes all the particles in
    hydro, then uses the stars particles array to
    put all the particles back into hydro (which
    will now all have new and proper tags correctly
    numbered for the current processors). Finally
    it removes all particles from grav and stars. Finally
    it reinitializes the particles in grav and stars
    with the new ones in hydro.
    """

    # First lets deal with the sink particles
    hydro.set_particle_pointers('sink')
    # Record all the relevant sink info into a particle set.
    ns             = hydro.get_number_of_particles()
    sinks          = Particles(ns)
    st             = hydro.get_particle_tags(ns)
    sinks.tag      = st
    sinks.mass     = hydro.get_particle_mass(st)
    sinks.position = hydro.get_particle_position(st)
    sinks.velocity = hydro.get_particle_velocity(st)
    # Remove all the sinks from hydro.
    hydro.remove_all_particles()
    # Now readd them.
    for i, sink in enumerate(sinks):
        new_tag = hydro.add_particles(sink.x, sink.y, sink.z)
        hydro.set_particle_mass(new_tag, sink.mass)
        hydro.set_particle_velocity(new_tag, sink.vx, sink.vy, sink.vz)
        state.all_masses[new_tag] = state.all_masses[sink.tag]

    # Now lets do the regular particles.
    # First remove all particles from hydro.
    hydro.set_particle_pointers('mass')
    hydro.remove_all_particles()

    # Now add all the particles back to hydro. Note I'm going to do
    # this with a loop (one at a time) so I definitely get the proper
    # new tag for each star that is in stars in correct order.
    t_hy = hydro.get_time()
    for star in stars:

        new_tag = hydro.add_particles(star.x, star.y, star.z)
        hydro.set_particle_mass(new_tag, star.mass)
        hydro.set_particle_oldmass(new_tag, star.initial_mass)
        hydro.set_particle_velocity(new_tag, star.vx, star.vy, star.vz)
        hydro.set_particle_creation_time(new_tag, star.age + t_hy)
        # Update the tag for this star in stars array.
        star.tag = new_tag

    # Sort the stars to be in the same order as the tags.
    stars = stars.sorted_by_attribute('tag')

    # Remove the improperly sorted particles and readd
    # the properly sorted ones.
    grav.particles.remove_particles(grav.particles)
    grav.particles.add_particles(stars)

    return


def remove_particles_outside_bndbox(hydro, stars, grav, mult_grav, with_multiples, bndbox, debug=True):

    ### Remove any particles that have left the simulation.

    stars_removed = False
    rem_index = 0
    num_particles = len(stars)

    ### x-direction check ###

    if (num_particles > 0):
        rem_index = np.where(
            np.abs(grav.particles.x.value_in(units.cm)) >= bndbox.value_in(units.cm)
        )[0]
        rem_key   = grav.particles.key[rem_index]]
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    rem_size = rem_index.size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        grav_rem_part.add_particles(grav.particles[grav_rem_index])

        for st in grav_rem_part:

            if st in stars:
                st_index = np.where(stars.key == st.key)[0]
                stars_rem_part.add_particles(stars[st_index])

            if (with_multiples):
                print "Trying to look in multiples for this star."
                # Cycle through the leaves and check if this particular star
                # is a multiple root particle.
                if st in mult_grav.root_to_tree:
                    print "Found it in multiples as a root particle."
                    tree = mult_grav.root_to_tree[st]
                    leaves = tree.get_leafs_subset()
                    print "Removing", len(leaves), "particles that are multiples by deleting the root of the tree."
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]

        rem_tag.append(stars_rem_part.tag)

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)

        num_particles = len(stars)

    ### y-direction check ###

    if (num_particles > 0):
        rem_index = np.where(np.abs(grav.particles.y.value_in(units.cm)) >= bndbox.value_in(units.cm))[0]
        rem_key   = grav.particles.key[rem_index]
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    rem_size = rem_index.size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        grav_rem_part.add_particles(grav.particles[grav_rem_index])

        for st in grav_rem_part:

            if st in stars:
                st_index = np.where(stars.key == st.key)[0]
                stars_rem_part.add_particles(stars[st_index])

            if (with_multiples):
                print "Trying to look in multiples for this star."
                # Cycle through the leaves and check if this particular star
                # is a multiple root particle.
                if st in mult_grav.root_to_tree:
                    print "Found it in multiples as a root particle."
                    tree = mult_grav.root_to_tree[st]
                    leaves = tree.get_leafs_subset()
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]

        rem_tag.append(stars_rem_part.tag)

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)

        num_particles = len(stars)

    ### z-direction check ###

    if (num_particles > 0):
        rem_index = np.where(np.abs(grav.particles.z.value_in(units.cm)) >= bndbox.value_in(units.cm))[0]
        rem_key   = grav.particles.key[rem_index]
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    rem_size = rem_index.size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        grav_rem_part.add_particles(grav.particles[grav_rem_index])

        for st in grav_rem_part:

            if st in stars:
                st_index = np.where(stars.key == st.key)[0]
                stars_rem_part.add_particles(stars[st_index])

            if (with_multiples):
                print "Trying to look in multiples for this star."
                # Cycle through the leaves and check if this particular star
                # is a multiple root particle.
                #st_grav = st.as_particle_in_set(grav.particles)
                #for root, tree in mult_grav.root_to_tree.iteritems():
                    #leaves = tree.get_leafs_subset()
                    #if st in leaves:
                if st in mult_grav.root_to_tree:
                    print "Found it in multiples as a root particle."
                    tree = mult_grav.root_to_tree[st]
                    leaves = tree.get_leafs_subset()
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]

        rem_tag.append(stars_rem_part.tag)

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)

        num_particles = len(stars)

    # Sync stars to the particles in grav.
    if (with_multiples):
        grav.particles.synchronize_to(mult_grav._inmemory_particles)
    else:
        grav.particles.synchronize_to(stars)

    return stars_removed


###############################################
###             Multiples module            ###
###############################################

SMALLN = None
def init_smalln(converter):
    global SMALLN
    SMALLN = SmallN(convert_nbody=converter)
    SMALLN.initialize_code()
def new_smalln():
    global SMALLN
    SMALLN.reset()
    return SMALLN
def stop_smalln():
    global SMALLN
    SMALLN.stop()


def initialize_multiples(stars, grav, conv, mult_debug_level=1, kep=None, new_smalln=None):

    print "Starting multiples codes."
    grav.parameters.epsilon_squared = 0.0 | units.cm**2.0
    print "Gravity epsilon set to zero."
    stopping_condition = grav.stopping_conditions.collision_detection
    stopping_condition.enable()
    print "Stopping condition enabled."

    if (kep is None):
        print "Starting Kepler."
        kep = Kepler(unit_converter=conv)
        print "Initializing Kepler."
        kep.initialize_code()

    print "Starting multiples."
    multiples_code = multiples.Multiples(grav, new_smalln, kep, constants.G)
    multiples_code.global_debug                = mult_debug_level
    multiples_code.neighbor_veto               = True
    multiples_code.check_tidal_perturbation    = True
    multiples_code.neighbor_perturbation_limit = 0.05
    multiples_code.wide_perturbation_limit     = 0.08
    mult_to_stars = multiples_code.stars.new_channel_to(stars)
    stars_to_mult = stars.new_channel_to(multiples_code.stars)
    print "Multiples initialized."
    return multiples_code, mult_to_stars, stars_to_mult


def initialize_gravity_codes(convert, stars = None,
                             num_grav_workers = 1,
                             eps = 15.0 | units.RSun,
                             with_ph4 = True, tree_exists = False,
                             with_multiples = False):
    print "Starting gravity code."

### NOTE! Trying to use the generic unit converter here after I saw this on the AMUSE website
### amusecode.org/doc/reference/quantities_and_units.html

    if with_ph4:
        grav = ph4(convert, number_of_workers=num_grav_workers, mode='cpu', redirection="none")
        grav.parameters.set_defaults()
        grav.parameters.timestep_parameter=0.14  # Timestep accuracy for PH4
        grav.parameters.force_sync=1 # Force the code to end exactly at the specified time.
    else:
        grav = Hermite(convert, number_of_workers=num_grav_workers)
        grav.parameters.end_time_accuracy_factor=0.0 #1e-8
        grav.parameters.dt_param=0.02 #0.14**2  # Timestep accuracy for Hermite

    # N-body softening radius is the actual radius of a large massive star here.
    if with_multiples:
        grav.parameters.epsilon_squared = 0.0 | units.cm**2.0
    else:
        grav.parameters.epsilon_squared = eps**2.0

    # Use particle set to pass information back and forth
    # between Flash and the gravity code as well as do the bridge kicks.
    stars_to_grav = stars.new_channel_to(grav.particles)
    grav_to_stars = grav.particles.new_channel_to(stars)
    return grav, stars_to_grav, grav_to_stars


def update_roots_from_leaves(mult_grav, grav):

    '''
    Update the center of mass particles from
    the leaves properties (in all codes!).
    '''
    for root, tree in mult_grav.root_to_tree.iteritems():
        root_particle = root.as_particle_in_set(mult_grav._inmemory_particles)
        leaves = tree.get_leafs_subset()
        com     = 0.0
        com_vel = 0.0

        for leaf in leaves:
            com += (leaf.mass.value_in(units.g)*leaf.position.value_in(units.cm))
            com_vel += (leaf.mass.value_in(units.g)*leaf.velocity.value_in(units.cm/units.s))
        com     = (com/(leaves.mass.sum()).value_in(units.g)) | units.cm
        com_vel = (com_vel/(leaves.mass.sum()).value_in(units.g)) | units.cm/units.s
        msum    = leaves.mass.sum()
        # Update the root particle in multiples.
        root_particle.mass     = msum
        root_particle.position = com
        root_particle.velocity = com_vel
        # Also update the tree.particle thing (top of the tree in the dictionary in multiples).
        tree.particle.mass     = msum
        tree.particle.position = com
        tree.particle.velocity = com_vel
        # To be consistent the same particle in the N body code must also be updated!
        grav_particle = root.as_particle_in_set(grav.particles)
        grav_particle.mass     = msum
        grav_particle.position = com
        grav_particle.velocity = com_vel

    return



if __name__ == '__main__':
    pass
