### Gravity bridge implementation for
### the Flash MHD code and a N-body solver.

### Joshua Wall
### Drexel University

print "At the top"

from amuse.lab import *
import numpy as np
from amuse.community.flash.interface import Flash
#from amuse.rfi.channel import AsyncRequestsPool
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.community.flash import josh_multiples as multiples
#from amuse.community.flash import steve_multiples as multiples
import sys
import time
import datetime
from amuse.community.sse.interface import SSE
#import main_sequence as ms
import ionizingflux as ion
from scipy.integrate import *
import glob
import pickle
import os

print "After imports."

np.set_printoptions(precision=3)

# Global variables

all_masses = {}
first_call_for_stars = False
old_sink_tags = []
mult_grav      = None

# Timer object from http://www.huyng.com/posts/python-performance-analysis
# Allows easy timing of chunks of python.

class Timer(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print 'elapsed time: %f ms' % self.msecs


class logfile:
    def __init__(self, filename="profiler.log"):
        self.filename=filename
        self.new_run()
    def new_run(self):
        self.write(120*"-")
        self.write("New run initiated.")
        self.write(120*"-")
        self.write(datetime.datetime.now().ctime())
        self.write("\n")
    def write(self, output):
        self.out = open(self.filename, 'a')
        self.out.write(output+'\n')
        self.out.close()
    def __del__(self):
        self.write(120*"-")
        self.write("End of run.")
        self.write(120*"-")
        self.write("\n")

# Function to determine the stellar luminosity at a particular wavelength, temperature and cross section.
# Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
# Note I left out sig0 here b/c we divide this by lum_wl_cs_per_ph that would
# also have sig0 in it.
def lum_wl_cs(l, l_max, T):

    h = 6.6261e-27 # Plank's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzman constant

    L = (2.0*h*(c**2.0)/(l**5.0)) * (l/l_max)**3.0 / (np.exp(h*c/(l*k*T)) - 1.0)

    return L
# Function to determine the number count of photons at a particular wavelength, temperature and cross section.
# Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
def lum_wl_cs_per_ph(l, l_max, T):

    h = 6.6261e-27 # Plank's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzman constant

    L = (2.0*h*(c**2.0)/(l**5.0)) * (l/l_max)**3.0 / (np.exp(h*c/(l*k*T)) - 1.0) / (h*c/l)

    return L


# Function to determine the stellar luminosity at a particular wavelength and temp.
# Uses the standard blackbody curve.
def lum_wl(l, l_max, T):

    h = 6.6261e-27 # Plank's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzman constant

    L = (2.0*h*(c**2.0)/(l**5.0)) / (np.exp(h*c/(l*k*T)) - 1.0)

    return L
# Function to determine the number count of photons at a particular wavelength and temp.
# Uses the standard blackbody curve.
def lum_wl_per_ph(l, l_max, T):

    h = 6.6261e-27 # Plank's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzman constant

    L = (2.0*h*(c**2.0)/(l**5.0)) / (np.exp(h*c/(l*k*T)) - 1.0) / (h*c/l)

    return L

### A function to move particles from Flash to AMUSE.

def load_rnd_state_files(restart, chknum, refresh_rand_seed_on_restart):

    global all_masses
    global first_call_for_stars
    global old_sink_tags
    # If this is a restart, shouldn't we load the random number state?

    rstatefile =  output_dir+'/rnd_state{:04d}.pickle'.format(chknum)
    massesfile = output_dir+'/all_masses{:04d}.pickle'.format(chknum)

    refresh_rand_seed_on_restart = (refresh_rand_seed_on_restart or
                                    not os.path.isfile(rstatefile) or
                                    not os.path.isfile(massesfile))

    if (restart and not refresh_rand_seed_on_restart):

            print rstatefile
            with open(rstatefile, 'r') as f:
                rnd_state = pickle.load(f)
            np.random.set_state(rnd_state)
            print "Random state set with file # "+rstatefile

            print massesfile
            with open(massesfile, 'r') as f:
                all_masses = pickle.load(f)
            print "Loaded all_masses dictionary from file # "+massesfile
            print "all_masses =", all_masses

            hydro.set_particle_pointers('sink')
            num_sinks     = hydro.get_number_of_particles()
            if (num_sinks > 0): old_sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))
            hydro.set_particle_pointers('mass')

            first_call_for_stars = False


    else: print "WARNING: Refreshing random state with a new seed. Pray to the R-N-Gesus."



def add_particles_to_grav(tags_keys, stars, tree_exists, newtags=None):

    global mult_grav

    # How to add stars.id from multiples test?!
    global stars_current_id_num

    add_parts_restart = False

    if (newtags is None):

        num_new_parts = hydro.get_number_of_new_tags()

        print "[add_particles_to_grav]: Number of new tags =", num_new_parts

        if (num_new_parts == 0):

            print "Assuming this is a restart, since Flash reports no new particles!"
            print "Fetching number of particles total as number of new particles for grav."
            num_new_parts = hydro.get_number_of_particles()
            newtags = hydro.get_particle_tags(range(1,num_new_parts+1))
            add_parts_restart = True

        else:

            newtags = hydro.get_new_tags(range(1,num_new_parts+1))


    # Sort the new tags here because sorting is currently broken in Fortran... (fix later?)
    newtags.sort()
    # Check the newtags to make sure they don't already exist in the tags array.
    if (np.shape(tags_keys)[0] > 0 and (not add_parts_restart)): test_tags(tags_keys, check_tags=newtags)

    position = hydro.get_particle_position(newtags)          # Get properties from hydro.
    velocity = hydro.get_particle_velocity(newtags)
    mass     = hydro.get_particle_mass(newtags)
    initMass = hydro.get_particle_oldmass(newtags)
    age      = hydro.get_time() - hydro.get_particle_creation_time(newtags)
    add_star = Particles(num_new_parts)             # Make new particles for grav code.
    keys     = add_star.key                         # Get the new keys for tags_keys.

    #print "tags =", newtags
    #print "pos =", position
    #print "vel =", velocity
    #print "mass =", mass
    #print "initMass =", initMass
    #print "age =", age

    for kk in range(num_new_parts):       # Add the properties from hydro to
                                          # the new particles. This will be added to
        add_star[kk].mass= mass[kk]       # stars and the grav code.
        add_star[kk].age = age[kk]
        add_star[kk].x   = position[kk][0]
        add_star[kk].y   = position[kk][1]
        add_star[kk].z   = position[kk][2]
        add_star[kk].vx  = velocity[kk][0]
        add_star[kk].vy  = velocity[kk][1]
        add_star[kk].vz  = velocity[kk][2]
        add_star[kk].tag = newtags[kk]
        add_star[kk].stellar_type = 1 |units.stellar_type # Start life as a ZAMS star.
        add_star[kk].radius = 100 | units.AU # set all the stars initial collision radius.
        if (initMass[kk].value_in(units.MSun) <= 0.001):
            add_star[kk].initial_mass = mass[kk] # record the initial mass of the star for SE/SN uses.
            hydro.set_particle_oldmass(newtags[kk], mass[kk])
        else:
            add_star[kk].initial_mass = initMass[kk] # record the initial mass of the star for SE/SN uses.
        stars_current_id_num += 1
        add_star[kk].id = stars_current_id_num

        if (debug_aptg):
            print "New star:"
            print "pos  =", add_star[kk].position.in_(units.cm)
            print "vel  =", add_star[kk].velocity.in_(units.km/units.s)
            print "mass =", add_star[kk].mass.in_(units.MSun)
            print "age  =", add_star[kk].age.in_(units.Myr)
            print "rad  =", add_star[kk].radius.in_(units.AU)
            print "type =", add_star[kk].stellar_type
            print "tag  =", add_star[kk].tag
            print "id   =", add_star[kk].id
    # Initialize radiation parameters to zero if this
    # is a radiation run.

        if (use_radiation):

            add_star[kk].nion = 0.0 # ionizing flux
            add_star[kk].eion = 0.0 #eion #0.0 # ionizing energy *OVER* 13.6 eV
            add_star[kk].sigh = 0.0 #sigh #0.0 # ionizing cross section.

    #print "Hydro code pos/vel/mass:"
    #print position
    #print velocity
    #print mass

    #print "Does the new tag match the new entry in tags?"
    #print newtags, tags[ind]

    if (np.shape(tags_keys)[0] == 0):


        tags_keys = np.zeros((num_new_parts, 2))
        tags_keys[:,0] = newtags
        tags_keys[:,1] = keys

    else:
        #print tags_keys
        #print [[newtags, keys]]
        jj = 0

        for ii in range(num_new_parts):

            tags_keys = np.append(tags_keys, [[newtags[ii], keys[ii]]], axis=0)

    # Now that we've modified the tags_keys array, check for repeated tags.
    test_tags(tags_keys)

    #grav.particles.remove_particles(stars)
    stars.add_particles(add_star)

    # Sort tags_keys so that we always have the tags in order.
    tags_keys = tags_keys[tags_keys[:,0].argsort()]
    # Sort the stars to be in the same order as the tags.
    stars = stars.sorted_by_attribute('tag')

    #print "Stars.id =", stars.id

    # Okay, I think I finally sorted this out. stars.sorted_by_attribute
    # returns a particle set, so it is okay to use this to set things like
    # grav.particles = stars.sorted_by_attribute('tag') even if you cant
    # sort grav.particles itself by tag (I think it just ignores things
    # like this in the copying over to grav.particles).

    if (tree_exists):
        #tree.particles.add_particles(add_star)
        # Sort these (or maybe copy from stars?)
        # Note neither of these methods gets the order right,
        # and we can't just set tree.particles = stars as that
        # changes tree.particles from a code to a pointer to stars.
        # so we just remove and then readd all the stars (after sorting them).
        #stars_to_tree.copy()
        #tree.particles = stars
        #tree.particles = tree.particles.sorted_by_attribute('tag')
        tree.particles = stars.sorted_by_attribute('tag')
    #stars_to_grav.copy()

    # Sort these (or maybe copy from stars?)
    # Note neither of these methods gets the order right,
    # and we can't just set grav.particles = stars because this
    # makes grav.particles a pointer to stars (and no longer a code
    # pointer to Hermite). So for now we just remove and then
    # readd all the stars particles (after sorting stars).
    #stars_to_grav.copy()
    #grav.particles = grav.particles.sorted_by_attribute('tag')
    #grav.particles = stars
    ### I wonder if we can get away without removing all the particles
    ### in the gravity code now that we've switched to ph4?? It would
    ### make multiples integration easier if we don't first remove all particles
    ### then immediately readd them all...
    #grav.particles.remove_particles(grav.particles)
    #grav.particles.add_particles(stars)
    grav.particles.add_particles(add_star)
    if (mult_grav is not None and not add_parts_restart):
        mult_grav._inmemory_particles.add_particles(add_star)
        mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
    #grav.particles = stars.sorted_by_attribute('tag')
    if (with_ph4 and add_parts_restart): grav.commit_particles()
    #if (with_ph4 and add_parts_restart): grav.commit_particles()
    #if (with_se):
        ##se.particles.add_particles(add_star)
        ##stars_to_se.copy()
        ##se.particles = se.particles.sorted_by_attribute('tag')
        ##se.particles = stars
        #se.particles.remove_particles(stars)
        #se.particles.add_particles(stars)
        #se.particles = stars.sorted_by_attribute('tag')

    # For now, set the same ionizing energy and cross section
    # for all stars in the hydro code.
    #hydro.set_particle_nion(newtags, 0.0 | units.s**-1.0)
    #hydro.set_particle_eion(newtags, eion)
    #hydro.set_particle_sigh(newtags, sigh)
    #if (with_winds):
    #    hydro.set_particle_wind_mass(newtags, 0.0 | units.g/units.s)
    #print "Does the new tag match the new entry in tags_keys?"
    #print newtags, tags_keys

    #print "Grav code pos/vel/mass:"
    #print grav.particles.position
    #print grav.particles.velocity
    #print grav.particles.mass

    # If this is a restart, clear all the hydro tags and reinitialize.
    # This fixes any problems with tags from using a different number
    # of processors on a restart.

    # Actually lets just try and sort out the proper local_tag_number
    # on each processor (the integer that keeps up with how many particles
    # each proc has made so far).

    if (add_parts_restart): hydro.set_starting_local_tag_numbers()

    if (add_parts_restart and clear_particles_on_restart): reinitialize_all_particles_from_stars(stars, hydro, grav, tags_keys)

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

    print "hydro.clear_new_tags called, num new particles is ", hydro.get_number_of_new_tags()

    if (add_star.mass.value_in(units.MSun).any() > min_mass.value_in(units.MSun)):
        hydro.new_source_flag(True)

    return tags_keys, stars

def reinitialize_all_particles_from_stars(stars, hydro, grav, tags_keys):

    # This function first deletes all the particles in
    # hydro, then uses the stars particles array to
    # put all the particles back into hydro (which
    # will now all have new and proper tags correctly
    # numbered for the current processors). Finally
    # it removes all particles from grav and stars,
    # deletes the current tags_keys array. Finally
    # it reinitializes the particles in grav and stars
    # hydro with the new ones in hydro.


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
        all_masses[new_tag] = all_masses[sink.tag]

    # Now lets do the regular particles.
    # First remove all particles from hydro.

    hydro.set_particle_pointers('mass')
    print "Before reinitializing tags_keys and stars.tags, they look like:"
    print tags_keys[:,0].astype(int)
    print stars.tag.astype(int)

    print "Hydro now reporting", hydro.get_number_of_particles(), "number of particles", \
          "before removal."

    hydro.remove_all_particles()

    # Did it work?

    print "Hydro now reporting", hydro.get_number_of_particles(), "number of particles", \
          "after removal."

    # Now add all the particles back to hydro. Note I'm going to do
    # this with a loop (one at a time) so I definitely get the proper
    # new tag for each star that is in stars in correct order.

    print len(tags_keys[:,0])

    t_hy = hydro.get_time()
    for i, star in enumerate(stars):

        new_tag = hydro.add_particles(star.x, star.y, star.z)
        hydro.set_particle_mass(new_tag, star.mass)
        hydro.set_particle_oldmass(new_tag, star.initial_mass)
        hydro.set_particle_velocity(new_tag, star.vx, star.vy, star.vz)
        hydro.set_particle_creation_time(new_tag, star.age + t_hy)
        # Update the tag for this star in stars array.
        star.tag   = new_tag
        # Update the tag in tags_keys array.
        tags_keys[i,0] = new_tag

    # Now that we sorted out the hydro stuff, lets sort out the particles
    # in the N body code and stars so both are sorted by tag.
    print len(tags_keys[:,0])
    # Sort tags_keys so that we always have the tags in order.
    #tags_keys = tags_keys[tags_keys[:,0].argsort()]
    tags_keys.sort(axis=0)

    # Did we break anything?
    test_tags(tags_keys)

    # Sort the stars to be in the same order as the tags.
    stars = stars.sorted_by_attribute('tag')

    # Remove the improperly sorted particles and readd
    # the properly sorted ones.
    grav.particles.remove_particles(grav.particles)
    grav.particles.add_particles(stars)
    #if (with_ph4): grav.commit_particles() #ph4 is tempermental.

    print "After reinitializing tags_keys and stars.tags, they look like:"
    print tags_keys[:,0].astype(int)
    print stars.tag.astype(int)

    # Done!
    return

def kick_grid(hydro, lim, tree, eps, dt):

    lim3 = lim**3

    # axis variables to pass for cell coords.
    a = np.empty(lim)
    a.fill(1)
    b = np.empty_like(a)
    b.fill(2)
    c = np.empty_like(a)
    c.fill(3)

    # variables for locations of each cell
    x = np.zeros(lim) | units.cm
    y = np.zeros(lim) | units.cm
    z = np.zeros(lim) | units.cm

    # variables to hold gravity for each cell in each direction

    #grav_x = np.zeros((8,8,8)) | units.cm * units.s**(-2)
    #grav_y = np.zeros((8,8,8)) | units.cm * units.s**(-2)
    #grav_z = np.zeros((8,8,8)) | units.cm * units.s**(-2)

    # Get the number of leaf grids and there indices on each processor.
    # Note that block indices are not unique across processors, only intra-
    # processor. But get_leaf_indices returns them in order and with an
    # array of the number on each processor (block_array).
    all_grids = hydro.get_number_of_grids()
    leaf_grids = np.zeros((all_grids*lim3))
    block_array = np.zeros((all_grids*lim3))
    [leaf_grids, block_array, num_leafs] = hydro.get_leaf_indices(range(all_grids))
    #print num_leafs
    numblks=num_leafs[0]
    # Note we can only pass 1-d flattened arrays using the AMUSE interface,
    # and all of them must be the same length within a function. Therefore
    # we also pass how long an array really is with it.
    leaf_grids = np.resize(leaf_grids,numblks*lim3)
    block_array = np.resize(block_array, numblks*lim3)
    #num_leafs=num_leafs[0]
    #numblks=2
    #leaf_grids=np.append(leaf_grids,np.zeros(num_leafs*numblks-all_grids))
    #print len(leaf_grids)
    #print leaf_grids[:numblks]
    #print len(leaf_grids)
    #print numblks
    #print block_array
    #sys.exit()

    grav_x = np.zeros((numblks,lim,lim,lim)) | units.cm * units.s**(-2)
    grav_y = np.zeros((numblks,lim,lim,lim)) | units.cm * units.s**(-2)
    grav_z = np.zeros((numblks,lim,lim,lim)) | units.cm * units.s**(-2)

    grav_x2 = np.zeros((numblks,lim3)) | units.cm * units.s**(-2)
    grav_y2 = np.zeros((numblks,lim3)) | units.cm * units.s**(-2)
    grav_z2 = np.zeros((numblks,lim3)) | units.cm * units.s**(-2)

    ii = 0


    # Loop over each block and get the gravity at each cell coordinate
    # from the tree code (much faster than getting it by direct summation).
    for leaf in leaf_grids[:numblks]:

    #def get_grav_block(leaf, grav_x, grav_y, grav_z):

        limits = hydro.get_grid_range(leaf)

        #print limits
        #print leaf
        # Cell coordinates (x,y,z).
        x = hydro.get_1blk_cell_coords(a,leaf,lim)
        y = hydro.get_1blk_cell_coords(b,leaf,lim)
        z = hydro.get_1blk_cell_coords(c,leaf,lim)

        #print len(x), len(y), len(z)
        # Meshgrid them to make a cube in the same shape as the block
        # we are looping over.
        [xx,yy,zz] = np.meshgrid(x.value_in(units.cm),y.value_in(units.cm),z.value_in(units.cm),indexing='ij') | units.cm
        # Now flatten this because again AMUSE can only take 1-d arrays as
        # arguments. Ravel=flatten. C-ordering is used here explicitly.
        xx = np.ravel(xx.value_in(units.cm),'C') | units.cm
        yy = np.ravel(yy.value_in(units.cm),'C') | units.cm
        zz = np.ravel(zz.value_in(units.cm),'C') | units.cm

        #print len(xx)
        # Actually get the gravity accel from the tree code now.
        [grav_x2[ii,:],grav_y2[ii,:],grav_z2[ii,:]] = tree.get_gravity_at_point(eps, xx, yy, zz)
        # Reshape the gravity arrays to be the same shape we want to reconstruct in Flash.
        # We used C ordering before because that's what reshape uses by default.
        grav_x[ii,:,:,:] = np.reshape(grav_x2[ii].value_in(units.cm * units.s**-2),(8,8,8)) | units.cm * units.s**-2
        grav_y[ii,:,:,:] = np.reshape(grav_y2[ii].value_in(units.cm * units.s**-2),(8,8,8)) | units.cm * units.s**-2
        grav_z[ii,:,:,:] = np.reshape(grav_z2[ii].value_in(units.cm * units.s**-2),(8,8,8)) | units.cm * units.s**-2

        #print grav_x[ii,:,:,:]

        ii += 1
    #print grav_x[1,4,1,2]
    #print grav_x[:,:,:,:]
    # Now flatten again to pass it to Flash. We'll have to reconstruct in Flash using
    # C ordering (last index changes fastest).
    grav_x = (grav_x.value_in(units.cm * units.s**-2)).flatten('C') | units.cm * units.s**-2
    grav_y = (grav_y.value_in(units.cm * units.s**-2)).flatten('C') | units.cm * units.s**-2
    grav_z = (grav_z.value_in(units.cm * units.s**-2)).flatten('C') | units.cm * units.s**-2

    #print len(grav_x)
    # Pass the gravity to Flash.
    hydro.kick_block(grav_x,grav_y,grav_z,leaf_grids, block_array, lim, dt)
    #print "Kick grid done!"

def bridge_kick(hydro, stars, eps, dt, time_in_hydro):

    grid_limits = 8 # num of cells in 1-d in one block.

    loc = grav.particles.position

    #gaccel = hydro.get_gravity_at_point(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])
    gaccel = hydro.get_accel_gas_on_particles(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])

    print "Got gravity at star locations."
    #print gaccel
    #print "Before kick star vel = "
    #print stars.velocity

    for k in range(len(grav.particles)):

        stars[k].vx = stars[k].vx + 0.5*dt*gaccel[0][k]
        stars[k].vy = stars[k].vy + 0.5*dt*gaccel[1][k]
        stars[k].vz = stars[k].vz + 0.5*dt*gaccel[2][k]


    ### Update gravity code with kicked velocity.

    print "Stars kicked."
    #print stars.velocity

    #print "Grav vel before update."
    #print grav.particles.velocity
    stars_to_grav.copy()

    print "Grav updated."
    #print grav.particles.velocity

    hydro.set_particle_velocity(tags_keys[:,0], grav.particles.vx, grav.particles.vy, grav.particles.vz)

    ### Second kick to hydro gas grid.
    with Timer(verbose=True) as grid_kick_timer:
    #kick_grid(hydro, grid_limits, tree, eps, 0.5*dt)
        hydro.kick_grid(0.5*dt)
    time_in_hydro = time_in_hydro + grid_kick_timer.secs
    #hydro_timer = hydro_timer + grid_kick_timer.secs
    print "Grid kicked."

def bridge_kick2(hydro, stars, eps, dt, time_in_hydro, kick_number):

    ### Second kick to hydro gas grid.
    #with Timer(verbose=True) as grid_kick_timer:
    ##kick_grid(hydro, grid_limits, tree, eps, 0.5*dt)
        #hydro.kick_grid(0.5*dt)
    #time_in_hydro = time_in_hydro + grid_kick_timer.secs
    ##hydro_timer = hydro_timer + grid_kick_timer.secs

    # Lets try kicking the grid first, since it updates the potential
    # everywhere on the grid during the kick, which is later used to
    # calculate the gas_on_particles gravity.

    # Note this function gets the gravity AND kicks the velocities of the gas.
    hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)

    print "Grid kicked."

    #loc = grav.particles.position

    #gaccel = hydro.get_gravity_at_point(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])
    #gaccel = hydro.get_accel_gas_on_particles(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])

    #print "Before kick star vel = "
    #print stars.velocity

    # Calculate the gravity at the star locations.

    hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)

    # Get the gravity at these points (stored in ACCX_PART_PROP, etc)

    #gaccel = hydro.get_particle_acceleration(tags_keys[:,0])

    #print "Got gravity at star locations."
    #print gaccel


    #for k in range(len(grav.particles)):

        #stars[k].vx = stars[k].vx + 0.5*dt*gaccel[k][0]
        #stars[k].vy = stars[k].vy + 0.5*dt*gaccel[k][1]
        #stars[k].vz = stars[k].vz + 0.5*dt*gaccel[k][2]

    ### Update gravity code with kicked velocity.

    vel = hydro.get_particle_velocity(tags_keys[:,0])
    stars.velocity = vel

    print "Stars kicked."
    #print stars.velocity
    return


def stellar_properties(stellar_mass = 1.0|units.MSun, metallicity = 0.02):
    print "Starting SSE."
    #stellar_evolution = SSE()
    stellar_evolution = SeBa(redirection = "none")
    #print "Refresh memory."
    #stellar_evolution.refresh_memory()
    print "Setting metallicity."
    stellar_evolution.parameters.metallicity = metallicity
    print "Making particle."
    star = Particle()
    print "Setting mass."
    star.mass = stellar_mass
    print "Adding star to stellar evo code."
    star = stellar_evolution.particles.add_particle(star)
    print "Doing stellar evolution."
    stellar_evolution.evolve_model(1.0|units.yr)
    print "Saving values."
    radius = stellar_evolution.particles.radius[0] #star.radius
    temperature = stellar_evolution.particles.temperature[0] #star.temperature
    luminosity = stellar_evolution.particles.luminosity[0] #star.luminosity
    print "Stopping SSE."
    stellar_evolution.stop()
    #print "Clean up code."
    #stellar_evolution.cleanup_code()
    return radius, temperature, luminosity

### Make a fractal cluster that is contained within a bounding box
### (generally the hydro box size) with an AMUSE particle set.

def make_cluster(conv_cluster, nm_part, bndbox, fractal=False, equal_mass=False, eq_mass=1.0 | units.MSun):

    from amuse.community.fractalcluster.interface import new_fractal_cluster_model

    stars_out=True
    n=0

    while (stars_out==True):

        n=n+1

        if (fractal):
            print "Making star cluster using fractal model."
            cluster = new_kroupa_mass_distribution(nm_part, mass_max = (100.0 |units.MSun))
            cluster = new_fractal_cluster_model(masses=cluster,convert_nbody=convert, do_scale=False, virial_ratio=5.0)
        else:
            print "Making star cluster using plummer model."
            cluster = new_plummer_sphere(nm_part, convert_nbody=convert, do_scale=False)

        if (equal_mass):
            cluster.mass = eq_mass
        else:
            #cluster.mass = new_kroupa_mass_distribution(nm_part, mass_max = (100.0 |units.MSun))
            cluster.mass = new_salpeter_mass_distribution(nm_part, mass_min = (0.1 | units.MSun), mass_max = (100.0 |units.MSun))

        remove_stars = cluster.select(lambda r: bndbox < max(abs(r)), ["position"])

        if (len(remove_stars) > 0):

            stars_out=True

        else:

            stars_out=False

    print "This took", n, "runs."
    print cluster.mass[np.where(cluster.mass.value_in(units.MSun) > 7.0)].value_in(units.MSun)

    return cluster

### Place a single star into the hydro simulation.

def make_single_star_in_hydro(x, y, z, mass, initMass = 0.0 | units.MSun, age = 0.0 | units.Myr,
                              vx=0.0 | units.cm/units.s, vy=0.0 | units.cm/units.s, vz=0.0 | units.cm/units.s):
    if (initMass.value_in(units.MSun) == 0.0): initMass = mass
    creation_time = hydro.get_time() - age
    tag = hydro.add_particles(x,y,z)
    hydro.set_particle_velocity(tag, vx, vy, vz)
    hydro.set_particle_mass(tag, mass)
    hydro.set_particle_oldmass(tag, initMass)
    hydro.set_particle_creation_time(tag, creation_time)

    return tag

### Place an entire cluster into the hydro simulation (after it was
### made in AMUSE.

def make_cluster_in_hydro(cluster, initial_x=0.0 | units.cm, initial_y=0.0 | units.cm, initial_z=0.0 | units.cm):

    x = cluster.x + initial_x; y = cluster.y + initial_y; z = cluster.z + initial_z

    print "Are stars outside the bndbox?"

    print np.where(np.any(x>bndbox)), np.where(np.any(y>bndbox)), np.where(np.any(z>bndbox))

    tag = hydro.add_particles(x,y,z)
    print "Length of tag =", len(tag)
    print "tag =", tag
    hydro.set_particle_velocity(tag, cluster.vx, cluster.vy, cluster.vz)
    hydro.set_particle_mass(tag, cluster.mass)
    #print hydro.get_particle_mass(tag)

    #hydro.particles_gather()

    return tag

def test_tags(tk, check_tags=None):

    check_passed = True
    if (check_tags is None):
        print "Checking tags_keys for any repeated tags."
        tag_ind = 0

        for tag in tk[:,0]:

            masked_tag = np.ma.array(tk[:,0], mask=False)
            masked_tag.mask[tag_ind] = True

            if ((np.equal(tag, masked_tag)).any()):
                print int(tag), masked_tag[np.where(np.equal(tag,masked_tag))].astype(int)
                check_passed = False
                break
            tag_ind += 1

    else:
        print "Checking tags_keys for the following tags:"
        #print check_tags
        tags_from_tk = tk[:,0]
        for tag in check_tags:

            if ((np.equal(tag, tags_from_tk)).any()):
                print int(tag), tags_from_tk[np.where(np.equal(tag,tags_from_tk))].astype(int)
                check_passed = False
                break

    if (not check_passed):

        print "Repeated tag found!"
        print "Did you restart with a different number of processors? This breaks the tags \
               which are assigned by processor id."
        print tk.astype(int)
        sys.stdout.flush()
        sys.exit()

    return

def kroupa(m,a):

    if (0.001 <= m < 0.08):
        k = a*m**(-0.3)
    elif (0.08 <= m < 0.5):
        k = a*(0.08)*m**(-1.3)
    elif (0.5 <= m):
        k = a*(0.08*0.5)*m**(-2.3)
    else:
        print "Invalid mass range!"
        k=0
    return k

def mkroupa(m,a):

    if (0.001 <= m < 0.08):
        k = m*a*m**(-0.3)
    elif (0.08 <= m < 0.5):
        k = m*a*(0.08)*m**(-1.3)
    elif (0.5 <= m):
        k = m*a*(0.08*0.5)*m**(-2.3)
    else:
        print "Invalid mass range!"
        k=0
    return k

def m_max_star(m_max_clust):
    # The max stellar mass for sampling the "normal" IMF
    # calculated from Weidner et. al. 2013 eqn 1,
    # based on the integrated galatic IMF of Weidner and Kroupa 2004.

    # m_max_clust is the maximum cluster mass
    # and should figure in losses due to jets and other
    # feedback. Generally, I just assume a SFE of 0.5.

    a0 = -0.66
    a1 =  1.08
    a2 = -0.15
    a3 = 0.0084

    Lmclust = np.log10(m_max_clust)

    if (m_max_clust <= 2.5E5):
        m_max = a0 + a1*Lmclust + a2*Lmclust**2. + a3*Lmclust**3.0
    else:
        m_max = np.log10(150.0)

    return 10**m_max

#def sample_stars_poisson(sink_mass, M_min, M_max, num_bins):
    #'''Return a poisson random sampling from the Kroupa IMF of sink total
       #mass from M_min to M_max separated into num_bins in logspace.

       #Returns:
               #n_stars: Number of stars in each logarithmic bin
               #binsL:   The bin edges, including the right most bin edge
               #lam:     The average number of stars in each bin that the
                        #Poisson sample is centered around.'''

    #from scipy.integrate import quad

    #norm_inv = quad(kroupa,M_min,M_max,args=(1))[0]

    #norm = 1/norm_inv

    #binsL = np.logspace(np.log10(M_min),np.log10(M_max),num_bins+1)
    #mass_per_bin = []
    #frac_per_bin = []

    ##print binsL

    #avg_mass = quad(mkroupa, M_min, M_max, args=(norm))[0]

    ##print avg_mass

    #for i in range(num_bins):

        ## m_i = int_bin_i_low^bin_i_high(m*f*dm) / int(f*dm)
        #mass_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    #quad(kroupa, binsL[i], binsL[i+1], args=(norm))[0] )
        ## f_i = int_bin_i_low^bin_i_high(m*f*dm) / int_M_low^M_high(m*f*dm)
        #frac_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    #avg_mass )

    #mass_per_bin = np.array(mass_per_bin)
    #frac_per_bin = np.array(frac_per_bin)

    ##print mass_per_bin
    ##print sum(mass_per_bin)

    ##print frac_per_bin
    ##print sum(frac_per_bin)

    #lam = sink_mass*frac_per_bin/mass_per_bin

    #n_stars = np.random.poisson(lam=lam)

    #return n_stars, binsL, lam #frac_per_bin

## remaining_mass needs to survive from call to call
## in this function, so we make it global.

#remaining_mass = 0.0


#def get_stellar_mass_sampling(sample_imf_mass, num_bins=10, min_samp_mass=1.0, max_samp_mass=150.0, eff=1.0):

    #'''Return a random sampling of an IMF from a Mass_min to Mass_max
    #using a Poisson method to sample the number of stars in each bin and a
    #Salpeter distribution to choose the exact stellar masses in each bin
    #of the individual stars.'''

    ## Get some number of stars given a total mass and some number of bins.
    #num_bins = 10
    #eff = 1.0
    #[n_stars, bins, lam] = sample_stars_poisson(eff*sample_imf_mass.value_in(units.MSun),
                                                #min_samp_mass, max_samp_mass, num_bins)
    ## Now use that to sample the IMF.

    ## Now fill out the masses of the stars in each bin.

    #mass_in_each_bin = np.zeros((num_bins-1, np.max(n_stars)))

    #for b in range(num_bins-1):

        #if (n_stars[b] > 0):

            #mass_in_each_bin[b,0:n_stars[b]] = new_salpeter_mass_distribution(n_stars[b],
                           #mass_min= (bins[b]   | units.MSun),
                           #mass_max= (bins[b+1] | units.MSun), alpha=-2.3).value_in(units.MSun)

    #all_samp_masses = np.ravel(mass_in_each_bin)
    #all_samp_masses = all_samp_masses[all_samp_masses!=0.0]

    ## Here we move all the stars smaller than 1 MSun into particles
    ## that are at least 1 MSun. To do this we do a bit of fancy
    ## footwork with the arrays.

    #small_masses = all_samp_masses[np.where(all_samp_masses < 1.0)] # Smaller than 1.0 MSun.
    #all_samp_masses = all_samp_masses[np.where(all_samp_masses >= 1.0)]  # Everyone else.

    #b = 0
    ## If there are any left smaller than 1.0 MSun, sum with others
    ## that are smaller than 1.0 MSun until there are none left.

    #if (len(small_masses) > 1):
        #while(small_masses[b] < 1.0 and len(small_masses[b:])>=1):

            #small_masses[b] = small_masses[b]+small_masses[b+1]
            #b = np.delete(small_masses, b+1)
            #if (len(small_masses[b:]) > 1):
                #if(small_masses[b] >= 1.0 and small_masses[i+1] < 1.0):
                    #b += 1

    ## If the last one is smaller than 1.0 MSun, lump that bit into
    ## the remaining mass after star formation.
    ##if (len(small_masses) >= 1):
        ##if (small_masses[-1] < 1.0):
            ##remaining_mass += small_masses[-1]
            ##small_masses = np.delete(small_masses, -1)

    #all_samp_masses = np.append(all_samp_masses, small_masses)

    ## Now randomly shuffle all the masses in the array.
    #np.random.shuffle(all_samp_masses)


    ##print all_samp_masses

    #all_stars_mass = np.sum(all_samp_masses)

    #print "Stars made are:", all_samp_masses
    #print "Total stellar mass is:", all_stars_mass

    #return all_samp_masses

def sample_stars_poisson(sink_mass, M_min, M_max, num_bins):
    '''Return a poisson random sampling from the Kroupa IMF of sink total
       mass from M_min to M_max separated into num_bins in logspace.

       Returns:
               n_stars: Number of stars in each logarithmic bin
               binsL:   The bin edges, including the right most bin edge
               lam:     The average number of stars in each bin that the
                        Poisson sample is centered around.
               norm:    Norm to be used to sample the Kroupa IMF
                        using the n_stars array.'''

    from scipy.integrate import quad

    norm_inv = quad(kroupa,M_min,M_max,args=(1))[0]

    norm = 1/norm_inv

    binsL = np.logspace(np.log10(M_min),np.log10(M_max),num_bins+1)
    mass_per_bin = []
    frac_per_bin = []

    #print binsL

    avg_mass = quad(mkroupa, M_min, M_max, args=(norm))[0]

    #print avg_mass

    for i in range(num_bins):

        # m_i = int_bin_i_low^bin_i_high(m*f*dm) / int(f*dm)
        mass_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    quad(kroupa, binsL[i], binsL[i+1], args=(norm))[0] )
        # f_i = int_bin_i_low^bin_i_high(m*f*dm) / int_M_low^M_high(m*f*dm)
        frac_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    avg_mass )

    mass_per_bin = np.array(mass_per_bin)
    frac_per_bin = np.array(frac_per_bin)

    #print mass_per_bin
    #print sum(mass_per_bin)

    #print frac_per_bin
    #print sum(frac_per_bin)

    lam = sink_mass*frac_per_bin/mass_per_bin

    n_stars = np.random.poisson(lam=lam)

    print "Mass from N stars ~", np.sum(n_stars*mass_per_bin)


    return n_stars, binsL, lam, norm #frac_per_bin

def collect_small_stars_mass(all_samp_masses):

    # Here we move all the stars smaller than 1 MSun into particles
    # that are at least 1 MSun. To do this we do a bit of fancy
    # footwork with the arrays.

    small_masses = all_samp_masses[np.where(all_samp_masses < 1.0)] # Smaller than 1.0 MSun.
    all_samp_masses = all_samp_masses[np.where(all_samp_masses >= 1.0)]  # Everyone else.

    b = 0
    # If there are any left smaller than 1.0 MSun, sum with others
    # that are smaller than 1.0 MSun until there are none left.

    #print "small_masses"
    #print small_masses

    if (len(small_masses) > 1):
        while(small_masses[-1] < 1.0 and len(small_masses[b:])>1):

            small_masses[b] = small_masses[b]+small_masses[-1]
            small_masses = np.delete(small_masses, -1)
            if (len(small_masses[b:]) > 1):
                if(small_masses[b] >= 1.0):
                    b += 1

        # If the last one is smaller than 1.0 MSun, lump that bit into
        # the last star.
        if (small_masses[-1] < 1.0):
            small_masses[-2] = small_masses[-2] + small_masses[-1]
            small_masses = np.delete(small_masses, -1)

    #print "small_masses"
    #print small_masses

    all_samp_masses = np.append(all_samp_masses, small_masses)

    return all_samp_masses

def get_stellar_mass_sampling(sample_imf_mass, num_bins=10, min_samp_mass=1.0, max_samp_mass=150.0, eff=1.0, sum_small=False):


    [n_stars, bins, lam, norm] = sample_stars_poisson(eff*sample_imf_mass.value_in(units.MSun),
                                                min_samp_mass, max_samp_mass, num_bins)

    print bins

    # Now use that to sample the IMF.
    masses = np.zeros(n_stars.sum())
    k = 0
    counter = 0
    for i,n in enumerate(n_stars):

        print "Pulling ", n, "stars from ranges ", bins[i], "to ", bins[i+1]
        print "counter =", counter
        counter = 0
        for j in range(n):

            while (masses[k] == 0):

                m = np.random.uniform(low=bins[i], high=bins[i+1])
                r = np.random.uniform()
                p = mkroupa(m, norm)

                if (p/r > 1.0): masses[k] = m

            k+=1
            counter+=1

    print "Just got a new sampling of the IMF from", min_samp_mass, "to", max_samp_mass, "."
    print "Masses          =", masses
    print "Number of stars =", len(masses)
    print "Total mass      =", masses.sum()
    print "Max mass        =", masses.max()

    # Sum all stars < 1 MSun into stars > 1 MSun.
    if (sum_small): masses = collect_small_stars_mass(masses)


    print "masses before shuffle."
    print masses

    np.random.shuffle(masses)

    print "masses after shuffle."
    print masses

    return masses

def make_stars_from_sinks(hydro, min_imf_mass):

    global remaining_mass

    formed_stars = False
    use_ang_mom  = False
    # Get the total mass in sink particles in the simulation

    # Note this requires us to move the particles pointer over to
    # the sink array, then MOVE IT BACK to the particles array.

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()

    print "Num sinks =", num_sinks

    if (num_sinks < 1):
        hydro.set_particle_pointers('mass')
        return

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))
    sink_tags.sort()
    print "sink_tags = ", sink_tags
    sink_masses = hydro.get_particle_mass(sink_tags)
    total_sink_mass = sink_masses.sum()

    # Efficiency (effective mass left over) assuming jets remove some mass.

    eff = 1.0

    eff_total_sink_mass = eff*total_sink_mass

    print "sink masses", sink_masses.as_quantity_in(units.MSun)

    print "total sink mass", total_sink_mass.as_quantity_in(units.MSun)

    if (eff_total_sink_mass.value_in(units.MSun) > min_imf_mass):

        #sink_mass_fractions = np.divide(sink_masses,total_sink_mass)
        #sink_mean_vel       = hydro.get_sink_gas_mean_velocity(sink_tags)
        #sink_var_vel        = hydro.get_sink_gas_var_velocity(sink_tags)
        sink_positions      = hydro.get_particle_position(sink_tags)
        sink_vel            = hydro.get_particle_velocity(sink_tags)

        if (use_ang_mom):
            sink_ang_mom        = hydro.get_sink_ang_mom(sink_tags)

        # This makes the loop work even in the case of a single sink.
        if (num_sinks == 1):

            #sink_mean_vel  = [sink_mean_vel]
            #sink_var_vel   = [sink_var_vel]
            sink_positions = [sink_positions]
            sink_vel       = [sink_vel]
            sink_tags      = [sink_tags]

            if (use_ang_mom):
                sink_ang_mom   = [sink_ang_mom]

        # Now lets randomize the sink tags indices so that we can start with
        # a different, random sink each time.

        sink_ind = np.arange(len(sink_tags))
        #sink_ind = np.where(sink_tags > 0.0)[0] # Get the indices for the tags array.

        np.random.shuffle(sink_ind) # Randomize them.

        print sink_ind
        print sink_masses[sink_ind[0]]
        print (sink_masses[sink_ind[0]]).value_in(units.MSun)
        print sink_masses[sink_ind[0]].value_in(units.MSun)


        num_bins = 10 #50

        all_stars_mass = 0.0

        for s in range(num_sinks):
            # If this sink has less than needed mass, cycle to the next sink.
            if (sink_masses[sink_ind[s]].value_in(units.MSun) < 1.0): continue
        # Poisson sample over the sink mass, obtaining the number of stars
        # in each mass bin.
            #print "sink position =", sink_positions
            print "sink position s =", sink_positions[sink_ind[s]]
            [n_stars, bins, lam] = sample_stars_poisson(eff*sink_masses[sink_ind[s]].value_in(units.MSun), 1.0, 150.0, num_bins)

            if (np.sum(n_stars) == 0): continue

        # Now fill out the masses of the stars in each bin.

            mass_in_each_bin = np.zeros((num_bins-1, np.max(n_stars)))

            for b in range(num_bins-1):

                if (n_stars[b] > 0):

                    mass_in_each_bin[b,0:n_stars[b]] = new_salpeter_mass_distribution(n_stars[b],
                                   mass_min= (bins[b]   | units.MSun),
                                   mass_max= (bins[b+1] | units.MSun), alpha=-2.3).value_in(units.MSun)

            all_masses = np.ravel(mass_in_each_bin)
            all_masses = all_masses[all_masses!=0.0]

            # Here we move all the stars smaller than 1 MSun into particles
            # that are at least 1 MSun. To do this we do a bit of fancy
            # footwork with the arrays.

            small_masses = all_masses[np.where(all_masses < 1.0)] # Smaller than 1.0 MSun.
            all_masses = all_masses[np.where(all_masses >= 1.0)]  # Everyone else.

            b = 0
            # If there are any left smaller than 1.0 MSun, sum with others
            # that are smaller than 1.0 MSun until there is at most one left.

            if (len(small_masses) > 1):
                while(small_masses[b] < 1.0 and len(small_masses[b:])>1):

                    small_masses[b] = small_masses[b]+small_masses[b+1]
                    b = np.delete(small_masses, b+1)
                    if (len(small_masses[b:]) > 1):
                        if(small_masses[b] >= 1.0 and small_masses[i+1] < 1.0):
                            b += 1

            # If the last one is smaller than 1.0 MSun, lump that bit into
            # the remaining mass after star formation.
            if (len(small_masses) >= 1):
                if (small_masses[-1] < 1.0):
                    remaining_mass += small_masses[-1]
                    small_masses = np.delete(small_masses, -1)

            all_masses = np.append(all_masses, small_masses)

            #print all_masses

            all_stars_mass = np.sum(all_masses)
            remaining_mass = sink_masses[sink_ind[s]].value_in(units.MSun) - all_stars_mass

            print "Sum of all masses", all_stars_mass
            print "remaining mass", remaining_mass

            # If we pulled too much mass, we'll delete stars until its correct.
            #if (remaining_mass < 0.0):
                # Randomly shuffle all the stellar masses.
            #    np.random.shuffle(all_masses)
                # Now remove stars until we are under the mass limit.
            #    while (remaining_mass < 0.0):
            #        all_masses = np.delete(all_masses, 0)
            #        all_stars_mass = np.sum(all_masses)
            #        remaining_mass = sink_masses[sink_ind[s]].value_in(units.MSun) - all_stars_mass

            # Update this sink mass to remove what we used to make stars.
            if (remaining_mass < 0.0):
                hydro.set_particle_mass(sink_tags[sink_ind[s]], (0.0 | units.MSun))
            else:
                hydro.set_particle_mass(sink_tags[sink_ind[s]], (remaining_mass | units.MSun))

            s1 = s

            # If we used more mass than was available from this sink, start borrowing forward from other sinks.
            # If we make it through all the sinks, we'll borrow forward from the next accretion onto sinks.
            # Note this is why we randomize the order that sinks are sampled from.
            while (remaining_mass < 0.0 and s1 < (num_sinks-1)):

                # Take mass from each sink we have left until either we meet the mass requirement or we run out of sinks.
                s1 = s1 + 1
                next_sink_mass = sink_masses[sink_ind[s1]].value_in(units.MSun)
                remaining_mass = remaining_mass + next_sink_mass

                if (remaining_mass < 0.0):
                    # If we dropped this sink's mass to zero, zero out the ang mom and mass appropriately.
                    sink_masses[sink_ind[s1]] = 0.0 | units.MSun
                    if (use_ang_mom):
                        sink_ang_mom[sink_ind[s1]][:] = 0.0 | units.cm**2.0 * units.g / units.s

                else:
                    # Reduce angular momentum by the mass fraction removed.
                    if (use_ang_mom):
                        sink_ang_mom[sink_ind[s1]][:] = remaining_mass / sink_masses[sink_ind[s1]].value_in(units.MSun) * sink_ang_mom[sink_ind[s1]][:]
                    # Reduce mass to what is left over.
                    sink_masses[sink_ind[s1]] = remaining_mass | units.MSun
                    # No need to take out anymore mass from other sinks, so reset remaining mass to zero.
                    remaining_mass = 0.0
                # Now actually set the sink mass and ang mom in the hydro code.
                hydro.set_particle_mass(sink_tags[sink_ind[s1]], sink_masses[sink_ind[s1]])
                if (use_ang_mom):
                    hydro.set_particle_ang_mom(sink_tags[sink_ind[s1]], sink_ang_mom[sink_ind[s1]][0], sink_ang_mom[sink_ind[s1]][1], sink_ang_mom[sink_ind[s1]][2])

            print "After removal total stars mass", all_stars_mass

            if (all_stars_mass <= 0.0): continue

            n_stars = len(all_masses)
            sum_mr2 = 0.0 | units.g*units.cm**2.0
            new_stars = Particles(n_stars)
            new_stars.mass = all_masses | units.MSun

            print "Total number of stars made =", n_stars
            print "Star masses =", new_stars.mass
            #for st in range(len(all_masses)):

                #print (random_three_vector()*(np.random.rand())**(1.0/3.0)*sink_rad).in_(units.cm) + sink_positions[sink_ind[s]]

            # Uniform distribution in the sink radius.
            #new_stars[:].position = (random_three_vector(n_stars)[:,:]*(np.random.rand(n_stars)**(1.0/3.0))[:,None]*sink_rad).in_(units.cm) + sink_positions[sink_ind[s]]

            # Singular isothermal spherical distribution.
            stars_rvec = (random_three_vector(n_stars)[:,:]*(np.random.rand(n_stars))[:,None]*sink_rad)
            print "stars_rvec=", stars_rvec
            rx = stars_rvec[:,0]
            ry = stars_rvec[:,1]
            rz = stars_rvec[:,2]
            print "rx, ry, rz =", rx, ry, rz
            r2           = (rx**2.0 + ry**2.0 + rz**2.0).in_(units.cm**2.0)
            new_stars[:].position = np.add(stars_rvec.value_in(units.cm),sink_positions[sink_ind[s]].value_in(units.cm)) | units.cm

            #print "star positions=",new_stars[:].position
            #print "star 1 position =",new_stars[0].position
            #print "star 1 x, y, z =", new_stars[0].x, new_stars[0].y, new_stars[0].z

            #print "sink position=", sink_positions[sink_ind[s]]


                # Plummer / BE sphere distribution in the sink radius.
                #u_vec = (np.array(random_three_vector())*(np.random.rand(3))*sink_rad) | units.cm + sink_positions[sink_ind[s]]

                #print sink_var_vel[sink_ind[s]]

                #v_vec = np.random.normal(loc=sink_mean_vel[sink_ind[s]], scale=sink_var_vel[sink_ind[s]])
            # There is an implicit assumption here that the velocity of the gas that is collapsing resembles the velocity of the gas
            # that is in the surrounding cells. I'm not too sure about this assumption actually...

            #print np.random.normal(loc=sink_vel[sink_ind[s]].value_in(units.cm/units.s), scale=np.sqrt(sink_var_vel[sink_ind[s]]).value_in(units.cm/units.s))

            new_stars[:].vx = sink_vel[sink_ind[s]][0].as_quantity_in(units.cm / units.s)
            new_stars[:].vy = sink_vel[sink_ind[s]][1].as_quantity_in(units.cm / units.s)
            new_stars[:].vz = sink_vel[sink_ind[s]][2].as_quantity_in(units.cm / units.s)

            print "new stars vx = ", new_stars[:].vx
            print "new stars vy = ", new_stars[:].vy
            print "new stars vz = ", new_stars[:].vz

            #new_stars[:].vx = 0.0 | units.cm / units.s
            #new_stars[:].vy = 0.0 | units.cm / units.s
            #new_stars[:].vz = 0.0 | units.cm / units.s

            #print "new stars vx = ", new_stars[:].vx
            #print "new stars vy = ", new_stars[:].vy
            #print "new stars vz = ", new_stars[:].vz
            #new_stars[:].vx = np.random.normal(loc=sink_vel[sink_ind[s]][0].value_in(units.cm/units.s),
            #                           scale=np.sqrt(sink_var_vel[sink_ind[s]][0].value_in(units.cm**2.0/units.s**2.0)),
            #                           size=n_stars) | units.cm / units.s
            #new_stars[:].vy = np.random.normal(loc=sink_vel[sink_ind[s]][1].value_in(units.cm/units.s),
            #                           scale=np.sqrt(sink_var_vel[sink_ind[s]][1].value_in(units.cm**2.0/units.s**2.0)),
            #                           size=n_stars) | units.cm / units.s
            #new_stars[:].vz = np.random.normal(loc=sink_vel[sink_ind[s]][2].value_in(units.cm/units.s),
            #                           scale=np.sqrt(sink_var_vel[sink_ind[s]][2].value_in(units.cm**2.0/units.s**2.0)),
            #                           size=n_stars) | units.cm / units.s


            # Calculate the rotational velocity omega from the original angular momentum of the sink (which came from the
            # infalling gas), reduced by the mass ratio of the new stars and the original sink.

            # Sum up the moments of inertia.

            if (use_ang_mom):
                sum_mr2 = (np.sum(new_stars[:].mass.in_(units.g) * r2)).as_quantity_in(units.g*units.cm**2.0)

                print "Stellar inertia =", sum_mr2
                mass_ratio = (np.sum(new_stars.mass.value_in(units.g)) / sink_masses[sink_ind[s]].value_in(units.g))

                print "sink ang mom =", sink_ang_mom[sink_ind[s]] * mass_ratio
                sink_ang_mag = (sink_ang_mom[sink_ind[s]].norm()).as_quantity_in(units.cm**2.0*units.g/units.s)
                print "mag of ang mom =", sink_ang_mag

                omega = (sink_ang_mom[sink_ind[s]] * mass_ratio / sum_mr2).as_quantity_in(units.s**-1)

                print "omega =", omega

            # Now add the velocity from r (x) omega to the velocity of each star.

                new_stars[:].vx = new_stars[:].vx + ry*omega[2] - rz*omega[1]
                new_stars[:].vy = new_stars[:].vy + rz*omega[0] - rx*omega[2]
                new_stars[:].vz = new_stars[:].vz + rx*omega[1] - ry*omega[0]

            #print "Star 0 has mass, position and velocity="

            #print new_stars[0].mass
            #print new_stars[0].position
            #print new_stars[0].velocity

            #print new_stars[:].y-sink_positions[sink_ind[s]][1]
            #print new_stars[:].vz - sink_vel[sink_ind[s]][2]

                lx = (np.sum(new_stars[:].mass.in_(units.g)*((ry).in_(units.cm)*(-new_stars[:].vz + sink_vel[sink_ind[s]][2]).in_(units.cm/units.s)
                                         - (rz).in_(units.cm)
                                         *(-new_stars[:].vy + sink_vel[sink_ind[s]][1]).in_(units.cm/units.s)))).as_quantity_in(units.cm**2.0*units.g/units.s)
                ly = (np.sum(new_stars[:].mass.in_(units.g)*((rz).in_(units.cm)*(-new_stars[:].vx + sink_vel[sink_ind[s]][0]).in_(units.cm/units.s)
                                         - (rx).in_(units.cm)
                                         *(-new_stars[:].vz + sink_vel[sink_ind[s]][2]).in_(units.cm/units.s)))).as_quantity_in(units.cm**2.0*units.g/units.s)
                lz = (np.sum(new_stars[:].mass.in_(units.g)*((rx).in_(units.cm)*(-new_stars[:].vy + sink_vel[sink_ind[s]][1]).in_(units.cm/units.s)
                                         - (ry).in_(units.cm)
                                         *(-new_stars[:].vx + sink_vel[sink_ind[s]][0]).in_(units.cm/units.s)))).as_quantity_in(units.cm**2.0*units.g/units.s)

                print "lx, ly, lz =", lx, ly, lz

                star_ang_mag = (np.sqrt(lx**2.0 + ly**2.0 + lz**2.0)).as_quantity_in(units.cm**2.0*units.g/units.s)

                print "Star total ang momentum from sink = ", star_ang_mag

                #ang_norm_factor = star_ang_mag / sink_ang_mag
                #lx = lx*ang_norm_factor
                #ly = ly*ang_norm_factor
                #lz = lz*ang_norm_factor
                #star_ang_mag = (np.sqrt(lx**2.0 + ly**2.0 + lz**2.0)).as_quantity_in(units.cm**2.0*units.g/units.s)

                print "Star total ang momentum from sink = ", star_ang_mag

                print "Before removing ang mom, sink ang mom =", hydro.get_sink_ang_mom(sink_tags[sink_ind[s]])

            # Don't forget to remove this angular momentum from the sink particle!

                hydro.set_particle_ang_mom(sink_tags[sink_ind[s]], sink_ang_mom[sink_ind[s]][0]-lx, sink_ang_mom[sink_ind[s]][1]-ly, sink_ang_mom[sink_ind[s]][2]-lz)

                print "Now sink ang mom =", hydro.get_sink_ang_mom(sink_tags[sink_ind[s]])


                #print "After adding angular momentum star 0 has velocity="

                #print new_stars[0].velocity

            print "new stars vx = ", new_stars[:].vx
            print "new stars vy = ", new_stars[:].vy
            print "new stars vz = ", new_stars[:].vz

            # Switch to massive particles to set properties.
            hydro.set_particle_pointers('mass')

            # Now all the stars have mass, position and velocity. So add them to the code.

            new_star_tags = hydro.add_particles(new_stars.x,new_stars.y,new_stars.z)

            print "new_star_tags =", new_star_tags

            new_star_tags.sort()

            print "new_star_tags =", new_star_tags

            hydro.set_particle_mass(new_star_tags, new_stars.mass)
            hydro.set_particle_velocity(new_star_tags,new_stars.vx,new_stars.vy,new_stars.vz)

            #print "New star tags are ", new_star_tags

            # Switch back to sinks to remove them if needed.
            #hydro.set_particle_pointers('sink')

            # If this sink is empty, get rid of it.
            #if (remaining_mass < 1.0):
            #    print "[make_stars_from_sink]: Removing sink", sink_tags[sink_ind[s]]
            #    hydro.remove_particles(sink_tags[sink_ind[s]])

            formed_stars = True

            # Save random number state.

            rnd_state = np.random.get_state()
            with open('rnd_state.pickle', 'wb') as f:
                pickle.dump(rnd_state, f)

    # We're done, so now switch back to massive star particles.
    hydro.set_particle_pointers('mass')

    return formed_stars


def make_stars_from_sinks2(hydro, min_imf_mass, max_imf_mass, sample_imf_mass=10000 | units.MSun,
                           local_sfe=1.0, sum_small=False):

    # Given an initial sampling of the IMF, distribute the stars randomly
    # as sinks accrete the required mass to form them.

    global first_call_for_stars
    global all_masses
    global old_sink_tags

    formed_stars = False
    formed_massive_star = False

    if (first_call_for_stars == True):

        print "Initializing all_masses and old_sink_tags arrays."
        all_masses = {} # emtpy dict. get_stellar_mass_sampling(sample_imf_mass)

        #star_ind = {} #0

        old_sink_tags = []

        first_call_for_stars = False

    # Now that we have a set of stars to pull from
    # we can check to see if a star should be formed.

    # Note this requires us to move the particles pointer over to
    # the sink array, then MOVE IT BACK to the particles array.

    hydro.set_particle_pointers('sink')
    num_sinks = hydro.get_number_of_particles()

    print "Num sinks =", num_sinks

    if (num_sinks < 1):
        hydro.set_particle_pointers('mass')
        return (formed_stars, formed_massive_star)

    sink_tags = hydro.get_particle_tags(range(1,num_sinks+1))
    sink_tags.sort()
    print "sink_tags = ", sink_tags
    print "len sink tags =", len(sink_tags)
    print "len old_sink_tags =", len(old_sink_tags)

    if (len(sink_tags) > len(old_sink_tags)):
        # Then we need to make a new list of star masses.
        #print "Inside setup for getting masses."

        if (len(old_sink_tags) > 0):
            tags_mask  = np.ones_like(sink_tags, dtype=bool)
            search_ind = np.searchsorted(sink_tags, old_sink_tags)
            tags_mask[search_ind] = False
            new_tags   = sink_tags[tags_mask]
        else:
            new_tags = sink_tags

        for new_tag in new_tags:
            #print "inside new_tag loop"
            # Make a new list of star masses for each using a dictionary.
            # Sample down to the low side of the Kroupa IMF.
            all_masses.update({new_tag:get_stellar_mass_sampling(sample_imf_mass,
                                num_bins=10,
                                min_samp_mass=min_imf_mass.value_in(units.MSun),
                                max_samp_mass=max_imf_mass.value_in(units.MSun),
                                eff=local_sfe, sum_small=sum_small)})
            # Instead of keeping up with the index, we could just pop / delete
            # the star from the list once we use it...
            #star_ind[new_tag]   = 0

        # make both dicts sorted
        #all_masses = collections.OrderedDict(sorted(all_masses.items()))
        #star_ind   = collections.OrderedDict(sorted(star_ind.items()))

        old_sink_tags = sink_tags

    sink_masses = hydro.get_particle_mass(sink_tags)
    print "sink_masses", sink_masses.as_quantity_in(units.MSun)
    #print "star_ind", star_ind
    print "all_masses", all_masses

    # Check the mass in each sink and if enough to make stars
    # then make some.

    for s,sink_mass in enumerate(sink_masses.value_in(units.MSun)):
        print "current star up for assignment:", \
              all_masses[sink_tags[s]][0], \
              "for sink with mass:", sink_mass
        # Does this sink have enough mass to make a star? If so, party on.
        while(sink_mass > all_masses[sink_tags[s]][0]):

            sink_position      = hydro.get_particle_position(sink_tags[s])
            sink_vel           = hydro.get_particle_velocity(sink_tags[s]).value_in(units.cm/units.s)
            sink_cs            = hydro.get_sink_mean_cs(sink_tags[s]).value_in(units.cm/units.s)

            print "Sink vel =", sink_vel
            print "Sink mean gas sound speed =", sink_cs

            # Check before pop.
            #print "all_masses[sink_tags[s],[0]]=", all_masses[sink_tags[s]][0]
            new_star_mass     = all_masses[sink_tags[s]][0]
            all_masses[sink_tags[s]] = np.delete(all_masses[sink_tags[s]],0)
            # Check new_star_mass
            #print "new_star_mass =", new_star_mass
            # Check that it deleted properly.
            #print "all_masses[sink_tags[s],[0]]=", all_masses[sink_tags[s]][0]

            new_star          = Particles(1)
            new_star.mass     = new_star_mass  | units.MSun
            new_star.velocity = np.random.uniform(sink_vel, np.ones(3)*sink_cs) | units.cm/units.s

            sink_mass         = sink_mass  - new_star_mass

            # Singular isothermal spherical distribution.
            stars_rvec = (random_three_vector(1)[:,:]*(np.random.rand(1))[:,None]*sink_rad)
            #print "stars_rvec=", stars_rvec
            rx = stars_rvec[:,0]
            ry = stars_rvec[:,1]
            rz = stars_rvec[:,2]
            #print "rx, ry, rz =", rx, ry, rz
            r2 = (rx**2.0 + ry**2.0 + rz**2.0).in_(units.cm**2.0)
            new_star.position = np.add(stars_rvec.value_in(units.cm),sink_position.value_in(units.cm)) | units.cm

            #print "new star position =", new_star.position.as_quantity_in(units.parsec)
            #print "new star velocity =", new_star.velocity.as_quantity_in(units.km/units.s)
            #print "new star mass =" , new_star.mass.as_quantity_in(units.MSun)

            # Remove the mass from the sink.
            hydro.set_particle_mass(sink_tags[s], (sink_mass | units.MSun))

            # Make the new star particle.
            hydro.set_particle_pointers('mass')
            new_star_tag = hydro.add_particles(new_star.x, new_star.y, new_star.z)
            hydro.set_particle_mass(new_star_tag, new_star.mass)
            hydro.set_particle_velocity(new_star_tag, new_star.vx, new_star.vy, new_star.vz)
            hydro.set_particle_oldmass(new_star_tag, new_star.mass) # Save initial stellar mass for SE code.
            # Switch back to sinks to continue the loop.
            hydro.set_particle_pointers('sink')

            # Tell the main code we made a star.
            if (formed_stars == False): formed_stars = True
            if (new_star_mass > min_mass.value_in(units.MSun)): formed_massive_star = True

    # Last thing is to ensure we are pointing back at massive particles.
    hydro.set_particle_pointers('mass')

    return (formed_stars, formed_massive_star)

def random_three_vector(n=1):
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """

    three_vector = np.zeros((n,3))

    phi = np.random.uniform(0,np.pi*2,n)
    costheta = np.random.uniform(-1,1,n)

    theta = np.arccos( costheta )
    three_vector[:,0] = np.sin( theta) * np.cos( phi )
    three_vector[:,1] = np.sin( theta) * np.sin( phi )
    three_vector[:,2] = np.cos( theta )
    return three_vector

def remove_particles_outside_bndbox(hydro, stars, grav, mult_grav, with_multiples, tags_keys, bndbox, debug=True):

    ### Remove any particles that have left the simulation.

    #print "Grav pos", grav.particles.x.value_in(units.cm), \
    #                  grav.particles.y.value_in(units.cm), \
    #                  grav.particles.z.value_in(units.cm), \
    #                  bndbox.value_in(units.cm), \
    #                  bndbox.value_in(units.cm)

    stars_removed = False
    rem_index = 0
    num_particles = len(stars)

    if (num_particles > 0):
       # rem_index = np.where(np.abs(grav.particles.x.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        if (debug_remove):
            print np.abs(grav.particles.x.value_in(units.cm)).max()
            print "grav x remove", np.where(np.greater_equal(np.abs(grav.particles.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            if (with_multiples): print "multiples.stars x remove", np.where(np.greater_equal(np.abs(mult_grav.stars.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "stars x remove", np.where(np.greater_equal(np.abs(stars.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "hydro x remove", np.where(np.greater_equal(np.abs((hydro.get_particle_position(tags_keys[:,0])[:,0]).value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
        #rem_index = np.where(np.greater_equal(np.abs(stars.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        if (debug_remove): print "x rem index=", rem_index
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
        if (debug_remove): print "grav x rem index=", grav_rem_index
        if (debug_remove): print "star x rem index=", star_rem_index
        if (debug_remove): print "x key =", rem_key
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    if (debug_remove): print "Remove index = ", rem_index
    if (debug_remove): print "Remove position x =", grav.particles.x[rem_index]
    #if (debug_remove): print "Remove position x =", stars.x[rem_index]
    rem_size = rem_index.size
    if (debug_remove): print "Remove index size = ", rem_size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):
        if (debug_remove): print "Remove position x from key =", grav.particles.x[grav_rem_index]

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        #for iii in np.nditer(rem_index):
        #    #rem_part.add_particle(stars[iii])
        #    grav_rem_part.add_particle(grav.particles[iii])
        #    if grav.particles[iii] in stars:
        #        if (debug_remove):
        #            print "Found this grav particle in stars."
        #            #print "Its", grav.particles.get_all_indices_in_store()[iii], "in grav and", \
        #            # stars.get_indices_of_keys([grav.particles.key[iii]]), "in stars."
        #            print "Its", grav.particles[np.where(grav.particles==grav.particles[iii])], "in grav and", \
        #             stars[np.where(stars==grav.particles[iii])], "in stars."
        #        stars_rem_part.add_particles(stars[np.where(stars==grav.particles[iii])])
        #        #stars_rem_part.add_particle(stars[stars.get_indices_of_keys([grav.particles.key[iii]])])
        #if (debug_remove): print "Added", stars[np.where(stars==stars_rem_part)], "to stars_rem_part."
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
                    #if (debug_multiples):
                    print "Removing", len(leaves), "particles that are multiples by deleting the root of the tree."
                    #stars_rem_part.add_particles(stars[stars.get_indices_of_keys(leaves.key)])
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    #rem_tag.append(stars[stars.get_indices_of_keys(leaves.key)].tag)
                    #rem_tag.append(stars[np.where(stars==leaves)].tag)
                    if (debug_remove):
                        print "Added", len(stars_rem_part.get_all_indices_in_store()), "to stars_rem_part."
                        #print "rem_tag is now", rem_tag
                    sys.stdout.flush()
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    #hydro.remove_particles(leaves.as_particle_in_set(stars).tag)
                    #stars.remove_particles(leaves.as_particle_in_set(stars))
                    #grav.particles.remove_particles(st)
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]
                    #rem_part.remove_particles(st)
                    #rem_size -= 1

        rem_tag.append(stars_rem_part.tag)
        if (debug_remove): print "Tags for removal are now", rem_tag

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        if (debug_remove): print "Tags for removal are now", rem_tag
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)
        #if (with_se):
            #se.particles.remove_particles(rem_part)
        if (debug_remove): print "Now ", len(grav.particles), "particles in grav."
        if (debug_remove): print "Now",  len(stars), "particles in stars."
        if (debug_remove): print "and", hydro.get_number_of_particles(), "particles in hydro."
        if (with_multiples and debug_remove):
            print "Now ", len(mult_grav.root_to_tree), "multiples in multiples."
            print "Now ", len(mult_grav.stars), "leaves in multiples."


    if (num_particles > 0):
        #rem_index = np.where(np.abs(grav.particles.y.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        if (debug_remove):
            print np.abs(grav.particles.y.value_in(units.cm)).max()
            print "grav y remove", np.where(np.greater_equal(np.abs(grav.particles.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            if (with_multiples): print "multiples.stars y remove", np.where(np.greater_equal(np.abs(mult_grav.stars.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "stars y remove", np.where(np.greater_equal(np.abs(stars.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "hydro y remove", np.where(np.greater_equal(np.abs((hydro.get_particle_position(tags_keys[:,0])[:,1]).value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        #print "y remove",  np.where(np.abs(grav.particles.y.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        #rem_index = np.where(np.greater_equal(np.abs(stars.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
        #rem_index = np.where(np.greater_equal(np.abs(stars.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        if (debug_remove): print "y rem index=", rem_index
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
        if (debug_remove): print "grav y rem index=", grav_rem_index
        if (debug_remove): print "star y rem index=", star_rem_index
        if (debug_remove): print "y key =", rem_key
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    if (debug_remove): print "Remove index = ", rem_index
    if (debug_remove): print "Remove position y =", grav.particles.y[rem_index]
    #if (debug_remove): print "Remove position y =", stars.y[rem_index]
    rem_size = rem_index.size
    if (debug_remove): print "Remove index size = ", rem_size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):
        if (debug_remove): print "Remove position y from key =", grav.particles.y[grav_rem_index]

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        #for iii in np.nditer(rem_index):
        #    #rem_part.add_particle(stars[iii])
        #    grav_rem_part.add_particle(grav.particles[iii])
        #    if grav.particles[iii] in stars:
        #        if (debug_remove):
        #            print "Found this grav particle in stars."
        #            #print "Its", grav.particles.get_all_indices_in_store()[iii], "in grav and", \
        #            # stars.get_indices_of_keys([grav.particles.key[iii]]), "in stars."
        #            print "Its", grav.particles[np.where(grav.particles==grav.particles[iii])], "in grav and", \
        #             stars[np.where(stars==grav.particles[iii])], "in stars."
        #        stars_rem_part.add_particles(stars[np.where(stars==grav.particles[iii])])
        #        #stars_rem_part.add_particle(stars[stars.get_indices_of_keys([grav.particles.key[iii]])])
        #if (debug_remove): print "Added", stars[np.where(stars==stars_rem_part)], "to stars_rem_part."
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
                    #if (debug_multiples):
                    print "Removing", len(leaves), "particles that are multiples by deleting the root of the tree."
                    #stars_rem_part.add_particles(stars[stars.get_indices_of_keys(leaves.key)])
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    #rem_tag.append(stars[stars.get_indices_of_keys(leaves.key)].tag)
                    #rem_tag.append(stars[np.where(stars==leaves)].tag)
                    if (debug_remove):
                        print "Added", len(stars_rem_part.get_all_indices_in_store()), "to stars_rem_part."
                        #print "rem_tag is now", rem_tag
                    sys.stdout.flush()
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    #hydro.remove_particles(leaves.as_particle_in_set(stars).tag)
                    #stars.remove_particles(leaves.as_particle_in_set(stars))
                    #grav.particles.remove_particles(st)
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]
                    #rem_part.remove_particles(st)
                    #rem_size -= 1

        rem_tag.append(stars_rem_part.tag)
        if (debug_remove): print "Tags for removal are now", rem_tag

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        if (debug_remove): print "Tags for removal are now", rem_tag
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)
        #if (with_se):
            #se.particles.remove_particles(rem_part)
        if (debug_remove): print "Now ", len(grav.particles), "particles in grav."
        if (debug_remove): print "Now",  len(stars), "particles in stars."
        if (debug_remove): print "and", hydro.get_number_of_particles(), "particles in hydro."
        if (with_multiples and debug_remove):
            print "Now ", len(mult_grav.root_to_tree), "multiples in multiples.root_to_tree."
            print "Now ", len(mult_grav.stars), "leaves in multiples."


    if (num_particles > 0):
        #rem_index = np.where(np.abs(grav.particles.z.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        if (debug_remove):
            print np.abs(grav.particles.z.value_in(units.cm)).max()
            print "grav z remove", np.where(np.greater_equal(np.abs(grav.particles.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            if (with_multiples): print "multiples.stars z remove", np.where(np.greater_equal(np.abs(mult_grav.stars.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "stars z remove", np.where(np.greater_equal(np.abs(stars.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
            print "hydro z remove", np.where(np.greater_equal(np.abs((hydro.get_particle_position(tags_keys[:,0])[:,2]).value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        #print "z remove",  np.where(np.abs(grav.particles.z.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        #rem_index = np.where(np.greater_equal(np.abs(stars.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
        #rem_index = np.where(np.greater_equal(np.abs(stars.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        if (debug_remove): print "z rem index=", rem_index
        grav_rem_index = np.array([np.where(grav.particles.key == x)[0] for x in rem_key]).flatten()
        star_rem_index = np.array([np.where(stars.key == x)[0] for x in rem_key]).flatten()
        if (debug_remove): print "grav z rem index=", grav_rem_index
        if (debug_remove): print "star z rem index=", star_rem_index
        if (debug_remove): print "z key =", rem_key
    else:
        rem_index = np.empty(0)
        rem_key   = np.empty(0)
    if (debug_remove): print "Remove index = ", rem_index
    if (debug_remove): print "Remove position z =", grav.particles.z[rem_index]
    #if (debug_remove): print "Remove position z =", stars.z[rem_index]
    rem_size = rem_index.size
    if (debug_remove): print "Remove index size = ", rem_size
    rem_tag = []
    grav_rem_part = Particles()
    stars_rem_part = Particles()

    if (rem_size > 0):
        if (debug_remove): print "Remove position z from key =", grav.particles.z[grav_rem_index]

        stars_removed = True
        print "About to try and remove", len(rem_index), "from grav."
        print "Currently", len(grav.particles), "particles in grav."
        #for iii in np.nditer(rem_index):
        #    #rem_part.add_particle(stars[iii])
        #    grav_rem_part.add_particle(grav.particles[iii])
        #    if grav.particles[iii] in stars:
        #        if (debug_remove):
        #            print "Found this grav particle in stars."
        #            #print "Its", grav.particles.get_all_indices_in_store()[iii], "in grav and", \
        #            # stars.get_indices_of_keys([grav.particles.key[iii]]), "in stars."
        #            print "Its", grav.particles[np.where(grav.particles==grav.particles[iii])], "in grav and", \
        #             stars[np.where(stars==grav.particles[iii])], "in stars."
        #        stars_rem_part.add_particles(stars[np.where(stars==grav.particles[iii])])
        #        #stars_rem_part.add_particle(stars[stars.get_indices_of_keys([grav.particles.key[iii]])])
        #if (debug_remove): print "Added", stars[np.where(stars==stars_rem_part)], "to stars_rem_part."
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
                    #if (debug_multiples):
                    print "Removing", len(leaves), "particles that are multiples by deleting the root of the tree."
                    #stars_rem_part.add_particles(stars[stars.get_indices_of_keys(leaves.key)])
                    for leaf in leaves:
                        stars_rem_part.add_particles(stars[np.where(stars.key==leaf.key)[0]])
                    #rem_tag.append(stars[stars.get_indices_of_keys(leaves.key)].tag)
                    #rem_tag.append(stars[np.where(stars==leaves)].tag)
                    if (debug_remove):
                        print "Added", len(stars_rem_part.get_all_indices_in_store()), "to stars_rem_part."
                        #print "rem_tag is now", rem_tag
                    sys.stdout.flush()
                    # Note only the stars particle set has a tag attribute.
                    # Also all leaves exist in stars and hydro while all roots
                    # exist in grav.
                    #hydro.remove_particles(leaves.as_particle_in_set(stars).tag)
                    #stars.remove_particles(leaves.as_particle_in_set(stars))
                    #grav.particles.remove_particles(st)
                    del mult_grav.root_to_tree[st.as_particle_in_set(mult_grav._inmemory_particles)]
                    #rem_part.remove_particles(st)
                    #rem_size -= 1

        rem_tag.append(stars_rem_part.tag)
        if (debug_remove): print "Tags for removal are now", rem_tag

        rem_tag = np.array(rem_tag).flatten()
        rem_tag.sort()
        if (debug_remove): print "Tags for removal are now", rem_tag
        grav.particles.remove_particles(grav_rem_part)
        mult_grav._inmemory_particles.remove_particles(grav_rem_part)
        hydro.remove_particles(rem_tag)
        stars.remove_particles(stars_rem_part)
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)
        #if (with_se):
            #se.particles.remove_particles(rem_part)
        if (debug_remove): print "Now ", len(grav.particles), "particles in grav."
        if (debug_remove): print "Now",  len(stars), "particles in stars."
        if (debug_remove): print "and", hydro.get_number_of_particles(), "particles in hydro."
        if (with_multiples and debug_remove):
            print "Now ", len(mult_grav.root_to_tree), "multiples in multiples."
            print "Now ", len(mult_grav.stars), "leaves in multiples."

        # Sync stars to the particles in grav.
        if (with_multiples):
            grav.particles.synchronize_to(mult_grav._inmemory_particles)
            #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
        else:
            grav.particles.synchronize_to(stars)

    return stars_removed, num_particles


def check_sanity(hydro, grav, stars = None, min_pos_diff = 1.0e-4 | units.AU ,
                              min_vel_diff = 1.0e-2 | units.km / units.s,
                              min_mass_diff = 0.01 | units.MSun,
                              kill=True, with_multiples=False):

    print "Num in grav = ", len(grav.particles)
    if (with_multiples):
         print "Num in multiples._inmemory_particles = ", len(mult_grav._inmemory_particles)
         print "Num in multiples.stars = ", len(mult_grav.stars)
         print "Num in multiples.root_to_tree = ", len(mult_grav.root_to_tree)
    print "Num in hydro = ", hydro.get_number_of_particles()
    print "Num in stars = ", len(stars)

    if (with_multiples):
        print "Grav id =", grav.particles.index_in_code
        print "Mult mem id =", mult_grav._inmemory_particles.id
        print "Multiples stars id =", mult_grav.stars.id
        print "Stars id =", stars.id
    print "Stars tags = ", stars.tag

    if (with_multiples):
        pos_diff     = (hydro.get_particle_position(tags_keys[:,0]).as_quantity_in(units.m)
                     - stars.position.as_quantity_in(units.m))
        vel_diff     = (hydro.get_particle_velocity(tags_keys[:,0]).as_quantity_in(units.m*units.s**(-1))
                     - stars.velocity.as_quantity_in(units.m*units.s**(-1)))
        mass_diff    = (hydro.get_particle_mass(tags_keys[:,0]).in_(units.MSun)
                     - stars.mass.in_(units.MSun))
    else:
        pos_diff     = hydro.get_particle_position(tags_keys[:,0]).as_quantity_in(units.m) - grav.particles.position
        vel_diff     = hydro.get_particle_velocity(tags_keys[:,0]).as_quantity_in(units.m*units.s**(-1)) - grav.particles.velocity
        mass_diff    = hydro.get_particle_mass(tags_keys[:,0]).in_(units.MSun) - grav.particles.mass.in_(units.MSun)

    pos_diff_max = np.abs(pos_diff.value_in(units.cm)).max() | units.cm
    print "Max Position diff = "
    print pos_diff_max.in_(units.AU)

    vel_diff_max = np.abs(vel_diff.value_in(units.cm*(units.s**-1))).max() | units.cm/units.s
    print "Max Velocity diff = "
    print vel_diff_max.in_(units.km/units.s)

    mass_diff_max =  np.abs(mass_diff.value_in(units.MSun)).max() | units.MSun
    print "Max Mass diff = "
    print mass_diff_max.in_(units.MSun)

    repeat_pos = False
    if (with_multiples):
        print "Checking stars for repeat position."
        repeat_pos = check_repeat_position(stars)
        print "Checking mult.stars for repeat position."
        repeat_pos = repeat_pos or check_repeat_position(mult_grav.stars)
        print "Grav id =", grav.particles.index_in_code
        print "Mult mem id =", mult_grav._inmemory_particles.id
        print "Multiples stars id =", mult_grav.stars.id
    if (repeat_pos): sys.exit()

    if (kill):
        if (pos_diff_max.value_in(units.cm)  > min_pos_diff.value_in(units.cm)
            or vel_diff_max.value_in(units.cm*(units.s**-1)) > min_vel_diff.value_in(units.cm*(units.s**-1))
            or mass_diff_max.value_in(units.MSun) > min_mass_diff.value_in(units.MSun)):

            print "hydro pos"
            print hydro.get_particle_position(tags_keys[:,0]).as_quantity_in(units.AU)
            print "grav pos"
            if (with_multiples):
                print stars.position.as_quantity_in(units.AU)
            else:
                print grav.particles.position.as_quantity_in(units.AU)
            print "hydro vel"
            print hydro.get_particle_velocity(tags_keys[:,0]).as_quantity_in(units.km*units.s**(-1))
            print "grav vel"
            if (with_multiples):
                print stars.velocity.as_quantity_in(units.km*units.s**(-1))
            else:
                print grav.particles.velocity.as_quantity_in(units.km*units.s**(-1))
            print "hydro mass"
            print hydro.get_particle_mass(tags_keys[:,0]).in_(units.MSun)
            print "grav mass"
            if (with_multiples):
                print stars.mass.in_(units.MSun)
            else:
                print grav.particles.mass.in_(units.MSun)


            print "Difference in position or velocity greater than tolerance. Stopping code now."
            print pos_diff.in_(units.AU)
            print vel_diff.in_(units.km/units.s)
            print mass_diff.in_(units.MSun)
            sys.exit()
    return

def check_repeat_position(stars):
    repeat_pos     = False
    x,c= np.unique(stars.x.number, return_counts=True)
    if (np.any(c>1)):
        print c[np.where(c>1)[0]], "repeated x position found in stars."
        print "Pos =", stars.x[np.where(stars.x.number==x[np.where(c>1)[0]])[0]]
        repeat_pos=True
    x,c= np.unique(stars.y.number, return_counts=True)
    if (np.any(c>1)):
        print c[np.where(c>1)[0]], "repeated y position found in stars."
        print "Pos =", stars.y[np.where(stars.y.number==x[np.where(c>1)[0]])[0]]
        repeat_pos=True
    x,c= np.unique(stars.z.number, return_counts=True)
    if (np.any(c>1)):
        print c[np.where(c>1)[0]], "repeated z position found in stars."
        print "Pos =", stars.z[np.where(stars.z.number==x[np.where(c>1)[0]])[0]]
        repeat_pos=True

    return repeat_pos

def oldvnew_position(stars, oldstars):
    repeat_pos     = False
    repeat_pos_any = False

    x_check = np.equal(oldstars.x.number, stars.x.number)
    y_check = np.equal(oldstars.y.number, stars.y.number)
    z_check = np.equal(oldstars.z.number, stars.z.number)

    repeat_pos = x_check.any()
    if (repeat_pos):
        print "Repeated x position found at:", stars.x[x_check]
        repeat_pos_any=True
    repeat_pos = y_check.any()
    if (repeat_pos):
        print "Repeated y position found at:", stars.y[y_check]
        repeat_pos_any=True
    repeat_pos = z_check.any()
    if (repeat_pos):
        print "Repeated z position found at:", stars.z[z_check]
        repeat_pos_any=True

    return repeat_pos_any

def check_stellar_type(stars):

    for star in stars:
        print "Star mass, type=", star.mass.in_(units.MSun), star.stellar_type, star.stellar_type.value_in(units.stellar_type)
    return

class stellar_wind(object):
    """Implementation of stellar winds based on Kudritzki and Puls ARAA 2000 and Vink A&A 2000."""

    def __init__(self, teff, mass, lum, radius): #, thom_sig=0.32, thom_Gam=0.0,
        #         vesc  = 0.0 | units.cm / units.s,
        #         vterm = 0.0 | units.cm / units.s,
        #         dm_dt = 0.0 | units.g / units.s):

        self.mass     = mass
        self.lum      = lum
        self.teff     = teff
        self.radius   = radius
        self.thom_sig()
        self.thom_Gam()
        self.vesc()
        self.vterm()
        self.dm_dt()

        return


    def thom_sig(self):

        if (self.teff.value_in(units.K) < 3e4):
            self.thom_sig = 0.31 # | units.cm**2.0 / units.g
        elif (3e4 <= self.teff.value_in(units.K) < 3.5e4):
            self.thom_sig = 0.32 # | units.cm**2.0 / units.g
        else:
            self.thom_sig = 0.33 # | units.cm**2.0 / units.g
        return

    def thom_Gam(self):

        self.thom_Gam = 7.66e-5*self.thom_sig/self.mass.value_in(units.MSun)*self.lum.value_in(units.LSun)

        #print "thom_Gam=", self.thom_Gam

        return

    def vesc(self):

        self.vesc = np.sqrt(2.0*units.constants.G*self.mass*(1-self.thom_Gam)
                            /(self.radius)).as_quantity_in(units.km / units.s)

        #print self.vesc

        return

    def vterm(self):

        if (self.teff.value_in(units.K) <= 1.0e4):
            self.vterm = (self.vesc)
        elif (1.0e4 < self.teff.value_in(units.K) < 2.1e4):
            self.vterm = (1.4*self.vesc)
        else:
            self.vterm = (2.65*self.vesc)

        #print self.vterm.value_in(units.km / units.s)

        return

    def dm_dt(self):
        # Above the bi-stability jump (larger than B1).
        if (self.teff.value_in(units.K) > 2.75e4):

            self.dm_dt = 10**(self.mass_loss1()) | units.MSun / units.yr
        # Below the bi-stability jump (smaller than B1).
        elif (self.teff.value_in(units.K) < 2.25e4):

            self.dm_dt = 10**(self.mass_loss2()) | units.MSun / units.yr

        # Linear interpolation between the two.
        else:

            xp = np.array([2.25e4, 2.75e4])
            fp = np.array([self.mass_loss2(2.25e4), self.mass_loss1(2.75e4)])

            self.dm_dt = 10**(np.interp(self.teff.value_in(units.K), xp, fp))  | units.MSun / units.yr

        return

    # Note we make the temp passable so we can interpolate if we need to
    # and we return a value here for the same reason.

    # Above the bi-stability jump (larger than B1).
    def mass_loss1(self, teff=None):

        if (teff is None):
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.697 + 2.194*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.313*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.226*np.log10(self.vterm/self.vesc/2.0) \
                            + 0.933*np.log10(teff/4e4) \
                            - 10.92*np.log10(teff/4e4)**2.0
        return log_dm_dt

    # Below the bi-stability jump (smaller than B1).
    def mass_loss2(self, teff=None):

        if (teff is None):
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.688 + 2.210*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.339*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.601*np.log10(self.vterm/self.vesc/2.0) \
                            + 1.07*np.log10(teff/2e4)
        return log_dm_dt


def get_np_from_run_script():
    np = -1
    with open("run.sh") as f:
        for i, line in enumerate(f):
            words = line.split()
            if (len(words) >= 3 and words[0] == '#SBATCH' and words[1] == '-n'):
                np = int(words[2])
                break
    if (np == -1): print "WARNING! No proper number of procs found!"
    return np

def get_restart_and_chk_num():
    restart = False
    chknum = None
    with open("flash.par") as f:
        for i, line in enumerate(f):
            words = line.split()
            if (len(words) > 1):
                if (words[0].lower() == 'restart'):
                    if (words[2].lower() == '.true.'): restart = True
                if (words[0].lower() == 'checkpointfilenumber'):
                    chknum = words[2]
                if (words[0].lower() == 'plotfilenumber'):
                    pltnum = words[2]
                if (words[0].lower() == 'basenm'):
                    basename = words[2].strip("\"")

    return restart, int(chknum), int(pltnum), basename

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
    #id = np.arange(len(stars)) # Moved to add_particles_to_grav
    #stars.id = id+1
    grav.parameters.epsilon_squared = 0.0 | units.cm**2.0
    print "Gravity epsilon set to zero."
    stopping_condition = grav.stopping_conditions.collision_detection
    stopping_condition.enable()
    print "Stopping condition enabled."

    #init_smalln(conv)

    if (kep == None):
        print "Starting Kepler."
        kep = Kepler(unit_converter=conv)
        print "Initializing Kepler."
        kep.initialize_code()

    print "Starting multiples."
    multiples_code = multiples.Multiples(grav, new_smalln, kep,
                                         constants.G)

    multiples_code.global_debug = mult_debug_level

    multiples_code.neighbor_veto               = True
    multiples_code.check_tidal_perturbation    = True
    multiples_code.neighbor_perturbation_limit = 0.05
    multiples_code.wide_perturbation_limit     = 0.08
    print "Setting up channels between stars and multiples."
    mult_to_stars = multiples_code.stars.new_channel_to(stars)
    stars_to_mult = stars.new_channel_to(multiples_code.stars)
    print "Multiples initialized."
    return multiples_code, mult_to_stars, stars_to_mult

def cleanup_multiples(multiples_code, kep_code):

    kep_code.stop()
    stop_smalln()
    del multiples
    return

def initialize_gravity_codes(convert, stars = None, start_time = None,
                             num_grav_workers = 1,
                             eps = 15.0 | units.RSun,
                             with_ph4 = True, tree_exists = False,
                             with_multiples = False):
    print "Starting gravity code."

### NOTE! Trying to use the generic unit converter here after I saw this on the AMUSE website
### amusecode.org/doc/reference/quantities_and_units.html

    if (with_ph4):
        grav = ph4(convert, number_of_workers=num_grav_workers, mode='cpu', redirection="none")
        grav.parameters.set_defaults()
        grav.parameters.timestep_parameter=0.14  # Timestep accuracy for PH4
        grav.parameters.force_sync=1 # Force the code to end exactly at the specified time.
#grav.parameters.epsilon_squared = (6.8359e14 | units.cm)**2.0
    else:
        grav = Hermite(convert, number_of_workers=num_grav_workers)
        #, debugger="gdb-remote", debugger_port=4343)
        grav.parameters.end_time_accuracy_factor=0.0 #1e-8
        grav.parameters.dt_param=0.02 #0.14**2  # Timestep accuracy for Hermite
        #grav.parameters.stopping_conditions_timeout = 60 | units.s
        #grav.parameters.stopping_conditions_number_of_steps = 999999999
        #grav.parameters.stopping_conditions_minimum_internal_energy = -1e99 | units.m**2 * units.s**-2
        #grav.parameters.stopping_conditions_maximum_internal_energy =  1e99 | units.m**2 * units.s**-2

# N-body softening radius is the actual radius of a large massive star here.
    if (with_multiples):
        grav.parameters.epsilon_squared = 0.0 | units.cm**2.0
    else:
        grav.parameters.epsilon_squared = eps**2.0
    #if (tree_exists):
    #    tree = Fi(convert)
    #    tree.parameters.epsilon_squared = (eps/2.5)**2.0

    if (start_time is not None):
        if (with_ph4):
            grav.parameters.begin_time=start_time
            grav.parameters.sync_time=start_time
            grav.parameters.force_sync=1
        else:
            grav.parameters.begin_time=start_time
            grav.evolve_model(start_time)


    # Make a particle set to pass information back and forth
    # between Flash and the gravity code as well as do the bridge kicks.
    if (stars is None):
        stars = Particles(0)
    else:
        stars = stars
    stars_to_grav = stars.new_channel_to(grav.particles)
    grav_to_stars = grav.particles.new_channel_to(stars)

    #if (tree_exists):
    #    stars_to_tree = stars.new_channel_to(tree.particles)

    #if (with_multiples):
        #multiples_code = initialize_multiples(stars, grav, convert)
    #else:
    multiples_code = None
    #if (with_multiples):
    return stars, multiples_code, grav, stars_to_grav, grav_to_stars
    #else:
    #    return stars, grav, stars_to_grav, grav_to_stars

def cleanup_gravity_codes(grav, tree_exists=False, with_multiples=False):

    grav.stop()

    if (with_multiples):
        kep.stop()
        stop_smalln()
    return


def update_psetA_from_psetB(setA, setB, debug=False):
    '''
    Updates the properties of particle set setA from the properties of setB,
    assuming they share common keys and are NOT in the same order.
    It updates mass, position and velocity.
    '''
    if (debug): print "In call to update_psetA_from_psetB."
    numBcycled = 0
    numAupdated = 0
    for p in setB:
        numBcycled += 1
        A_ind = np.where(setA.key == p.key)[0]
        if   (len(A_ind)>1): print "Warning, found more than one key in set A that matches this key from set B!"
        elif (len(A_ind)<1): print "Warning, did not find this key from set B in set A!"
        elif (len(A_ind)==1):
            numAupdated += 1
            setA[A_ind].mass = p.mass
            setA[A_ind].position = p.position
            setA[A_ind].velocity = p.velocity
            #setA[A_ind].x = p.x
            #setA[A_ind].y = p.y
            #setA[A_ind].z = p.z
            #setA[A_ind].vx = p.vx
            #setA[A_ind].vy = p.vy
            #setA[A_ind].vz = p.vz

    if (debug): print "Cycled", numBcycled, "in setB."
    if (debug): print "Found", numAupdated, "in setA"
    return

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

def check_root_and_leaves(mult_grav, grav, stars, kill=True):

    print "In check_root_and_leaves."
    for star in stars:
        for root, tree in mult_grav.root_to_tree.iteritems():
            root_particle = root.as_particle_in_set(mult_grav._inmemory_particles)
            grav_particle = root.as_particle_in_set(grav.particles)
            leaves = tree.get_leafs_subset()
            if star in leaves:
                mult_ind         = np.where(mult_grav.stars.id == star.id)[0]
                mult_star        = mult_grav.stars[mult_ind]
                same_spot = mult_star.x.value_in(units.cm) == star.x.value_in(units.cm)

                if (not same_spot):
                    print "Star id   =", star.id
                    print "Mult id   =", mult_star.id
                    print "Root id   =", root_particle.id
                    print "Grav id   =", grav_particle.index_in_code
                    print "Leaf id   =", leaves.id
                    print "Star age  =", star.age
                    print "Star mass =", star.mass.in_(units.MSun)
                    print "Mult mass =", mult_star.mass.in_(units.MSun)
                    print "Leaf mass =", leaves.mass.in_(units.MSun)
                    print "Root mass =", root_particle.mass.in_(units.MSun)
                    print "Tree mass =", tree.particle.mass.in_(units.MSun)
                    print "Grav mass =", grav_particle.mass.in_(units.MSun)
                    print "Star pos  =", star.position.in_(units.cm)
                    print "Mult pos  =", mult_star.position.in_(units.cm)
                    print "Leaf pos  =", leaves.position.in_(units.cm)
                    print "Root pos  =", root_particle.position.in_(units.cm)
                    print "Tree pos  =", tree.particle.position.in_(units.cm)
                    print "Grav pos  =", grav_particle.position.in_(units.cm)
                    print "Star vel  =", star.velocity.in_(units.km/units.s)
                    print "Mult vel  =", mult_star.velocity.in_(units.km/units.s)
                    print "Leaf vel  =", leaves.velocity.in_(units.km/units.s)
                    print "Root vel  =", root_particle.velocity.in_(units.km/units.s)
                    print "Tree vel  =", tree.particle.velocity.in_(units.km/units.s)
                    print "Grav vel  =", grav_particle.velocity.in_(units.km/units.s)
                    sys.stdout.flush()

                    if (kill):
                        print "Mult_star and star are not in the same spot. Exiting!"
                        sys.exit()

    return

def update_stars_from_leaves(mult_grav, grav, stars):

    local_debug = False
    print "In update_stars_from_leaves."
    for star in stars:
        for root, tree in mult_grav.root_to_tree.iteritems():
            root_particle = root.as_particle_in_set(mult_grav._inmemory_particles)
            grav_particle = root.as_particle_in_set(grav.particles)
            leaves = tree.get_leafs_subset()
            for leaf in leaves:
                if (star == leaf):
                    star.position = leaf.position
                    star.velocity = leaf.velocity
                    star.mass     = leaf.mass
                    if (local_debug):
                        print "Star id   =", star.id
                        print "Root id   =", root_particle.id
                        print "Grav id   =", grav_particle.index_in_code
                        print "Leaf id   =", leaves.id
                        print "Star age  =", star.age
                        print "Star mass =", star.mass.in_(units.MSun)
                        print "Leaf mass =", leaves.mass.in_(units.MSun)
                        print "Root mass =", root_particle.mass.in_(units.MSun)
                        print "Grav mass =", grav_particle.mass.in_(units.MSun)
                        print "Star pos  =", star.position.in_(units.cm)
                        print "Leaf pos  =", leaves.position.in_(units.cm)
                        print "Root pos  =", root_particle.position.in_(units.cm)
                        print "Grav pos  =", grav_particle.position.in_(units.cm)
                        print "Star vel  =", star.velocity.in_(units.km/units.s)
                        print "Leaf vel  =", leaves.velocity.in_(units.km/units.s)
                        print "Root vel  =", root_particle.velocity.in_(units.km/units.s)
                        print "Grav vel  =", grav_particle.velocity.in_(units.km/units.s)
                        sys.stdout.flush()

    return

# Write out pickles with the current random number state and the
# dictionary with all the future stars in it.
def write_rnd_and_mass_pickles(all_masses, output_dir, pklnum):
    print "Writing all_masses{:04d}.pickle".format(pklnum)
    with open(output_dir+'/all_masses{:04d}.pickle'.format(pklnum), 'wb') as f:
        pickle.dump(all_masses, f)
    print "Writing rnd_state{:04d}.pickle".format(pklnum)
    rnd_state = np.random.get_state()
    with open(output_dir+'/rnd_state{:04d}.pickle'.format(pklnum), 'wb') as f:
        pickle.dump(rnd_state, f)
    return


#def main():

### Set up an async pool so that we can evolve the two codes
### simultaneously.

#pool = AsyncRequestsPool()

### Define a function to handle requests to evolve the models which
### returns a result when both codes are finished evolving.

#def handle_result(request, index):
#    print "finished:", request.is_result_available(), index

#def handle_result(request, index):
            #self.assertTrue(request.is_result_available())
            #finished_requests.append(index)
            #print "finished:", request.is_result_available(), index


#########################################################################
# Converters

#convert = nbody.nbody_to_si(2.0e15 | units.cm, 3.0*1.0e-10 | units.MSun)
#convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
#                   1.0| units.cm,1.0 | units.g, 1 | units.s)
#cloud_dens = 3.82e-18 | units.g*units.cm**-3.0
#cloud_dens = 1.0e-21 | units.g*units.cm**-3
cloud_rad = 3.0*3.086e18 | units.cm # 5 pc  # 0.04029 | units.parsec
cloud_dens = 1.0e-18 | units.g*units.cm**-3
#cloud_rad  = 5e16 | units.cm # Cloud collapse test problem.
#cloud_mass = 4.0/3.0*np.pi*cloud_dens*cloud_rad**3.0
#cloud_mass = 1.0e3 | units.MSun
cloud_mass = 1000.0 | units.MSun
print cloud_mass
#total_mass = 1.1 | units.MSun
#bndbox = 6.0047e18 | units.cm # For Peters et al.
#bndbox = 1.5428388e+19 | units.cm  # For Richard's cloud model.
#bndbox = 7.0e18 | units.cm # For cloud collapse testing.#1.29093e+17 | units.cm # for SMT #
#bndbox = 4e19 | units.cm # For rad test.
#bndbox = 1.29093e+17 | units.cm # For SinkMomTest

# Converter for the N-body code.
convert = nbody.nbody_to_si(1.0 | units.parsec, 1000.0 | units.MSun)   # Note previously we had to multiply this radius by 0.1

#Converter for the hydro code.
convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
                   1.0| units.cm,1.0 | units.g, 1 | units.s)


# Radiation parameters
####################################
# Ionizing energy over 13.6 eV that heats the gas.
eion = 2.0 | units.eV #6.0 | units.eV
# Lyman cross section.
sigh = 6.3e-18 | units.cm**2.0
# Metallicity
z = 0.02
# Minimum mass of stars to turn on feedback.
min_mass = 7.0 | units.MSun #4.0 | units.MSun

# Min and max star formation masses. Max from IGIMF restrictions figuring 0.5 SFE.
cluster_formation_eff = 1.0
min_sf_mass = 0.08 | units.MSun
max_sf_mass = m_max_star(cluster_formation_eff*cloud_mass.value_in(units.MSun)) | units.MSun

print "Min IMF star mass=", min_sf_mass
print "Max IMF star mass=", max_sf_mass

h = 6.6261e-27 # Plank's constant
c = 2.9979e10  # Speed of light
k = 1.3807e-16 # Boltzman constant
sig0 = 6.304e-18 # Photoionization cross section at threshold for hydrogen

l_ev = 1.2398e-4 # wavelength for 1 eV.
l_min = 1e-7 #0.0 #np.finfo(float).eps # Something really small.
E_ev = 1.60222497096e-12 # energy of 1 eV.

E_min = 13.6*E_ev  # 13.6 eV
l_max = h*c/E_min  # wavelength of 13.6 eV.

nu_min = E_min / h # freq of 13.6 eV.
nu_max = np.inf    # Max freq.

sigDust = 1e-21 | units.cm**2.0 # Cross section for dust from Draine 2011
####################################


# Logging parameters for timing code.
##########################


# Number of processors for each code.
##########################
default_hy_procs = 2
num_grav_workers = 1

try:
    num_hy_workers= get_np_from_run_script() - num_grav_workers - 4
    # hydro + ph4 + 4 threads for amuse, smallN, kepler, seba
except:
    print "WARNING: Setting num_hy_workers from script failed. Defaulting to {} procs.".format(default_hy_procs)
    num_hy_workers = default_hy_procs

print "Number of hydro procs = ", num_hy_workers
print "Number of Nbody procs = ", num_grav_workers

##########################

# Debugging flags.
##########################
insane      = False
type_insane = False
debug_aptg  = False
debug_se    = False
debug_remove = False
debug_multiples = False
mult_debug_level = 1

# Check if any tags are repeats in the tags_keys array.
# All tags should be unique!
test_unique_tags = False

# Record particle sets to file.
write_psets = False
pdir = "./psets"

# Print timer summary every step?
profile = True

# Start with a star or cluster (either from file or made here).
##########################
start_with_cluster     = False #True
read_cluster_from_file = True
nm_part=10
start_with_star = False #True #
stars_current_id_num = 0

# Runtime parameters.
##########################
with_bridge = True
no_sinks    = False
tree_exists = False
with_ph4    = True #False
with_multiples = True #False

use_radiation   = True
pe_heat         = True

with_se         = True
with_sn         = True
do_sn_once      = False #True

with_winds      = True
with_massloss   = True
massloss_method = 'puls'

min_pos_diff = 0.01*3.086e16 #| units.m #((0.01*3.085e16) | units.m)
min_vel_diff = 1e3 # | units.m/units.s


# N-body softening radius is the actual radius of a large massive star here.
eps = 15.0 | units.RSun

stars = None
hydro_time = None

if (not with_bridge): print "WARNING NO GRAVITY BRIDGE."
if (not use_radiation): print "WARNING NO RADIATION."
if (not with_winds): print "WARNING NO WINDS."

sys.stdout.flush()
print "Starting smallN."
init_smalln(convert)
sys.stdout.flush()

print "Starting Kepler."
sys.stdout.flush()
kep = Kepler(unit_converter=convert)
print "Initializing Kepler."
sys.stdout.flush()
kep.initialize_code()
sys.stdout.flush()

print "Starting stellar evolution code."
sys.stdout.flush()
if (with_se):
#    #se = SSE()
    se = SeBa()
    se.initialize_code()
#else:
#    se = SeBa()
sys.stdout.flush()

stars, mult_grav, grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
         convert, stars = stars, start_time = hydro_time,
         num_grav_workers = num_grav_workers,
         eps = eps,
         with_ph4 = with_ph4,
         with_multiples = with_multiples)
print "Gravity code initialized."
sys.stdout.flush()

time.sleep(10)

print "Starting hydro code."
sys.stdout.flush()
hydro = Flash(unit_converter = convert2, number_of_workers=num_hy_workers, redirection='none')
print "Hydro code initialized."
sys.stdout.flush()

hydro.initialize_code()
print "hydro.initialize_code() called."
sys.stdout.flush()

### Get the simulation end time, current AMUSE step and
### current simulation time from Flash (in case of restart).

tmax = hydro.get_end_time()
t = hydro.get_time()
hydro_time = t

bndbox = hydro.get_runtime_parameter('xmax') | units.cm
print 'bndbox = ', bndbox

if (with_ph4):
    grav.parameters.begin_time=hydro_time
    grav.parameters.sync_time=hydro_time
    grav.parameters.force_sync=1
else:
    grav.parameters.begin_time=hydro_time
    grav.evolve_model(hydro_time)

logname = "profiler"+`num_hy_workers`+".log"
### Add the star particles to both codes.

#grav.particles.add_particles(stars)

###############################################
### Initialize either clusters or stars.
###############################################


if (start_with_cluster):

    from amuse.io import write_set_to_file, read_set_from_file

    if (read_cluster_from_file):

        print "Reading initial cluster from file."

        initial_cluster = read_set_from_file('stars.hdf5', 'hdf5')

    else:

    #make_single_star_in_hydro(x = -6.999e17 | units.cm, y = 0.0 | units.cm, z = 0.0 | units.cm,
                              #mass = 30.0 | units.MSun,
                              #vx = -2e8 | units.cm/units.s, vy = 0.0 | units.cm/units.s, vz = 0.0 | units.cm/units.s)

    #Converter for the initial cluster distribution, if there is one.
        conv_cluster = nbody_system.nbody_to_si(3.0 | units.parsec, 300.0 | units.MSun)

    # Build a plummer sphere with a Kroupa IMF
    #initial_cluster = make_cluster(conv_cluster, bndbox, fractal=False,  equal_mass=True, eq_mass=50.0 | units.MSun)
        initial_cluster = make_cluster(conv_cluster, nm_part, bndbox,
                                       fractal=True ,equal_mass=False,
                                       eq_mass=20.0 | units.MSun)

        print "Writing initial cluster to file."
        write_set_to_file(initial_cluster, 'starting_cluster.hdf5', 'hdf5')
        print "Done."

    #num_sub_clusters = 3

    #num_sub_stars = len(initial_cluster) // num_sub_clusters

    #posx1 = -3.086e18 | units.cm
    #posx2 =  3.086e18 | units.cm



    # Now lets distribute the stars as subclusters at several locations.

    print "Setting up initial cluster in hydro."
    make_cluster_in_hydro(initial_cluster)
    #make_cluster_in_hydro(initial_cluster[num_sub_stars:num_sub_stars*2], initial_x=posx2)
    #make_cluster_in_hydro(initial_cluster[num_sub_stars*2:])



if (start_with_star):

    x  = 0.0 | units.cm #(0.5*smallest_dx)
    y  = 0.0 | units.cm #(0.5*smallest_dx)
    z  = 0.0 | units.cm #(0.5*smallest_dx)
    vx = 0.0 | units.km/units.s
    vy = 0.0 | units.km/units.s
    vz = 0.0 | units.km/units.s
    m  = 30. | units.MSun

    make_single_star_in_hydro(x = x, y = y, z = z,
                              mass = m, initMass = m, age = 0.0 | units.Myr,
                              vx = vx, vy = vy, vz = vz)

    #make_single_star_in_hydro(x = 25.0 | units.parsec, y = 0.0 | units.cm, z = 0.0 | units.cm,
                              #mass = 60.0 | units.MSun,
                              #vx = 0.0 | units.cm/units.s, vy = 0.0 | units.cm/units.s, vz = 0.0 | units.cm/units.s)

### Initialize some parameters.

#np.random.seed(103180)     # Set initial random seed for testing/debugging.

num_particles = 0      # Total number of sinks in Flash
num_particles_old = 0  # Old number to compare if new sinks formed.
check_particles = 0    # Temp storage for num particles checks.
first_particle = False # Logical that makes the gravity code wait for a particle.
i = 0                  # Step number for bridge.
tags_keys = np.empty(0)# Tags (ID in Flash) and keys (ID in AMUSE) for all particles.
first_loop = True
made_stars = False
gridChanged = True
old_chk_file = '-1'
old_plt_file = '-1'

max_ref  = hydro.get_max_refinement()
sink_rad = 2.5*(2.0*bndbox)/(8.0*2.0**(float(max_ref)-1.0))
smallest_dx = (2.0*bndbox)/(8.0*2.0**(float(max_ref)-1.0))

time.sleep(5)

'''
stars = None
stars, mult_grav, grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
         convert, stars = stars, start_time = hydro_time,
         num_grav_workers = num_grav_workers,
         eps = eps,
         with_ph4 = with_ph4,
         with_multiples = with_multiples)
print "Gravity code initialized."
grav.initialize_code()
print "grav.initialize_code() called."
'''
# This should be automated!
restart,chknum,pltnum,basename = get_restart_and_chk_num()
refresh_rand_seed_on_restart = False #True

if (refresh_rand_seed_on_restart):
    print "WARNING WARNING WARNING! Resetting the rand seed!!!!!!!!!!"

output_dir = hydro.get_output_dir()

print "Writing output files to: ", output_dir

print "Is this a restart?", restart
print "We are assuming the checkpointfilenumber =", chknum
load_rnd_state_files(restart, chknum, refresh_rand_seed_on_restart)

if (not restart):
    write_rnd_and_mass_pickles(all_masses, output_dir, chknum)

# After FLASH inits (whether restart or not), it increments
# io_checkpointFileNumber for the /next/ file write.
chknum = hydro.IO_num('chk')

# This is used if you change the number of processors on restart!
# Or if you need to clear out tracers or some other thing.
# NOTE this is still bugged. 2/7/18 - JW
clear_particles_on_restart = False #True

### Set dt based on the current cloud free fall time.

tff = (np.sqrt(3.0*np.pi/(32.0*units.constants.G*cloud_dens))).as_quantity_in(units.s)

print "Hydro dynamical time in secs: %3.3e s" % (tff).value_in(units.s)
print (tff/100.0).as_quantity_in(units.yr)
print "Currently setting initial dt artbitrarily high so that dt = hydro dt."
dtinit = 1.0e10 | units.s #(tcross/100.0).as_quantity_in(units.s)
dtmax = 1.0e13 | units.s
dt = dtinit
print "dt in secs: %3.3e s" % dt.value_in(units.s)
print "End time in secs: %3.3e s" %tmax.value_in(units.s)

first_step = True

i = int(t/dt)

max_hy_steps = hydro.get_max_num_steps() # max number of iterations.
curr_hy_step = hydro.get_current_step()
total_time    = 0.0
time_in_hydro = 0.0
time_in_grav  = 0.0

log = logfile(filename=logname)
log.write("Number of hydro threads = "+`num_hy_workers`)
log.write("Step\t Tot time\t \t Hyd time\t \t Grav time\t \t Script time")

hydro.set_particle_pointers('mass')

try:

        while ((t < tmax) and (curr_hy_step < max_hy_steps)):
        #for i in range(1,2): ### For testing one iteration!


            i = i + 1

        ### Check for proper bridge timestep based on hydro timestep and crossing time.
        ### Have to write a proper routine to get timestep from hydro.
            #print "Checking dt. Current dt is :", dt
            hy_dt = hydro.get_timestep()
            print "Current hydro dt is: {:.2E}".format(hy_dt.value_in(units.s))
            #dt = min(dtinit, (2.5*hy_dt),(tmax-t))
            print "dt is = ", "{:.2e}".format(dt.value_in(units.s))
            if (first_step):
                dt = dtinit
                # Check for existing particles in case of restart.
                check_particles = hydro.get_number_of_particles()
                if (check_particles > 0):
                    first_particle = True

                    # If this is a restart, FLASH may still have all the
                    # particles mis-sorted in the particles array. Lets check
                    # for this.
                    print "Entering check_particles loop, check_particles=", check_particles
                    hydro.particles_sort()
                    check_particles = hydro.get_number_of_particles()
                    print "After sorting check_particles=", check_particles

                    print "Evolving grav to current hydro time."
                    hydro_time = hydro.get_time()

                    #stars, mult_grav, grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
                             #convert, stars = stars, start_time = hydro_time,
                             #num_grav_workers = num_grav_workers,
                             #eps = eps,
                             #with_ph4 = with_ph4,
                             #with_multiples = False)

                    #with Timer(verbose=True) as grav_timer:
                    #    if (with_ph4):
                    #        grav.parameters.begin_time=hydro_time
                    #        grav.parameters.sync_time=hydro_time
                    #        grav.parameters.force_sync=1
                    #    else:
                    #        grav.evolve_model(hydro_time)
                    grav_time = grav.get_time()
                    print "Hydro time:", hydro_time
                    print "Grav time:", grav_time
                    tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
                    num_particles = check_particles
                    if (with_multiples):
                        mult_grav, mult_to_stars, stars_to_mult = initialize_multiples(stars, grav, convert,
                                                                  mult_debug_level=mult_debug_level, kep=kep, new_smalln=new_smalln)
                        print "Gravity softening radius**2 =", grav.parameters.epsilon_squared
                    #if (with_se): se = SeBa()
                    if (not refresh_rand_seed_on_restart): first_call_for_stars = False
                first_step = False
            else:
                dt = min(dtmax, 1.5*hy_dt, (tmax-t), 2.0*dt_old)
            #dt = 2.0*hy_dt

            print "dt is now :", "{:.2e}".format(dt.value_in(units.s))
            dt_old = dt
            t_old  = t
            #t = t + dt

            print "Starting step ", i
            print "Number of stars = ", num_particles
            print "Num in grav = ", len(grav.particles)
            if (mult_grav is not None):
                 print "Num in multiples.stars = ", len(mult_grav.stars)
                 print "Num in multiples.root_to_tree = ", len(mult_grav.root_to_tree)
            print "Num in hydro = ", hydro.get_number_of_particles()
            sys.stdout.flush()
            #sys.exit()

            #if (tree_exists):
            #    print "Num in tree = ", len(tree.particles)

            made_stars,made_massive_star = make_stars_from_sinks2(hydro, min_sf_mass, max_sf_mass)
            sys.stdout.flush()
            print "Did we make stars?", made_stars

            if (made_stars):

                gridChanged = True
                if (first_particle == True):
                #    first_particle = True
                #    print "Evolving grav to current hydro time."
                #    hydro_time = hydro.get_time()
                #    with Timer(verbose=True) as grav_timer:
                #        grav.evolve_model(hydro_time)
                    tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
                    num_particles = hydro.get_number_of_particles()

        ### Wait for the first star to form.

            if (first_particle == False):

                check_particles = hydro.get_number_of_particles()

                print "Number of stars =", check_particles

                if check_particles != 0:

                    with Timer(verbose=True) as loop_timer:
                        first_particle = True
                        print "We got our first star!"
                        #print "Setting its mass = 90 MSun FOR THIS TEST CASE!!!"

                        print "Current simulation time:", t
                        hydro_time = hydro.get_time()

                        #stars, mult_grav, grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
                             #convert, stars = stars, start_time = hydro_time,
                             #num_grav_workers = num_grav_workers,
                             #eps = eps,
                             #with_ph4 = with_ph4,
                             #with_multiples = False)

                        grav_time  = grav.get_time()
                        print "Hydro time:", hydro_time
                        print "Grav time:", grav_time

                    # Try to evolve grav to the current time.
                        #print "Evolving grav to current hydro time."
                        #with Timer(verbose=True) as grav_timer:
                        #    #hydro_time = hydro.get_time()
                        #    if (with_ph4):
                        #        grav.parameters.begin_time=hydro_time
                        #        grav.parameters.sync_time=hydro_time
                        #        grav.parameters.force_sync=1
                        #    else:
                        #        grav.evolve_model(hydro_time)

                        print "Num particles in grav:", len(grav.particles)

                        tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
                        if (with_multiples):
                            mult_grav, mult_to_stars, stars_to_mult = initialize_multiples(stars, grav, convert,
                                                                      mult_debug_level=mult_debug_level, kep=kep, new_smalln=new_smalln)
                            print "Gravity softening radius**2 =", grav.parameters.epsilon_squared
                        #if (with_se): se = SeBa()
                        if (test_unique_tags):

                            test_tags(tags_keys)

                        print "Num particles in grav:", len(grav.particles)

                        if (tree_exists):
                            print "Num particles in tree:", len(tree.particles)

                        num_particles = check_particles

                        #hydro.set_particle_mass(tags_keys[0,0], 90 | units.MSun)
                        #stars.mass[0] = 90 | units.MSun
                        #stars_to_grav.copy()
                        #print "Hydro mass = ", hydro.get_particle_mass(tags_keys[:,0])
                        #print "Grav mass = ", grav.particles.mass
                    # Lets look for differences in the particle positions or velocities
                    # in the two codes.
                        if (insane):
                            if (with_multiples):
                                check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                            else:
                                check_sanity(hydro, grav)

                        #if (with_ph4):
                        #    grav.parameters.begin_time=hydro_time
                        #    grav.parameters.sync_time=hydro_time
                        #    grav.parameters.force_sync=1     # For PH4.
                        #else:
################    Not currently working for some reason!!!! (At least in Hermite) #########################
                            #grav.parameters.begin_time=hydro_time # For Hermite.
                            #grav.set_begin_time(hydro_time)
                        #    pass
                        print "Current simulation time:", t
                        hydro_time = hydro.get_time()
                        grav_time  = grav.get_time()
                        print "Hydro time:", hydro_time
                        print "Grav time:", grav_time


                    # Lets just try to evolve it before we add any particles above ^^^^

                    # Note this is important b/c if the gas isn't set up on the
                    # grid yet, the stars don't get kicked properly on the very first
                    # bridge step, since there is no gas gravity until Grid_solvePoisson is
                    # run at least once.
                        #if (first_loop):
                            #hydro.evolve_model(1.0e5 | units.s)  # Just set up, don't evolve.
                            #hydro.set_timestep(1.0e8 | units.s)
                            #first_loop = False


                else:

                    ### Global timer for one loop

                    with Timer(verbose=True) as loop_timer:
                        t = t + dt

                        if (do_sn_once and not first_loop):

                            inj_x = 0.0 | units.cm
                            inj_y = 0.0 | units.cm
                            inj_z = 0.0 | units.cm
                            tot_e = 1e51 | units.erg
                            fracKin = -1.0 #0.22734 #1.0
                            inj_mass = (5.*1.989e33) | units.g

                            print "BOOOOOOOMMMMMMM!"
                            dt = hydro.energy_injection(tot_e, fracKin, inj_mass, inj_x, inj_y, inj_z)
                            print "Timestep after SN is =", dt
                            t = t_old + dt
                            print "Now evolving until t =", t
                            do_sn_once = False
                            hydro.write_chpt()
                        print "I'm about to evolve hydro without evolving grav for :" , dt, "to evolve to t =", t

                        with Timer(verbose=True) as hydro_timer:
                            hydro.evolve_model(t) # Can't use i*dt if dt can get smaller each timestep.

                        if (first_step): first_step = False


                    # Note that if you are trying to output on nstep, Flash
                    # assumes this is checked during the Driver_evolveFlash loop
                    # and this won't output files properly. It will work
                    # properly if you use output time parameters, however.
                        hy_pltnum = hydro.IO_out('pltpart')
                        hy_chknum = hydro.IO_out('chk')

                        # If a checkpoint file was written,
                        # then store the list of stars up for creation.
                        if (hy_chknum != chknum):  # allow for possibility of rolling chk
                            write_rnd_and_mass_pickles(all_masses, output_dir, chknum)
                            chknum = hy_chknum

                        # wrote plt file
                        if (hy_pltnum > pltnum):
                            pltnum = hy_pltnum
                        elif (hy_pltnum < pltnum):
                            raise Exception("Error: hy_pltnum={} < pltnum={}".format(hy_pltnum, pltnum))

                    print "Current simulation time:", t
                    hydro_time = hydro.get_time()
                    grav_time  = grav.get_time()
                    print "Hydro time:", hydro_time
                    print "Grav time:", grav_time

                    #hydro.make_stars(made_stars)

                    if (first_loop): first_loop=False
                    time_in_hydro = time_in_hydro + hydro_timer.secs
                    #time_in_grav  = time_in_grav + grav_timer.secs
                    total_time    = total_time + loop_timer.secs

                    print "Total time in Flash = %f s" %time_in_hydro
                    print "Total time in N-body = %f s" %time_in_grav
                    print "Total time in AMUSE = %f s" %(total_time - time_in_grav - time_in_hydro)

                    log.write(`i`+"\t \t"+"{0:.2e}".format(total_time)+
                      "\t \t"+"{0:.2e}".format(time_in_hydro)+
                      "\t \t" +"{0:.2e}".format(time_in_grav)
                      + "\t \t" + "{0:.2e}".format(total_time - time_in_grav - time_in_hydro))

                    continue
                    #if (first_loop): first_loop=False


            ### Start the gravity bridge.

        ### Global timer for one loop

            with Timer(verbose=True) as loop_timer:

                print "Starting the gravity bridge."

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                #if (write_psets):
                #    stars.dt = dt
                #    write_set_to_file(stars, pdir+'/stars'+'BeforeSE'+'.amuse')

                if (with_se):

                    star_mass = np.zeros(num_particles) | units.MSun
                    star_age  = np.zeros(num_particles) | units.s
                    star_type = stars.stellar_type
                    dm_dt     = np.zeros(num_particles) | units.g / units.s
                    vterm     = np.zeros(num_particles) | units.cm / units.s
                    nphot     = np.zeros(num_particles) | units.s**-1.0
                    eion      = np.zeros(num_particles) | units.erg
                    sigh      = np.zeros(num_particles) | units.cm**2.0
                    npe       = np.zeros(num_particles) | units.s**-1.0
                    epe       = np.zeros(num_particles) | units.erg
                    sigpe     = np.zeros(num_particles) | units.cm**2.0

                    part_inds = []

                    #timing_1 = time.time()
                    star_mass = hydro.get_particle_mass(tags_keys[:,0])
                    #stars.mass = hydro.get_particle_mass(tags_keys[:,0])
                    #timing_2 = time.time()
                    #print "max star mass = ", star_mass.max().in_(units.MSun)
                    #print "getting mass took ", timing_2 - timing_1, "secs."

                    # Note creation time will be wrong if you do anything wacky like take a cluster from
                    # a different run or something weird.

                    # NOTE FURTHER: You must add in the dt here, otherwise newly formed stars will try to
                    # evolve for 0.0 seconds and things will get ugly.
                    star_age  = hydro.get_time() + dt - hydro.get_particle_creation_time(tags_keys[:,0])
                    #if ((stars.age).number.min() < star_age.number.min()):
                    #    print "Stars.age max =",(stars.age).number.max()
                    #    print "Star_age max =",(star_age).number.max()
                    #    stars.age = star_age
                    #    print "Using star_age instead of stars.age array."
                    #stars.age = star_age

                    print "Doing stellar evolution."
                    if (debug_se): print "min star age = ", star_age.min()

                    #print "Number of particles reported before SE stuff =", num_particles

                    for part in range(num_particles):

                        #print "Star age is", star_age.in_(units.Myr)
                        if (star_age[part].value_in(units.yr) < 1.0): star_age[part] = 1.0 | units.yr
                        #if (stars.age[part].value_in(units.yr) < 1.0): stars.age[part] = 1.0 | units.yr
                        #star_radius, star_temp, star_lum = stellar_properties(star_mass.in_(units.MSun),z)
                        if (debug_se):
                            print "this star initial mass before seba =", stars.initial_mass[part].in_(units.MSun)
                            print "this star mass before seba =", stars.mass[part].in_(units.MSun)
                            print "this star type before seba =", stars.stellar_type[part], stars.stellar_type[part].value_in(units.stellar_type)
                            print "this star age =", star_age[part].value_in(units.Myr)

                        # Do stellar evolution unless I already went SN, in which case skip me.
                        if (13 <=  stars.stellar_type[part].value_in(units.stellar_type) <= 15):
                            print "Skipping this star that already went SN, current stellar type =", stars.stellar_type[part]
                            continue
                        else:
                            # Note here we want to use the SE code on the initial mass of the star, not the current mass.
                            #st_time, st_mass, star_radius, star_lum, star_temp, st_evol_time, st_type = se.evolve_star(stars.initial_mass[part], star_age[part], 0.02)
                            st_time, st_mass, star_radius, star_lum, star_temp, st_evol_time, st_type = se.evolve_star(stars.initial_mass[part], star_age[part], 0.02)
                            star_type[part] = st_type
                            #stars.stellar_type[part] = st_type
                        if (debug_se):
                            print "this star mass after seba =", st_mass
                            print "this star type after seba =", st_type, st_type.value_in(units.stellar_type)
                            print "this star evolve time in seba =", star_age[part].value_in(units.Myr), st_evol_time
                            print "this star temp after seba =", star_temp
                            print "this star lum after seba =", star_lum

                        if (with_massloss and (massloss_method == 'seba' or st_mass.value_in(units.MSun) < min_mass.value_in(units.MSun))):
                            if (debug_se): print "using seba method"
                            dm_dt[part] = ((stars.mass[part]-st_mass)/st_time).in_(units.g / units.s)
                            # If below the mass cutoff for feedback, no wind present.
                            if (st_mass.value_in(units.MSun) < min_mass.value_in(units.MSun)):
                                vterm[part] = 0.0 | units.km/units.s
                            # Since we are using less certain mass loss rates anyway, just use velocity from Leitherer et. al. 1992.
                            else:
                                vterm[part] = 10**(1.23 - 0.30 * np.log10(star_lum.value_in(units.LSun))
                                        + 0.55*np.log10(st_mass.value_in(units.MSun))
                                        + 0.64*np.log10(star_temp.value_in(units.K))) | units.km/units.s
                        # Shouldn't Lietherer and Puls calculations use the old mass (stars[part].mass)?
                        elif (with_massloss and massloss_method == 'leit'):
                            # From Leitherer et. al. 1992.
                            if (debug_se): print "using Leither method"
                            dm_dt[part] = 10**(-24.06 + 2.45 * np.log10(star_lum.value_in(units.LSun))
                                        -1.10*np.log10(stars[part].mass.value_in(units.MSun))
                                        + 1.31*np.log10(star_temp.value_in(units.K))) | units.MSun/units.yr
                            vterm[part] = 10**(1.23 - 0.30 * np.log10(star_lum.value_in(units.LSun))
                                        + 0.55*np.log10(stars[part].mass.value_in(units.MSun))
                                        + 0.64*np.log10(star_temp.value_in(units.K))) | units.km/units.s
                        elif (with_massloss and massloss_method == 'puls'):
                            # Kudritzki and Puls winds, see Kudritzki & Puls 2000, Markova & Puls 2004, 2008 and Vink 2000
                            if (debug_se): print "using Kudritzki method"
                            star_wind   = stellar_wind(star_temp, stars[part].mass, star_lum, star_radius)
                            dm_dt[part] = star_wind.dm_dt.as_quantity_in(units.g / units.s)
                            vterm[part] = star_wind.vterm.as_quantity_in(units.cm / units.s)
                        elif (with_massloss and (massloss_method == 'test')):
                            if (debug_se): print "using constant wind mass loss and velocity from Weaver."
                            dm_dt[part] = 1e-6 | units.MSun / units.yr
                            vterm[part] = 2e3 | units.km / units.s
                        else:
                            if (debug_se): print "error: no method selected"
                            dm_dt[part] = 0.0 | units.g / units.s
                            vterm[part] = 0.0 | units.cm / units.s

                        # If with energy injection, check to see if anything went supernova. If so, inject 10^51 ergs of
                        # energy into the grid

                        if (13 <=  st_type.value_in(units.stellar_type) <= 15):
                            print "A star just went SN on you. Should be calling that SN code now!"

                        if (with_sn):

                            if (do_sn_once and not first_loop):

                                inj_x = 0.0 | units.cm
                                inj_y = 0.0 | units.cm
                                inj_z = 0.0 | units.cm
                                tot_e = 1e51 | units.erg
                                fracKin = -1.0 #0.22734 #1.0
                                inj_mass = (5.*1.989e33) | units.g

                                print "BOOOOOOOMMMMMMM!"
                                dt = hydro.energy_injection(tot_e, fracKin, inj_mass, inj_x, inj_y, inj_z)
                                print "Timestep after SN is =", dt
                                #t = t_old + dt
                                print "Now evolving until t =", t+dt
                                do_sn_once = False
                                #hydro.IO_out('pltpart')

                            else:

                                #for part in range(num_particles):

                                #print "Stellar type =", stars.stellar_type[part]

                                if (13 <=  st_type.value_in(units.stellar_type) <= 15):

                                    print "Going supernova at", stars.x[part], stars.y[part], stars.z[part]

                                    inj_x = stars.x[part]
                                    inj_y = stars.y[part]
                                    inj_z = stars.z[part]
                                    tot_e = 1e51 | units.erg
                                    fracKin = -1.0 #0.22734 #1.0
                                    # Here the injected mass is the current mass minus the remnant mass.
                                    inj_mass = (star_mass[part] - st_mass).in_(units.g)
                                    #inj_mass = (stars.mass[part] - st_mass).in_(units.g)
                                    # This should never be more than 10 solar masses though.
                                    if (inj_mass.value_in(units.MSun) > 10.0):
                                        print "[bridge:SN]: WARNING! SN MASS > 10 MSun!"
                                        inj_mass = 10.0 | units.MSun
                                    print "Injection mass is", inj_mass.in_(units.MSun)

                                    dt =  min(dt.value_in(units.s), hydro.energy_injection(tot_e, fracKin, inj_mass, inj_x, inj_y, inj_z).value_in(units.s)) | units.s
                                    print "Timestep after SN is =", dt
                                    #t = t_old + dt
                                    print "Now evolving until t =", t + dt

                                    # Set proper mass for remnant (SeBa does fine for this).
                                    star_mass[part] = st_mass
                                    # Set proper remnant stellar type so that we don't get any feedback from remnants.
                                    stars.stellar_type[part] = st_type
                                    #stars.mass[part] = st_mass
                                    # Switch off all feedback for this star.
                                    nphot[part] = 0.0 | units.s**-1
                                    eion[part]  = 0.0 | units.erg
                                    if (use_radiation):
                                        hydro.set_particle_nion(tags_keys[part,0], nphot[part])
                                        hydro.set_particle_eion(tags_keys[part,0], eion[part])
                                    npe[part] = 0.0 | units.s**-1
                                    epe[part] = 0.0 | units.eV
                                    if (pe_heat):
                                        hydro.set_particle_npep(tags_keys[part,0], npe[part])
                                        hydro.set_particle_epep(tags_keys[part,0], epe[part])
                                    # Winds don't matter for a SN star and we don't want to lose any more mass.
                                    dm_dt[part] = 0.0 | units.g / units.s
                                    vterm[part] = 0.0 | units.cm / units.s
                                    if (with_winds): # and (star_mass.in_(units.MSun) >= 10.0 | units.MSun)):
                                        #print "Setting wind with dm_dt =", dm_dt[np.where(dm_dt.value_in(units.MSun / units.yr) > 0.0)].in_(units.MSun / units.yr)
                                        #print "Setting wind velocity =", vterm[np.where(vterm.value_in(units.km / units.s) > 0.0)].in_(units.km / units.s)
                                        hydro.set_particle_wind_mass(tags_keys[part,0], dm_dt[part])
                                        hydro.set_particle_wind_vel(tags_keys[part,0], vterm[part])
                                        #print "mass loss should be ", dm_dt[part_inds].in_(units.MSun/units.s)*dt.in_(units.s), "for star of mass ", stars.mass[part_inds]
                                    # Do nothing else with this star, just jump to the next one.
                                    continue

                        if ((star_mass[part].in_(units.MSun) >= min_mass) and not (13 <=  st_type.value_in(units.stellar_type) <= 15)):

                            print "Found massive star, star mass =", star_mass[part].in_(units.MSun)

                            part_inds.append(part)


                            if (use_radiation):

                                #print "Entering radiation calculation."

                                flux = ion.ionizing_photon_flux(st_mass, star_radius, star_temp)

                                # Calculate the average ionizing photon energy based on the blackbody curve.

                                # First integrate the power from the BB curve at this stars temp.
                                # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                                [power, err] = quad(lum_wl_cs, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                                #print "The ionizing energy flux and error for this star is:"
                                #print power, err
                                # Now integrate to find the number of photons.
                                [per_ph, err] = quad(lum_wl_cs_per_ph, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                                #print "The average ionizing photon number flux and error for this star is:"
                                #print per_ph, err
                                #print "The average energy per ionizing photon for this star is:"
                                #print power/per_ph
                                avg_E = power/per_ph / E_ev
                                #print avg_E
                                # Calculate the average frequency of an ionizing photon for this star
                                avg_nu = avg_E*E_ev/h
                                #print "The average frequency of an ionizing photon for this star and the min ionizing frequency is:"
                                #print avg_nu, nu_min
                                # Cross section calculation
                                # Make sure you convert energy back to ergs if you
                                # use it to calculate the frequency!
                                sig = sig0*(avg_nu/nu_min)**(-3.0)
                                #print "The cross section for these photons is:"
                                #print sig
                                eion[part] = (avg_E | units.eV) - (13.6 |units.eV) #2.0 | units.eV #6.0 | units.eV
                                sigh[part] = sig | units.cm**2.0 #6.3e-18 | units.cm**2.0
                                nphot[part] = (flux*4*np.pi*star_radius**2.0).as_quantity_in(units.s**-1.0) #5e48 | units.s**(-1.0)


                                if (pe_heat):

                                    # First integrate the power from the BB curve at this stars temp.
                                    # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                                    l_min_dust = h*c / E_min # wavelength at 13.6 eV
                                    l_max_dust = h*c / (5.6*E_ev) # wavelength at 5.6 eV
                                    [power, err] = quad(lum_wl, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))
                                    #print "The ionizing energy flux and error for this star is:"
                                    #print power, err
                                    # Now integrate to find the number of photons.
                                    [per_ph, err] = quad(lum_wl_per_ph, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))
                                    #print "The average ionizing photon number flux and error for this star is:"
                                    #print per_ph, err
                                    #print "The average energy per ionizing photon for this star is:"

                                    avg_E = power/per_ph / E_ev
                                    #print "Avg PE energy (eV) =", avg_E
                                    # Calculate the average frequency of an ionizing photon for this star
                                    #avg_nu = avg_E*E_ev/h
                                    #print "The average frequency of an ionizing photon for this star and the min ionizing frequency is:"
                                    #print avg_nu, nu_min
                                    # Cross section calculation
                                    # We assume constant cross section for dust per hydrogen atom.
                                    # Value = tau / N_H where tau = gamma * Av (Draine and Bertoli 96)
                                    # Av = N_H,tot / (1.87e21 cm^2) (Bohlin et al 78)
                                    # gamma = 2.5 (Bergin et al 2004)
                                    sigpe[part] = sigDust
                                    #print "The cross section for these photons is:"
                                    #print sig
                                    # Eion is the actual average energy of the photons WITH the ionizing potential still in there!
                                    epe[part] = avg_E | units.eV # should be around 8 eV
                                    #sigh = sig | units.cm**2.0 #6.3e-18 | units.cm**2.0
                                    # Calculate total number of photons from stellar surface with stellar radius.
                                    npe[part] = ((per_ph | units.cm**-2*units.s**-1)*4*np.pi*star_radius**2.0).as_quantity_in(units.s**-1.0) #5e48 | units.s**(-1.0)

                        if ((dm_dt[part]*dt).value_in(units.MSun) > 0.0):
                            if (st_type.value_in(units.stellar_type) == 1):
                                star_mass[part] = ((stars.mass[part] - dm_dt[part]*dt).value_in(units.MSun)) | units.MSun
                            # Note other evolutionary things besides winds could have reduced the stars mass.
                            else:
                                star_mass[part] = min(st_mass.value_in(units.MSun), (stars.mass[part] - dm_dt[part]*dt).value_in(units.MSun)) | units.MSun
                        else:
                            star_mass[part] = st_mass

                        if (with_multiples):
                            # Cycle through the leaves and check if this particular star
                            # is in the leaves in multiples.
                            for root, tree in mult_grav.root_to_tree.iteritems():
                                leaves = tree.get_leafs_subset()
                                if stars[part] in leaves:
                                    if (debug_multiples and debug_se):
                                        print "Star evolve time =", stars.age[part]
                                        print "Before update from stars to multiples leaves."
                                        print "Stars mass =", stars[part].mass.in_(units.MSun)
                                        print "SE returned mass =", star_mass[part].in_(units.MSun)
                                        print "Leaf mass =", leaves.mass.in_(units.MSun)
                                        print "Root mass =", root.mass.in_(units.MSun)
                                        print "Star pos =", stars[part].position.in_(units.cm)
                                        print "Leaf pos =", leaves.position.in_(units.cm)
                                        print "Root pos =", root.position.in_(units.cm)
                                        print "Star vel =", stars[part].velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                    # If yes, then set the mass of the star to the updated to the new mass from
                                    # the SE code. This should be possible by using the as_particle_in_set(leaves)
                                    # method on the stars particle. We'll check with print statements.
                                    stars[part].as_particle_in_set(leaves).mass = star_mass[part]

                                    if (debug_multiples and debug_se):
                                        print "After update from stars to multiples leaves."
                                        print "Stars mass =", stars[part].mass.in_(units.MSun)
                                        print "SE returned mass =", star_mass[part].in_(units.MSun)
                                        print "Leaf mass =", leaves.mass.in_(units.MSun)
                                        print "Root mass =", root.mass.in_(units.MSun)
                                        print "Star pos =", stars[part].position.in_(units.cm)
                                        print "Leaf pos =", leaves.position.in_(units.cm)
                                        print "Root pos =", root.position.in_(units.cm)
                                        print "Star vel =", stars[part].velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                        # Don't forget at the end of this we need to update the masses
                                        # in multiples and the gravity code using the proper root mass
                                        # if that didn't happen automatically.

                    # If the leaf mass changed, so does the com particle's properties.
                    if (with_multiples):
                        update_roots_from_leaves(mult_grav, grav)
                        if (debug_multiples):
                            # Now check the changed particles.
                            check_root_and_leaves(mult_grav, grav, stars)
                            if (debug_se):
                                for st, st_mass in zip(stars, star_mass):  # reusing st_mass...
                                    for root, tree in mult_grav.root_to_tree.iteritems():
                                        leaves = tree.get_leafs_subset()
                                        if st in leaves:
                                            print "Star evolve time =", st.age
                                            print "After update from multiples leaves to root."
                                            print "Stars mass =", st.mass.in_(units.MSun)  # beware st.mass != st_mass
                                            print "SE returned mass =", st_mass.in_(units.MSun)
                                            print "Leaf mass =", leaves.mass.in_(units.MSun)
                                            print "Root mass =", root.mass.in_(units.MSun)
                                            print "Star pos =", st.position.in_(units.cm)
                                            print "Leaf pos =", leaves.position.in_(units.cm)
                                            print "Root pos =", root.position.in_(units.cm)
                                            print "Star vel =", st.velocity.in_(units.km/units.s)
                                            print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                            print "Root vel =", root.velocity.in_(units.km/units.s)
                                            sys.stdout.flush()

                    # Are there any massive stars?

                    if not (part_inds==[]):

                        if (use_radiation):

                            print  tags_keys[part_inds,0]
                            print "Stellar Mass and N photons=", star_mass[part_inds], nphot[part_inds]
                            print "Eion (eV), SigH=", eion[part_inds], sigh[part_inds]
                            hydro.set_particle_nion(tags_keys[part_inds,0], nphot[part_inds])
                            hydro.set_particle_eion(tags_keys[part_inds,0], eion[part_inds].as_quantity_in(units.erg))
                            hydro.set_particle_sigh(tags_keys[part_inds,0], sigh[part_inds])


                        if (pe_heat):
                            print "Npe photons=", npe[part_inds]
                            print "Eion PE (eV), SigD=", epe[part_inds], sigpe[part_inds]
                            hydro.set_particle_npep(tags_keys[part_inds,0], npe[part_inds])

                            # Set average energy of PE photon
                            hydro.set_particle_epep(tags_keys[part_inds,0], epe[part_inds].as_quantity_in(units.erg))
                            # Set cross section of dust to PE photons.
                            hydro.set_particle_sigd(tags_keys[part_inds,0], sigpe[part_inds])


                        if (with_winds): # and (star_mass.in_(units.MSun) >= 10.0 | units.MSun)):

                            #print "Setting wind with dm_dt =", dm_dt[np.where(dm_dt.value_in(units.MSun / units.yr) > 0.0)].in_(units.MSun / units.yr)
                            #print "Setting wind velocity =", vterm[np.where(vterm.value_in(units.km / units.s) > 0.0)].in_(units.km / units.s)
                            hydro.set_particle_wind_mass(tags_keys[part_inds,0], dm_dt[part_inds])
                            hydro.set_particle_wind_vel(tags_keys[part_inds,0], vterm[part_inds])
                            #print "mass loss should be ", dm_dt[part_inds].in_(units.MSun/units.s)*dt.in_(units.s), "for star of mass ", stars.mass[part_inds]

                #if (made_massive_star):
                    #print "Massive star made, checking if dt should be reduced. Current dt is", dt
                    #dt = min(dt.value_in(units.s), 0.3*smallest_dx.value_in(units.cm)/np.max(vterm.value_in(units.cm/units.s))) | units.s
                    #print "dt now is", dt

                # Remove any mass loss due to winds and update to this
                # mass. Note this assumes steps are relatively small
                # in the mass loss rate of stars, so that graviationally
                # we can use the mass after all the wind mass loss
                # has occcured. Otherwise we'd have to average
                # mass loss and keep up with old and new masses and
                # it just gets ugly.
                stars.mass = star_mass # - dm_dt*dt
                stars.age  = stars.age + dt
                stars.stellar_type = star_type
                #print "new star mass = ", stars.mass

                # Note when using the multiples module we can't copy stars to grav
                # because the stars are all the particles (including leaves) and only
                # the roots are in grav.
                # Instead we should modify the values in multiples directly such that:
                # 1. Leaves have updated mass
                # 2. Roots have updated mass and velocity.
                # 3. We then need to update the radii of all the particles (to update the root radii).
                # 4. We then copy from multiples roots to grav code.

                if (with_multiples):
                    pass
                    #mult_grav.channel_from_memory_to_code.copy()
                    #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                else:
                    stars_to_grav.copy()
                #stars_to_grav.copy_attributes(["mass"])
                hydro.set_particle_mass(tags_keys[:,0], stars.mass)

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)


                # Set the bridge timestep.
                t = t + dt
                print "I'm about to evolve hydro and grav for :" , dt, "to evolve to t =", t

                #if (np.abs(hydro.get_particle_mass(tags_keys[:,0]).value_in(units.MSun) - grav.particles.mass.value_in(units.MSun)).any() > 1.0):

                    #print "Masses not equal at end of SE, stopping code."
                    #sys.exit()

            ### First kick.

                # Before the kick, lets look for differences in the particle positions or velocities
                # in the two codes.

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                #print "Saving hydro positions"
                #old_pos = hydro.get_particle_position(tags_keys[:,0]).in_(units.m)

                #print "Before first kick, max star velocity = ", np.sqrt((stars.vx**2.0 + stars.vy**2.0 + stars.vz**2.0).value_in(units.km**2.0/units.s**2.0)).max() | units.km / units.s

                if (with_bridge):

                    print "First kick."
                    step = 1
                    if (gridChanged):
                        step = 2
                        gridChanged = False

                    bridge_kick2(hydro, stars, eps, dt, time_in_hydro, step)

                    #print "Grav vel before update."
                    #print grav.particles.velocity

                    if (with_multiples):
                        # First copy updated velocities from stars to grav for
                        # all particles in grav that are also in stars
                        # (i.e. things that are leaves in grav and stars).
                        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
                        mult_grav.channel_from_code_to_memory.copy()
                        # Now check for any particle that is not in grav (that is
                        # a leaf in multiples but replaced by a root particle in grav).
                        # Update the velocity of the leaves with that in stars and copy
                        # this over to grav.
                        for st in stars:
                            # Cycle through the leaves and check if this particular star
                            # is in the leaves in multiples.
                            for root, tree in mult_grav.root_to_tree.iteritems():
                                leaves = tree.get_leafs_subset()
                                if st in leaves:
                                    if (debug_multiples):
                                        print "After kick 1."
                                        print "Star vel =", st.velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                    # If yes, then set the mass of the star to the updated to the new mass from
                                    # the SE code. This should be possible by using the as_particle_in_set(leaves)
                                    # method on the stars particle. We'll check with print statements.
                                    st.as_particle_in_set(leaves).velocity = st.velocity
                                    #com_vel = 0.0
                                    #for leaf in leaves:
                                    #    com_vel += (leaf.mass.value_in(units.MSun)*leaf.velocity.value_in(units.km/units.s))
                                    #com_vel = (com_vel/(leaves.mass.sum()).value_in(units.MSun)) | units.km/units.s
                                    #root.velocity = com_vel
                                    #root.velocity = (np.array([l.mass*l.velocity for l in leaves]).sum()/leaves.mass.sum()).in_(units.km/units.s)
                                    #root.velocity = (leaves.mass[:,np.newaxis]*leaves.velocity)/leaves.mass.sum()
                                    if (debug_multiples):
                                        print "After kick 1."
                                        print "Star vel =", st.velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                    # Don't forget at the end of this we need to update the velocity
                                    # in multiples and the gravity code using the proper root mass
                                    # if that didn't happen automatically.

                        # Update pset shouldn't be needed now that we use update_roots_from_leaves and the above code.
                        #if (debug_multiples):
                        #    "Trying using update_pset."
                        #update_psetA_from_psetB(mult_grav.stars, stars)

                        update_roots_from_leaves(mult_grav, grav)
                        if (debug_multiples):
                            print "First kick after update_roots_from_leaves."
                            check_root_and_leaves(mult_grav, grav, stars)
                            for st in stars:
                                # Cycle through the leaves and check if this particular star
                                # is in the leaves in multiples.
                                for root, tree in mult_grav.root_to_tree.iteritems():
                                    leaves = tree.get_leafs_subset()
                                    if st in leaves:
                                        if (debug_multiples):
                                            print "Star vel =", st.velocity.in_(units.km/units.s)
                                            print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                            print "Root vel =", root.velocity.in_(units.km/units.s)
                                            sys.stdout.flush()
                        #mult_grav.channel_from_memory_to_code.copy() # Now copy from the memory of multiples to grav.
                        #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                    else:
                        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
                        #stars_to_grav.copy()

                    print "Grav updated."
                    #print grav.particles.velocity

                    #print "After first kick, max star velocity = ", np.sqrt((stars.vx**2.0 + stars.vy**2.0 + stars.vz**2.0).value_in(units.km**2.0/units.s**2.0)).max() | units.km / units.s

                    #print "Checking hydro positions"
                    #print old_pos - hydro.get_particle_position(tags_keys[:,0]).in_(units.m)

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

            ### Get the location of the particles so that we can get the
            ### gravitational acceleration from the gas at these locations from
            ### Flash.

            #    loc = grav.particles.position

            ### Record the potential at the particle locations in GPOT_PART_PROP.
            ### Note that AMUSE's original bridge implementation requires access to this
            ### information.

                #gpot = hydro.get_potential_at_point(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])

            ### Set the gravitational potential at the sink locations in Flash
            ### so we have access to this data for analysis later.

                #hydro.set_particle_gpot(tags_keys[:,0], gpot)

            ### Get the gravitational acceleration from the gas at the particle
            ### locations.

                #gaccel = hydro.get_gravity_at_point(0.0 | units.m, loc[:,0], loc[:,1], loc[:,2])

                #print "Got gravity at star locations."

                #for k in range(len(grav.particles)):

                    #stars[k].vx = stars[k].vx + 0.5*dt*gaccel[0][k]
                    #stars[k].vy = stars[k].vy + 0.5*dt*gaccel[1][k]
                    #stars[k].vz = stars[k].vz + 0.5*dt*gaccel[2][k]

            #### Update the gravity code with the kicked velocities.

                #stars_to_grav.copy()
            ##   stars_to_grav.copy_attributes(["vx", "vy", "vz"])

            #### First kick on the gas grid.

                #hydro.kick_grid(0.5*dt)




            #   Update FLASH's sink velocities (from the kick), because it makes plot files during evolve.
            #   NOTE: Plotting should probably be moved to outside of the FLASH evolve and into here.

            # Something in here is not right... Flash is reporting the correct velocities to AMUSE but the plot
            # file doesn't have the correct updated velocities from the particles... have to actually print this during
            # a Flash evolve step to check it!

                #hydro.set_particle_velocity(tags_keys[:,0], grav.particles.vx, grav.particles.vy, grav.particles.vz)

                #print hydro.get_particle_velocity(tags_keys[:,0]).as_quantity_in(units.m*units.s**(-1))


            #   Record the old pos/vel to compare against hydro evol AFTER we
            #   bridge kick but BEFORE evol, so that we can properly find the
            #   changes from accretion in hydro. This will ensure momentum conservation.

                #old_r = stars.position
                #old_v = stars.velocity


                #if (np.abs(hydro.get_particle_mass(tags_keys[:,0]).value_in(units.MSun) - grav.particles.mass.value_in(units.MSun)).any() > 1.0):

                    #print "Masses not equal before evolve, stopping code."
                    #sys.exit()

                #if (with_multiples and insane):
                #    print "Storing the original stars and mult_grav.star."
                #    oldstars = stars.copy_to_new_particles()
                #    oldmultstars = mult_grav.stars.copy_to_new_particles()

                #    print "Checking originals vs originals for repeats (should ALL be repeats!)"
                #    oldvnew_position(stars, oldstars)
                #    oldvnew_position(mult_grav.stars, oldmultstars)


            ### Evolve models.

                print "Evolving models."

            #    request1 = hydro.evolve_model.async(t)
            #    request2 = grav.evolve_model.async(t)

            ### Submit the request for async evolution.

            #    pool.add_request(request1, handle_result, [1])
            #    pool.wait()

            #    pool.add_request(request2, handle_result, [2])

            ### Wait until all requests are finished.


            #    pool.wait()

            #    pool.waitall()
                print "Before evolve, max star velocity = ", stars.velocity.norm().max().in_(units.km/units.s)
                print "Before evolve, max grav velocity = ", grav.particles.velocity.norm().max().in_(units.km/units.s)
                if (with_multiples and insane):
                    print "Before evolve, max multiples.stars velocity = ", mult_grav.stars.velocity.norm().max().in_(units.km/units.s)

                print "Calling grav."
                with Timer(verbose=True) as grav_timer:
                    if (with_multiples):
                        mult_grav.evolve_model(t)
                    else:
                        grav.evolve_model(t)

                grav_evolve_time = grav.get_time()

                if (t > grav_evolve_time):
                    print "WARNING: grav didn't evolve properly. Try again?"
                    print "Calling grav."
                    with Timer(verbose=True) as grav_timer:
                        if (with_multiples):
                            mult_grav.evolve_model(t)
                        else:
                            grav.evolve_model(t)

                print "Calling hydro."
                with Timer(verbose=True) as hydro_timer:
                    hydro.evolve_model(t)

                if (first_step): first_step = False

                hydro_time = hydro.get_time()
                grav_time  = grav.get_time()

                time_diff = (hydro_time - grav_time).value_in(units.s)

                if (time_diff > 1e4):

                    print "Warning: Grav is behind Hydro by:", time_diff
                    print "Trying again to evolve Grav to Hydro time."
                    with Timer(verbose=True) as grav_timer:
                        if (with_multiples):
                            mult_grav.evolve_model(hydro_time)
                        else:
                            grav.evolve_model(hydro_time)

                # Update the current hydro step.
                curr_hy_step = hydro.get_current_step()

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, kill=False, with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars, kill=False)
                    else:
                        check_sanity(hydro, grav, kill=False)

                print "After evolve / before update, max star velocity = ", stars.velocity.norm().max().in_(units.km/units.s)
                print "After evolve / before update, max grav velocity = ", grav.particles.velocity.norm().max().in_(units.km/units.s)
                if (with_multiples):
                    print "After evolve / before update, max multiples.stars velocity = ", mult_grav.stars.velocity.norm().max().in_(units.km/units.s)


                if (with_multiples):
                    if (debug_multiples):
                        print "###########################"
                        print "Before grav_to_stars"
                        print "###########################"
                        check_root_and_leaves(mult_grav, grav, stars, kill=False)
                # Always call this now, followed by updating any stars that
                # from leaves in multiples if they are present.
                grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])

                ### Update the positions and velocities in stars from gravity.
                #grav_to_stars.copy()

                if (with_multiples):
                    ### Lets force updating (again) from the grav to the
                    ### _inmemory_particles just to be sure.
                    #mult_grav.channel_from_code_to_memory.copy()

                    if (debug_multiples):
                        print "###########################"
                        print "Before update_leaves_pos_vel"
                        print "###########################"
                        check_root_and_leaves(mult_grav, grav, stars, kill=False)
                    mult_grav.update_leaves_pos_vel()

                    if (debug_multiples):
                        print "###########################"
                        print "Before any update to stars."
                        print "###########################"
                        check_root_and_leaves(mult_grav, grav, stars, kill=False)

                    mult_grav.stars.copy_values_of_attributes_to(["x", "y", "z", "vx", "vy", "vz"], stars)
                    if (debug_multiples):
                        print "###########################"
                        print "After copy_values_of_attributes."
                        print "###########################"
                        check_root_and_leaves(mult_grav, grav, stars, kill=True)


                    #update_stars_from_leaves(mult_grav, grav, stars)
                    #if (debug_multiples):
                    #    print "###########################"
                    #    print "After update_stars_from_leaves"
                    #    print "###########################"
                    #    check_root_and_leaves(mult_grav, grav, stars, kill=False)

                    #mult_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])
                    #update_psetA_from_psetB(stars,mult_grav.stars,debug=True)
                    #if (debug_multiples):
                    #    print "###########################"
                    #    print "After update_psetA_from_psetB."
                    #    print "###########################"
                    #    check_root_and_leaves(mult_grav, grav, stars, kill=False)


                    #mult_grav.channel_from_memory_to_code.copy() # Now copy from the memory of multiples to grav.
                    #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                else:
                    grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])

                # Lets look for differences in the particle positions or velocities
                # in the two codes.

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, kill=False,  with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                if (test_unique_tags):
                    test_tags(tags_keys)

                print "After evolve / after update, max star velocity = ", stars.velocity.norm().max().in_(units.km/units.s)
                print "After evolve / after update, max grav velocity = ", grav.particles.velocity.norm().max().in_(units.km/units.s)
                if (with_multiples):
                    print "After evolve / after update, max multiples.stars velocity = ", mult_grav.stars.velocity.norm().max().in_(units.km/units.s)
                    #if (insane):
                    #    #check_root_and_leaves(mult_grav, grav, stars)
                    #    #print "Checking new stars against originals (should be different)."
                    #    #oldvnew_position(stars, oldstars)
                    #    #oldvnew_position(mult_grav.stars, oldmultstars)

            #   Calculate the change in position and velocity of sinks in hydro
            #   that was due to accretion. Note the comparison must be made
            #   against the old positions and velocities after the bridge kick but
            #   before hydro evolve.

                #hy_r  = hydro.get_particle_position(tags_keys[:,0]).in_(units.m)
                #hy_v  = hydro.get_particle_velocity(tags_keys[:,0]).in_(units.m*units.s**(-1))

                #dx    = np.subtract(hy_r, old_r)
                #dv    = np.subtract(hy_v, old_v)


            ### Accretion occured in Flash, so get the new mass of the stars to the gravity code.
            ### Note it is important to do this now, because the mass has been removed from the gas
            ### already and won't be there for the next bridge kick on the gravity code.

                #print "Before accretion mass is ", stars.mass.as_quantity_in(units.MSun)
                #print "Hydro reports mass = ", hydro.get_particle_mass(tags_keys[:,0])
                #stars.mass = hydro.get_particle_mass(tags_keys[:,0])

                #print "After accretion mass is ", stars.mass.as_quantity_in(units.MSun)

                #print "Tags keys =", tags_keys
            ### Copy the updated positions and velocities from the gravity to the stars.

                #stars.x  = np.add(stars.x, dx[:,0])
                #stars.y  = np.add(stars.y, dx[:,1])
                #stars.z  = np.add(stars.z, dx[:,2])
                #stars.vx = np.add(stars.vx, dv[:,0])
                #stars.vy = np.add(stars.vy, dv[:,1])
                #stars.vz = np.add(stars.vz, dv[:,2])

            ### Shouldn't we also update the position and velocity of the stars
            ### in hydro here before we kick again??? (Note PM bridge uses star
            ### positions to map the mass of the stars to the grid).
            ### DONE A FEW MORE LINES BELOW.

                # Can't set grav.particles = stars, it just makes
                # grav.particles a pointer to stars (and not a Hermite code instance).
                # since copy() doesn't preserve sorting by tag.
                #stars_to_grav.copy()
                # grav.particles = stars


                #if (with_se):
                    #stars_to_se.copy()
                    ##se.particles.mass = stars.mass.as_quantity_in(units.MSun)
                    #print "Star mass", stars.mass.in_(units.MSun)
                    #print "SE mass", se.particles.mass.in_(units.MSun)

                #if (tree_exists):
                #    stars_to_tree.copy()


                #if (np.abs(hydro.get_particle_mass(tags_keys[:,0]).value_in(units.MSun) - grav.particles.mass.value_in(units.MSun)).any() > 1.0):

                    #print "Masses not equal after correcting momentum, stopping code."
                    #sys.exit()


                # If not using sinks, call star creation routine.
                #if (no_sinks):
                #    hydro.make_stars(dt)

            ### Update the star locations in FLASH so that the grid is kicked at the proper place by
            ### the stars. This is required because the grid kick all happens internally
            ### inside Flash since that is much faster. (Do this BEFORE hydro.kick_grid!)

            ### Also note, we need hydro to see the latest particle positions so that it can remove anything
            ### that has left the computational domain when we call hydro.sort_particles below.

                if (with_multiples):
                    #mult_grav.update_leaves_pos_vel()
                    hydro.set_particle_position(tags_keys[:,0], stars.x, stars.y, stars.z)
                    hydro.set_particle_velocity(tags_keys[:,0], stars.vx, stars.vy, stars.vz)
                else:
                    hydro.set_particle_position(tags_keys[:,0], grav.particles.x, grav.particles.y, grav.particles.z)
                    hydro.set_particle_velocity(tags_keys[:,0], grav.particles.vx, grav.particles.vy, grav.particles.vz)

                if (test_unique_tags):

                    test_tags(tags_keys)


                    print "Num particles in grav:", len(grav.particles)

                    #num_particles = check_particles

                #else:

                    #print "No particles created during step ", i


            ### Check: Are the particles in sync across the two codes?

                if (tags_keys[:,0].all() == hydro.get_particle_tags(range(1,num_particles+1)).all() and
                    tags_keys[:,1].all() == grav.particles.key.all()):

                    pass
                    #print "Particles in sync."

                else:

                    print "Particles out of sync! Stopping!"
                    print "hydro tags"
                    print hydro.get_particle_tags(range(1,num_particles+1))
                    print "tags tags"
                    print tags_keys[:,0]
                    sys.stdout.flush()
                    sys.exit()


                # Lets look for differences in the particle positions or velocities
                # in the two codes.

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                if (test_unique_tags):
                    test_tags(tags_keys)

                # Remove any particles that have left the simulation domain.
                stars_removed, num_particles  = remove_particles_outside_bndbox(hydro, stars, grav, mult_grav, with_multiples, tags_keys, bndbox, debug=True)

                # Sync stars to the particles in grav.
                if (with_multiples):
                    grav.particles.synchronize_to(mult_grav._inmemory_particles)
                    #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
                else:
                    grav.particles.synchronize_to(stars)


                if (test_unique_tags):
                    test_tags(tags_keys)

                # If stars were removed, resort particles and check if we need to return to hydro only.
                # NOTE: I think this is important to include here and not in the actual particle removal
                # from hydro. This is because the following steps can occur:
                # 1. Star A leaves out of the x-dir boundary and star B leaves out of the y dir boundary.
                # 2. The x check is done and star A is found outside and removed from FLASH. But when FLASH
                #    sorts the particles, it also removes star B since the code is smart enough to do this.
                # 3. The bridge module now tries to remove star B from FLASH, and FLASH complains that the
                #    block for this star is block 0 (because removing it from FLASH just means zeroing
                #    its properties.
                # But also note, I do think we should still sort at the end, because it sets all the particle
                # info in pt_typeInfo array correctly. - JW

                hydro.particles_sort()
                print "Hydro particles sort called."
                print "Now ", len(grav.particles), "particles in grav."
                print "and", hydro.get_number_of_particles(), "particles in hydro."
                print "and", num_particles, "particles in stars."
                if (stars_removed):
                    if (num_particles == 0):
                        print "No particles left! Going back to hydro only!"
                        first_particle = False
                        continue


                if (test_unique_tags):

                    test_tags(tags_keys)
                # Lets look for differences in the particle positions or velocities
                # in the two codes.

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                #if (np.abs(pos_diff.value_in(units.m)).any() > min_pos_diff or np.abs(vel_diff.value_in(units.m*(units.s**-1))).any() > min_vel_diff):

                #    print "Difference in position or velocity greater than tolerance. Stopping code now."
                #    sys.exit()

                #print "Before second kick, max star velocity = ", np.sqrt((stars.vx**2.0 + stars.vy**2.0 + stars.vz**2.0).value_in(units.km**2.0/units.s**2.0)).max() | units.km / units.s

                if (with_bridge):

                    print "Second kick."
                    step = 2
                    bridge_kick2(hydro, stars, eps, dt, time_in_hydro, step)

                    #print "Grav vel before update."
                    #print grav.particles.velocity

                    if (with_multiples):
                        # First copy updated velocities from stars to grav for
                        # all particles in grav that are also in stars.
                        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
                        mult_grav.channel_from_code_to_memory.copy()
                        # Now check for any particle that is not in grav (that is
                        # a leaf in multiples but replaced by a root particle in grav).
                        # Update the velocity of the leaves with that in stars and copy
                        # this over to grav.
                        for st in stars:
                            # Cycle through the leaves and check if this particular star
                            # is in the leaves in multiples.
                            for root, tree in mult_grav.root_to_tree.iteritems():
                                leaves = tree.get_leafs_subset()
                                if st in leaves:
                                    if (debug_multiples):
                                        print "After kick 2."
                                        print "Star vel =", st.velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                    # If yes, then set the mass of the star to the updated to the new mass from
                                    # the SE code. This should be possible by using the as_particle_in_set(leaves)
                                    # method on the stars particle. We'll check with print statements.
                                    st.as_particle_in_set(leaves).velocity = st.velocity
                                    #com_vel = 0.0
                                    #for leaf in leaves:
                                    #    com_vel += (leaf.mass.value_in(units.MSun)*leaf.velocity.value_in(units.km/units.s))
                                    #com_vel = (com_vel/(leaves.mass.sum()).value_in(units.MSun)) | units.km/units.s
                                    #root.velocity = com_vel
                                    #root.velocity = (np.array([l.mass*l.velocity for l in leaves]).sum()/leaves.mass.sum()).in_(units.km/units.s)
                                    #root.velocity = (leaves.mass[:,np.newaxis]*leaves.velocity)/leaves.mass.sum()
                                    if (debug_multiples):
                                        print "After kick 2."
                                        print "Star vel =", st.velocity.in_(units.km/units.s)
                                        print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                        print "Root vel =", root.velocity.in_(units.km/units.s)
                                        sys.stdout.flush()
                                    # Don't forget at the end of this we need to update the velocity
                                    # in multiples and the gravity code using the proper root mass
                                    # if that didn't happen automatically.

                        # Update pset shouldn't be needed now that we use update_roots_from_leaves and the above code.
                        #if (debug_multiples):
                        #    "Trying using update_pset."
                        #update_psetA_from_psetB(mult_grav.stars, stars)

                        update_roots_from_leaves(mult_grav, grav)
                        if (debug_multiples):
                            check_root_and_leaves(mult_grav, grav, stars)
                            print "Second kick after update_roots_from_leaves."
                            for st in stars:
                                # Cycle through the leaves and check if this particular star
                                # is in the leaves in multiples.
                                for root, tree in mult_grav.root_to_tree.iteritems():
                                    leaves = tree.get_leafs_subset()
                                    if st in leaves:
                                        if (debug_multiples):
                                            print "Star vel =", st.velocity.in_(units.km/units.s)
                                            print "Leaf vel =", leaves.velocity.in_(units.km/units.s)
                                            print "Root vel =", root.velocity.in_(units.km/units.s)
                                            sys.stdout.flush()
                        #mult_grav.channel_from_memory_to_code.copy() # Now copy from the memory of multiples to grav.
                    else:
                        stars_to_grav.copy_attributes(["vx", "vy", "vz"])
                        #stars_to_grav.copy()

                    print "Grav updated."
                    #print grav.particles.velocity

                #print "After second kick, max star velocity = ", np.sqrt((stars.vx**2.0 + stars.vy**2.0 + stars.vz**2.0).value_in(units.km**2.0/units.s**2.0)).max() | units.km / units.s

                #print "Checking hydro positions"
                #print old_pos - hydro.get_particle_position(tags_keys[:,0]).in_(units.m)

                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                        check_root_and_leaves(mult_grav, grav, stars)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)


                ### Check if output files need to be written.

                print "Checking for plot."

                hy_pltnum = hydro.IO_out('pltpart')
                hy_chknum = hydro.IO_out('chk')

                # If a checkpoint file was written,
                # then store the list of stars up for creation.
                if (hy_chknum != chknum):  # allow for possibility of rolling chk
                    write_rnd_and_mass_pickles(all_masses, output_dir, chknum)
                    chknum = hy_chknum  # must increment after write

                # wrote plt file
                if (hy_pltnum > pltnum):
                    if (write_psets):
                        stars.dt = dt
                        write_set_to_file(stars, pdir+'/stars{:04d}.amuse'.format(pltnum))
                        multstars = mult_grav.stars.copy_to_new_particles()
                        write_set_to_file(multstars, pdir+'/mult{:04d}.amuse'.format(pltnum))
                    pltnum = hy_pltnum  # must increment after write
                elif (hy_pltnum < pltnum):
                    raise Exception("Error: hy_pltnum={} < pltnum={}".format(hy_pltnum, pltnum))


                if (insane):
                    if (with_multiples):
                        check_sanity(hydro, mult_grav, stars, with_multiples=with_multiples)
                    else:
                        check_sanity(hydro, grav)
                if (type_insane): check_stellar_type(stars)

                print "Current simulation time:", t
                hydro_time = hydro.get_time()
                grav_time  = grav.get_time()
                print "Hydro time", hydro_time
                print "Grav time:", grav_time

                time_diff = (hydro_time - grav_time).value_in(units.s)

                print "Hydro time between initial hydro and grav bridge:", time_diff

                if (first_loop): first_loop=False


            ### Update timer information for the simulation.

            time_in_hydro = time_in_hydro + hydro_timer.secs
            time_in_grav  = time_in_grav + grav_timer.secs
            total_time    = total_time + loop_timer.secs

            print "Total time in Flash = %f s" %time_in_hydro
            print "Total time in N-body = %f s" %time_in_grav
            print "Total time in AMUSE = %f s" %(total_time - time_in_grav - time_in_hydro)


        #    if (tree_exists):
        #        if (hydro.get_number_of_particles() != len(grav.particles) != len(tree.particles)):
        #            print "Number of particles doesn't match in the three codes!!!"
        #            print "Aborting!"
        #            print "Num in grav = ", len(grav.particles)
        #            print "Num in hydro = ", hydro.get_number_of_particles()
        #            print "Num in tree = ", len(tree.particles)
        #            sys.exit()
        #    else:

        #        if (hydro.get_number_of_particles() != len(grav.particles)):
        #            print "Number of particles doesn't match in hydro and grav!"
        #            print "Aborting!"
        #            print "Num in grav = ", len(grav.particles)
        #            print "Num in hydro = ", hydro.get_number_of_particles()
        #            sys.exit()

            log.write(`i`+"\t \t"+"{0:.2e}".format(total_time)+
                      "\t \t"+"{0:.2e}".format(time_in_hydro)+
                      "\t \t" +"{0:.2e}".format(time_in_grav)
                      + "\t \t" + "{0:.2e}".format(total_time - time_in_grav - time_in_hydro))

            if (np.abs(time_diff) > 1e4):

                print "Time difference b/t hydro and n-body larger than 1e4, exiting."
                break

                hydro.timer_summary()

            if (profile):
                hydro.timer_summary()

except:
        raise
        log.write("Time spent in hydro: "+"{0:.2e}".format(time_in_hydro/total_time))
        log.write("Time spent in grav: "+"{0:.2e}".format(time_in_grav/total_time))
        log.write("Time spent in script: "+
                  "{0:.2e}".format((total_time - time_in_grav - time_in_hydro)/total_time))
        hydro.timer_summary()


log.write("Time spent in hydro: "+"{0:.2e}".format(time_in_hydro/total_time))
log.write("Time spent in grav: "+"{0:.2e}".format(time_in_grav/total_time))
log.write("Time spent in script: "+
          "{0:.2e}".format((total_time - time_in_grav - time_in_hydro)/total_time))

### Clean up the codes.
del log
#grav.stop()
#tree.stop()
#hydro.cleanup_code()

#if __name__ == "__main__":
    #main()
