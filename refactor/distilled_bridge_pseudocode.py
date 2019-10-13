#!/usr/bin/env python

### Gravity bridge implementation for
### the Flash MHD code and a N-body solver.

### Joshua Wall
### Drexel University

### IDEA: strip out debugging code and lots of special case handling
### like restart code, print statements, etc...
### Remove some of method abstraction to see the exact AMUSE calls.
### Remove checks for successful file write and such

print "At the top"

import datetime
import glob
import numpy as np
import os
import pickle
from scipy.integrate import *
import sys
import time

from amuse.lab import *
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.community.flash import josh_multiples as multiples
#from amuse.community.flash import steve_multiples as multiples
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
#from amuse.rfi.channel import AsyncRequestsPool
from amuse.community.smalln.interface import SmallN
from amuse.community.sse.interface import SSE

#import main_sequence as ms
import ionizingflux as ion

print "After imports."

np.set_printoptions(precision=3)

# Global variables

all_masses = {}
first_call_for_stars = False
old_sink_tags = []
mult_grav      = None

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

    rstatefile =  output_dir+'/rnd_state'+chknum+'.pickle'
    massesfile = output_dir+'/all_masses'+chknum+'.pickle'

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

    position = hydro.get_particle_position(newtags)
    velocity = hydro.get_particle_velocity(newtags)
    mass     = hydro.get_particle_mass(newtags)
    initMass = hydro.get_particle_oldmass(newtags)
    age      = hydro.get_time() - hydro.get_particle_creation_time(newtags)

    add_star = Particles(num_new_parts)             # Make new particles for grav code.
    keys     = add_star.key                         # Get the new keys for tags_keys.

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
        add_star[kk].radius = 100 | units.AU # initial collision radius.
        if (initMass[kk].value_in(units.MSun) <= 0.001):
            add_star[kk].initial_mass = mass[kk] # record the initial mass of the star for SE/SN uses.
            hydro.set_particle_oldmass(newtags[kk], mass[kk])
        else:
            add_star[kk].initial_mass = initMass[kk] # record the initial mass of the star for SE/SN uses.
        stars_current_id_num += 1
        add_star[kk].id = stars_current_id_num

        if (use_radiation):
            add_star[kk].nion = 0.0 # ionizing flux
            add_star[kk].eion = 0.0 #eion #0.0 # ionizing energy *OVER* 13.6 eV
            add_star[kk].sigh = 0.0 #sigh #0.0 # ionizing cross section.

    if (np.shape(tags_keys)[0] == 0):
        tags_keys = np.zeros((num_new_parts, 2))
        tags_keys[:,0] = newtags
        tags_keys[:,1] = keys
    else:
        jj = 0
        for ii in range(num_new_parts):
            tags_keys = np.append(tags_keys, [[newtags[ii], keys[ii]]], axis=0)

    stars.add_particles(add_star)

    # Sort tags_keys so that we always have the tags in order.
    tags_keys = tags_keys[tags_keys[:,0].argsort()]
    # Sort the stars to be in the same order as the tags.
    stars = stars.sorted_by_attribute('tag')

    if (tree_exists):
        tree.particles = stars.sorted_by_attribute('tag')

    grav.particles.add_particles(add_star)
    if (mult_grav is not None and not add_parts_restart):
        mult_grav._inmemory_particles.add_particles(add_star)
        mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
    if (with_ph4 and add_parts_restart):
        grav.commit_particles()

    if (add_parts_restart):
        hydro.set_starting_local_tag_numbers()

    if (add_parts_restart and clear_particles_on_restart):
        reinitialize_all_particles_from_stars(stars, hydro, grav, tags_keys)

    # Clear any stored new tags in FLASH now that we've successfully added the particles
    # to the gravity code.
    hydro.clear_new_tags()

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

    tags_keys.sort(axis=0)

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

    avg_mass = quad(mkroupa, M_min, M_max, args=(norm))[0]

    for i in range(num_bins):
        # m_i = int_bin_i_low^bin_i_high(m*f*dm) / int(f*dm)
        mass_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    quad(kroupa, binsL[i], binsL[i+1], args=(norm))[0] )
        # f_i = int_bin_i_low^bin_i_high(m*f*dm) / int_M_low^M_high(m*f*dm)
        frac_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0] /
                    avg_mass )

    mass_per_bin = np.array(mass_per_bin)
    frac_per_bin = np.array(frac_per_bin)

    lam = sink_mass*frac_per_bin/mass_per_bin

    n_stars = np.random.poisson(lam=lam)

    print "Mass from N stars ~", np.sum(n_stars*mass_per_bin)

    return n_stars, binsL, lam, norm

def collect_small_stars_mass(all_samp_masses):

    # Here we move all the stars smaller than 1 MSun into particles
    # that are at least 1 MSun. To do this we do a bit of fancy
    # footwork with the arrays.

    small_masses = all_samp_masses[np.where(all_samp_masses < 1.0)] # Smaller than 1.0 MSun.
    all_samp_masses = all_samp_masses[np.where(all_samp_masses >= 1.0)]  # Everyone else.

    b = 0
    # If there are any left smaller than 1.0 MSun, sum with others
    # that are smaller than 1.0 MSun until there are none left.
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
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
        #rem_index = np.where(np.greater_equal(np.abs(stars.x.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
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
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)


    if (num_particles > 0):
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.y.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
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
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)
        #if (with_se):
            #se.particles.remove_particles(rem_part)

    if (num_particles > 0):
        #rem_index = np.where(np.abs(grav.particles.z.value_in(units.cm)) > bndbox.value_in(units.cm))[0]
        rem_index = np.where(np.greater_equal(np.abs(grav.particles.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]
        rem_key   = grav.particles.key[np.where(np.greater_equal(np.abs(grav.particles.z.value_in(units.cm)), bndbox.value_in(units.cm)))[0]]
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
        for rt in rem_tag:
            tags_keys = tags_keys[~(tags_keys[:,0]==rt),:]

        num_particles = len(stars)

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
                np = words[2]
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

    return restart, chknum, pltnum, basename

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

    multiples_code = None
    return stars, multiples_code, grav, stars_to_grav, grav_to_stars

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

# Check to see if we wrote a plot or checkpoint file.
def poll_for_new_output_files(file_name, file_modified_time):
    import os
    # There's two scenarios.
    # 1) We wrote a new file and
    # 2) We overwrote an existing file.
    new_file = False
    # If file didn't exist before, file_modified_time < 0
    if (file_modified_time < 0.0):
        if (os.path.isfile(file_name)):
            new_file = True
            file_modified_time = os.stat(file_name).st_mtime
    else: # we pass the last known modified time.
        if (os.path.isfile(file_name)):
            curr_file_modified_time = os.stat(file_name).st_mtime
            if (curr_file_modified_time > file_modified_time):
                new_file = True
                file_modified_time = curr_file_modified_time
        else:
            file_modified_time = -1.0 # There isn't a file with this name.
    return new_file, file_modified_time

# Write out pickles with the current random number state and the
# dictionary with all the future stars in it.
def write_rnd_and_mass_pickles(all_masses, output_dir, new_chpt_file):
    print "Writing out all masses to all_masses.pickle"
    with open(output_dir+'/all_masses'+new_chpt_file+'.pickle', 'wb') as f:
        pickle.dump(all_masses, f)
    old_chpt_file = new_chpt_file
    print "Writing out random state to rnd_state.pickle"
    rnd_state = np.random.get_state()
    with open(output_dir+'/rnd_state'+new_chpt_file+'.pickle', 'wb') as f:
        pickle.dump(rnd_state, f)
    return


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
    num_hy_workers= int(get_np_from_run_script()) - num_grav_workers - 4
except:
    print "WARNING: Setting num_hy_workers from script failed. Defaulting to {} procs.".format(default_hy_procs)
    num_hy_workers = default_hy_procs

print "Number of hydro procs = ", num_hy_workers
print "Number of Nbody procs = ", num_grav_workers

##########################

# Debugging flags.
##########################
mult_debug_level = 1

# Record particle sets to file.
write_psets = False
pdir = "./psets"

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

print "Starting smallN."
init_smalln(convert)

print "Starting Kepler."
kep = Kepler(unit_converter=convert)
kep.initialize_code()

print "Starting stellar evolution code."
#se = SSE()
se = SeBa()
se.initialize_code()

# (AT) stars is an AMUSE particle set
#      mult_grav is AMUSE's multiples module
#      grav is ph4 or Hermite
#      stars_to_grav, grav_to_stars are "channels" between gravity code and prtl set
print "Starting gravity code."
stars, mult_grav, grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
         convert, stars = stars, start_time = hydro_time,
         num_grav_workers = num_grav_workers,
         eps = eps,
         with_ph4 = with_ph4,
         with_multiples = with_multiples)

print "Starting hydro code."
hydro = Flash(unit_converter = convert2, number_of_workers=num_hy_workers, redirection='none')
hydro.initialize_code()

### Get the simulation end time, current AMUSE step and
### current simulation time from Flash (in case of restart).

tmax = hydro.get_end_time()
t = hydro.get_time()
hydro_time = t
bndbox = hydro.get_runtime_parameter('xmax') | units.cm

if (with_ph4):
    grav.parameters.begin_time=hydro_time
    grav.parameters.sync_time=hydro_time
    grav.parameters.force_sync=1
else:
    grav.parameters.begin_time=hydro_time
    grav.evolve_model(hydro_time)

###############################################
### Initialize either clusters or stars.
###############################################

if (start_with_cluster):

    # Make fractal cluster as AMUSE prtl set within hydro bounding box
    stars_out = True
    while (stars_out):
        init_cluster = new_kroupa_mass_distribution(nm_part, mass_max = (100.0 |units.MSun))
        init_cluster = new_fractal_cluster_model(masses=init_cluster, convert_nbody=convert, do_scale=False, virial_ratio=5.0)
        init_cluster.mass = new_salpeter_mass_distribution(nm_part, mass_min = (0.1 | units.MSun), mass_max = (100.0 |units.MSun))  # ??
        remove_stars = cluster.select(lambda r: bndbox < max(abs(r)), ["position"])
        stars_out = len(remove_stars) > 0

    # Push cluster to hydro
    tag = hydro.add_particles(init_cluster.x + 0.0 | units.cm,
                              init_cluster.y + 0.0 | units.cm,
                              init_cluster.z + 0.0 | units.cm,)
    hydro.set_particle_velocity(tag, init_cluster.vx, init_cluster.vy, init_cluster.vz)
    hydro.set_particle_mass(tag, init_cluster.mass)
    del tag


if (start_with_star):

    x  = 0.0 | units.cm #(0.5*smallest_dx)
    y  = 0.0 | units.cm #(0.5*smallest_dx)
    z  = 0.0 | units.cm #(0.5*smallest_dx)
    vx = 0.0 | units.km/units.s
    vy = 0.0 | units.km/units.s
    vz = 0.0 | units.km/units.s
    m  = 30. | units.MSun

    tag = hydro.add_particles(0.0|units.cm,
                              0.0|units.cm,
                              0.0|units.cm,)
    hydro.set_particle_velocity(tag, 0.0|units.km/units.s,
                                     0.0|units.km/units.s,
                                     0.0|units.km/units.s,)
    hydro.set_particle_mass(tag, 30.|units.MSun)
    hydro.set_particle_oldmass(tag, 30.|units.MSun)
    hydro.set_particle_creation_time(tag, hydro.get_time() - 0.0|units.Myr)

# for my code study, assume no prtl to start; I threw out some special case
# handling of restart
assert not start_with_cluster
assert not start_with_star

###############################################
### Initialize some parameters.
###############################################

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
new_chk_file = '-1'
old_chk_file = '-1'
new_plt_file = '-1'
old_plt_file = '-1'

max_ref  = hydro.get_max_refinement()
sink_rad = 2.5*(2.0*bndbox)/(8.0*2.0**(float(max_ref)-1.0))
smallest_dx = (2.0*bndbox)/(8.0*2.0**(float(max_ref)-1.0))

time.sleep(5)

# This should be automated!
restart,chknum,pltnum,basename = get_restart_and_chk_num()
chknum = str(chknum).zfill(4)
pltnum = str(pltnum).zfill(4)
refresh_rand_seed_on_restart = False #True

if (refresh_rand_seed_on_restart):
    print "WARNING WARNING WARNING! Resetting the rand seed!!!!!!!!!!"

output_dir = hydro.get_output_dir()
chk_pre = output_dir+'/'+basename+'hdf5_chk_'
plt_pre = output_dir+'/'+basename+'hdf5_plt_cnt_'
current_chk_file = chk_pre+chknum
new_chk_file = chk_pre+str(int(chknum)+1).zfill(4)
current_plt_file = plt_pre+pltnum
new_plt_file = plt_pre+str(int(pltnum)+1).zfill(4)
chk_mtime = -1.0
plt_mtime = -1.0

print "Writing output files to: ", output_dir

print "Is this a restart?", restart
print "We are assuming the checkpointfilenumber =", chknum
load_rnd_state_files(restart, chknum, refresh_rand_seed_on_restart)

wrote_chk_file = False
wrote_plt_file = False
# Get file modified times for current chk and plt files.
wrote_chk_file, chk_mtime = poll_for_new_output_files(current_chk_file, chk_mtime)
wrote_plt_file, plt_mtime = poll_for_new_output_files(current_plt_file, plt_mtime)

print "Curr chk_file=", current_chk_file
print "Curr chknum=", chknum
print "Curr chk_mtime=", chk_mtime
print "New chk_file=", new_chk_file

print "Curr plt_file=", current_plt_file
print "Curr pltnum=", pltnum
print "Curr plt_mtime=", plt_mtime
print "New plt_file=", new_plt_file

if (wrote_chk_file and not restart):
    write_rnd_and_mass_pickles(all_masses, output_dir, chknum)

# This is used if you change the number of processors on restart!
# Or if you need to clear out tracers or some other thing.
# NOTE this is still bugged. 2/7/18 - JW
clear_particles_on_restart = False #True

### Set dt based on the current cloud free fall time.
tff = (np.sqrt(3.0*np.pi/(32.0*units.constants.G*cloud_dens))).as_quantity_in(units.s)

print "Hydro dynamical time in secs: %3.3e s" % (tff).value_in(units.s)
print (tff/100.0).as_quantity_in(units.yr)
print "Currently setting initial dt artbitrarily high so that dt = hydro dt."
dtinit = 1.0e10 | units.s
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

hydro.set_particle_pointers('mass')


while ((t < tmax) and (curr_hy_step < max_hy_steps)):
    i = i + 1

    ### Check for proper bridge timestep based on hydro timestep and crossing time.
    ### Have to write a proper routine to get timestep from hydro.
    hy_dt = hydro.get_timestep()

    dt = min(dtmax, 1.5*hy_dt, (tmax-t), 2.0*dt_old)
    if (first_step):
        dt = dtinit
        first_step = False

    dt_old = dt
    t_old  = t

    # (AT) if FLASH created a new sink prtl that's "massive enough",
    # allocate a new massive prtl.  Josh uses hydro.set_particle_pointers
    # to toggle between sink and massive particles in FLASH
    made_stars, made_massive_star = make_stars_from_sinks2(hydro, min_sf_mass, max_sf_mass)

    # (AT) grav can be ph4 or Hermite
    if (made_stars):
        gridChanged = True
        if (first_particle):
            tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
            num_particles = hydro.get_number_of_particles()

    ### Wait for the first star to form.
    if (not first_particle):

        check_particles = hydro.get_number_of_particles()

        if check_particles != 0:

            first_particle = True
            hydro_time = hydro.get_time()
            grav_time  = grav.get_time()

            tags_keys, stars = add_particles_to_grav(tags_keys, stars, tree_exists)
            if (with_multiples):
                # (AT) mult_grav is an instance of multiples.Multiples(...)
                # which might be Kepler or SmallN, I'm not sure.
                # See AMUSE paper: https://ui.adsabs.harvard.edu/#abs/2013A&A...557A..84P/abstract
                _ = initialize_multiples(stars, grav, convert,
                                         mult_debug_level=mult_debug_level,
                                         kep=kep, new_smalln=new_smalln)
                mult_grav, mult_to_stars, stars_to_mult = _
                print "Gravity softening radius**2 =", grav.parameters.epsilon_squared

            print "Num particles in grav:", len(grav.particles)
            if (tree_exists):
                print "Num particles in tree:", len(tree.particles)

            num_particles = check_particles
            hydro_time = hydro.get_time()
            grav_time  = grav.get_time()

        else:

            t = t + dt

            if (do_sn_once and not first_loop):
                dt = hydro.energy_injection(1e51|units.erg, -1.0,
                                            (5.*1.989e33)|units.g,
                                            0.0|units.cm,
                                            0.0|units.cm,
                                            0.0|units.cm)
                t = t_old + dt
                do_sn_once = False
                hydro.write_chpt()

            print "I'm about to evolve hydro without evolving grav for :" , dt, "to evolve to t =", t
            hydro.evolve_model(t) # Can't use i*dt if dt can get smaller each timestep.

            # Note that if you are trying to output on nstep, Flash
            # assumes this is checked during the Driver_evolveFlash loop
            # and this won't output files properly. It will work
            # properly if you use output time parameters, however.
            #hydro.IO_out('all')
            hydro.IO_out('pltpart')
            hydro.IO_out('chk')

            chknum = str(int(chknum)+1).zfill(4)
            pltnum = str(int(pltnum)+1).zfill(4)
            write_rnd_and_mass_pickles(all_masses, output_dir, chknum)

            hydro_time = hydro.get_time()
            grav_time  = grav.get_time()

            first_loop = False

            continue  # main while loop

    print "Starting the gravity bridge."

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

        star_mass = hydro.get_particle_mass(tags_keys[:,0])

        # NOTE: You must add in the dt here, otherwise newly formed stars will try to
        # evolve for 0.0 seconds and things will get ugly.
        star_age  = hydro.get_time() + dt - hydro.get_particle_creation_time(tags_keys[:,0])

        for part in range(num_particles):

            if (star_age[part].value_in(units.yr) < 1.0):
                star_age[part] = 1.0 | units.yr

            if (13 <=  stars.stellar_type[part].value_in(units.stellar_type) <= 15):
                print "Skipping this star that already went SN, current stellar type =", stars.stellar_type[part]
                continue

            # Use SE code on star's initial mass, not current mass.
            _ = se.evolve_star(stars.initial_mass[part], star_age[part], 0.02)
            st_time, st_mass, star_radius, star_lum, star_temp, st_evol_time, st_type = _
            star_type[part] = st_type

            if (with_massloss and (massloss_method == 'seba' or st_mass.value_in(units.MSun) < min_mass.value_in(units.MSun))):
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
                dm_dt[part] = 10**(-24.06 + 2.45 * np.log10(star_lum.value_in(units.LSun))
                            -1.10*np.log10(stars[part].mass.value_in(units.MSun))
                            + 1.31*np.log10(star_temp.value_in(units.K))) | units.MSun/units.yr
                vterm[part] = 10**(1.23 - 0.30 * np.log10(star_lum.value_in(units.LSun))
                            + 0.55*np.log10(stars[part].mass.value_in(units.MSun))
                            + 0.64*np.log10(star_temp.value_in(units.K))) | units.km/units.s
            elif (with_massloss and (massloss_method == 'puls' or st_mass.value_in(units.MSun) >= min_mass.value_in(units.MSun))):
                # Kudritzki and Puls winds, see Kudritzki & Puls 2000, Markova & Puls 2004, 2008 and Vink 2000
                star_wind   = stellar_wind(star_temp, stars[part].mass, star_lum, star_radius)
                dm_dt[part] = star_wind.dm_dt.as_quantity_in(units.g / units.s)
                vterm[part] = star_wind.vterm.as_quantity_in(units.cm / units.s)
            else:
                dm_dt[part] = 0.0 | units.g / units.s
                vterm[part] = 0.0 | units.cm / units.s

            # If with energy injection, check to see if anything went supernova. If so, inject 10^51 ergs of
            # energy into the grid

            if (with_sn):

                # (AT) I don't understand this logic yet
                # Here initial SN does not depend on particles; is it
                # supposed to preclude multiple supernovae in same dt?
                if (do_sn_once and not first_loop):
                    dt = hydro.energy_injection(1e51|units.erg, -1.0,
                                                (5.*1.989e33)|units.g,
                                                0.0|units.cm,
                                                0.0|units.cm,
                                                0.0|units.cm)
                    do_sn_once = False

                elif (13 <=  st_type.value_in(units.stellar_type) <= 15):

                    # injected mass = current mass minus the remnant mass.
                    inj_mass = (star_mass[part] - st_mass).in_(units.g)
                    if (inj_mass.value_in(units.MSun) > 10.0):  # should never be more than 10 Msun
                        inj_mass = 10.0 | units.MSun

                    dt_sn = hydro.energy_injection(1e51|units.erg, -1.0,
                                                   inj_mass,
                                                   stars.x[part],
                                                   stars.y[part],
                                                   stars.z[part])

                    dt =  min(dt.value_in(units.s), dt_sn.value_in(units.s)) | units.s
                    del dt_sn

                    # Set proper mass for remnant (SeBa does fine for this).
                    star_mass[part] = st_mass
                    # Set proper remnant stellar type so that we don't get any feedback from remnants.
                    stars.stellar_type[part] = st_type
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
                    if (with_winds):
                        hydro.set_particle_wind_mass(tags_keys[part,0], dm_dt[part])
                        hydro.set_particle_wind_vel(tags_keys[part,0], vterm[part])
                    # Do nothing else with this star, just jump to the next one.
                    continue

            if ((star_mass[part].in_(units.MSun) >= min_mass) and not (13 <=  st_type.value_in(units.stellar_type) <= 15)):

                print "Found massive star, star mass =", star_mass[part].in_(units.MSun)

                part_inds.append(part)

                if (use_radiation):

                    flux = ion.ionizing_photon_flux(st_mass, star_radius, star_temp)

                    # Calculate the average ionizing photon energy based on the blackbody curve.

                    # First integrate the power from the BB curve at this stars temp.
                    # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                    [power, err] = quad(lum_wl_cs, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                    # Now integrate to find the number of photons.
                    [per_ph, err] = quad(lum_wl_cs_per_ph, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                    avg_E = power/per_ph / E_ev
                    # Calculate the average frequency of an ionizing photon for this star
                    avg_nu = avg_E*E_ev/h
                    # Cross section calculation
                    # Make sure you convert energy back to ergs if you
                    # use it to calculate the frequency!
                    sig = sig0*(avg_nu/nu_min)**(-3.0)
                    eion[part] = (avg_E | units.eV) - (13.6 |units.eV) #2.0 | units.eV #6.0 | units.eV
                    sigh[part] = sig | units.cm**2.0 #6.3e-18 | units.cm**2.0
                    nphot[part] = (flux*4*np.pi*star_radius**2.0).as_quantity_in(units.s**-1.0) #5e48 | units.s**(-1.0)

                    if (pe_heat):

                        # First integrate the power from the BB curve at this stars temp.
                        # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                        l_min_dust = h*c / E_min # wavelength at 13.6 eV
                        l_max_dust = h*c / (5.6*E_ev) # wavelength at 5.6 eV
                        [power, err] = quad(lum_wl, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))
                        # Now integrate to find the number of photons.
                        [per_ph, err] = quad(lum_wl_per_ph, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))
                        avg_E = power/per_ph / E_ev
                        # Calculate the average frequency of an ionizing photon for this star
                        #avg_nu = avg_E*E_ev/h
                        # Cross section calculation
                        # We assume constant cross section for dust per hydrogen atom.
                        # Value = tau / N_H where tau = gamma * Av (Draine and Bertoli 96)
                        # Av = N_H,tot / (1.87e21 cm^2) (Bohlin et al 78)
                        # gamma = 2.5 (Bergin et al 2004)
                        sigpe[part] = sigDust
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
                        # If yes, then set the mass of the star to the updated to the new mass from
                        # the SE code. This should be possible by using the as_particle_in_set(leaves)
                        # method on the stars particle. We'll check with print statements.
                        stars[part].as_particle_in_set(leaves).mass = star_mass[part]

                        # Don't forget at the end of this we need to update the masses
                        # in multiples and the gravity code using the proper root mass
                        # if that didn't happen automatically.

                        # Update the new root com and root vel.
                        update_roots_from_leaves(mult_grav, grav) # If works, move outside se loop!

        # If the leaf mass changed, so does the com particle's properties.
        update_roots_from_leaves(mult_grav, grav) # If works, move outside se loop!

        # Are there any massive stars?

        if not (part_inds==[]):

            if (use_radiation):
                hydro.set_particle_nion(tags_keys[part_inds,0], nphot[part_inds])
                hydro.set_particle_eion(tags_keys[part_inds,0], eion[part_inds].as_quantity_in(units.erg))
                hydro.set_particle_sigh(tags_keys[part_inds,0], sigh[part_inds])

            if (pe_heat):
                hydro.set_particle_npep(tags_keys[part_inds,0], npe[part_inds])
                # Set average energy of PE photon
                hydro.set_particle_epep(tags_keys[part_inds,0], epe[part_inds].as_quantity_in(units.erg))
                # Set cross section of dust to PE photons.
                hydro.set_particle_sigd(tags_keys[part_inds,0], sigpe[part_inds])

            if (with_winds):
                hydro.set_particle_wind_mass(tags_keys[part_inds,0], dm_dt[part_inds])
                hydro.set_particle_wind_vel(tags_keys[part_inds,0], vterm[part_inds])

    # Remove any mass loss due to winds and update to this
    # mass. Note this assumes steps are relatively small
    # in the mass loss rate of stars, so that gravitationally
    # we can use the mass after all the wind mass loss
    # has occcured. Otherwise we'd have to average
    # mass loss and keep up with old and new masses and
    # it just gets ugly.
    stars.mass = star_mass # - dm_dt*dt
    stars.age  = stars.age + dt
    stars.stellar_type = star_type

    # Note when using the multiples module we can't copy stars to grav
    # because the stars are all the particles (including leaves) and only
    # the roots are in grav.
    # Instead we should modify the values in multiples directly such that:
    # 1. Leaves have updated mass
    # 2. Roots have updated mass and velocity.
    # 3. We then need to update the radii of all the particles (to update the root radii).
    # 4. We then copy from multiples roots to grav code.
    if (not with_multiples):
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
            gridChanged = False

        hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # this gets gravity AND kicks gas velocities
        print "Grid kicked."
        hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # calculate gravity at star locations
        stars.velocity = hydro.get_particle_velocity(tags_keys[:,0])  # update gravity code with kicked velocity
        print "Stars kicked."

        stars_to_grav.copy_attributes(["vx", "vy", "vz"])

        if (with_multiples):
            # First copy updated velocities from stars to grav for
            # all particles in grav that are also in stars
            # (i.e. things that are leaves in grav and stars).
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
                        # If yes, then set the mass of the star to the updated to the new mass from
                        # the SE code. This should be possible by using the as_particle_in_set(leaves)
                        # method on the stars particle. We'll check with print statements.
                        st.as_particle_in_set(leaves).velocity = st.velocity
                        # Don't forget at the end of this we need to update the velocity
                        # in multiples and the gravity code using the proper root mass
                        # if that didn't happen automatically.

            # Update pset shouldn't be needed now that we use update_roots_from_leaves and the above code.
            update_roots_from_leaves(mult_grav, grav)

        print "Grav updated."

    print "Evolving models."

    if (with_multiples):
        mult_grav.evolve_model(t)
    else:
        grav.evolve_model(t)

    hydro.evolve_model(t)

    hydro_time = hydro.get_time()
    grav_time  = grav.get_time()
    time_diff = (hydro_time - grav_time).value_in(units.s)
    if (time_diff > 1e4):
        print "Warning: Grav behind Hydro by:", time_diff, "retry Grav evolve"
        if (with_multiples):
            mult_grav.evolve_model(hydro_time)
        else:
            grav.evolve_model(hydro_time)

    curr_hy_step = hydro.get_current_step()

    grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])
    if (with_multiples):
        mult_grav.update_leaves_pos_vel()
        mult_grav.stars.copy_values_of_attributes_to(["x", "y", "z", "vx", "vy", "vz"], stars)

### Update the star locations in FLASH so that the grid is kicked at the proper place by
### the stars. This is required because the grid kick all happens internally
### inside Flash since that is much faster. (Do this BEFORE hydro.kick_grid!)

### Also note, we need hydro to see the latest particle positions so that it can remove anything
### that has left the computational domain when we call hydro.sort_particles below.

    if (with_multiples):
        hydro.set_particle_position(tags_keys[:,0], stars.x, stars.y, stars.z)
        hydro.set_particle_velocity(tags_keys[:,0], stars.vx, stars.vy, stars.vz)
    else:
        hydro.set_particle_position(tags_keys[:,0], grav.particles.x, grav.particles.y, grav.particles.z)
        hydro.set_particle_velocity(tags_keys[:,0], grav.particles.vx, grav.particles.vy, grav.particles.vz)

    stars_removed, num_particles  = remove_particles_outside_bndbox(hydro, stars, grav, mult_grav, with_multiples, tags_keys, bndbox, debug=True)

    # Sync stars to the particles in grav.
    if (with_multiples):
        grav.particles.synchronize_to(mult_grav._inmemory_particles)
    else:
        grav.particles.synchronize_to(stars)

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

    if (stars_removed and num_particles == 0):
        print "No particles left! Going back to hydro only!"
        first_particle = False
        continue

    if (with_bridge):

        print "Second kick."
        kick_number = 2

        hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)  # this gets gravity AND kicks gas velocities
        print "Grid kicked."
        hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)  # calculate gravity at star locations
        stars.velocity = hydro.get_particle_velocity(tags_keys[:,0])  # update gravity code with kicked velocity
        print "Stars kicked."

        stars_to_grav.copy_attributes(["vx", "vy", "vz"])

        if (with_multiples):
            # First copy updated velocities from stars to grav for
            # all particles in grav that are also in stars.
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
                        # If yes, then set the mass of the star to the updated to the new mass from
                        # the SE code. This should be possible by using the as_particle_in_set(leaves)
                        # method on the stars particle. We'll check with print statements.
                        st.as_particle_in_set(leaves).velocity = st.velocity
                        # Don't forget at the end of this we need to update the velocity
                        # in multiples and the gravity code using the proper root mass
                        # if that didn't happen automatically.

            # Update pset shouldn't be needed now that we use update_roots_from_leaves and the above code.

            update_roots_from_leaves(mult_grav, grav)

        print "Grav updated."

    # (AT) write output and checkpoint

    hydro.IO_out('pltpart')
    hydro.IO_out('chk')

    # (AT) for learning purposes: don't check if files are written,
    # just assume it worked.  Maybe we could get success/fail from
    # FLASH's own IO unit?
    chknum = str(int(chknum)+1).zfill(4)
    pltnum = str(int(pltnum)+1).zfill(4)
    write_rnd_and_mass_pickles(all_masses, output_dir, chknum)

    # store the list of stars up for creation.
    if (write_psets):
        stars.dt = dt
        write_set_to_file(stars, pdir+'/stars'+pltnum+'.amuse')
        multstars = mult_grav.stars.copy_to_new_particles()
        write_set_to_file(multstars, pdir+'/mult'+pltnum+'.amuse')

    hydro_time = hydro.get_time()
    grav_time  = grav.get_time()

    first_loop = False

### Clean up the codes.
#grav.stop()
#tree.stop()
#hydro.cleanup_code()
