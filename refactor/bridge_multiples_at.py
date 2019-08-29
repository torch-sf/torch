### Gravity bridge implementation for
### the Flash MHD code and a N-body solver.

### Joshua Wall
### Drexel University

print "At the top"

import glob
import os
import numpy as np
import pickle
import sys
import time

from scipy.integrate import *

from amuse.lab import *
from amuse.community.flash.interface import Flash
from amuse.community.kepler.interface import Kepler
from amuse.community.flash import josh_multiples as multiples
from amuse.community.smalln.interface import SmallN
from amuse.community.sse.interface import SSE

from amuse.io import write_set_to_file, read_set_from_file
#from amuse.rfi.channel import AsyncRequestsPool

from torch_ic import make_single_star_in_hydro, make_cluster_in_hydro, make_cluster
from torch_state import TorchIOState, add_particles_to_grav,  \
                        remove_particles_outside_bndbox, \
                        initialize_gravity_codes, initialize_multiples, \
                        init_smalln, new_smalln, stop_smalln, \
                        update_roots_from_leaves
from torch_se import do_stellar_evolution

from imf_sample import m_max_star
import ionizingflux as ion

from torch_state import add_particles_to_grav

print "After imports."

np.set_printoptions(precision=3)


#########################################################################
# Converters


# Converter for the N-body code.
convert = nbody.nbody_to_si(1.0 | units.parsec, 1000.0 | units.MSun)

#Converter for the hydro code.
convert2 = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
                   1.0| units.cm,1.0 | units.g, 1 | units.s)


####################################
# Radiation parameters
eion = 2.0 | units.eV  # ionizing photon energy in excess of 13.6 eV that's available to heat gas
sigh = 6.3e-18 | units.cm**2.0  # Lyman cross section
#z = 0.02  # Metallicity, don't think this is used anywhere
min_mass = 7.0 | units.MSun  # minimum stellar mass to turn on feedback

# Min and max star formation masses. Max from IGIMF restrictions figuring 0.5 SFE.
cluster_formation_eff = 1.0
cloud_mass = 1000.0 | units.MSun
min_sf_mass = 0.08 | units.MSun
max_sf_mass = m_max_star(cluster_formation_eff*cloud_mass.value_in(units.MSun)) | units.MSun

####################################
# Radiation constants
h = 6.6261e-27 # Planck's constant
c = 2.9979e10  # Speed of light
k = 1.3807e-16 # Boltzmann constant
sig0 = 6.304e-18 # Photoionization cross section at threshold for hydrogen

#l_ev = 1.2398e-4 # wavelength for 1 eV.
l_min = 1e-7 # Something really small.
E_ev = 1.60222497096e-12 # energy of 1 eV.

E_min = 13.6*E_ev  # 13.6 eV
l_max = h*c/E_min  # wavelength of 13.6 eV.

nu_min = E_min / h # freq of 13.6 eV.
nu_max = np.inf    # Max freq.

sigDust = 1e-21 | units.cm**2.0 # Cross section for dust from Draine 2011


# Number of processors for each code.
##########################

nproc = 42
num_grav_workers = 1
num_hy_workers = nproc - num_grav_workers - 4

##########################

# Debugging flags.
mult_debug_level = 1

# Record particle sets to file.
write_psets = False
pdir = "./psets"

# Start with a star or cluster (either from file or made here).
start_with_cluster     = False #True
read_cluster_from_file = True
nm_part=10
start_with_star = False #True #

# Runtime parameters.
with_bridge = True
no_sinks    = False
tree_exists = False
with_ph4    = True #False
with_multiples = True #False

use_radiation   = True
pe_heat         = True

with_se         = True
with_sn         = True

with_winds      = True
with_massloss   = True
massloss_method = 'puls'

min_pos_diff = 0.01*3.086e16 #| units.m #((0.01*3.085e16) | units.m)
min_vel_diff = 1e3 # | units.m/units.s


# N-body softening radius is the actual radius of a large massive star here.
eps = 15.0 | units.RSun

stars = Particles(0) # AMUSE Particles(...) set

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
    #se = SSE()
    se = SeBa()
    se.initialize_code()
sys.stdout.flush()

mult_grav = None
grav, stars_to_grav, grav_to_stars = initialize_gravity_codes(
         convert,
         stars,
         num_grav_workers = num_grav_workers,
         eps = eps,
         with_ph4 = with_ph4,
         with_multiples = with_multiples
)
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

if with_ph4:
    grav.parameters.begin_time = t
    grav.parameters.sync_time  = t
    grav.parameters.force_sync = 1
else:
    grav.parameters.begin_time = t
    grav.evolve_model(t)

### ======================================================================
### USER CONFIGURATION - INITIAL CONDITIONS

if (start_with_cluster):

    if (read_cluster_from_file):
        print "Reading initial cluster from file."
        initial_cluster = read_set_from_file('stars.hdf5', 'hdf5')
    else:
        initial_cluster = make_cluster(convert, nm_part, bndbox,
                                       fractal=True ,equal_mass=False,
                                       eq_mass=20.0 | units.MSun)
        print "Writing initial cluster to file."
        write_set_to_file(initial_cluster, 'starting_cluster.hdf5', 'hdf5')
        print "Done."

    print "Setting up initial cluster in hydro."
    make_cluster_in_hydro(hydro, initial_cluster)

if (start_with_star):

    make_single_star_in_hydro(
            hydro,
            x = 0.0 | units.cm,
            y = 0.0 | units.cm,
            z = 0.0 | units.cm,
            mass  = 30. | units.MSun,
            initMass  = 30. | units.MSun,
            age = 0.0 | units.Myr,
            vx = 0.0 | units.km/units.s,
            vy = 0.0 | units.km/units.s,
            vz = 0.0 | units.km/units.s,
    )

### ======================================================================
### USER CONFIGURATION

#np.random.seed(103180) # Set initial random seed for testing/debugging.

max_ref  = hydro.get_max_refinement()
SINK_RAD = 2.5*(2.0*bndbox)/(8.0*2.0**(float(max_ref)-1.0))

clear_particles_on_restart = False  # set True if nproc changed on restart.  # NOTE this is still bugged. 2/7/18 - JW
refresh_rand_seed_on_restart = False #True
dtinit = 1.0e10 | units.s
dtmax = 1.0e13 | units.s

### ======================================================================
### CODE INITIALIZATION, CONTINUED

state = TorchIOState(hydro, stars, mult_grav,
                     pdir=pdir,
                     refresh=refresh_rand_seed_on_restart)
state.prepare()  # handles I/O and restart loading

hydro.set_particle_pointers('mass')

# Check for existing particles in case of restart.
# FLASH particles array may be unsorted if this is a restart
hydro.particles_sort()  # (AT) apparently this can remove prtl?? - 2019 Jul 01
num_particles = hydro.get_number_of_particles()
if num_particles > 0:
    # Propagate existing (restart) stars to gravity codes
    add_particles_to_grav(state, mult_grav, stars, tree_exists)
    if with_multiples and mult_grav is None:
        mult_grav, mult_to_stars, stars_to_mult = initialize_multiples(
                stars, grav, convert,
                mult_debug_level=mult_debug_level,
                kep=kep,
                new_smalln=new_smalln
        )
        state.mult_grav = mult_grav


### ======================================================================
### MAIN EVOLUTION LOOP

first_step = True
max_hy_steps = hydro.get_max_num_steps()
curr_hy_step = hydro.get_current_step()

while (t < tmax and curr_hy_step < max_hy_steps):

    if first_step:
        first_step = False
        dt = dtinit
    else:
        hy_dt = hydro.get_timestep()
        dt = min(dtmax, 1.5*hy_dt, (tmax-t), 2.0*dt_old)

    dt_old = dt
    t_old  = t

    # Create star particles in hydro; update sink masses and sink
    # star-formation queues.
    made_stars = make_stars_from_sinks2(state, hydro, min_sf_mass, max_sf_mass)

    if made_stars:
        # Propagate new stars to gravity codes
        add_particles_to_grav(state, mult_grav, stars, tree_exists)
        if with_multiples and mult_grav is None:
            mult_grav, mult_to_stars, stars_to_mult = initialize_multiples(
                    stars, grav, convert,
                    mult_debug_level=mult_debug_level,
                    kep=kep,
                    new_smalln=new_smalln
            )
            state.mult_grav = mult_grav

    num_particles = hydro.get_number_of_particles()

    if num_particles == 0:

        t = t + dt
        print "I'm about to evolve hydro without evolving grav for :" , dt, "to evolve to t =", t

        ### --------------
        ### Evolve models.
        ### --------------
        hydro.evolve_model(t)

        ### -------------
        ### Output files.
        ### -------------
        torch_io.out(write_psets=False)

        continue

    else:

        print "Starting the gravity bridge."

        if with_se:
            do_stellar_evolution(stars, hydro)

        # when using multiples, we can't copy stars to grav because the
        # stars are all the particles (including leaves) but only the
        # roots are in grav.
        # Instead we should modify the values in multiples directly such that:
        # 1. Leaves have updated mass
        # 2. Roots have updated mass and velocity.
        # 3. We then need to update the radii of all the particles (to update the root radii).
        # 4. We then copy from multiples roots to grav code.
        if with_multiples:
            pass
            #mult_grav.channel_from_memory_to_code.copy()
            #mult_grav.channel_from_code_to_memory.copy_attribute("index_in_code", "id")
        else:
            stars_to_grav.copy()
        hydro.set_particle_mass(stars.tag, stars.mass)

        # Set the bridge timestep.
        t = t + dt
        print "I'm about to evolve hydro and grav for :" , dt, "to evolve to t =", t

        # key principle for understanding below code:
        # without multiples, stars <-> grav have a 1-to-1 prtl
        # correspondence.
        # with multiples, the tight binaries are separated out.
        #   grav tracks a center-of-mass (root) particle,
        #   multiples tracks the leaf particles
        # the stars data-structure only cares about leaf particles,
        # but changes to leaves also affect roots.
        # need to sync changes between:
        #     mult_grav, grav, stars, hydro
        # throughout the bridge in a careful manner.

        ### -----------
        ### First kick.
        ### -----------
        if (with_bridge):
            kick_number = 1
            if made_stars:  # grid changed, do kick#2 instead?...
                kick_number = 2

            ### Kick stuff in FLASH

            # gets gravity AND kicks the velocities of the gas.
            hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)
            # Calculate the gravity at the star locations.
            hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)

            ### Push veloc update from FLASH to stars, grav workers

            # pull hydro veloc -> stars
            stars.velocity = hydro.get_particle_velocity(stars.tag)
            # push stars veloc -> grav leaf
            stars_to_grav.copy_attributes(["vx", "vy", "vz"])
            # push stars veloc -> multiples leaf, grav root
            if (with_multiples):
                mult_grav.channel_from_code_to_memory.copy()
                for st in stars:
                    for root, tree in mult_grav.root_to_tree.iteritems():
                        leaves = tree.get_leafs_subset()
                        if st in leaves:
                            st.as_particle_in_set(leaves).velocity = st.velocity
                update_roots_from_leaves(mult_grav, grav)

        ### --------------
        ### Evolve models.
        ### --------------

        if (with_multiples):
            mult_grav.evolve_model(t)
        else:
            grav.evolve_model(t)

        hydro.evolve_model(t)

        hydro_time = hydro.get_time()
        grav_time  = grav.get_time()
        time_diff = (hydro_time - grav_time).value_in(units.s)
        if (time_diff > 1e4):
            print "Warning: Grav is behind Hydro by:", time_diff
            print "Trying again to evolve Grav to Hydro time."
            if (with_multiples):
                mult_grav.evolve_model(hydro_time)
            else:
                grav.evolve_model(hydro_time)

        curr_hy_step = hydro.get_current_step()

        # pull grav (leaf) pos,veloc -> stars pos,veloc
        grav_to_stars.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])
        # pull multiples (leaf) pos,veloc -> stars pos,veloc
        if (with_multiples):
            mult_grav.update_leaves_pos_vel()
            mult_grav.stars.copy_values_of_attributes_to(["x", "y", "z", "vx", "vy", "vz"], stars)

        # push stars pos,veloc -> hydro
        if (with_multiples):
            hydro.set_particle_position(stars.tag, stars.x, stars.y, stars.z)
            hydro.set_particle_velocity(stars.tag, stars.vx, stars.vy, stars.vz)
        else:
            hydro.set_particle_position(stars.tag, grav.particles.x, grav.particles.y, grav.particles.z)
            hydro.set_particle_velocity(stars.tag, grav.particles.vx, grav.particles.vy, grav.particles.vz)

        # removes from: hydro, stars, grav, and mult_grav
        remove_particles_outside_bndbox(hydro, stars, grav, mult_grav,
                with_multiples, bndbox, debug=True)

        # sort and also remove stars outside computational domain,
        # though remove_particles_outside_bndbox(...) should have us
        # covered
        hydro.particles_sort()

        ### ------------
        ### Second kick.
        ### ------------
        if (with_bridge):
            kick_number = 2

            ### Kick stuff in FLASH

            # gets gravity AND kicks the velocities of the gas.
            hydro.get_gravity_particles_on_gas(0.5*dt, kick_number)
            # Calculate the gravity at the star locations.
            hydro.get_gravity_gas_on_particles(0.5*dt, kick_number)

            ### Push veloc update from FLASH to stars, grav workers

            # pull hydro veloc -> stars
            stars.velocity = hydro.get_particle_velocity(stars.tag)
            # push stars veloc -> grav leaf
            stars_to_grav.copy_attributes(["vx", "vy", "vz"])
            # push stars veloc -> multiples leaf, grav root
            if (with_multiples):
                mult_grav.channel_from_code_to_memory.copy()
                for st in stars:
                    for root, tree in mult_grav.root_to_tree.iteritems():
                        leaves = tree.get_leafs_subset()
                        if st in leaves:
                            st.as_particle_in_set(leaves).velocity = st.velocity
                update_roots_from_leaves(mult_grav, grav)

        ### -------------
        ### Output files.
        ### -------------

        # is it OK to update stars.dt even if not doing write_psets?
        # why do we need to update stars.dt?
        # - 2019 Jul 01
        if write_psets:
            stars.dt = dt
        torch_io.out(write_psets=write_psets)

        print "Current simulation time:", t
        hydro_time = hydro.get_time()
        grav_time  = grav.get_time()
        print "Hydro time", hydro_time
        print "Grav time:", grav_time
        time_diff = (hydro_time - grav_time).value_in(units.s)
        print "Hydro time between initial hydro and grav bridge:", time_diff
        if (np.abs(time_diff) > 1e4):
            print "Time difference b/t hydro and n-body larger than 1e4, exiting."
            break

### Clean up the codes.

hydro.timer_summary()

#grav.stop()
#if (with_multiples):
#    kep.stop()
#    stop_smalln()
#tree.stop()
#hydro.cleanup_code()
