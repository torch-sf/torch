#!/usr/bin/env python



import numpy as np
from os import path
import pickle

from amuse.datamodel import Particles
from amuse.io import write_set_to_file
from amuse.units import units

from torch_param import FlashPar
from torch_stdout import tprint

class TorchState(object):
    """
    Container for AMUSE framework global objects, to help
    (1) hold things, (2) perform I/O for all torch workers.
    """

    def __init__(self, hydro, grav, mult, se):

        self.hydro = hydro
        self.grav  = grav
        self.mult  = mult
        self.se    = se     #CCC 26/04/2024 to match above
        
        # "Global" AMUSE-level data structures
        self.all_masses = {}
        self.loop = {}
        self.stars = Particles(0)
        self.stars_next_id = 0  # to supply ID attribute for ph4

        # Adding these to permit primordial binaries -CCC, May 1, 2020
        self.system_masses = {}
        self.all_positions = {}
        self.all_velocities = {}

        self.stars_to_grav = self.stars.new_channel_to(grav.particles)
        self.grav_to_stars = grav.particles.new_channel_to(self.stars)

        # Stellar evolution to stars, CCC 26/04/2024
        self.stars_to_se = self.stars.new_channel_to(se.particles)
        self.se_to_stars = se.particles.new_channel_to(self.stars)

        # TODO enhancement - read from FLASH's own RuntimeParameter interface,
        # instead of duplicating the flash.par file parsing and default case
        # behavior.  Needs new interface code, see hydro.get_runtime_parameter(...)
        # which only works for doubles presently. -AT, 2019jul01, 2019oct14
        #self.restart = False  # flash defaults, in case not specified in flash.par
        #self.chknum = 0
        #self.pltnum = 0

        p = FlashPar("flash.par")
        assert 'checkpointFileNumber' in p  # guard against, e.g., wrong-case typos
        assert 'plotFileNumber' in p
        assert 'restart' in p
        self.chknum = int(p['checkpointFileNumber'])
        self.pltnum = int(p['plotFileNumber'])
        self.restart = p['restart']

        self.output_dir = hydro.get_output_dir()

    def initial_io(self, overwrite, refresh=False):
        """Load restart files or write starting Torch state
        Kwargs:
            refresh (False): for restarts,
            ... use new random number generator state.
            ... force code to draw new list of stars from IMF for all sinks.
            ... if run is not restart, no effect.
        """

        # flash numbering for file I/O:
        #   restart:
        #     load {chknum}
        #     ...
        #     write {chknum+1}
        #     write {pltnum}
        #     ...
        #   not restart:
        #     write {chknum}
        #     write {pltnum}
        #     ...
        #     write {chknum+1}
        #     write {pltnum+1}
        #     ...
        # mirror flash logic in our writing/loading of RNG state.
        # and tracking of output mtimes

        if self.restart:

            loopfile = path.join(self.output_dir,
                'torch_loop{:04d}.pickle'.format(self.chknum))

            with open(loopfile, 'rb') as f:                          #Added b to avoid unicode errors, CCC 03/2021
                self.loop = pickle.load(f, encoding='latin1')        #Added encoding to restart python2 with python3, CCC 27/01/2022
            tprint("Loaded torch loop state from "+loopfile)

            if not refresh:

                rstatefile = path.join(self.output_dir,
                    'rnd_state{:04d}.pickle'.format(self.chknum))
                massesfile = path.join(self.output_dir,
                    'all_masses{:04d}.pickle'.format(self.chknum))

                # Addind these to permit binaries - CCC, May 3, 2020
                systemsfile = path.join(self.output_dir,
                    'system_masses{:04d}.pickle'.format(self.chknum))
                positionsfile = path.join(self.output_dir,
                    'all_positions{:04d}.pickle'.format(self.chknum))
                velocitiesfile = path.join(self.output_dir,
                    'all_velocities{:04d}.pickle'.format(self.chknum))

                with open(rstatefile, 'rb') as f:                 #Added b, CCC 03/2021
                    rnd_state = pickle.load(f, encoding='latin1') #Added encoding to restart python2 with python3, CCC 27/01/2022
                np.random.set_state(rnd_state)
                tprint("Random state set from "+rstatefile)

                with open(massesfile, 'rb') as f:                       #Added b, CCC 03/2021
                    self.all_masses = pickle.load(f, encoding='latin1') #Added encoding to restart python2 with python3, CCC 27/01/2022
                tprint("Loaded all_masses dictionary from "+massesfile)
                
                # Adding these to permit primordial binaries -CCC, May 3, 2020
                with open(systemsfile, 'rb') as f:                          #Added b, CCC 03/2021
                    self.system_masses = pickle.load(f, encoding='latin1')  #Added encoding to restart python2 with python3, CCC 27/01/2022 
                tprint("Loaded system_masses dictionary from "+systemsfile)
                
                with open(positionsfile, 'rb') as f:                        #Added b, CCC 03/2021
                    self.all_positions = pickle.load(f, encoding='latin1')  #Added encoding to restart python2 with python3, CCC 27/01/2022
                tprint("Loaded all_positions dictionary from "+positionsfile)
                
                with open(velocitiesfile, 'rb') as f:                       #Added b, CCC 03/2021
                    self.all_velocities = pickle.load(f, encoding='latin1') #Added encoding to restart python2 with python3, CCC 27/01/2022
                tprint("Loaded all_velocities dictionary from "+velocitiesfile)

            else:

                tprint("WARNING: Refreshing random state with a new seed.")

        else:

            # state from "previous" loop, so gets incremented
            self.loop['it'] = 0
            self.loop['dt'] = 0.0 | units.s

            # This call will try to write chk/plt files.
            # FLASH shouldn't write chk/plt again, but hy_chknum != self.chknum
            # and hy_pltnum != self.pltnum, so we will write torch files.
            self.output(overwrite)

        # After FLASH inits (whether restart or not), it increments
        # io_checkpointFileNumber for the /next/ file write.
        # Must update our value of chknum, too.
        self.chknum = self.hydro.IO_num('chk')

    # Force a checkpoint (to use for stalls), CCC 09/03/2023
    # Also used to ensure a checkpoint is written immediately after sink formation, CCC 04/05/2023
    def force_output(self, overwrite):
        """
        Force write chk and full Torch state.
        """

        # Force a checkpoint to be written, CCC 09/03/2023
        # Torch state files
        self.out_loop()
        self.out_mass()
        self.out_system()
        self.out_position()
        self.out_velocity()
        self.out_rnd()

        hy_chknum = self.hydro.IO_out('chk')
        # Check if a checkpoint would have been written at the end of the step, and update nums accordingly
        if hy_chknum > self.chknum:
            self.chknum = hy_chknum
        else:
            self.hydro.write_chpt()
            self.chknum += 1

                
    def output(self, overwrite):
        """Try to write chk and plt files; if FLASH chk written,
        also dump full Torch state to disk.
        """
        hy_chknum = self.hydro.IO_out('chk')

        # If a checkpoint file was written, dump all Torch state files
        # conditional allows for possibility of rolling chk
        if hy_chknum != self.chknum:
            self.out_loop()
            self.out_mass()
        # Adding the three below to include primordial binaries -CCC, May 3, 2020
            self.out_system()
            self.out_position()
            self.out_velocity()
            self.out_rnd()
            self.chknum = hy_chknum

        hy_pltnum = self.hydro.IO_out('pltpart')

        # If a plt file was written, dump star properties
        if hy_pltnum > self.pltnum:
            self.out_stars(overwrite)
            self.pltnum = hy_pltnum
        elif hy_pltnum < self.pltnum:
            raise Exception("Error: hy_pltnum={} < pltnum={}".format(hy_pltnum, self.pltnum))

    def out_loop(self):
        """Write dict with bridge loop state to pickle"""
        fname = path.join(self.output_dir,
                          "torch_loop{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.loop, f)
        tprint("*** Wrote bridge loop state to {:s} ****".format(fname))

    def out_mass(self):
        """Write dict with all future stars to pickle"""
        fname = path.join(self.output_dir,
                          "all_masses{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.all_masses, f)
        tprint("*** Wrote queued stars to {:s} ****".format(fname))

    # Adding these to permit primordial binaries
    def out_system(self):
        """Write dict with all future system masses to pickle"""
        fname = path.join(self.output_dir,
                          "system_masses{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.system_masses, f)
        tprint("*** Wrote queued system masses to {:s} ****".format(fname))

    def out_position(self):
        """Write dict with all future positions wrt COM to pickle"""
        fname = path.join(self.output_dir,
                          "all_positions{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.all_positions, f)
        tprint("*** Wrote queued positions wrt COM to {:s} ****".format(fname))

    def out_velocity(self):
        """Write dict with all future velocities wrt COM to pickle"""
        fname = path.join(self.output_dir,
                          "all_velocities{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.all_velocities, f)
        tprint("*** Wrote queued velocities wrt COM to {:s} ****".format(fname))

    def out_rnd(self):
        """Write current random number state to pickle"""
        fname = path.join(self.output_dir,
                          "rnd_state{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(np.random.get_state(), f)
        tprint("*** Wrote numpy random state to {:s} ****".format(fname))

    def out_stars(self, overwrite):
        """Write star particles to AMUSE file"""
        stars_fname = path.join(self.output_dir,
                               "stars{:04d}.amuse".format(self.pltnum))
        write_set_to_file(self.stars, stars_fname, format='hdf5', append_to_file=False, overwrite_file=overwrite)  # hdf5 works with Particles(0), csv breaks
        #mult_file = path.join(self.output_dir,
        #                      "mult{:04d}.amuse".format(self.pltnum))
        #multstars = mult.stars.copy_to_new_particles(, format='hdf5')
        #write_set_to_file(multstars, mult_file)
        tprint("*** Wrote existing stars to {:s} ****".format(stars_fname))
        
    def out_merged_stars(self, removed_stars, overwrite):
        """Write merged star particles to AMUSE file"""
        stars_fname = path.join(self.output_dir,
                               "merged_stars.amuse")
        write_set_to_file(removed_stars, stars_fname, format='hdf5', append_to_file=True, overwrite_file=overwrite)  # hdf5 works with Particles(0), csv breaks
        tprint("*** Wrote merged stars to merged_stars.amuse")

    def stars_to_mult_grav_copy(self, attr):
        """
        Copy attribute from stars Particles() set to leaves AND center-of-mass
        particles tracked by both Multiples and gravity code.

        Stars NOT tracked by multiples are not updated.

        mult.root_to_tree is usually not updated.  A pure N-body simulation
        would update only COM particles, barring encounter.
        We must update mult.root_to_tree AND the COM particles because
        (1) gas+sinks can kick individual stars in binaries, and
        (2) stellar evolution modifies individual stars, not COM particles.
        """
        assert attr in ["mass", "velocity"]

        for s in self.stars:
            for root, tree in self.mult.root_to_tree.items():
                leaves = tree.get_leafs_subset()
                if s in leaves:
                    if attr == "mass":
                        s.as_particle_in_set(leaves).mass = s.mass
                    elif attr == "velocity":
                        s.as_particle_in_set(leaves).velocity = s.velocity

        update_roots_from_leaves(self.mult, self.grav)


def update_roots_from_leaves(mult, grav):
    """
    Update the center of mass particles from
    the leaves properties (in all codes!).
    """
    for root, tree in mult.root_to_tree.items():
        leaves = tree.get_leafs_subset()
        msum    = leaves.mass.sum()
        com = leaves.center_of_mass().as_quantity_in(units.cm)
        com_vel = leaves.center_of_mass_velocity().as_quantity_in(units.cm/units.s)

        # Update the root particle in multiples.
        root_particle = root.as_particle_in_set(mult._inmemory_particles)
        root_particle.mass     = msum
        root_particle.position = com
        root_particle.velocity = com_vel
        # Also update the tree.particle thing (top of the tree in the dictionary in multiples).
        tree.particle.mass     = msum
        tree.particle.position = com
        tree.particle.velocity = com_vel
        # The same particle in the N body code must also be updated!
        grav_particle = root.as_particle_in_set(grav.particles)
        grav_particle.mass     = msum
        grav_particle.position = com
        grav_particle.velocity = com_vel

    return
