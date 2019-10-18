#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np
from os import path
import pickle

from amuse.datamodel import Particles
from amuse.io import write_set_to_file

from torch_param import FlashPar
from torch_stdout import tprint

class TorchState(object):
    """
    Container for AMUSE framework global objects, to help
    (1) hold things, (2) perform I/O for all torch workers.
    """

    def __init__(self, hydro, grav, mult):

        self.hydro = hydro
        self.grav  = grav
        self.mult  = mult

        # "Global" AMUSE-level data structures
        self.all_masses = {}
        self.stars = Particles(0)

        self.stars_to_grav = self.stars.new_channel_to(grav.particles)
        self.grav_to_stars = grav.particles.new_channel_to(self.stars)
        if self.mult is not None:
            self.stars_to_mult = self.stars.new_channel_to(mult.stars)
            self.mult_to_stars = mult.stars.new_channel_to(self.stars)

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

    def initial_io(self, refresh=False):
        """Load restart files or write starting Torch state"""

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

            if not refresh:

                rstatefile = path.join(self.output_dir,
                    'rnd_state{:04d}.pickle'.format(self.chknum))
                massesfile = path.join(self.output_dir,
                    'all_masses{:04d}.pickle'.format(self.chknum))
                assert path.isfile(rstatefile)
                assert path.isfile(massesfile)

                with open(rstatefile, 'r') as f:
                    rnd_state = pickle.load(f)
                np.random.set_state(rnd_state)
                tprint("Random state set from "+rstatefile)

                with open(massesfile, 'r') as f:
                    self.all_masses = pickle.load(f)
                tprint("Loaded all_masses dictionary from "+massesfile)

            else:

                tprint("WARNING: Refreshing random state with a new seed.")

        else:
            # This call will try to write chk/plt files.
            # FLASH shouldn't write chk/plt again, but hy_chknum != self.chknum
            # and hy_pltnum != self.pltnum, so we will write torch files.
            self.output()

        # After FLASH inits (whether restart or not), it increments
        # io_checkpointFileNumber for the /next/ file write.
        # Must update our value of chknum, too.
        self.chknum = self.hydro.IO_num('chk')

    def output(self):
        """Try to write chk and plt files; if FLASH chk written,
        also dump full Torch state to disk.
        """
        hy_chknum = self.hydro.IO_out('chk')

        # If a checkpoint file was written, dump all Torch state files
        # conditional allows for possibility of rolling chk
        if hy_chknum != self.chknum:
            self.out_rnd()
            self.out_mass()
            self.chknum = hy_chknum

        hy_pltnum = self.hydro.IO_out('pltpart')

        # If a plt file was written, dump star properties
        if hy_pltnum > self.pltnum:
            self.out_stars()
            self.pltnum = hy_pltnum
        elif hy_pltnum < self.pltnum:
            raise Exception("Error: hy_pltnum={} < pltnum={}".format(hy_pltnum, self.pltnum))

    def out_mass(self):
        """Write dict with all future stars to pickle"""
        fname = path.join(self.output_dir,
                          "all_masses{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(self.all_masses, f)
        tprint("*** Wrote queued stars to {:s} ****".format(fname))

    def out_rnd(self):
        """Write current random number state to pickle"""
        fname = path.join(self.output_dir,
                          "rnd_state{:04d}.pickle".format(self.chknum))
        with open(fname, 'wb') as f:
            pickle.dump(np.random.get_state(), f)
        tprint("*** Wrote numpy random state to {:s} ****".format(fname))

    def out_stars(self):
        """Write star particles to AMUSE file"""
        stars_fname = path.join(self.output_dir,
                               "stars{:04d}.amuse".format(self.pltnum))
        write_set_to_file(self.stars, stars_fname, format='hdf5', append_to_file=False)  # hdf5 works with Particles(0), csv breaks
        #mult_file = path.join(self.output_dir,
        #                      "mult{:04d}.amuse".format(self.pltnum))
        #multstars = mult_grav.stars.copy_to_new_particles(, format='hdf5')
        #write_set_to_file(multstars, mult_file)
        tprint("*** Wrote existing stars to {:s} ****".format(stars_fname))
