#!/usr/bin/env python



import numpy as np
from os import path
import pickle

from amuse.datamodel import Particles
from amuse.io import write_set_to_file
from amuse.units import units
# Import for binaries, CCC 19/07/2023
from amuse.ext.orbital_elements import get_orbital_elements_from_arrays

from torch_param import FlashPar
from torch_stdout import tprint

class TorchState(object):
    """
    Container for AMUSE framework global objects, to help
    (1) hold things, (2) perform I/O for all torch workers.
    """

    def __init__(self, hydro, grav, mult, se): #Add se, CCC 04/11/2023

        self.hydro = hydro
        self.grav  = grav
        self.mult  = mult

        # "Global" AMUSE-level data structures
        self.all_masses = {}
        self.loop = {}
        self.stars = Particles(0)
        self.binaries = Particles(0) # CCC 19/07/2023
        self.stars_next_id = 0  # to supply ID attribute for ph4

        # For primordial binaries -CCC, 01/05/2020
        self.system_masses = {}
        self.all_positions = {}
        self.all_velocities = {}

        self.stars_to_grav = self.stars.new_channel_to(grav.particles)
        self.grav_to_stars = grav.particles.new_channel_to(self.stars)
        # Stellar evolution to stars, CCC 04/11/2023
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
    def force_output(self, overwrite):
        """
        Force write chk and full Torch state.
        """
        hy_chknum = self.hydro.IO_out('chk')
        #print("hy_chknum:", hy_chknum)
        #print("self.chknum:", self.chknum)

        # Force a checkpoint to be written, CCC 09/03/2023
        self.hydro.write_chpt()
        # Torch state files
        self.out_loop()
        self.out_mass()
        self.out_system()
        self.out_position()
        self.out_velocity()
        self.out_rnd()
        self.chknum = hy_chknum

                
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
            self.out_binaries(overwrite) #CCC 19/07/2023
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
        
    def out_binaries(self, overwrite):
        """Write binary particles to AMUSE file"""
        binaries_fname = path.join(self.output_dir,
                               "binaries{:04d}.amuse".format(self.pltnum))
        write_set_to_file(self.binaries, binaries_fname, format='hdf5', append_to_file=False, overwrite_file=overwrite)  # hdf5 works with Particles(0), csv breaks
        #mult_file = path.join(self.output_dir,
        #                      "mult{:04d}.amuse".format(self.pltnum))
        #multstars = mult.stars.copy_to_new_particles(, format='hdf5')
        #write_set_to_file(multstars, mult_file)
        tprint("*** Wrote existing binaries to {:s} ****".format(binaries_fname))

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
        
    def binaries_from_stars(self):
        """
        Identifies binaries from the stars particle set, and saves
        binaries with semi-major axis and eccentricity, as well as
        the particle information.
        """
        
        def relative_v2(vx1, vy1, vz1, vx2, vy2, vz2):
            v2 = (vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2
            return v2

        def relative_r(px1, py1, pz1, px2, py2, pz2):
            r2 = (px1-px2)**2 + (py1-py2)**2 + (pz1-pz2)**2
            r = np.zeros(len(r2)) | units.cm
            for i in range(len(r)):
                r[i] = np.sqrt(r2[i])
            return r

        def relative_r_scalar(px1, py1, pz1, px2, py2, pz2):
            r2 = (px1-px2)**2 + (py1-py2)**2 + (pz1-pz2)**2
            r = np.sqrt(r2)
            return r

        def rel_v(vx1, vy1, vz1, vx2, vy2, vz2):
            v = np.array([vx1-vx2, vy1-vy2, vz1-vz2])
            return v

        def rel_r(px1, py1, pz1, px2, py2, pz2):
            r = np.array([px1-px2, py1-py2, pz1-pz2])
            return r

        def binding_energy(m1, m2, r, v2):
            E = m1 * m2 * v2 / (2 * (m1 + m2)) - ((units.constants.G * m1 * m2) / r)
            return E

        def semi_major(m1, m2, r, v2):
            a = abs(1 / (2/r - v2/(units.constants.G*(m1+m2))))
            return a

        def perturbation(m1, m2, m3, a, d):
            gamma = abs((m1*m3/(d-a)**2) - (m2*m3/(d+a)**2)) * 4 * a**2/(m1*m2)
            return gamma

        def E_bind(m1, m2, a):
            E = 0.5 * units.constants.G * m1 * m2 / a
            return E
        
        stars = self.stars[np.argsort(self.stars.mass.value_in(units.MSun))][::-1]
    
        arg1 = 0
    
        tags_primaries  = []
        tags_companions = []
        singles    = []
        primaries  = []
        companions = []
        rel_pos    = []
        rel_vel    = []
        E_bind     = []
    
        for star1 in stars:

            arg1 += 1
            stars_ = stars[arg1:]
            if len(stars_) == 0:
                break
            arg2 = 0
        
            v2 = relative_v2(star1.vx, star1.vy, star1.vz, stars_.vx, stars_.vy, stars_.vz)
            r  = relative_r(star1.x, star1.y, star1.z, stars_.x, stars_.y, stars_.z)
            E  = binding_energy(star1.mass, stars_.mass, r, v2)
            a  = semi_major(star1.mass, stars_.mass, r, v2)
            # Keep checking for stars in the same location
            if len(r) > len(np.nonzero(r)[0]):
                print('Same loc') 
            # Save all bound companions, then order
            bound = np.where(E < 0 | units.erg)[0]
            for b in bound:
                tags_primaries.append(star1.tag)
                tags_companions.append(stars_[b].tag)
                primaries.append(star1.mass.value_in(units.MSun))
                companions.append(stars_[b].mass.value_in(units.MSun))
                rel_pos.append(rel_r(star1.x.value_in(units.cm), star1.y.value_in(units.cm), star1.z.value_in(units.cm),
                                stars_[b].x.value_in(units.cm), stars_[b].y.value_in(units.cm), stars_[b].z.value_in(units.cm)))
                rel_vel.append(rel_v(star1.vx.value_in(units.cm/units.s), star1.vy.value_in(units.cm/units.s), star1.vz.value_in(units.cm/units.s),
                                stars_[b].vx.value_in(units.cm/units.s), stars_[b].vy.value_in(units.cm/units.s), 
                                stars_[b].vz.value_in(units.cm/units.s)))
                E_bind.append(abs(E[b].value_in(units.erg)))
        
    
        sort = np.argsort(E_bind)[::-1]
        # Keep the most bound for each star
        args = []
        tags_p = []
        tags_c = []
        for s in sort:
            # Skip pair if primary already saved
            if (tags_primaries[s] in tags_p) or (tags_primaries[s] in tags_c):
                pass
            # Skip pair if companion already saved
            elif (tags_companions[s] in tags_p) or (tags_companions[s] in tags_c):
                pass
            # Save if most bound combination
            else:
                args.append(s)
                tags_p.append(tags_primaries[s])
                tags_c.append(tags_companions[s])
    
        args = np.array(args)
    
        tags_primaries  = np.array(tags_primaries)[args]
        tags_companions = np.array(tags_companions)[args]
        primaries  = np.array(primaries)[args] | units.MSun
        companions = np.array(companions)[args] | units.MSun
        rel_pos = np.array(rel_pos)[args] | units.cm
        rel_vel = np.array(rel_vel)[args] | units.cm / units.s
        #E_bind  = np.array(E_bind)[args] | units.erg
    
        semi_major_axes, eccentricities, _, _, _, _ = get_orbital_elements_from_arrays(rel_pos, rel_vel, primaries+companions, G=units.constants.G)
        
        # Data structure
        binaries = Particles(len(tags_primaries))
        binaries.semi_major_axis = semi_major_axes
        binaries.eccentricity = eccentricities
        binaries.initial_semi_major_axis = semi_major_axes
        binaries.initial_eccentricity = eccentricities
        for i in range(len(binaries)):
            binaries[i].child1, binaries[i].child2 = Particles(2)
            j = np.where(stars.tag == tags_primaries[i])[0]
            binaries[i].child1.initial_mass = stars[j].initial_mass
            binaries[i].child1.mass = stars[j].mass
            binaries[i].child1.tag  = stars[j].tag
            k = np.where(stars.tag == tags_companions[i])[0]
            binaries[i].child2.initial_mass = stars[k].initial_mass
            binaries[i].child2.mass = stars[k].mass
            binaries[i].child2.tag  = stars[k].tag

        return binaries


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
