import unittest

from amuse.lab import *

from torch_state import TorchState
from torch_se import (
    stellar_evolution,
    remove_merged_stars,
)
from torch_sf import (
    add_particles_to_grav,
    remove_particles_outside_bndbox,
    make_stars_from_sinks,
    queue_stars,
    random_three_vector,
)

class Test_SE(unittest.TestCase):
    def setUp(self):
        """Runs before each test, set up a new Torch instance"""
        self.grav = None
        self.hydro = None
        self.mult = None

        self.stars = Particles(0)

        se = SeBa()
        se.initialize_code()
        self.se = se 
        self.stars_to_se = self.stars.new_channel_to(self.se.particles)
        self.se_to_stars = self.se.particles.new_channel_to(self.stars)

    def tearDown(self):
        """Runs after each test, cleanup Torch instance"""
        self.se.stop()
        del self.stars


    def test_supernova(self):
        """
        Two part test
        Part 1: checks that then went_sn code block is triggered
        Part 2: checks that the went_sn block is not triggered a second time 
                for already exploded star
        """
        dt = 1e99

        add_star = Particles(2)
        # 50 Msun star should go SN in < 5 Myr
        add_star.mass = [50,50] | units.MSun
        add_star.initial_mass = [50,50] | units.MSun
        add_star.tag = [1, 2]
        add_star.stellar_type = [1, 1] | units.stellar_type # ZAMS star
        add_star.age = [0.0, 0.0] | units.Myr
        time = 10.0 | units.Myr

        self.stars.add_particles(add_star)
        self.se.particles.add_particles(add_star)

        result = stellar_evolution(time, dt, self, self.hydro, self.se,
                                   with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                                   massloss_method="puls", min_feedback_mass=20|units.MSun, unit_test="SUPERNOVA_PART1")
        self.assertEqual(result, 1, msg="SN not triggered.")
        time += 1.0 | units.Myr
        result = stellar_evolution(time, dt, self, self.hydro, self.se,
                                   with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                                   massloss_method="puls", min_feedback_mass=20|units.MSun, unit_test="SUPERNOVA_PART2")
        self.assertEqual(result, 1, msg="Remnant triggered further SNe.")




if __name__ == "__main__":
    unittest.main()
