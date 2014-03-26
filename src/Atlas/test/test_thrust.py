import unittest
import numpy as np

from openmdao.util.testutil import assert_rel_error

from Atlas import thrust

class ThrustTestCase(unittest.TestCase):
    def setup(self):
        pass
    
    def test_induced_velocity(self):
        comp = thrust.InducedVelocity()
        
        comp.vc = 0.
        comp.b = 2.
        comp.rho = 1.1800
        comp.R = 10.
        comp.h = 1.5000
        
        comp.Ns  = 10

        comp.dr = np.array([
                1., 1., 1., 1.,
                1., 1., 1., 1.,
                1., 1.,
                ])

        comp.r  = np.array([
                0.5, 1.5, 2.5, 3.5,
                4.5, 6.5, 7.5, 8.5,
                9.5, 10.5
                ])

        comp.dT = np.array([
                0.06634278, 1.58865797, 6.80466217, 11.25513959, 
                15.78396483, 27.78374808, 30.97045427, 33.04010879,
                34.03753803,  22.77205544,
                ])

        vi = np.array([
                -1.0343, -0.7853, 0.0044, 0.1058,
                0.1622, 0.1846, 0.1890, 0.2333,
                0.3749, 0.3463,
                ])
        
        comp.execute()
        print comp.vi
        print vi
        assert_rel_error(self, comp.vi, vi, 1e-4)

        
    def test_thrust(self):
        comp = thrust.Thrust()
        comp.yN = np.array([
                0, 1, 2, 3, 
                4, 5, 6, 7,
                8, 9, 10,
                ])

        comp.rho = 1.1800
        comp.ycmax = 1.4656
        comp.Omega = 1.0367

        comp.c = np.array([
                0.279, 1.3903, 1.1757,
                1.0176, 0.8818, 0.7602,
                0.6507, 0.5528, 0.4666,
                0.3925,
                ])

        comp.Cl = np.array([
                1.5000, 1.4987, 1.4604,
                1.4239, 1.3940, 1.3642,
                1.3344, 1.3046, 1.2747,
                0.8299,
                ])

        comp.Ns = 10

        comp.dr = np.array([
                1., 1., 1., 1.,
                1., 1., 1., 1.,
                1., 1.,
                ])

        comp.r  = np.array([
                0.5, 1.5, 2.5, 3.5,
                4.5, 6.5, 7.5, 8.5,
                9.5, 10.5
                ])
     
        chordFrac = np.array([
            1., 0.5344, 1., 1.,
            1., 1., 1., 1.,
            1., 1.,
            ])
        
        dT = np.array([
                0.06634278, 1.58865797, 6.80466217, 11.25513959, 
                15.78396483, 27.78374808, 30.97045427, 33.04010879,
                34.03753803,  22.77205544,
                ])
        
        comp.execute()
        print comp.chordFrac
        assert_rel_error(self, comp.chordFrac, chordFrac, 1e-4)
        assert_rel_error(self, comp.dT, dT, 1e-4)
    
if __name__ == "__main__":
    unittest.main()
