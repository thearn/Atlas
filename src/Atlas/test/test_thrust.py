import unittest
import numpy as np

from openmdao.util.testutil import assert_rel_error

from Atlas import thrust


class ThrustTestCase(unittest.TestCase):
    def setup(self):
        pass

    def test_induced_velocity(self):
        comp = thrust.ActuatorDiskInducedVelocity()

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
            4.5, 5.5, 6.5, 7.5,
            8.5, 9.5
        ])

        comp.dT = np.array([
            0.064902, 1.588709, 6.804887, 11.255622,
            15.785445, 19.894067, 23.261656, 25.721612,
            27.253708, 18.644574,
        ])

        vi = np.array([
            0.035025, 0.100048, 0.160388, 0.174335,
            0.182077, 0.184890, 0.183906, 0.180033,
            0.174075, 0.136191,
        ])

        comp.execute()
        assert_rel_error(self, comp.vi, vi, 1e-4)

    def test_thrust(self):
        comp = thrust.Thrust()
        comp.yN = np.array([
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        ])

        comp.rho = 1.1800
        comp.ycmax = 1.4656
        comp.Omega = 1.0367

        comp.c = np.array([
            0.27293, 1.39029, 1.17568,
            1.01762, 0.88181, 0.76021,
            0.65066, 0.55275, 0.46665,
            0.39253,
        ])

        comp.Cl = np.array([
            1.5000,  1.49868, 1.46041,
            1.42387, 1.39404, 1.36422,
            1.33439, 1.30456, 1.27474,
            0.82994,
        ])

        comp.Ns = 10

        comp.dr = np.array([
            1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        ])

        comp.r  = np.array([
            0.5, 1.5, 2.5, 3.5,
            4.5, 5.5, 6.5, 7.5,
            8.5, 9.5
        ])

        chordFrac = np.array([
            1., 0.5344, 1., 1.,
            1., 1., 1., 1.,
            1., 1.,
        ])

        dT = np.array([
            0.064902, 1.588709, 6.804887, 11.255622,
            15.785445, 19.894067, 23.261656, 25.721612,
            27.253708, 18.644574,
        ])

        comp.execute()
        assert_rel_error(self, comp.chordFrac, chordFrac, 1e-4)
        assert_rel_error(self, comp.dT, dT, 1e-4)


if __name__ == "__main__":
    unittest.main()
