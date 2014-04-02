# pylint: disable=line-too-long, invalid-name, bad-whitespace, trailing-whitespace, too-many-locals, line-too-long
from Atlas import LiftDrag
import numpy as np
import unittest


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


def absolute_err(x, y):
    return np.abs(x-y).max()


class AtlasTestLiftDrag(unittest.TestCase):

    def test_LiftDrag(self):
        comp = LiftDrag()

        # set inputs
        comp.Ns = 10
        comp.Omega = 1.0367
        comp.r = np.array([
            0.50000, 1.50000, 2.50000, 3.50000, 4.50000,
            5.50000, 6.50000, 7.50000, 8.50000, 9.50000
        ])
        comp.t = np.array([
            0.1400, 0.1400, 0.1400, 0.1400, 0.1400,
            0.1400, 0.1400, 0.1400, 0.1400, 0.1400
        ])
        comp.xtU = np.array([
            0.0500, 0.1500, 0.1500, 0.1500, 0.1500,
            0.1500, 0.1500, 0.1500, 0.1500, 0.1500
        ])
        comp.xtL = np.array([
            0.0500, 0.3000, 0.3000, 0.3000, 0.3000,
            0.3000, 0.3000, 0.3000, 0.3000, 0.3000
        ])
        comp.vw = 0.0
        comp.vc = 0.0
        comp.vi = np.array([
            0.035025, 0.100048, 0.160388, 0.174335, 0.182077,
            0.184890, 0.183906, 0.180033, 0.174075, 0.136191
        ])
        comp.c = np.array([
            0.2729, 1.3903, 1.1757, 1.0176, 0.8818,
            0.7602, 0.6507, 0.5528, 0.4666, 0.3925
        ])
        comp.rho = 1.18
        comp.Cl = np.array([
            1.5000, 1.4987, 1.4604, 1.4239, 1.3940,
            1.3642, 1.3344, 1.3046, 1.2747, 0.8299
        ])
        comp.dr = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        comp.visc = 1.7800e-05
        comp.d = np.array([
            0.0843, 0.0780, 0.0718, 0.0655, 0.0592,
            0.0530, 0.0477, 0.0431, 0.0384, 0.0338
        ])
        comp.yWire = np.array([5.8852, ])
        comp.zWire = 1.0
        comp.yN = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        comp.tWire = 0.0016
        comp.chordFrac = np.array([
            1.00000, 0.53440, 1.00000, 1.00000, 1.00000,
            1.00000, 1.00000, 1.00000, 1.00000, 1.00000
        ])
        comp.Cm = - np.array([
            0.1500, 0.1494, 0.1330, 0.1200, 0.1200,
            0.1200, 0.1200, 0.1200, 0.1200, 0.1200
        ])

        # run
        comp.run()

        # check outputs
        tol = 5e-4
        phi = np.array([
            0.067466, 0.064247, 0.061804, 0.048008, 0.039008,
            0.032414, 0.027284, 0.023150, 0.019751, 0.013827
        ])
        self.assertLess(relative_err(phi, comp.phi), tol)

        Fx = np.array([
            0.0067933, 0.1285898, 0.5270091, 0.7115166, 0.8550350,
            0.9471831, 0.9555218, 0.9603501, 0.9392176, 0.6869232
        ])
        self.assertLess(relative_err(Fx, comp.Fblade.Fx), tol)

        Fz = np.array([
            0.064888, 1.590310, 6.811400, 11.260433, 15.788155,
            19.894738, 23.261571, 25.720063, 27.251109, 18.640423
        ])
        self.assertLess(relative_err(Fz, comp.Fblade.Fz), tol)

        My = np.array([
            -0.0017794, -0.2211508, -0.7315582, -0.9675347, -1.2000506,
            -1.3317214, -1.3621147, -1.3085101, -1.1976834, -1.0583930
        ])
        self.assertLess(relative_err(My, comp.Fblade.My), tol)

        Q = np.array([
            0.0033967, 0.1928848, 1.3175227, 2.4903082, 3.8476573,
            5.2095070, 6.2108916, 7.2026258, 7.9833492, 6.5257705
        ])
        self.assertLess(relative_err(Q, comp.Fblade.Q), tol)

        P = np.array([
            0.0035214, 0.1999686, 1.3659095, 2.5817662, 3.9889647,
            5.4008291, 6.4389902, 7.4671464, 8.2765423, 6.7654332
        ])
        self.assertLess(relative_err(P, comp.Fblade.P), tol)

        Pi = np.array([
            0.0022783, 0.1592760, 1.0935131, 1.9645069, 2.8763529,
            3.6801445, 4.2795547, 4.6319694, 4.7451104, 2.5394582
        ])
        self.assertLess(relative_err(Pi, comp.Fblade.Pi), tol)

        Pp = np.array([
            0.0012431, 0.0406925, 0.2723964, 0.6172593, 1.1126119,
            1.7206846, 2.1594355, 2.8351770, 3.5314319, 4.2259749
        ])
        self.assertLess(relative_err(Pp, comp.Fblade.Pp), tol)

        Cd = np.array([
            0.045366, 0.022832, 0.020535, 0.019363, 0.018735,
            0.018439, 0.018376, 0.018489, 0.018740, 0.019098
        ])
        self.assertLess(relative_err(Cd, comp.Cd), tol)

        Re = np.array([
            9.4000e+03, 1.4362e+05, 2.0239e+05, 2.4506e+05, 2.7293e+05,
            2.8751e+05, 2.9077e+05, 2.8499e+05, 2.7266e+05, 2.5631e+05
        ])
        self.assertLess(relative_err(Re, comp.Re), tol)


if __name__ == "__main__":
    unittest.main()
