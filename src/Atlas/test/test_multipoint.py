import numpy as np
import unittest
from nose import SkipTest

from Atlas import Multipoint, HeliOptM

from openmdao.main.api import set_as_top


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


class Test_Multipoint(unittest.TestCase):
    """
    Test the MultiPoint assembly against the HeliOpt.m module from the MATLAB model
    """

    def test_Multipoint(self):
        """ Test the MultiPoint assembly
        """
        mp = set_as_top(Multipoint())

        # set inputs
        mp.alt_low   = 0.5       # low altitude
        mp.alt_high  = 3.5       # high altitude
        mp.alt_ratio = 35./60.   # proportion of time near ground

        mp.TWire_high = 900
        mp.TWire_wind = 2100
        mp.TWire_grav = 110

        mp.OmegaRatio = 2

        mp.vw = 0/3.6

        mp.Cl_max = [1.4, 1.35, 1.55]  # max control

        mp.Omega_low  = 1.0512
        mp.Omega_high = 1.0771
        mp.Cl0_high   = 1.4000
        mp.Cl1_high   = 1.3000

        # run
        mp.run()

        # check outputs against values from MATLAB code
        self.assertAlmostEquals(mp.low.Ptot,   421.3185, 0)
        self.assertAlmostEquals(mp.Mtot_low,   126.1670, 0)
        self.assertAlmostEquals(mp.Ttot_low,  1236.4369, 0)

        self.assertAlmostEquals(mp.high.Ptot,  846.6429, 0)
        self.assertAlmostEquals(mp.Mtot_high,  126.1670, 0)
        self.assertAlmostEquals(mp.Ttot_high, 1236.4369, 0)

        # Future Work:
        # self.assertAlmostEquals(mp.wind.Ptot, 1815.7617, 0)
        # self.assertAlmostEquals(mp.wind.Mtot,  126.1670, 0)
        # self.assertAlmostEquals(mp.wind.Ttot, 2239.7556, 0)

        # self.assertAlmostEquals(mp.grav.Ptot,   0.00000, 0)
        # self.assertAlmostEquals(mp.grav.Mtot,  126.1670, 0)
        # self.assertAlmostEquals(mp.grav.Ttot,   0.00000, 0)

        self.assertAlmostEquals(mp.P, 598.537, 0)

    def test_HeliOptM(self):
        """ Test the multipoint optimization (HeliOptM) assembly
        """
        # comment the following line to run this test
        raise SkipTest('Test skipped to save execution time')

        opt = set_as_top(HeliOptM())

        # run
        opt.run()

        # check results against values from MATLAB code
        self.assertAlmostEquals(opt.mp.P,         598.537, 0)
        self.assertAlmostEquals(opt.mp.Omega_low,  1.0512, 2)
        self.assertAlmostEquals(opt.mp.Omega_high, 1.0771, 2)

        self.assertAlmostEquals(0, opt.mp.Mtot_low*9.8-opt.mp.Ttot_low, 3)
        self.assertAlmostEquals(0, opt.mp.Mtot_high*9.8-opt.mp.Ttot_high, 3)


if __name__ == "__main__":
    unittest.main()
