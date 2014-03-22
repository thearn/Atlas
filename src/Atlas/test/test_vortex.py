from testvals import values
from Atlas import vortexRing, inducedVelocity
import numpy as np
import unittest

def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()

class Test_vortexRing(unittest.TestCase):
    """
    Tests the vortexRing component
    """

    @classmethod
    def setUpClass(cls):
        """
        This sets up a component and runs it once, using stored input values.
        Each test below then checks individual outputs that result from the run.
        """
        cls.comp = vortexRing()
        cls.comp.yN = values.yN
        cls.comp.rho = values.rho
        cls.comp.dT = values.dT
        cls.comp.vc = values.vc
        cls.comp.Omega = values.Omega
        cls.comp.b = values.b
        cls.comp.h = values.h
        cls.comp.Nw = values.Nw
        cls.comp.Ntt = values.Ntt
        cls.comp.Ntheta = values.Ntheta
        cls.comp.qh = values.qh

        cls.comp.run()

    def test_gamma(self):
        assert relative_err(self.comp.Gamma, values.Gamma) < 1e-7

    def test_z(self):
        assert relative_err(self.comp.z, values.z) < 1e-7

    def test_r(self):
        assert relative_err(self.comp.r, values.r) < 1e-7

    def test_vz(self):
        assert relative_err(self.comp.vz, values.vz) < 1e-7


    def test_vr(self):
        assert relative_err(self.comp.vr, values.vr) < 1e-7

class Test_inducedVelocity(unittest.TestCase):

    def test_vi(self):

        comp = inducedVelocity()
        comp.qh = values.qh
        comp.Gamma = values.Gamma
        comp.z = values.z
        comp.r = values.r
        comp.thetaArray = values.thetaArray
        comp.yE = values.yE
        comp.cr = values.cr
        comp.Ns = values.Ns
        comp.Nw = values.Nw
        comp.dtheta = values.dtheta

        comp.run()
        assert relative_err(comp.vi, values.vi) < 1e-7




if __name__ == "__main__":
    unittest.main()