import os
import numpy as np
import unittest

from scipy.io import loadmat

from Atlas import VortexRing, VortexRingC


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


class Test_VortexRing(unittest.TestCase):
    """
    Tests the VortexRing component
    """

    def test_Vortex(self):
        """ Test the vortex ring component
        """
        comp = VortexRing()

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'vortex.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        comp.b        = int(data['b'][0][0])
        comp.yN       = data['yN'].flatten()
        comp.Ns       = max(comp.yN.shape) - 1

        comp.rho      = data['rho'][0][0]
        comp.vc       = data['vc'][0][0]
        comp.Omega    = data['Omega'][0][0]
        comp.h        = data['h'][0][0]
        comp.dT       = data['dT']
        comp.q        = data['q']
        comp.anhedral = data['anhedral'][0][0]

        # run
        comp.run()

        # check outputs
        assert relative_err(comp.vi, data['vi']) < 1e-7

        for i, val in enumerate(data['vi']):
            self.assertAlmostEquals(comp.vi[i], val, 4,
                msg='vi[%d] is %f, should be %f' % (i, comp.vi[i], val))

        for i, row in enumerate(data['ring']['Gamma'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.Gamma[i, j], val, 4,
                    msg='Gamma[%d, %d] is %f, should be %f' % (i, j, comp.Gamma[i, j], val))

        for i, row in enumerate(data['ring']['z'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.z[i, j], val, 4,
                    msg='z[%d, %d] is %f, should be %f' % (i, j, comp.z[i, j], val))

        for i, row in enumerate(data['ring']['r'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.r[i, j], val, 4,
                    msg='r[%d, %d] is %f, should be %f' % (i, j, comp.r[i, j], val))

        for i, row in enumerate(data['ring']['vz'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.vz[i, j], val, 4,
                    msg='vz[%d, %d] is %f, should be %f' % (i, j, comp.vz[i, j], val))

        for i, row in enumerate(data['ring']['vr'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.vr[i, j], val, 4,
                    msg='vr[%d, %d] is %f, should be %f' % (i, j, comp.vr[i, j], val))


class Test_VortexRingC(unittest.TestCase):
    """
    Tests the VortexRingC component
    """

    def test_Vortex(self):
        """ Test the vortex ring component
        """
        comp = VortexRingC()

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'vortex.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        comp.b        = int(data['b'][0][0])
        comp.yN       = data['yN'].flatten()
        comp.Ns       = max(comp.yN.shape) - 1

        comp.rho      = data['rho'][0][0]
        comp.vc       = data['vc'][0][0]
        comp.Omega    = data['Omega'][0][0]
        comp.h        = data['h'][0][0]
        comp.dT       = data['dT']
        comp.q        = data['q']
        comp.anhedral = data['anhedral'][0][0]

        # run
        comp.run()

        # check outputs
        assert relative_err(comp.vi, data['vi']) < 1e-7

        for i, val in enumerate(data['vi']):
            self.assertAlmostEquals(comp.vi[i], val, 4,
                msg='vi[%d] is %f, should be %f' % (i, comp.vi[i], val))

        for i, row in enumerate(data['ring']['Gamma'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.Gamma[i, j], val, 4,
                    msg='Gamma[%d, %d] is %f, should be %f' % (i, j, comp.Gamma[i, j], val))

        for i, row in enumerate(data['ring']['z'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.z[i, j], val, 4,
                    msg='z[%d, %d] is %f, should be %f' % (i, j, comp.z[i, j], val))

        for i, row in enumerate(data['ring']['r'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.r[i, j], val, 4,
                    msg='r[%d, %d] is %f, should be %f' % (i, j, comp.r[i, j], val))

        for i, row in enumerate(data['ring']['vz'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.vz[i, j], val, 4,
                    msg='vz[%d, %d] is %f, should be %f' % (i, j, comp.vz[i, j], val))

        for i, row in enumerate(data['ring']['vr'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.vr[i, j], val, 4,
                    msg='vr[%d, %d] is %f, should be %f' % (i, j, comp.vr[i, j], val))


if __name__ == "__main__":
    unittest.main()
