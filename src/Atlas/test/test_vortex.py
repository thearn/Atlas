import unittest

import os
import numpy as np
from scipy.io import loadmat

import Atlas


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


class Test_Vortex(unittest.TestCase):
    """
    Tests the VortexRing component
    """

    def initialize(self, classname):
        """ create and initialize VortexRing component """

        # create VortexRing or VortexRingC component
        comp = getattr(Atlas, classname)()

        # populate inputs from MATLAB data
        path = os.path.join(os.path.dirname(Atlas.__file__), 'test/vortex.mat')
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

        return comp, data

    def check_outputs(self, comp, data):
        """ check VortexRing outputs against MATLAB data """

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

    def test_VortexRing(self):
        """ Test the python version of the vortex ring component
        """
        comp, data = self.initialize('VortexRing')

        comp.run()

        self.check_outputs(comp, data)

    def test_VortexRingC(self):
        """ Test the cython version of the vortex ring component
        """
        comp, data = self.initialize('VortexRingC')

        comp.run()

        self.check_outputs(comp, data)


if __name__ == "__main__":
    unittest.main()
