import os
import unittest

from scipy.io import loadmat

from Atlas.aero import Aero, Aero2


class AeroTestCase(unittest.TestCase):
    def setup(self):
        pass

    def test_aero(self):
        """ test Aero
            (with simple actuator disk method for induced velocity calculation)
        """
        comp = Aero()

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'aero.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        comp.b        = int(data['b'][0][0])
        comp.yN       = data['yN'].flatten()
        comp.Ns       = max(comp.yN.shape) - 1
        comp.R        = 10.0
        comp.dr       = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
        comp.r        = [ 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
        comp.h        = data['h'][0][0]
        comp.ycmax    = data['ycmax'][0][0]

        comp.rho      = data['rho'][0][0]
        comp.visc     = data['visc'][0][0]
        comp.vw       = data['vw'][0][0]
        comp.vc       = data['vc'][0][0]
        comp.Omega    = data['Omega'][0][0]

        comp.c        = data['cE']
        comp.Cl       = data['Cl']
        comp.d        = data['d']

        comp.yWire    = data['yWire'][0]
        comp.zWire    = data['zWire'][0][0]
        comp.tWire    = data['tWire'][0][0]

        comp.Cm       = data['Cm']
        comp.xtL      = data['xtL']
        comp.xtU      = data['xtU']

        # run
        comp.run()

        # check outputs
        for i, val in enumerate(data['vi']):
            self.assertAlmostEquals(comp.vi[i], val, 4,
                msg='vi[%d] is %f, should be %f' % (i, comp.vi[i], val))

        for i, val in enumerate(data['phi']):
            self.assertAlmostEquals(comp.phi[i], val, 4,
                msg='phi[%d] is %f, should be %f' % (i, comp.phi[i], val))

        for i, val in enumerate(data['Re']):
            self.assertAlmostEquals(comp.Re[i], val, 4,
                msg='Re[%d] is %f, should be %f' % (i, comp.Re[i], val))

        for i, val in enumerate(data['Cd']):
            self.assertAlmostEquals(comp.Cd[i], val, 4,
                msg='Cd[%d] is %f, should be %f' % (i, comp.Cd[i], val))

        for i, val in enumerate(data['Fblade']['Fx'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Fx[i], val, 4,
                msg='Fx[%d] is %f, should be %f' % (i, comp.Fblade.Fx[i], val))

        for i, val in enumerate(data['Fblade']['Fz'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Fz[i], val, 4,
                msg='Fz[%d] is %f, should be %f' % (i, comp.Fblade.Fz[i], val))

        for i, val in enumerate(data['Fblade']['My'][0][0]):
            self.assertAlmostEquals(comp.Fblade.My[i], val, 4,
                msg='My[%d] is %f, should be %f' % (i, comp.Fblade.My[i], val))

        for i, val in enumerate(data['Fblade']['Q'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Q[i], val, 4,
                msg='Q[%d] is %f, should be %f' % (i, comp.Fblade.Q[i], val))

        for i, val in enumerate(data['Fblade']['P'][0][0]):
            self.assertAlmostEquals(comp.Fblade.P[i], val, 4,
                msg='P[%d] is %f, should be %f' % (i, comp.Fblade.P[i], val))

        for i, val in enumerate(data['Fblade']['Pi'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Pi[i], val, 4,
                msg='Pi[%d] is %f, should be %f' % (i, comp.Fblade.Pi[i], val))

        for i, val in enumerate(data['Fblade']['Pp'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Pp[i], val, 4,
                msg='Pp[%d] is %f, should be %f' % (i, comp.Fblade.Pp[i], val))

    def test_aero2(self):
        """ test Aero2
            (with vortex method for induced velocity calculation)
        """
        comp = Aero2()

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'aero2.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        comp.b        = int(data['b'][0][0])
        comp.yN       = data['yN'].flatten()
        comp.Ns       = max(comp.yN.shape) - 1
        comp.R        = 10.0
        comp.dr       = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
        comp.r        = [ 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
        comp.h        = data['h'][0][0]
        comp.ycmax    = data['ycmax'][0][0]

        comp.rho      = data['rho'][0][0]
        comp.visc     = data['visc'][0][0]
        comp.vw       = data['vw'][0][0]
        comp.vc       = data['vc'][0][0]
        comp.Omega    = data['Omega'][0][0]

        comp.c        = data['cE']
        comp.Cl       = data['Cl']
        comp.d        = data['d']

        comp.yWire    = data['yWire'][0]
        comp.zWire    = data['zWire'][0][0]
        comp.tWire    = data['tWire'][0][0]

        comp.Cm       = data['Cm']
        comp.xtL      = data['xtL']
        comp.xtU      = data['xtU']

        comp.q       = data['q']

        comp.anhedral = data['anhedral'][0][0]

        # run
        comp.run()

        # check outputs
        for i, val in enumerate(data['vi']):
            self.assertAlmostEquals(comp.vi[i], val, 4,
                msg='vi[%d] is %f, should be %f' % (i, comp.vi[i], val))

        for i, val in enumerate(data['phi']):
            self.assertAlmostEquals(comp.phi[i], val, 4,
                msg='phi[%d] is %f, should be %f' % (i, comp.phi[i], val))

        for i, val in enumerate(data['Re']):
            self.assertAlmostEquals(comp.Re[i], val, 4,
                msg='Re[%d] is %f, should be %f' % (i, comp.Re[i], val))

        for i, val in enumerate(data['Cd']):
            self.assertAlmostEquals(comp.Cd[i], val, 4,
                msg='Cd[%d] is %f, should be %f' % (i, comp.Cd[i], val))

        for i, val in enumerate(data['Fblade']['Fx'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Fx[i], val, 4,
                msg='Fx[%d] is %f, should be %f' % (i, comp.Fblade.Fx[i], val))

        for i, val in enumerate(data['Fblade']['Fz'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Fz[i], val, 4,
                msg='Fz[%d] is %f, should be %f' % (i, comp.Fblade.Fz[i], val))

        for i, val in enumerate(data['Fblade']['My'][0][0]):
            self.assertAlmostEquals(comp.Fblade.My[i], val, 4,
                msg='My[%d] is %f, should be %f' % (i, comp.Fblade.My[i], val))

        for i, val in enumerate(data['Fblade']['Q'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Q[i], val, 4,
                msg='Q[%d] is %f, should be %f' % (i, comp.Fblade.Q[i], val))

        for i, val in enumerate(data['Fblade']['P'][0][0]):
            self.assertAlmostEquals(comp.Fblade.P[i], val, 4,
                msg='P[%d] is %f, should be %f' % (i, comp.Fblade.P[i], val))

        for i, val in enumerate(data['Fblade']['Pi'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Pi[i], val, 4,
                msg='Pi[%d] is %f, should be %f' % (i, comp.Fblade.Pi[i], val))

        for i, val in enumerate(data['Fblade']['Pp'][0][0]):
            self.assertAlmostEquals(comp.Fblade.Pp[i], val, 4,
                msg='Pp[%d] is %f, should be %f' % (i, comp.Fblade.Pp[i], val))


if __name__ == "__main__":
    unittest.main()
