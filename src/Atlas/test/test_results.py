import os
import numpy as np
import unittest

from scipy.io import loadmat

from Atlas import Results, Fblade


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


class Test_Results(unittest.TestCase):
    """
    Tests the Results component
    """

    def test_Results(self):
        """ Test the Results component
        """
        comp = Results()

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'results.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        comp.b          = int(data['b'][0][0])
        comp.yN         = data['yN']
        comp.Ns         = max(comp.yN.shape) - 1
        comp.yE         = data['yE']
        comp.cE         = data['cE']
        comp.Cl         = data['Cl']
        comp.q          = data['q']
        comp.phi        = data['phi']
        comp.collective = int(data['collective'][0][0])

        comp.fblade = Fblade()
        comp.fblade.Fx = data['Fblade']['Fx'][0][0].flatten()
        comp.fblade.Fz = data['Fblade']['Fz'][0][0].flatten()
        comp.fblade.My = data['Fblade']['My'][0][0].flatten()
        comp.fblade.Q  = data['Fblade']['Q'][0][0].flatten()
        comp.fblade.P  = data['Fblade']['P'][0][0].flatten()
        comp.fblade.Pi = data['Fblade']['Pi'][0][0].flatten()
        comp.fblade.Pp = data['Fblade']['Pp'][0][0].flatten()

        # run
        comp.run()

        # check outputs
        for i, val in enumerate(data['di']):
            self.assertAlmostEquals(comp.di[i], val, 4,
                msg='di[%d] is %f, should be %f' % (i, comp.di[i], val))

        for i, val in enumerate(data['alphaJig']):
            self.assertAlmostEquals(comp.alphaJig[i], val, 4,
                msg='alphaJig[%d] is %f, should be %f' % (i, comp.alphaJig[i], val))

        val = data['Ttot'][0][0]
        self.assertAlmostEquals(comp.Ttot, val, 4,
                msg='Ttot is %f, should be %f' % (comp.Ttot, val))

        val = data['Qtot'][0][0]
        self.assertAlmostEquals(comp.Qtot, val, 4,
                msg='Qtot is %f, should be %f' % (comp.Qtot, val))

        val = data['Ptot'][0][0]
        self.assertAlmostEquals(comp.Ptot, val, 4,
                msg='Ptot is %f, should be %f' % (comp.Ptot, val))

        val = data['MomRot'][0][0]
        self.assertAlmostEquals(comp.MomRot, val, 4,
                msg='MomRot is %f, should be %f' % (comp.MomRot, val))


if __name__ == "__main__":
    unittest.main()
