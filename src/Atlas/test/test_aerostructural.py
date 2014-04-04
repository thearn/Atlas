import os
import numpy as np
import unittest

from scipy.io import loadmat

from Atlas import AeroStructural, Flags, PrescribedLoad

from openmdao.main.api import set_as_top


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


class Test_AeroStructural(unittest.TestCase):
    """
    Test the AeroStructural assembly against the HeliCalc.m module from the MATLAB model
    """

    def test_AeroStructural(self):
        """ Test the AeroStructural assembly
        """
        asm = set_as_top(AeroStructural())

        # populate inputs
        path = os.path.join(os.path.dirname(__file__), 'HeliCalc.mat')
        data = loadmat(path, struct_as_record=True, mat_dtype=True)

        config = asm.config

        config.Ns           = int(data['Ns'][0][0])
        config.ycmax        = data['ycmax'][0]
        config.rho          = data['rho'][0][0]
        config.visc         = data['visc'][0][0]
        config.vw           = data['vw'][0][0]
        config.vc           = data['vc'][0][0]
        config.R            = data['R'][0][0]
        config.b            = int(data['b'][0][0])
        config.h            = data['h'][0][0]
        config.Omega        = data['Omega'][0][0]
        config.c            = data['c_'][0]
        config.Cl           = data['Cl_'][0]
        config.Cm           = data['Cm_'][0]
        config.t            = data['t_'][0]
        config.xtU          = data['xtU_'][0]
        config.xtL          = data['xtL_'][0]
        config.yWire        = data['yWire'][0]
        config.zWire        = data['zWire'][0][0]
        config.tWire        = data['tWire'][0][0]
        config.TWire        = data['TWire'].flatten()
        config.TEtension    = data['TEtension'][0][0]
        config.xEA          = data['xEA_'][0]
        config.d            = data['d_'][0]
        config.theta        = data['theta_'][0]
        config.nTube        = data['nTube_'][0]
        config.nCap         = data['nCap_'][0]
        config.lBiscuit     = data['lBiscuit_'][0]
        config.anhedral     = data['anhedral'][0][0]
        config.collective   = data['collective'][0][0]
        config.dQuad        = data['dQuad'][0][0]
        config.thetaQuad    = data['thetaQuad'][0][0]
        config.nTubeQuad    = int(data['nTubeQuad'][0][0])
        config.lBiscuitQuad = data['lBiscuitQuad'][0][0]
        config.hQuad        = data['hQuad'][0][0]
        config.mElseRotor   = data['mElseRotor'][0][0]
        config.mElseCentre  = data['mElseCentre'][0][0]
        config.mElseR       = data['mElseR'][0][0]
        config.mPilot       = data['mPilot'][0][0]

        config.flags = Flags()
        config.flags.Cover    = int(data['flags']['Cover'][0][0][0][0])
        config.flags.Quad     = int(data['flags']['Quad'][0][0][0][0])
        config.flags.Load     = int(data['flags']['Load'][0][0][0][0])
        config.flags.wingWarp = int(data['flags']['wingWarp'][0][0][0][0])

        CFRPType = int(data['flags']['CFRPType'][0][0][0][0])
        if CFRPType == 1:
            config.flags.CFRPType = 'NCT301-1X HS40 G150 33 +/-2%RW'
        elif CFRPType == 2:
            config.flags.CFRPType = 'HexPly 6376 HTS-12K 268gsm 35%RW'
        elif CFRPType == 3:
            config.flags.CFRPType = 'MTM28-1/IM7-GP-145gsm 32 +/-2%Rw'
        elif CFRPType == 4:
            config.flags.CFRPType = 'HexPly 8552 IM7 160gsm 35%RW'
        elif CFRPType == 5:
            config.flags.CFRPType = 'NCT304-1 HR40 G80 40 +/-2%RW'
        elif CFRPType == 6:
            config.flags.CFRPType = 'MTM28-1B/M46J-140-37%RW'
        else:
            raise Exception('Unable to decode CFRPType from MATLAB data')

        WireType = int(data['flags']['WireType'][0][0][0][0])
        if WireType == 1:
            config.flags.WireType = 'Pianowire'
        elif WireType == 2:
            config.flags.WireType = 'Vectran'
        else:
            raise Exception('Unable to decode WireType from MATLAB data')

        config.presLoad = PrescribedLoad()
        config.presLoad.y = data['presLoad']['y'][0][0][0][0]
        config.presLoad.pointZ = data['presLoad']['pointZ'][0][0][0][0]
        config.presLoad.pointM = data['presLoad']['pointM'][0][0][0][0]
        config.presLoad.distributedX = data['presLoad']['distributedX'][0][0][0][0]
        config.presLoad.distributedZ = data['presLoad']['distributedZ'][0][0][0][0]
        config.presLoad.distributedM = data['presLoad']['distributedM'][0][0][0][0]

        # run
        asm.run()

        # compare outputs to MATLAB data
        for i, val in enumerate(data['out']['alphaJig'][0][0]):
            msg = 'alphaJig[%d] is %f, compared to %f' % (i, asm.results.alphaJig[i], val)
            print msg

        val = data['out']['Ttot'][0][0]
        msg = 'Ttot is %f, compared to %f' % (asm.results.Ttot, val)
        print msg

        val = data['out']['Ptot'][0][0]
        msg = 'Ptot is %f, compared to %f' % (asm.results.Ptot, val)
        print msg

        val = data['out']['MomRot'][0][0]
        msg = 'MomRot is %f, compared to %f' % (asm.results.MomRot, val)
        print msg

        val = data['out']['Mtot'][0][0]
        msg = 'Mtot is %f, compared to %f' % (asm.struc.Mtot, val)
        print msg

        val = data['out']['mQuad'][0][0]
        msg = 'mQuad is %f, compared to %f' % (asm.struc.mass.mQuad, val)
        print msg

        val = data['out']['mCover'][0][0]
        msg = 'mCover is %f, compared to %f' % (asm.struc.mass.mCover, val)
        print msg

        val = data['out']['mWire'][0][0]
        msg = 'mWire is %f, compared to %f' % (asm.struc.mass.mWire, val)
        print msg

        # for i, val in enumerate(data['out']['mSpar'][0][0]):
        #     msg = 'mSpar[%d] is %f, compared to %f' % (i, asm.struc.mass.mSpar[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['mChord'][0][0]):
        #     msg = 'mChord[%d] is %f, compared to %f' % (i, asm.struc.mass.mChord[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['Cl'][0][0]):
        #     msg = 'Cl[%d] is %f, compared to %f' % (i, asm.results.Cl[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['Cm'][0][0]):
        #     msg = 'Cm[%d] is %f, compared to %f' % (i, asm.aero2.Cm[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['xtU'][0][0]):
        #     msg = 'xtU[%d] is %f, compared to %f' % (i, asm.aero2.xtU[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['xtL'][0][0]):
        #     msg = 'xtL[%d] is %f, compared to %f' % (i, asm.aero2.xtL[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['xEA'][0][0]):
        #     msg = 'xEA[%d] is %f, compared to %f' % (i, asm.struc.xEA[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['d'][0][0]):
        #     msg = 'd[%d] is %f, compared to %f' % (i, asm.struc.d[i], val)
        #     print msg

        # for i, val in enumerate(data['out']['q'][0][0]):
        #     msg = 'q[%d] is %f, compared to %f' % (i, asm.struc.q[i], val)
        #     print msg


if __name__ == "__main__":
    unittest.main()
