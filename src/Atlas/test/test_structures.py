from Atlas import Flags, JointProperties, PrescribedLoad, FBlade, Strain, \
                  MassProperties, FEM, Strains, Failures, Structures

from scipy.io import loadmat

import unittest


class TestStructures(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Mass(self):
        """ test of mass properties calculations """
        comp = MassProperties()

        data = loadmat('mass.mat', struct_as_record=True, mat_dtype=True)

        # populate inputs
        comp.flags = Flags()
        comp.flags.Load = int(data['flags']['Cover'][0][0][0][0])
        comp.flags.wingWarp = int(data['flags']['Quad'][0][0][0][0])

        comp.b        = int(data['b'][0][0])

        comp.mSpar    = data['mSpar'].flatten()
        comp.mChord   = data['mChord'].flatten()
        comp.xCGChord = data['xCGChord'].flatten()
        comp.xEA      = data['xEA'].flatten()
        comp.mQuad    = data['mQuad'][0][0]

        comp.ycmax    = data['ycmax'][0][0]

        comp.yWire    = data['yWire'][0]
        comp.zWire    = data['zWire'][0][0]
        comp.tWire    = data['tWire'][0][0]
        comp.RHOWire  = data['RHOWire'][0][0]

        comp.mElseRotor  = data['mElseRotor'][0][0]
        comp.mElseCentre = data['mElseCentre'][0][0]
        comp.mElseR      = data['mElseR'][0][0]
        comp.R           = data['R'][0][0]
        comp.mPilot      = data['mPilot'][0][0]

        # run
        comp.run()

        # check outputs
        self.assertAlmostEquals(comp.Mtot, data['Mtot'][0][0], 10)

        for i, val in enumerate(data['xCG']):
            self.assertAlmostEquals(comp.xCG[i], val, 4,
                msg='xCG[%d] mismatch (%f vs %f)' % (i, comp.xCG[i], val))

    def check_FEM(self, comp, data):
        """ check component FEM results against MATLAB data  """
        for h, plane in enumerate(data['k']):
            for i, row in enumerate(plane):
                for j, val in enumerate(row):
                    self.assertAlmostEquals(comp.k[h, i, j], val, 4,
                        msg='k[%d, %d, %d] mismatch (%f vs %f)' % (h, i, j, comp.k[h, i, j], val))

        for i, row in enumerate(data['K']):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.K[i, j], val, 4,
                    msg='K[%d, %d] mismatch (%f vs %f)' % (i, j, comp.K[i, j], val))

        for i, val in enumerate(data['F']):
            self.assertAlmostEquals(comp.F[i], val, 4,
                msg='F[%d] mismatch (%f vs %f)' % (i, comp.F[i], val))

        for i, val in enumerate(data['q']):
            self.assertAlmostEquals(comp.q[i], val, 4,
                msg='q[%d] mismatch (%f vs %f)' % (i, comp.q[i], val))

    def test_FEM(self):
        """ test of FEM calculations """
        comp = FEM()

        data = loadmat('FEM.mat', struct_as_record=True, mat_dtype=True)

        # populate inputs
        comp.flags = Flags()
        comp.flags.Load = int(data['flags']['Load'][0][0][0][0])
        comp.flags.wingWarp = int(data['flags']['wingWarp'][0][0][0][0])

        comp.yN  = data['yN']

        comp.EIx = data['EIx']
        comp.EIz = data['EIz']
        comp.EA  = data['EA']
        comp.GJ  = data['GJ']

        comp.cE  = data['cE']
        comp.xEA = data['xEA']

        comp.Fblade = FBlade()
        comp.Fblade.Fx = data['Fblade']['Fx'][0][0].flatten()
        comp.Fblade.Fz = data['Fblade']['Fz'][0][0].flatten()
        comp.Fblade.My = data['Fblade']['My'][0][0].flatten()
        comp.Fblade.Q  = data['Fblade']['Q'][0][0].flatten()
        comp.Fblade.P  = data['Fblade']['P'][0][0].flatten()
        comp.Fblade.Pi = data['Fblade']['Pi'][0][0].flatten()
        comp.Fblade.Pp = data['Fblade']['Pp'][0][0].flatten()

        comp.mSpar  = data['mSpar']
        comp.mChord = data['mChord']
        comp.xCG    = data['xCG']

        comp.yWire = data['yWire'].flatten()
        comp.zWire = data['zWire'][0][0]
        comp.TWire = data['TWire'].flatten()

        comp.presLoad = PrescribedLoad()
        comp.presLoad.y = data['presLoad']['y'][0][0][0][0]
        comp.presLoad.pointZ = data['presLoad']['pointZ'][0][0][0][0]
        comp.presLoad.pointM = data['presLoad']['pointM'][0][0][0][0]
        comp.presLoad.distributedX = data['presLoad']['distributedX'][0][0][0][0]
        comp.presLoad.distributedZ = data['presLoad']['distributedZ'][0][0][0][0]
        comp.presLoad.distributedM = data['presLoad']['distributedM'][0][0][0][0]

        # run
        comp.run()

        # check outputs
        self.check_FEM(comp, data)

    def check_strains(self, comp, data):
        """ check component internal force and strain results against MATLAB data  """

        for i, row in enumerate(data['Finternal']):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.Finternal[i, j], val, 4,
                    msg='Finternal[%d, %d] mismatch (%f vs %f)' % (i, j, comp.Finternal[i, j], val))

        for i, row in enumerate(data['strain']['top'][0]):
            for j, val in enumerate(row[0]):
                self.assertAlmostEquals(comp.strain.top[i, j], val, 4,
                    msg='strain.top[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.top[i, j], val))

        for i, row in enumerate(data['strain']['bottom'][0]):
            for j, val in enumerate(row[0]):
                self.assertAlmostEquals(comp.strain.bottom[i, j], val, 4,
                    msg='strain.bottom[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.bottom[i, j], val))

        for i, row in enumerate(data['strain']['back'][0]):
            for j, val in enumerate(row[0]):
                self.assertAlmostEquals(comp.strain.back[i, j], val, 4,
                    msg='strain.back[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.back[i, j], val))

        for i, row in enumerate(data['strain']['front'][0]):
            for j, val in enumerate(row[0]):
                self.assertAlmostEquals(comp.strain.front[i, j], val, 4,
                    msg='strain.bottom[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.front[i, j], val))

        for i, row in enumerate(data['strain']['bending_x']):
            for j, val in enumerate(row[0][0]):
                self.assertAlmostEquals(comp.strain.bending_x[i][j], val, 4,
                    msg='strain.bending_x[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.bending_x[i, j], val))

        for i, row in enumerate(data['strain']['bending_z']):
            for j, val in enumerate(row[0][0]):
                self.assertAlmostEquals(comp.strain.bending_z[i][j], val, 4,
                    msg='strain.bending_z[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.bending_z[i, j], val))

        for i, row in enumerate(data['strain']['axial_y']):
            for j, val in enumerate(row[0][0]):
                self.assertAlmostEquals(comp.strain.axial_y[i][j], val, 4,
                    msg='strain.axial_y[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.axial_y[i, j], val))

        for i, row in enumerate(data['strain']['torsion_y']):
            for j, val in enumerate(row[0][0]):
                self.assertAlmostEquals(comp.strain.torsion_y[i][j], val, 4,
                    msg='strain.torsion_y[%d, %d] mismatch (%f vs %f)' % (i, j, comp.strain.torsion_y[i, j], val))

    def test_Strains(self):
        """ test of internal force and strain calculations  """
        comp = Strains()

        # populate inputs from MATLAB test data
        data = loadmat('strains.mat', struct_as_record=True)
        comp.yN = data['yN']
        comp.d  = data['d']
        comp.k  = data['k']
        comp.F  = data['F']
        comp.q  = data['q']

        # run
        comp.run()

        # check outputs
        self.check_strains(comp, data)

    def check_failures(self, comp, data):
        """ check component failure results against MATLAB data  """

        print
        # print "data['fail']\n", data['fail'].keys()
        print "data['fail']['top']",      data['fail']['top'].shape
        print "data['fail']['top']",      data['fail']['top'][0, 0]['cap'].shape
        print "data['fail']['top']",      data['fail']['top'][0, 0]['plus'].shape
        print "data['fail']['top']",      data['fail']['top'][0, 0]['minus'].shape

        print "data['fail']['top']",      data['fail']['bottom'][0, 0]['cap'].shape
        print "data['fail']['top']",      data['fail']['bottom'][0, 0]['plus'].shape
        print "data['fail']['top']",      data['fail']['bottom'][0, 0]['minus'].shape

        print "data['fail']['top']",      data['fail']['back'][0, 0]['cap'].shape
        print "data['fail']['top']",      data['fail']['back'][0, 0]['plus'].shape
        print "data['fail']['top']",      data['fail']['back'][0, 0]['minus'].shape

        print "data['fail']['top']",      data['fail']['front'][0, 0]['cap'].shape
        print "data['fail']['top']",      data['fail']['front'][0, 0]['plus'].shape
        print "data['fail']['top']",      data['fail']['front'][0, 0]['minus'].shape

        print "data['fail']['buckling']", data['fail']['buckling'][0, 0]['x'].shape
        print "data['fail']['buckling']", data['fail']['buckling'][0, 0]['z'].shape
        print "data['fail']['buckling']", data['fail']['buckling'][0, 0]['torsion'].shape

        print "data['fail']['quad']",     data['fail']['quad'][0, 0]['buckling'].shape
        print "data['fail']['quad']",     data['fail']['quad'][0, 0]['bend'].shape
        print "data['fail']['quad']",     data['fail']['quad'][0, 0]['torsion'].shape
        print "data['fail']['quad']",     data['fail']['quad'][0, 0]['torbuck'].shape
        print

        # top
        for i, row in enumerate(data['fail']['top'][0, 0]['cap'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.top.cap[i, j], val, 4,
                    msg='fail.top.cap[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.top.cap[i, j], val))

        for i, row in enumerate(data['fail']['top'][0, 0]['plus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.top.plus[i, j], val, 4,
                    msg='fail.top.plus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.top.plus[i, j], val))

        for i, row in enumerate(data['fail']['top'][0, 0]['minus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.top.minus[i, j], val, 4,
                    msg='fail.top.minus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.top.minus[i, j], val))

        # bottom
        for i, row in enumerate(data['fail']['bottom'][0, 0]['cap'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.bottom.cap[i, j], val, 4,
                    msg='fail.bottom.cap[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.bottom.cap[i, j], val))

        for i, row in enumerate(data['fail']['bottom'][0, 0]['plus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.bottom.plus[i, j], val, 4,
                    msg='fail.bottom.plus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.bottom.plus[i, j], val))

        for i, row in enumerate(data['fail']['bottom'][0, 0]['minus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.bottom.minus[i, j], val, 4,
                    msg='fail.bottom.minus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.bottom.minus[i, j], val))

        # back
        for i, row in enumerate(data['fail']['back'][0, 0]['cap'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.back.cap[i, j], val, 4,
                    msg='fail.back.cap[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.back.cap[i, j], val))

        for i, row in enumerate(data['fail']['back'][0, 0]['plus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.back.plus[i, j], val, 4,
                    msg='fail.back.plus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.back.plus[i, j], val))

        for i, row in enumerate(data['fail']['back'][0, 0]['minus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.back.minus[i, j], val, 4,
                    msg='fail.back.minus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.back.minus[i, j], val))

        # front
        for i, row in enumerate(data['fail']['front'][0, 0]['cap'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.front.cap[i, j], val, 4,
                    msg='fail.front.cap[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.front.cap[i, j], val))

        for i, row in enumerate(data['fail']['front'][0, 0]['plus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.front.plus[i, j], val, 4,
                    msg='fail.front.plus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.front.plus[i, j], val))

        for i, row in enumerate(data['fail']['front'][0, 0]['minus'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.front.minus[i, j], val, 4,
                    msg='fail.front.minus[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.front.minus[i, j], val))

        # buckling
        for i, row in enumerate(data['buckling']['buckling'][0, 0]['x'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.buckling.x[i, j], val, 4,
                    msg='fail.buckling.x[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.buckling.x[i, j], val))

        for i, row in enumerate(data['buckling']['buckling'][0, 0]['z'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.buckling.z[i, j], val, 4,
                    msg='fail.buckling.z[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.buckling.z[i, j], val))

        for i, row in enumerate(data['buckling']['buckling'][0, 0]['torsion'][0][0]):
            for j, val in enumerate(row):
                self.assertAlmostEquals(comp.fail.buckling.torsion[i, j], val, 4,
                    msg='fail.buckling.torsion[%d, %d] mismatch (%f vs %f)' % (i, j, comp.fail.buckling.torsion[i, j], val))

        # quad  (TODO)

        # wire
        for i, val in enumerate(data['fail']['wire'][0][0][0]):
            self.assertAlmostEquals(comp.fail.wire[i], val, 4,
                msg='fail.wire[%d] mismatch (%f vs %f)' % (i, comp.fail.wire[i], val))

    def test_Failures(self):
        """ test of failure calculations  """
        comp = Failures()

        data = loadmat('failure.mat', struct_as_record=True, mat_dtype=True)

        # populate inputs
        comp.flags = Flags()
        CFRPType = int(data['flags']['CFRPType'][0][0][0][0])
        if CFRPType == 1:
            comp.flags.CFRPType = 'NCT301-1X HS40 G150 33 +/-2%RW'
        elif CFRPType == 2:
            comp.flags.CFRPType = 'HexPly 6376 HTS-12K 268gsm 35%RW'
        elif CFRPType == 3:
            comp.flags.CFRPType = 'MTM28-1/IM7-GP-145gsm 32 +/-2%Rw'
        elif CFRPType == 4:
            comp.flags.CFRPType = 'HexPly 8552 IM7 160gsm 35%RW'
        elif CFRPType == 5:
            comp.flags.CFRPType = 'NCT304-1 HR40 G80 40 +/-2%RW'
        elif CFRPType == 6:
            comp.flags.CFRPType = 'MTM28-1B/M46J-140-37%RW'
        else:
            raise Exception('Unable to decode CFRPType from MATLAB data')
        WireType = int(data['flags']['WireType'][0][0][0][0])
        if WireType == 1:
            comp.flags.WireType = 'Pianowire'
        elif WireType == 2:
            comp.flags.WireType = 'Vectran'
        else:
            raise Exception('Unable to decode WireType from MATLAB data')

        comp.yN       = data['yN']

        comp.Finternal = data['Finternal']

        comp.strain = Strain()
        comp.strain.top       = data['strain']['top'][0][0]
        comp.strain.bottom    = data['strain']['bottom'][0][0]
        comp.strain.back      = data['strain']['back'][0][0]
        comp.strain.front     = data['strain']['front'][0][0]
        comp.strain.bending_x = data['strain']['bending_x'][0][0]
        comp.strain.bending_z = data['strain']['bending_z'][0][0]
        comp.strain.axial_y   = data['strain']['axial_y'][0][0]
        comp.strain.torsion_y = data['strain']['torsion_y'][0][0]

        comp.d        = data['d']
        comp.theta    = data['theta']
        comp.nTube    = data['nTube'].flatten()
        comp.nCap     = data['nCap']

        comp.yWire    = data['yWire'][0]
        comp.zWire    = data['zWire'][0][0]
        comp.EIxJ     = data['EIxJ']
        comp.EIzJ     = data['EIzJ']

        comp.lBiscuit     = data['lBiscuit']
        comp.dQuad        = data['dQuad'][0][0]
        comp.thetaQuad    = data['thetaQuad'][0][0]
        comp.nTubeQuad    = int(data['nTubeQuad'][0][0])
        comp.lBiscuitQuad = data['lBiscuitQuad'][0][0]
        comp.RQuad        = data['RQuad'][0][0]
        comp.hQuad        = data['hQuad'][0][0]
        comp.EIQuad       = data['EIQuad'].flatten()
        comp.GJQuad       = data['GJQuad'].flatten()
        comp.tWire        = data['tWire'][0][0]
        comp.TWire        = data['TWire'].flatten()
        comp.TEtension    = data['TEtension'][0][0]

        comp.b          = int(data['b'][0][0])
        comp.Fblade     = FBlade()
        comp.Fblade.Fx  = data['Fblade']['Fx'][0][0].flatten()
        comp.Fblade.Fz  = data['Fblade']['Fz'][0][0].flatten()
        comp.Fblade.My  = data['Fblade']['My'][0][0].flatten()
        comp.Fblade.Q   = data['Fblade']['Q'][0][0].flatten()
        comp.Fblade.P   = data['Fblade']['P'][0][0].flatten()
        comp.Fblade.Pi  = data['Fblade']['Pi'][0][0].flatten()
        comp.Fblade.Pp  = data['Fblade']['Pp'][0][0].flatten()
        comp.mSpar      = data['mSpar']
        comp.mChord     = data['mChord']
        comp.mElseRotor = data['mElseRotor'][0][0]

        # run
        comp.run()

        # check outputs
        self.check_failures(comp, data)

    def test_Structures(self):
        """ full up test of integrated structures calculations """
        comp = Structures()

        data = loadmat('StrCalc.mat', struct_as_record=True, mat_dtype=True)

        # populate inputs
        comp.flags = Flags()
        comp.flags.Load = int(data['flags']['Cover'][0][0][0][0])
        comp.flags.wingWarp = int(data['flags']['Quad'][0][0][0][0])
        comp.flags.Load = int(data['flags']['Load'][0][0][0][0])
        comp.flags.wingWarp = int(data['flags']['wingWarp'][0][0][0][0])

        CFRPType = int(data['flags']['CFRPType'][0][0][0][0])
        if CFRPType == 1:
            comp.flags.CFRPType = 'NCT301-1X HS40 G150 33 +/-2%RW'
        elif CFRPType == 2:
            comp.flags.CFRPType = 'HexPly 6376 HTS-12K 268gsm 35%RW'
        elif CFRPType == 3:
            comp.flags.CFRPType = 'MTM28-1/IM7-GP-145gsm 32 +/-2%Rw'
        elif CFRPType == 4:
            comp.flags.CFRPType = 'HexPly 8552 IM7 160gsm 35%RW'
        elif CFRPType == 5:
            comp.flags.CFRPType = 'NCT304-1 HR40 G80 40 +/-2%RW'
        elif CFRPType == 6:
            comp.flags.CFRPType = 'MTM28-1B/M46J-140-37%RW'
        else:
            raise Exception('Unable to decode CFRPType from MATLAB data')

        WireType = int(data['flags']['WireType'][0][0][0][0])
        if WireType == 1:
            comp.flags.WireType = 'Pianowire'
        elif WireType == 2:
            comp.flags.WireType = 'Vectran'
        else:
            raise Exception('Unable to decode WireType from MATLAB data')

        comp.yN       = data['yN']
        comp.d        = data['d']
        comp.theta    = data['theta']
        comp.nTube    = data['nTube']
        comp.nCap     = data['nCap']
        comp.lBiscuit = data['lBiscuit']

        comp.Jprop = JointProperties()
        comp.Jprop.d        = data['Jprop']['d'][0][0][0][0]
        comp.Jprop.theta    = data['Jprop']['theta'][0][0][0][0]
        comp.Jprop.nTube    = int(data['Jprop']['nTube'][0][0][0][0])
        comp.Jprop.nCap     = int(data['Jprop']['nCap'][0][0][0][0])
        comp.Jprop.lBiscuit = data['Jprop']['lBiscuit'][0][0][0][0]

        comp.b        = int(data['b'][0][0])
        comp.cE       = data['cE'].flatten()
        comp.xEA      = data['xEA'].flatten()
        comp.xtU      = data['xtU'].flatten()

        comp.dQuad        = data['dQuad'][0][0]
        comp.thetaQuad    = data['thetaQuad'][0][0]
        comp.nTubeQuad    = int(data['nTubeQuad'][0][0])
        comp.lBiscuitQuad = data['lBiscuitQuad'][0][0]
        comp.RQuad        = data['RQuad'][0][0]
        comp.hQuad        = data['hQuad'][0][0]

        comp.ycmax    = data['ycmax'][0][0]

        comp.yWire    = data['yWire'][0]
        comp.zWire    = data['zWire'][0][0]
        comp.tWire    = data['tWire'][0][0]
        comp.TWire    = data['TWire'].flatten()

        comp.mElseRotor  = data['mElseRotor'][0][0]
        comp.mElseCentre = data['mElseCentre'][0][0]
        comp.mElseR      = data['mElseR'][0][0]
        comp.R           = data['R'][0][0]
        comp.mPilot      = data['mPilot'][0][0]

        comp.Fblade = FBlade()
        comp.Fblade.Fx = data['Fblade']['Fx'][0][0].flatten()
        comp.Fblade.Fz = data['Fblade']['Fz'][0][0].flatten()
        comp.Fblade.My = data['Fblade']['My'][0][0].flatten()
        comp.Fblade.Q  = data['Fblade']['Q'][0][0].flatten()
        comp.Fblade.P  = data['Fblade']['P'][0][0].flatten()
        comp.Fblade.Pi = data['Fblade']['Pi'][0][0].flatten()
        comp.Fblade.Pp = data['Fblade']['Pp'][0][0].flatten()

        comp.presLoad = PrescribedLoad()
        comp.presLoad.y = data['presLoad']['y'][0][0][0][0]
        comp.presLoad.pointZ = data['presLoad']['pointZ'][0][0][0][0]
        comp.presLoad.pointM = data['presLoad']['pointM'][0][0][0][0]
        comp.presLoad.distributedX = data['presLoad']['distributedX'][0][0][0][0]
        comp.presLoad.distributedZ = data['presLoad']['distributedZ'][0][0][0][0]
        comp.presLoad.distributedM = data['presLoad']['distributedM'][0][0][0][0]

        # run
        comp.run()

        # check outputs

        self.assertAlmostEquals(comp.Mtot, data['Mtot'], 10)

        for i, val in enumerate(data['q']):
            self.assertAlmostEquals(comp.q[i], val, 4,
                msg='q[%d] mismatch (%f vs %f)' % (i, comp.q[i], val))

        self.check_strains(comp, data)

if __name__ == "__main__":
    unittest.main()
