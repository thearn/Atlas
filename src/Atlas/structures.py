import numpy as np

from math import pi, sqrt, sin, cos, atan2

from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.lib.datatypes.api import Int, Float, Array, VarTree, Str, Enum

from properties import SparProperties, ChordProperties, WireProperties


# data structures used in structural calculations

class Flags(VariableTree):
    Opt          = Int(1, desc='0 - single run, 1 - optimization')
    ConFail      = Int(0, desc='1 - use structural failure as a constraint on optimization')
    ConWireCont  = Int(0, desc='1 - use wire length continuity as a constraint to set appropriate wire forces in multi-point optimizations')
    ConJigCont   = Int(0, desc='1 - use jig continuity')
    ConDef       = Int(0, desc='1 - constraints on maximum deformation of the rotor')
    MultiPoint   = Int(4, desc='0 - single point optimization, 1 - 4 point optimization (h=0.5, h=3, wind case, gravity load)')
    Quad         = Int(1, desc='0 - prop drive, 1 - quad rotor')
    FreeWake     = Int(1, desc='0 - momentum theory, 1 - free vortex ring wake')
    PlotWake     = Int(0, desc='0 - dont plot wake, 1 - plot wake ')
    DynamicClimb = Int(0, desc='0 - vc imposes downward velocity, 1 - vc represents climb (final altitude depends on Nw)')
    Cover        = Int(0, desc='0 - no cover over root rotor blades, 1 - cover')
    Load         = Int(0, desc='0 - normal run, 1 - gravity forces only, 2 - prescribed load from pLoad')
    Cdfit        = Int(1, desc='0 - analytic model for drag coefficient, 1 - curve fit on BE airfoils')
    GWing        = Int(1, desc='0 - Daedalus style wing, 1 - Gossamer style wing (changes amount of laminar flow)')
    AeroStr      = Int(1, desc='0 - Assume flat wing, 1 - take deformation into account')
    Movie        = Int(0, desc='0 - dont save animation, 1 - save animation')
    wingWarp     = Int(0, desc='0 - no twist constraint, >0 - twist constraint at wingWarp')

    CFRPType     = Str('NCT301-1X HS40 G150 33 +/-2%RW', desc='type of carbon fibre reinforced polymer')

    WireType     = Enum('Pianowire',
                        ('Pianowire', 'Vectran'),
                        desc='Material to be used for lift wire')


class JointProperties(VariableTree):
    """ Properties at joint location for buckling analysis """
    theta    = Float(desc='wrap angle')
    nTube    = Int(desc='number of tube layers')
    nCap     = Int(desc='number of cap strips')
    lBiscuit = Float(desc='unsupported biscuit length')


class PrescribedLoad(VariableTree):
    y = Float(9.9999, desc='Point load location')
    pointZ = Float(0.15*9.8, desc='N')
    pointM = Float(0, desc='Nm')
    distributedX = Float(0, desc='N/m')
    distributedZ = Float(0, desc='N/m')
    distributedM = Float(0, desc='Nm/m')


class JointSparProperties(SparProperties):
    """ subclass of SparProperties for the joints
        (needed to dynamically create yN and get other properties from Jprop)
    """

    # inputs
    Jprop = VarTree(JointProperties(), iotype='in')

    def execute(self):
        self.yN = np.array([0, 1]).reshape(1, -1)
        self.d = self.Jprop.d
        self.theta = self.Jprop.theta
        self.nTube = self.Jprop.nTube
        self.nCap = self.Jprop.nCap
        self.lBiscuit = self.Jprop.lBiscuit
        super(JointSparProperties, self).execute()


class QuadSparProperties(SparProperties):
    """ subclass of SparProperties for the QuadCopter-specific spars
        (needed to dynamically create yN and nCap and provide scalar I/O)
    """
    # inputs
    dQuad     = Float(iotype='in', desc='')
    thetaQuad = Float(iotype='in', desc='')
    nTubeQuad = Int(iotype='in', desc='number of tube layers')
    lBiscuitQuad = Float(iotype='in', desc='')

    RQuad  = Float(iotype='in', desc='distance from centre of helicopter to centre of quad rotors')
    hQuad  = Float(iotype='in', desc='height of quad-rotor truss')

    # outputs
    mQuad        = Float(iotype='out', desc='mass of Quad spar (scalar')

    def execute(self):
        lQuad = sqrt(self.RQuad**2 + self.hQuad**2)
        self.yN = np.array([0, lQuad]).reshape(1, -1)
        self.nCap = np.array([0, 0]).reshape(1, -1)
        self.d = [self.dQuad]
        self.theta = [self.thetaQuad]
        self.nTube = [self.nTubeQuad]
        self.lBiscuit = [self.lBiscuitQuad]

        super(QuadSparProperties, self).execute()

        self.mQuad = self.mSpar[0]


class FBlade(VariableTree):
    Fx = Array(desc='')
    Fz = Array(desc='')
    My = Array(desc='')
    Q  = Array(desc='')
    P  = Array(desc='')
    Pi = Array(desc='')
    Pp = Array(desc='')


class Strain(object):
    def __init__(self, Ns):
        self.top    = np.zeros((3, Ns + 1))
        self.bottom = np.zeros((3, Ns + 1))
        self.back   = np.zeros((3, Ns + 1))
        self.front  = np.zeros((3, Ns + 1))

        self.bending_x = np.zeros((1, Ns + 1))
        self.bending_z = np.zeros((1, Ns + 1))
        self.axial_y   = np.zeros((1, Ns + 1))
        self.torsion_y = np.zeros((1, Ns + 1))


# components that perform structural calculations

class MassProperties(Component):
    """
    Computes the total mass and CG of the helicopter
    """

    # inputs
    flags    = VarTree(Flags(), iotype='in')
    b        = Float(iotype='in', desc='number of blades')
    mSpar    = Array(iotype='in', desc='mass of spars')
    mChord   = Array(iotype='in', desc='mass of chords')
    xCGChord = Array(iotype='in', desc='xCG of chords')
    xEA      = Array(iotype='in', desc='')
    mQuad    = Float(iotype='in', desc='')

    # inputs for cover
    ycmax    = Float(iotype='in', desc='')

    # inputs for wire
    yWire   = Array(iotype='in', desc='location of wire attachment along span')
    zWire   = Float(iotype='in', desc='depth of wire attachement')
    tWire   = Float(iotype='in', desc='thickness of wire')
    RHOWire = Float(iotype='in', desc='RHO of wire')

    # inputs for other
    mElseRotor  = Float(iotype='in', desc='')
    mElseCentre = Float(iotype='in', desc='')
    mElseR      = Float(iotype='in', desc='')
    R           = Float(iotype='in', desc='')
    mPilot      = Float(iotype='in', desc='mass of pilot (kg)')

    # outputs
    xCG  = Array(iotype='out', desc='')
    Mtot = Float(0.0, iotype='out', desc='total mass')

    def execute(self):
        self.xCG = (self.xCGChord.dot(self.mChord) + self.xEA.dot(self.mSpar)) / (self.mChord + self.mSpar)

        if self.flags.Cover:
            mCover = (self.ycmax**2 * 0.0528 + self.ycmax * 0.605 / 4) * 1.15
        else:
            mCover = 0

        LWire = sqrt(self.zWire**2 + self.yWire**2)
        mWire = pi * (self.tWire / 2)**2 * self.RHOWire * LWire

        if self.flags.Quad:
            self.Mtot = (np.sum(self.mSpar)*self.b + np.sum(self.mChord)*self.b + mWire*self.b + self.mQuad + mCover) * 4 \
                      + self.mElseRotor + self.mElseCentre + self.mElseR * self.R + self.mPilot
        else:
            self.Mtot = np.sum(self.mSpar)*self.b + np.sum(self.mChord)*self.b + mWire*self.b + mCover \
                      + self.mElseRotor + self.mElseCentre + self.mElseR * self.R + self.mPilot


class FEM(Component):
    """
    Computes the deformation of the spar
    """

    # flags
    flags    = VarTree(Flags(), iotype='in')

    yN  = Array(iotype='in', desc='')
    d   = Array(iotype='in', desc='spar diameter')

    EIx = Array(iotype='in', desc='')
    EIz = Array(iotype='in', desc='')
    EA  = Array(iotype='in', desc='')
    GJ  = Array(iotype='in', desc='')

    cE  = Array(iotype='in', desc='chord of each element')
    xEA = Array(iotype='in', desc='')

    Fblade = VarTree(FBlade(), iotype='in')

    mSpar  = Array(iotype='in', desc='mass of spars')
    mChord = Array(iotype='in', desc='mass of chords')

    xCG = Array(iotype='in', desc='')

    # inputs for wire
    yWire = Array(iotype='in', desc='location of wire attachment along span')
    zWire = Float(iotype='in', desc='depth of wire attachement')
    TWire = Float(iotype='in', desc='')

    presLoad = VarTree(PrescribedLoad(), iotype='in')

    # outputs
    q      = Array(iotype='out', desc='deformation')
    EIQuad = Float(iotype='out', desc='')
    GJQuad = Float(iotype='out', desc='')

    # internal forces and strains
    Finternal = Array(iotype='out', desc='')
    strain    = Array(iotype='out', desc='')

    def execute(self):
        # short aliases
        yN = self.yN
        EIx = self.EIx
        EIz = self.EIz
        EA  = self.EA
        GJ  = self.GJ
        xEA = self.xEA
        cE  = self.cE
        mSpar  = self.mSpar
        mChord  = self.mChord
        xCG = self.xCG
        yWire = self.yWire
        zWire = self.zWire
        TWire = self.TWire
        Fblade = self.Fblade
        presLoad = self.presLoad
        d = self.d

        Ns = max(yN.shape) - 1  # number of elements
        dy = np.zeros((Ns, 1))
        for s in range(1, (Ns+1)):
            dy[(s-1)] = yN[(s)] - yN[(s-1)]  # length of each element

        # FEM computation for structural deformations
        # -------------------------------------------

        # Initialize global stiffness matrix
        K = np.zeros(((Ns+1) * 6, (Ns+1) * 6))  # global stiffnes
        F = np.zeros(((Ns+1) * 6, 1))           # global force vector

        # Create global stiffness maxtrix and force vector
        k = np.zeros((12, 12, Ns))

        for s in range(0, Ns-1):

            # Local elastic stiffness matrix
            k[0,   0, s] = 12 * EIx[s] / (dy[s] * dy[s] * dy[s])
            k[0,   5, s] = -6 * EIx[s] / (dy[s] * dy[s])
            k[5,   0, s] = k[0, 5, s]
            k[0,   6, s] = -12 * EIx[s] / (dy[s] * dy[s] * dy[s])
            k[6,   0, s] = k[0, 6, s]
            k[0,  11, s] = -6 * EIx[s] / (dy[s] * dy[s])
            k[11,  0, s] = k[0, 11, s]
            k[1,   1, s] = EA[s] / dy[s]
            k[1,   7, s] = -EA[s] / dy[s]
            k[7,   1, s] = k[1, 7, s]
            k[2,   2, s] = 12 * EIz[s] / (dy[s] * dy[s] * dy[s])
            k[2,   3, s] = 6 * EIz[s] / (dy[s] * dy[s])
            k[3,   2, s] = k[2, 3, s]
            k[2,   8, s] = -12 * EIz[s] / (dy[s] * dy[s] * dy[s])
            k[8,   2, s] = k[2, 8, s]
            k[2,   9, s] = 6 * EIz[s] / (dy[s] * dy[s])
            k[9,   2, s] = k[2, 9, s]
            k[3,   3, s] = 4 * EIz[s] / dy[s]
            k[3,   8, s] = -6 * EIz[s] / (dy[s] * dy[s])
            k[8,   3, s] = k[3, 8, s]
            k[3,   9, s] = 2 * EIz[s] / dy[s]
            k[9,   3, s] = k[3, 9, s]
            k[4,   4, s] = GJ[s] / dy[s]
            k[4,  10, s] = -GJ[s] / dy[s]
            k[10,  4, s] = k[4, 10, s]
            k[5,   5, s] = 4 * EIx[s] / dy[s]
            k[5,   6, s] = 6 * EIx[s] / (dy[s] * dy[s])
            k[6,   5, s] = k[5, 6, s]
            k[5,  11, s] = 2 * EIx[s] / dy[s]
            k[11,  5, s] = k[5, 11, s]
            k[6,   6, s] = 12 * EIx[s] / (dy[s] * dy[s] * dy[s])
            k[6,  11, s] = 6 * EIx[s] / (dy[s] * dy[s])
            k[11,  6, s] = k[6, 11, s]
            k[7,   7, s] = EA[s] / dy[s]
            k[8,   8, s] = 12 * EIz[s] / (dy[s] * dy[s] * dy[s])
            k[8,   9, s] = -6 * EIz[s] / (dy[s] * dy[s])
            k[9,   8, s] = k[8, 9, s]
            k[9,   9, s] = 4 * EIz[s] / dy[s]
            k[10, 10, s] = GJ[s] / dy[s]
            k[11, 11, s] = 4 * EIx[s] / dy[s]

            # Perform dihedral and sweep rotations here if needed

            # Assemble global stiffness matrix
            K[(6*s):(6*s + 12), (6*s):(6*s + 12)] = \
                K[(6*s):(6*s + 12), (6*s):(6*s + 12)] + k[:, :, s]

            Faero = np.zeros((6, 1))
            if self.flags.Load == 0:  # include aero forces
                # aerodynamic forces
                xAC = 0.25
                Faero[0] = Fblade.Fx[s] / 2
                Faero[1] = 0
                Faero[2] = Fblade.Fz[s] / 2
                Faero[3] = Fblade.Fz[s] * dy[s] / 12
                Faero[4] = Fblade.My[s] / 2 + (xEA[s] - xAC) * cE[s] * Fblade.Fz[s] / 2
                Faero[5] = -Fblade.Fx[s] * dy[s] / 12

            Fg = np.zeros((6, 1))
            Fwire = np.zeros((12, 1))

            if (self.flags.Load == 0) or (self.flags.Load == 1):
                # gravitational forces
                g = 9.81
                Fg[0] = 0
                Fg[1] = 0
                Fg[2] = -(mSpar[s] + mChord[s]) * g / 2
                Fg[3] = -(mSpar[s] + mChord[s]) * g * dy[s] / 12
                Fg[4] = (xCG[s] - xEA[s]) * cE[s] * (mSpar[s] + mChord[s]) * g / 2
                Fg[5] = 0

                # Wire forces (using consistent force vector)
                for w in range(0, len(yWire)-1):
                    if (yWire[(w)] >= yN[s]) and (yWire[(w)] < yN[s+1]):
                        thetaWire = atan2(zWire, yWire[(w)])
                        a = yWire[(w)] - yN[s]
                        L = dy[s]
                        FxWire = -cos(thetaWire) * TWire[(w)]
                        FzWire = -sin(thetaWire) * TWire[(w)]
                        Fwire[0] = 0
                        Fwire[1] = FxWire * (1 - a / L)
                        Fwire[2] = FzWire * (2 * (a / L)**3 - 3 * (a / L)**2 + 1)
                        Fwire[3] = FzWire * a * ((a / L)**2 - 2 * (a / L) + 1)
                        Fwire[4] = 0
                        Fwire[5] = 0
                        Fwire[6] = 0
                        Fwire[7] = FxWire * (a / L)
                        Fwire[8] = FzWire * (- 2 * (a / L)**3 + 3 * (a / L)**2)
                        Fwire[9] = FzWire * a * ((a / L)**2 - (a / L))
                        Fwire[10] = 0
                        Fwire[11] = 0
                    else:
                        Fwire = np.zeros((12, 1))

            Fpres = np.zeros((12, 1))

            if self.flags.Load == 2:
                # Prescribed point load (using consistent force vector)
                if (presLoad.y >= yN[s]) and (presLoad.y < yN[s+1]):
                    a = presLoad.y - yN[s]
                    L = dy[s]
                    Fpres[0]  = 0
                    Fpres[1]  = 0
                    Fpres[2]  = presLoad.pointZ * (2 * (a / L)**3 - 3 * (a / L)**2 + 1)
                    Fpres[3]  = presLoad.pointZ * a * ((a / L)**2 - 2 * (a / L) + 1)
                    Fpres[4]  = presLoad.pointM * (1 - a / L)
                    Fpres[5]  = 0
                    Fpres[6]  = 0
                    Fpres[7]  = 0
                    Fpres[8]  = presLoad.pointZ * (- 2 * (a / L)**3 + 3 * (a / L)**2)
                    Fpres[9]  = presLoad.pointZ * a * ((a / L)**2 - (a / L))
                    Fpres[10] = presLoad.pointM * (a / L)
                    Fpres[11] = 0

                # Prescribed distributed load
                Fpres[0]  = Fpres[0]  + presLoad.distributedX * dy[s] / 2
                Fpres[1]  = Fpres[1]  + 0
                Fpres[2]  = Fpres[2]  + presLoad.distributedZ * dy[s] / 2
                Fpres[3]  = Fpres[3]  + presLoad.distributedZ * dy[s] * dy[s] / 12
                Fpres[4]  = Fpres[4]  + presLoad.distributedM * dy[s] / 2
                Fpres[5]  = Fpres[5]  - presLoad.distributedX * dy[s] * dy[s] / 12
                Fpres[6]  = Fpres[6]  + presLoad.distributedX * dy[s] / 2
                Fpres[7]  = Fpres[7]  + 0
                Fpres[8]  = Fpres[8]  + presLoad.distributedZ * dy[s] / 2
                Fpres[9]  = Fpres[9]  - presLoad.distributedZ * dy[s] * dy[s] / 12
                Fpres[10] = Fpres[10] + presLoad.distributedM * dy[s] / 2
                Fpres[11] = Fpres[11] + presLoad.distributedX * dy[s] * dy[s] / 12

            # Assemble global force vector
            F[(s*6 + 0)]  = F[(s*6 + 0)]  + Fpres[0]  + Fwire[0]  + Fg[0] + Faero[0]  # x force
            F[(s*6 + 1)]  = F[(s*6 + 1)]  + Fpres[1]  + Fwire[1]  + Fg[1] + Faero[1]  # y force
            F[(s*6 + 2)]  = F[(s*6 + 2)]  + Fpres[2]  + Fwire[2]  + Fg[2] + Faero[2]  # z force
            F[(s*6 + 3)]  = F[(s*6 + 3)]  + Fpres[3]  + Fwire[3]  + Fg[3] + Faero[3]  # x moment
            F[(s*6 + 4)]  = F[(s*6 + 4)]  + Fpres[4]  + Fwire[4]  + Fg[4] + Faero[4]  # y moment
            F[(s*6 + 5)]  = F[(s*6 + 5)]  + Fpres[5]  + Fwire[5]  + Fg[5] + Faero[5]  # z moment
            F[(s*6 + 6)]  = F[(s*6 + 6)]  + Fpres[6]  + Fwire[6]  + Fg[0] + Faero[0]  # x force
            F[(s*6 + 7)]  = F[(s*6 + 7)]  + Fpres[7]  + Fwire[7]  + Fg[1] + Faero[1]  # y force
            F[(s*6 + 8)]  = F[(s*6 + 8)]  + Fpres[8]  + Fwire[8]  + Fg[2] + Faero[2]  # z force
            F[(s*6 + 9)]  = F[(s*6 + 9)]  + Fpres[9]  + Fwire[9]  - Fg[3] - Faero[3]  # x moment
            F[(s*6 + 10)] = F[(s*6 + 10)] + Fpres[10] + Fwire[10] + Fg[4] + Faero[4]  # y moment
            F[(s*6 + 11)] = F[(s*6 + 11)] + Fpres[11] + Fwire[11] - Fg[5] - Faero[5]  # z moment

        # Add constraints to all 6 dof at root

        if self.flags.wingWarp > 0:  # Also add wingWarping constraint
            raise Exception('FEM is untested and surely broken for wingWarp > 0')
            ii = np.array([])
            for ss in range(0, ((Ns + 1) * 6 - 1)):
                if (ss > 5) and (ss != self.flags.wingWarp*6 + 5):
                    ii = np.array([ii, ss]).reshape(1, -1)
            Fc = F[(ii-1)]
            Kc = K[(ii-1), (ii-1)]
        else:
            Fc = F[6:]
            Kc = K[6:, 6:]

        # Solve constrained system
        # qc = np.linalg.solve(Kc, Fc)
        qc, _, _, _ = np.linalg.lstsq(Kc, Fc)

        if self.flags.wingWarp > 0:
            self.q[ii, 1] = qc
        else:
            # self.q = np.array([0, 0, 0, 0, 0, 0, qc]).reshape(1, -1)
            self.q = np.append(np.array([0, 0, 0, 0, 0, 0]).reshape(-1, 1), qc).reshape(-1, 1)

        # Compute internal forces and strains
        # -----------------------------------

        Ftemp = np.zeros((12, Ns))
        Finternal = np.zeros((6, Ns + 1))

        strain = Strain(Ns)

        for s in range(0, Ns-1):
            # Determine internal forces acting at the nodes of each element
            Ftemp[:, s] = -(np.dot(k[:, :, s], self.q[s*6:s*6 + 12]) - F[s*6:s*6 + 12]).squeeze()
            # Finternal[:, s] = Ftemp[0:5, s]
            Finternal[0, s] = Ftemp[0, s]  # x-shear load
            Finternal[1, s] = Ftemp[1, s]  # y-axial load
            Finternal[2, s] = Ftemp[2, s]  # z-shear load
            Finternal[3, s] = Ftemp[3, s]  # x-bending moment
            Finternal[4, s] = Ftemp[4, s]  # y-torsional load
            Finternal[5, s] = Ftemp[5, s]  # z-bending moment

            # Determine strains at each node
            x_hat = d[s] / 2
            z_hat = d[s] / 2
            r_hat = d[s] / 2

            # Break out displacement vector for element
            qq = self.q[s*6:s*6 + 12]

            strain.bending_x[0, s] = np.dot(-np.array([(-(6*x_hat) / (dy[s]**2)), ((4*x_hat) / dy[s]), ((6*x_hat) / (dy[s]**2)),  ((2*x_hat) / dy[s])]).reshape(1, -1),
                                            np.array([qq[0], qq[5], qq[6], qq[11]]).reshape(1, -1).T)

            strain.bending_z[0, s] = np.dot(-np.array([(-(6*z_hat) / (dy[s]**2)), ((-4*z_hat) / dy[s]), ((+6*z_hat) / (dy[s]**2)), ((-2*z_hat) / dy[s])]).reshape(1, -1),
                                            np.array([qq[2], qq[3], qq[8], qq[9]]).reshape(1, -1).T)

            strain.axial_y[0, s]   = np.dot(np.array([(-1 / dy[s]), (1 / dy[s])]).reshape(1, -1),
                                            np.array([qq[1], qq[7]]).reshape(1, -1).T)

            strain.torsion_y[0, s] = np.dot(r_hat * np.array([(-1 / dy[s]), (1 / dy[s])]).reshape(1, -1),
                                            np.array([qq[4], qq[10]]).reshape(1, -1).T)

            strain.top   [0, s] = strain.bending_z[0, s] + strain.axial_y[0, s]
            strain.top   [1, s] = 0
            strain.top   [2, s] = strain.torsion_y[0, s]
            strain.bottom[0, s] = -strain.bending_z[0, s] + strain.axial_y[0, s]
            strain.bottom[1, s] = 0
            strain.bottom[2, s] = strain.torsion_y[0, s]
            strain.back  [0, s] = strain.bending_x[0, s] + strain.axial_y[0, s]
            strain.back  [1, s] = 0
            strain.back  [2, s] = strain.torsion_y[0, s]
            strain.front [0, s] = -strain.bending_x[0, s] + strain.axial_y[0, s]
            strain.front [1, s] = 0
            strain.front [2, s] = strain.torsion_y[0, s]

        # Loads at the tip are zero
        Finternal[0, Ns] = 0  # x-shear load
        Finternal[1, Ns] = 0  # y-axial load
        Finternal[2, Ns] = 0  # z-shear load
        Finternal[3, Ns] = 0  # x-bending moment
        Finternal[4, Ns] = 0  # y-torsional load
        Finternal[5, Ns] = 0  # z-bending moment

        # Strains at tip are zero
        strain.top   [0, Ns] = 0
        strain.top   [1, Ns] = 0
        strain.top   [2, Ns] = 0

        strain.bottom[0, Ns] = 0
        strain.bottom[1, Ns] = 0
        strain.bottom[2, Ns] = 0

        strain.back  [0, Ns] = 0
        strain.back  [1, Ns] = 0
        strain.back  [2, Ns] = 0

        strain.front [0, Ns] = 0
        strain.front [1, Ns] = 0
        strain.front [2, Ns] = 0

        # Strains at tip are zero
        strain.bending_x[0, Ns] = 0
        strain.bending_z[0, Ns] = 0
        strain.axial_y  [0, Ns] = 0
        strain.torsion_y[0, Ns] = 0


class Structural(Assembly):
    """
    structural computation, first computes the mass of the helicopter based on
    the structural description of the spars and chord lengths. It then
    computes the deformation of the spars, the strains, and the resulting
    factor of safety for each of the failure modes.

    3. Computes the factor of safety for each of the failure modes of the
       spar.
    """

    # flags
    flags    = VarTree(Flags(), iotype='in')

    # inputs for spars
    yN       = Array(iotype='in', desc='node locations for each element along the span')
    d        = Array(iotype='in', desc='spar diameter')
    theta    = Array(iotype='in', desc='wrap angle')
    nTube    = Array(iotype='in', desc='number of tube layers')
    nCap     = Array(iotype='in', desc='number of cap strips')
    lBiscuit = Array(iotype='in', desc='unsupported biscuit length')

    # joint properties
    Jprop    = VarTree(JointProperties(), iotype='in')

    # inputs for chord
    b   = Float(iotype='in', desc='number of blades')
    cE  = Array(iotype='in', desc='chord of each element')
    xEA = Array(iotype='in', desc='')
    xtU = Array(iotype='in', desc='')

    # inputs for quad
    dQuad        = Float(iotype='in', desc='diameter of quad rotor struts')
    thetaQuad    = Float(iotype='in', desc='wrap angle of quad rotor struts')
    nTubeQuad    = Int(iotype='in', desc='number of CFRP layers in quad rotor struts')
    lBiscuitQuad = Float(iotype='in', desc='')
    RQuad        = Float(iotype='in', desc='distance from centre of helicopter to centre of quad rotors')
    hQuad        = Float(iotype='in', desc='height of quad-rotor truss')

    # inputs for cover
    ycmax        = Float(iotype='in', desc='')

    # inputs for wire
    yWire        = Array(iotype='in', desc='location of wire attachment along span')
    zWire        = Float(iotype='in', desc='depth of wire attachement')
    tWire        = Float(iotype='in', desc='thickness of wire')
    TWire        = Float(iotype='in', desc='thickness of wire')

    # inputs for 'other stuff'
    mElseRotor   = Float(iotype='in', desc='')
    mElseCentre  = Float(iotype='in', desc='')
    mElseR       = Float(iotype='in', desc='')
    R            = Float(iotype='in', desc='rotor radius')
    mPilot       = Float(iotype='in', desc='mass of pilot')

    # inputs for FEM
    presLoad     = VarTree(PrescribedLoad(), iotype='in')

    # outputs
    Mtot   = Float(iotype='out', desc='total helicopter mass')
    mSpar  = Array(iotype='out', desc='mass of spars')
    mChord = Array(iotype='out', desc='mass ribs, skin, trailing edge, LES for each spanwise element')
    mQuad  = Float(iotype='out', desc='mass of quad rotor struts')
    mCover = Array(iotype='out', desc='mass of covering')
    mWire  = Array(iotype='out', desc='mass of wire')

    def configure(self):
        self.add('spar', SparProperties())
        self.connect('yN', 'spar.yN')
        self.connect('d', 'spar.d')
        self.connect('theta', 'spar.theta')
        self.connect('nTube', 'spar.nTube')
        self.connect('nCap', 'spar.nCap')
        self.connect('lBiscuit', 'spar.lBiscuit')
        self.connect('flags.CFRPType', 'spar.CFRPType')

        self.add('joint', JointSparProperties())
        self.connect('Jprop', 'joint.Jprop')
        self.connect('flags.CFRPType', 'joint.CFRPType')

        self.add('chord', ChordProperties())
        self.connect('yN', 'chord.yN')
        self.connect('cE', 'chord.c')
        self.connect('d',  'chord.d')
        self.connect('flags.GWing', 'chord.GWing')
        self.connect('xtU', 'chord.xtU')

        self.add('quad', QuadSparProperties())
        self.connect('dQuad', 'quad.dQuad')
        self.connect('thetaQuad', 'quad.thetaQuad')
        self.connect('nTubeQuad', 'quad.nTubeQuad')
        self.connect('lBiscuitQuad', 'quad.lBiscuitQuad')
        self.connect('flags.CFRPType', 'quad.CFRPType')
        self.connect('RQuad', 'quad.RQuad')
        self.connect('hQuad', 'quad.hQuad')

        self.add('wire', WireProperties())
        self.connect('flags.WireType', 'wire.material')

        self.add('mass', MassProperties())
        self.connect('flags', 'mass.flags')
        self.connect('spar.mSpar', 'mass.mSpar')
        self.connect('chord.mChord', 'mass.mChord')
        self.connect('chord.xCGChord', 'mass.xCGChord')
        self.connect('quad.mQuad', 'mass.mQuad')
        self.connect('wire.RHO', 'mass.RHOWire')

        self.connect('xEA', 'mass.xEA')
        self.connect('ycmax', 'mass.ycmax')
        self.connect('zWire', 'mass.zWire')
        self.connect('yWire', 'mass.yWire')
        self.connect('mElseRotor', 'mass.mElseRotor')
        self.connect('mElseCentre', 'mass.mElseCentre')
        self.connect('mElseR', 'mass.mElseR')
        self.connect('R', 'mass.R')
        self.connect('mPilot', 'mass.mPilot')

        self.add('fem', FEM())
        self.connect('yN', 'fem.yN')
        self.connect('d', 'fem.d')
        self.connect('spar.EIx', 'fem.EIx')
        self.connect('spar.EIz', 'fem.EIz')
        self.connect('spar.EA', 'fem.EA')
        self.connect('spar.GJ', 'fem.GJ')
        self.connect('cE', 'fem.cE')
        self.connect('xEA', 'fem.xEA')
        self.connect('spar.mSpar', 'fem.mSpar')
        self.connect('chord.mChord', 'fem.mChord')
        self.connect('mass.xCG', 'fem.xCG')
        self.connect('zWire', 'fem.zWire')
        self.connect('yWire', 'fem.yWire')
        self.connect('TWire', 'fem.TWire')
        self.connect('presLoad', 'fem.presLoad')



    # Compute factor of safety for each failure mode
    # ----------------------------------------------

    # factor of safety for each failure mode
    # fail = Array(iotype='out', desc='')

    # TQuad = np.sum(Fblade.Fz) * b - (np.sum(mSpar + mChord) * b + mElseRotor / 4) * 9.81

    # fail = FailureCalc(yN, Finternal, strain, d, theta, nTube, nCap, yWire, zWire, EIxJ, EIzJ, lBiscuit,
    #                    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, TQuad, EIQuad, GJQuad,
    #                    tWire, TWire, TEtension, flags)
    # return Mtot,mSpar,mChord,mQuad,mCover,mWire,EIx,EIz,EA,GJ,q,EIQuad,GJQuad,Finternal,strain,fail
