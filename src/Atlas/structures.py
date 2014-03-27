import numpy as np

from math import pi, sqrt, sin, cos, tan, atan2

from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.lib.datatypes.api import Int, Float, Array, VarTree, Str, Enum

from properties import SparProperties, ChordProperties, wireProperties, prepregProperties


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
    d        = Float(desc='diameter')
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
        self.yN = np.array([0, 1])
        self.d = [self.Jprop.d]
        self.theta = [self.Jprop.theta]
        self.nTube = [self.Jprop.nTube]
        self.nCap = [self.Jprop.nCap]
        self.lBiscuit = [self.Jprop.lBiscuit]
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
        self.yN = np.array([0, lQuad])
        self.nCap = np.array([0, 0])
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


class Strain(VariableTree):
    top    = Array(desc='')
    bottom = Array(desc='')
    back   = Array(desc='')
    front  = Array(desc='')

    bending_x = Array(desc='')
    bending_z = Array(desc='')
    axial_y   = Array(desc='')
    torsion_y = Array(desc='')


class MaterialFailure(VariableTree):
    cap    = Array(desc='')
    plus   = Array(desc='')
    minus  = Array(desc='')


class BucklingFailure(VariableTree):
    x = Array(desc='Euler Buckling failure in main spar from wire force')
    z = Array(desc='Euler Buckling failure in main spar from wire force')
    torsion = Array(desc='Torsional Buckling failure')


class Failure(VariableTree):
    top    = VarTree(MaterialFailure())
    bottom = VarTree(MaterialFailure())
    back   = VarTree(MaterialFailure())
    front  = VarTree(MaterialFailure())

    buckling = VarTree(BucklingFailure())

    quad_buckling = Float(desc='Quad Buckling failure')
    quad_bend     = Float(desc='Quad bending moment failure')
    quad_torsion  = Float(desc='Quad torsional material failure')
    quad_torbuck  = Float(desc='Quad torsional buckling failure')

    wire = Array(desc='Wire tensile failure')


# components that perform structural calculations

class MassProperties(Component):
    """
    Computes the total mass and CG of the helicopter
    """

    # inputs
    flags    = VarTree(Flags(), iotype='in')
    b        = Int(iotype='in', desc='number of blades')
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
        self.xCG = ((self.xCGChord * self.mChord) + (self.xEA * self.mSpar)) / (self.mChord + self.mSpar)

        if self.flags.Cover:
            mCover = (self.ycmax**2 * 0.0528 + self.ycmax * 0.605 / 4) * 1.15
        else:
            mCover = 0

        wire_props = wireProperties[self.flags.WireType]

        LWire = sqrt(self.zWire**2 + self.yWire**2)
        mWire = pi * (self.tWire / 2)**2 * wire_props['RHO'] * LWire

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

    EIx = Array(iotype='in', desc='')
    EIz = Array(iotype='in', desc='')
    EA  = Array(iotype='in', desc='')
    GJ  = Array(iotype='in', desc='')

    cE  = Array(iotype='in', desc='chord of each element')
    xEA = Array(iotype='in', desc='')

    Fblade = VarTree(FBlade(), iotype='in')

    mSpar  = Array(iotype='in', desc='mass of spars')
    mChord = Array(iotype='in', desc='mass of chords')
    xCG    = Array(iotype='in', desc='')

    # inputs for wire
    yWire = Array(iotype='in', desc='location of wire attachment along span')
    zWire = Float(iotype='in', desc='depth of wire attachement')
    TWire = Array(iotype='in', desc='')

    presLoad = VarTree(PrescribedLoad(), iotype='in')

    # outputs
    k = Array(iotype='out', desc='local elastic stiffness matrix')
    K = Array(iotype='out', desc='global stiffness matrix')
    F = Array(iotype='out', desc='global force vector')
    q = Array(iotype='out', desc='deformation')

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

        Ns = max(yN.shape) - 1  # number of elements
        dy = np.zeros((Ns, 1))
        for s in range(1, (Ns+1)):
            dy[(s-1)] = yN[(s)] - yN[(s-1)]  # length of each element

        # FEM computation for structural deformations
        # -------------------------------------------

        # Initialize global stiffness matrix
        K = np.zeros(((Ns+1)*6, (Ns+1)*6))  # global stiffness
        F = np.zeros(((Ns+1)*6, 1))         # global force vector

        # Create global stiffness maxtrix and force vector
        k = np.zeros((12, 12, Ns))

        for s in range(0, Ns):

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
                for w in range(0, len(yWire)):
                    if (yWire[w] >= yN[s]) and (yWire[w] < yN[s+1]):
                        thetaWire = atan2(zWire, yWire[w])
                        a = yWire[w] - yN[s]
                        L = dy[s]
                        FxWire = -cos(thetaWire) * TWire[w]
                        FzWire = -sin(thetaWire) * TWire[w]
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

        # output the stiffness and force arrays as well (for use in failure analysis)
        self.k = k
        self.K = K
        self.F = F


class Strains(Component):
    """
    Computes internal forces and strains
    """

    # inputs
    yN = Array(iotype='in', desc='')
    d  = Array(iotype='in', desc='')
    k  = Array(iotype='in', desc='Local elastic stiffness matrix')
    F  = Array(iotype='in', desc='global force vector')
    q  = Array(iotype='in', desc='deformation')

    # outputs
    Finternal = Array(iotype='out', desc='internal forces')
    strain    = VarTree(Strain(), iotype='out', desc='strains')

    def execute(self):
        # short alias
        d = self.d

        Ns = max(self.yN.shape) - 1  # number of elements
        dy = np.zeros((Ns, 1))
        for s in range(1, (Ns+1)):
            dy[(s-1)] = self.yN[(s)] - self.yN[(s-1)]  # length of each element

        Ftemp = np.zeros((12, Ns))
        Finternal = np.zeros((6, Ns + 1))

        strain = Strain()
        strain.top    = np.zeros((3, Ns + 1))
        strain.bottom = np.zeros((3, Ns + 1))
        strain.back   = np.zeros((3, Ns + 1))
        strain.front  = np.zeros((3, Ns + 1))

        strain.bending_x = np.zeros((1, Ns + 1))
        strain.bending_z = np.zeros((1, Ns + 1))
        strain.axial_y   = np.zeros((1, Ns + 1))
        strain.torsion_y = np.zeros((1, Ns + 1))

        for s in range(0, Ns):
            # Determine internal forces acting at the nodes of each element
            Ftemp[:, s] = -(np.dot(self.k[:, :, s], self.q[s*6:s*6 + 12]) - self.F[s*6:s*6 + 12]).squeeze()
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

        # set outputs
        self.Finternal = Finternal
        self.strain = strain


class Failures(Component):
    """
    Computes the factor of safety for each of the failure modes of the spar.
    """
    # inputs
    flags        = VarTree(Flags(), iotype='in')

    yN           = Array(iotype='in', desc='')

    Finternal    = Array(iotype='in', desc='')
    strain       = VarTree(Strain(), iotype='in')

    d            = Array(iotype='in', desc='')
    theta        = Array(iotype='in', desc='')
    nTube        = Array(iotype='in', desc='')
    nCap         = Array(iotype='in', desc='')

    yWire        = Array(iotype='in', desc='')
    zWire        = Float(iotype='in', desc='')
    EIxJ         = Array(iotype='in', desc='')
    EIzJ         = Array(iotype='in', desc='')

    lBiscuit     = Array(iotype='in', desc='')
    dQuad        = Float(iotype='in', desc='')
    thetaQuad    = Float(iotype='in', desc='')
    nTubeQuad    = Float(iotype='in', desc='')
    lBiscuitQuad = Float(iotype='in', desc='')
    RQuad        = Float(iotype='in', desc='')
    hQuad        = Float(iotype='in', desc='')
    EIQuad       = Array(iotype='in', desc='')
    GJQuad       = Array(iotype='in', desc='')
    tWire        = Float(iotype='in', desc='')
    TWire        = Array(iotype='in', desc='')
    TEtension    = Float(iotype='in', desc='')

    # all this to get TQuad... maybe should be split out
    b            = Int(iotype='in', desc='number of blades')
    Fblade       = VarTree(FBlade(), iotype='in')
    mSpar        = Array(iotype='in', desc='mass of spars')
    mChord       = Array(iotype='in', desc='mass of chords')
    mElseRotor   = Float(iotype='in', desc='')

    # outputs
    fail         = VarTree(Failure(), iotype='out')

    def execute(self):
        # Compute factor of safety for each failure mode
        # ----------------------------------------------
        # computes material failure, euler buckling failure and torsional
        # buckling failure. All failure modes are computed at the nodes.
        #
        # For material failure it looks at a sample laminate on the top,
        # bottom, back and front of the spar. For each side of the spar
        # it computes the material failure for the cap, the plus angle
        # plys and the minus angle plys. For each ply it computes the
        # fracture of failure in the fibre direction, matrix direction
        # and shear. Thus there are 4x3x3=36 (side x lamina x direction)
        # possible failure modes in the tube.
        # ex. fail.top.cap = 3x(Ns+1) vector
        #     fail.top.plus = 3x(Ns+1) vector
        #     fail.top.minus = 3x(Ns+1) vector
        #
        # Stresses and strain are given as a 3x(Ns+1) vector
        # ex. sigma_11, sigma_22, sigma_12
        # Positive sign denotes tensile stresses, negative denotes compressive

        # factor of safety for each failure mode

        # short aliases
        flags        = self.flags
        yN           = self.yN
        Finternal    = self.Finternal
        strain       = self.strain
        d            = self.d
        theta        = self.theta
        nTube        = self.nTube
        nCap         = self.nCap
        yWire        = self.yWire
        zWire        = self.zWire
        EIxJ         = self.EIxJ
        EIzJ         = self.EIzJ
        lBiscuit     = self.lBiscuit
        dQuad        = self.dQuad
        thetaQuad    = self.thetaQuad
        nTubeQuad    = self.nTubeQuad
        lBiscuitQuad = self.lBiscuitQuad
        RQuad        = self.RQuad
        hQuad        = self.hQuad
        EIQuad       = self.EIQuad[0]  # convert single element array to scalar
        GJQuad       = self.GJQuad[0]  # convert single element array to scalar
        tWire        = self.tWire
        TWire        = self.TWire
        TEtension    = self.TEtension
        b            = self.b
        Fblade       = self.Fblade
        mSpar        = self.mSpar
        mChord       = self.mChord
        mElseRotor   = self.mElseRotor

        TQuad = np.sum(Fblade.Fz)*b - (np.sum(mSpar + mChord)*b + mElseRotor/4) * 9.81

        Ns = max(yN.shape) - 1  # number of elements
        dy = np.zeros((Ns, 1))
        for s in range(0, Ns):
            dy[s] = yN[s] - yN[s]

        fail = Failure()

        # Material failure
        fail.top    = self.material_failure(Ns, strain.top,    theta, nCap, flags)
        fail.bottom = self.material_failure(Ns, strain.bottom, theta, nCap, flags)
        fail.back   = self.material_failure(Ns, strain.back,   theta, [],   flags)
        fail.front  = self.material_failure(Ns, strain.front,  theta, [],   flags)

        # Euler Buckling failure in main spar from wire force
        k  = 0.7    # pinned-pinned = 1, fixed-pinned = 0.7 with correction factor
        #kk = 1.42  # correction factor
        kk = 1      # correction factor was in error, never saw buckling failure
        thetaWire = atan2(zWire, yWire)
        L = yWire   # wire attachment provides pinned end
        F = TWire * cos(thetaWire) + TEtension

        fail.buckling.x = np.zeros(Ns+1)
        fail.buckling.z = np.zeros(Ns+1)

        for s in range(0, Ns):
            if yN[s] <= yWire:
                critical_load_x = pi**2 * EIxJ / (k * L)**2
                critical_load_z = pi**2 * EIzJ / (k * L)**2
                fail.buckling.x[s] = kk * F / critical_load_x
                fail.buckling.z[s] = kk * F / critical_load_z
            else:
                fail.buckling.x[s] = 0
                fail.buckling.z[s] = 0

        # no buckling at tip
        fail.buckling.x[Ns] = 0
        fail.buckling.z[Ns] = 0

        # Torsional Buckling failure
        fail.buckling.torsion = self.torsional_buckling_failure(Ns,
            Finternal, d, theta, nTube, nCap, lBiscuit, flags)

        # Quad Buckling failure
        if EIQuad != 0:
            k = 1
            L = sqrt(RQuad**2 + hQuad**2)
            alpha = atan2(hQuad, RQuad)
            P = TQuad / sin(alpha)
            critical_load = pi**2 * EIQuad / (k * L)**2
            fail.quad_buckling = P / critical_load
        else:
            fail.quad_buckling = 0

        # Quad bending moment failure (does not include torsion since they don't occur at the same time)
        RotorMoment = 1400
        if EIQuad != 0:
            TbottomWire = TQuad / tan(alpha)
            BM = TbottomWire * zWire + RotorMoment
            strainQuad = -np.array([BM * (dQuad / 2) / EIQuad, 0, 0]).reshape(1, -1).T  # strain on compression side
            mf = self.material_failure(1, strainQuad, [thetaQuad], [], flags)
            fail.quad_bend = abs(mf.plus[0, 0])  # use only compressive failure in fibre direction
        else:
            fail.quad_bend = 0

        # Quad torsional material failure
        if GJQuad != 0:
            strainQuad = np.array([0,  0, dQuad / 2 * RotorMoment / GJQuad]).reshape(1,  -1).T
            mf = self.material_failure(1, strainQuad, [thetaQuad], [], flags)
            fail.quad_torsion = abs(mf.plus[0, 0])
        else:
            fail.quad_torsion = 0

        # Quad torsional buckling failure
        FRotorMoment = np.array([0, 0, 0, 0, RotorMoment]).reshape(1,  -1).T
        tbf = self.torsional_buckling_failure(1,
            FRotorMoment, [dQuad], [thetaQuad], [nTubeQuad], [0], [lBiscuitQuad], flags)
        fail.quad_torbuck = tbf[0]

        # Wire tensile failure
        wire_props = wireProperties[flags.WireType]
        fail.wire = np.zeros(len(yWire))
        for i in range(0, len(yWire)):
            stress_wire = TWire[i] / (pi*(tWire/2)**2)
            fail.wire[i] = stress_wire / wire_props['ULTIMATE']

        self.fail = fail

    def material_failure(self,  Ns, strain, theta, nCap, flags):
        failure = MaterialFailure()
        failure.cap = np.zeros((3, Ns+1))
        failure.plus = np.zeros((3, Ns+1))
        failure.minus = np.zeros((3, Ns+1))

        # Material Properties
        tube_props = prepregProperties[flags.CFRPType]

        # Cap Prepreg Properties (MTM28-M46J 140 37 %RW 12")
        cap_props = prepregProperties[flags.CFRPType]

        # Populate Q matrix for tube
        Q_TUBE = np.zeros((3, 3))
        Q_TUBE[0, 0] = tube_props['E_11']
        Q_TUBE[1, 1] = tube_props['E_22']
        Q_TUBE[0, 1] = tube_props['E_22'] * tube_props['V_12']
        Q_TUBE[1, 0] = Q_TUBE[0, 1]
        Q_TUBE[2, 2] = tube_props['G_12']

        # Populate Q matrix for caps
        Q_CAP = np.zeros((3, 3))
        Q_CAP[0, 0] = cap_props['E_11']
        Q_CAP[1, 1] = cap_props['E_22']
        Q_CAP[0, 1] = cap_props['E_22'] * cap_props['V_12']
        Q_CAP[1, 0] = Q_CAP[0, 1]
        Q_CAP[2, 2] = cap_props['G_12']

        stress = MaterialFailure()
        stress.cap = np.zeros((3, Ns))
        stress.plus = np.zeros((3, Ns))
        stress.minus = np.zeros((3, Ns))

        for s in range(0, Ns):
            # Compute stresses in structural axes for each lamina angle
            # Q is the matrix of elastic constants in the material axis
            # Q_bar is the matrix of elastic constants in the structural axes for a
            # lamina at a ply angle theta.

            # Failure is computed at each node, but using the ply angle of the
            # element (this shouldn't cause a large discrepency). Failure at tip is
            # zero, since stresses/strains at tip are zero.

            # Transform the elastic constants for the plus angle ply
            x = theta[s]  # Composite angle (in radians)
            T_PLUS = np.array([
                [ cos(x)**2,     sin(x)**2,      2*sin(x)*cos(x)       ],
                [ sin(x)**2,     cos(x)**2,     -2*sin(x)*cos(x)       ],
                [-sin(x)*cos(x), sin(x)*cos(x), (cos(x)**2)-(sin(x)**2)]
            ])
            # MATLAB version: Qbar_TUBE_PLUS = (T_PLUS\Q_TUBE)/T_PLUS'
            tmp = np.linalg.solve(T_PLUS, Q_TUBE)
            Qbar_TUBE_PLUS = np.linalg.solve(T_PLUS, tmp.T)

            # Transform the elastic constants for the minus angle ply
            x = -theta[s]  # Composite angle (in radians)
            T_MINUS = np.array([
                [ cos(x)**2,     sin(x)**2,      2*sin(x)*cos(x)       ],
                [ sin(x)**2,     cos(x)**2,     -2*sin(x)*cos(x)       ],
                [-sin(x)*cos(x), sin(x)*cos(x), (cos(x)**2)-(sin(x)**2)]
            ])
            # MATLAB version: Qbar_TUBE_MINUS = (T_MINUS\Q_TUBE)/T_MINUS'
            tmp = np.linalg.solve(T_MINUS, Q_TUBE)
            Qbar_TUBE_MINUS = np.linalg.solve(T_MINUS, tmp.T)

            # Compute stresses in structural coordinates
            stress.cap[:, s]   = np.dot(Q_CAP, strain[:, s])            # using Q for the cap
            stress.plus[:, s]  = np.dot(Qbar_TUBE_PLUS, strain[:, s])   # using Q_bar for the + angle plys
            stress.minus[:, s] = np.dot(Qbar_TUBE_MINUS, strain[:, s])  # using Q_bar for the - angle plys

            # Transform stresses to material axes for each lamina angle
            stress.plus[:, s]  = np.dot(T_PLUS, stress.plus[:, s])
            stress.minus[:, s] = np.dot(T_MINUS, stress.minus[:, s])

            # Determine fraction of failure for each lamina angle

            # ULTIMATE_11_TENS and ULTIMATE_11_COMP are both positive values
            # indicating the maximum tensile and compressive stress before failure.
            # failure.cap[0,s] will be positive for tensile failures and negative for
            # compressive failures.

            # Cap failure
            if len(nCap) == 0 or (nCap == 0).all():
                failure.cap[0, s] = 0
                failure.cap[1, s] = 0
                failure.cap[2, s] = 0
            else:
                if stress.cap[0, s] > 0:  # tensile stress in fibre
                    failure.cap[0, s] = stress.cap[0, s] / cap_props['ULTIMATE_11_TENS']
                else:
                    failure.cap[0, s] = stress.cap[0, s] / cap_props['ULTIMATE_11_COMP']

                if stress.cap[1, s] > 0:  # tensile stress in matrix
                    failure.cap[1, s] = stress.cap[1, s] / cap_props['ULTIMATE_22_TENS']
                else:
                    failure.cap[1, s] = stress.cap[1, s] / cap_props['ULTIMATE_22_COMP']

                failure.cap[2, s] = stress.cap[2, s] / cap_props['ULTIMATE_12']

            # Plus angle ply failure
            if stress.plus[0, s] > 0:  # tensile stress in fibre
                failure.plus[0, s] = stress.plus[0, s] / tube_props['ULTIMATE_11_TENS']
            else:
                failure.plus[0, s] = stress.plus[0, s] / tube_props['ULTIMATE_11_COMP']

            if stress.plus[1, s] > 0:  # tensile stress in matrix
                failure.plus[1, s] = stress.plus[1, s] / tube_props['ULTIMATE_22_TENS']
            else:
                failure.plus[1, s] = stress.plus[1, s] / tube_props['ULTIMATE_22_COMP']

            failure.plus[2, s] = stress.plus[2, s] / tube_props['ULTIMATE_12']

            # Minus angle ply failure
            if stress.minus[0, s] > 0:  # tensile stress in fibre
                failure.minus[0, s] = stress.minus[0, s] / tube_props['ULTIMATE_11_TENS']
            else:
                failure.minus[0, s] = stress.minus[0, s] / tube_props['ULTIMATE_11_COMP']

            if stress.minus[1, s] > 0:  # tensile stress in matrix
                failure.minus[1, s] = stress.minus[1, s] / tube_props['ULTIMATE_22_TENS']
            else:
                failure.minus[1, s] = stress.minus[1, s] / tube_props['ULTIMATE_22_COMP']

            failure.minus[2, s] = stress.minus[2, s] / tube_props['ULTIMATE_12']

        return failure

    def torsional_buckling_failure(self, Ns, Finternal, d, theta, nTube, nCap, lBiscuit, flags):
        # Material Properties
        tube_props = prepregProperties[flags.CFRPType]
        V_21_TUBE = tube_props['V_12'] * (tube_props['E_22'] / tube_props['E_11'])

        # Coordinate system: x is axial direction, theta is circumferential direction
        mu_prime_x = tube_props['V_12']
        mu_prime_theta = V_21_TUBE

        Q = np.zeros((3, 3))  # Preallocate Q-matrix

        # Matrix of elastic constants (AER1401, Composite Lamina, slide 8)
        Q[0, 0] = tube_props['E_11'] / (1 - tube_props['V_12'] * V_21_TUBE)
        Q[1, 1] = tube_props['E_22'] / (1 - tube_props['V_12'] * V_21_TUBE)
        Q[0, 1] = tube_props['E_22'] * tube_props['V_12'] / (1 - tube_props['V_12'] * V_21_TUBE)
        Q[1, 0] = Q[0, 1]
        Q[2, 2] = tube_props['G_12']

        failure = np.zeros(Ns+1)

        for s in range(0, Ns):

            if nCap[s] != 0:
                AF_torsional_buckling = 1.25  # See "Validation - Torsional Buckling.xlsx"
            else:
                AF_torsional_buckling = 1

            # Determine elastic properties of rotated tube laminate
            x = theta[s]  # Composite angle (in radians)

            # Transformation matrix
            T = np.array([
                [cos(x)**2,      sin(x)**2,      2*sin(x)*cos(x)],
                [sin(x)**2,      cos(x)**2,     -2*sin(x)*cos(x)],
                [-sin(x)*cos(x), sin(x)*cos(x),  cos(x)**2-sin(x)**2]
            ])

            # Transform the elastic constants using the transformation matrix to obtain the
            # elastic constants at the composite angle.
            # MATLAB version: Qbar = (T\Q)/T'
            tmp = np.linalg.solve(T, Q)
            Qbar = np.linalg.solve(T, tmp.T)

            # Breakout tube elastic constants at the transformed angle
            E_x = Qbar[0, 0]
            E_theta = Qbar[1, 1]

            # Calculate tube geometric properties
            t_tube = nTube[s] * tube_props['T_PLY']  # Shell thickness, CR-912 p.xxxvi
            R = (d[s] + t_tube) / 2                  # Radius from axis of rotation to centroidal surface of cylinder wall, CR-912 p.xxxiv
            L = lBiscuit[s]                          # Unsupported length of cylinder, CR-912 p.xxx

            # Calculate tube elastic properties (CR-912 p.576)
            D_x     = E_x*(1./12)*(t_tube**3)
            D_theta = E_theta*(1./12)*(t_tube**3)
            B_x     = E_x*t_tube
            B_theta = E_theta*t_tube

            # Calculate Gamma
            rho = ((D_x*D_theta)/(B_x*B_theta))**(1./4)
            Gamma = 3.6125e-07  * ((R / (rho * 1000)) ** 6)  + \
                    -1.9724e-05 * ((R / (rho * 1000)) ** 5)  + \
                     0.0004283  * ((R / (rho * 1000)) ** 4)  + \
                    -0.0048315  * ((R / (rho * 1000)) ** 3)  + \
                     0.031801   * ((R / (rho * 1000)) ** 2)  + \
                    -0.12975    *  (R / (rho * 1000))        + \
                     0.88309

            # Calculate factors required in critical torque equation
            Z   = ((B_theta*(1-mu_prime_x*mu_prime_theta)*(L**4)) / (12*D_x*(R**2)))**(1./2)
            Z_s = ((D_theta/D_x)**(5./6))*((B_x/B_theta)**(1./2))*Z
            K_s = 0.89*(Z_s**(3./4))
            N_x_theta = (Gamma*K_s*(pi**2)*D_x)/(L**2)

            # Calculate critical torque
            critical_torque = AF_torsional_buckling * N_x_theta * 2 * pi * (R**2)

            failure[s] = abs(Finternal[4, s] / critical_torque)

        failure[Ns] = 0  # no torsion at tip

        return failure


class Structures(Assembly):
    """
    structural computation, first computes the mass of the helicopter based on
    the structural description of the spars and chord lengths. It then
    computes the deformation of the spars, the strains, and the resulting
    factor of safety for each of the failure modes.
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
    b   = Int(iotype='in', desc='number of blades')
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
    TWire        = Array(iotype='in', desc='')
    TEtension    = Float(iotype='in', desc='')

    # inputs for 'other stuff'
    mElseRotor   = Float(iotype='in', desc='')
    mElseCentre  = Float(iotype='in', desc='')
    mElseR       = Float(iotype='in', desc='')
    R            = Float(iotype='in', desc='rotor radius')
    mPilot       = Float(iotype='in', desc='mass of pilot')

    # inputs for FEM
    Fblade       = VarTree(FBlade(), iotype='in')
    presLoad     = VarTree(PrescribedLoad(), iotype='in')

    def configure(self):
        self.add('spar', SparProperties())
        self.connect('yN',             'spar.yN')
        self.connect('d',              'spar.d')
        self.connect('theta',          'spar.theta')
        self.connect('nTube',          'spar.nTube')
        self.connect('nCap',           'spar.nCap')
        self.connect('lBiscuit',       'spar.lBiscuit')
        self.connect('flags.CFRPType', 'spar.CFRPType')

        self.add('joint', JointSparProperties())
        self.connect('flags.CFRPType', 'joint.CFRPType')
        self.connect('Jprop',          'joint.Jprop')

        self.add('chord', ChordProperties())
        self.connect('yN',             'chord.yN')
        self.connect('cE',             'chord.c')
        self.connect('d',              'chord.d')
        self.connect('flags.GWing',    'chord.GWing')
        self.connect('xtU',            'chord.xtU')

        self.add('quad', QuadSparProperties())
        self.connect('dQuad',          'quad.dQuad')
        self.connect('thetaQuad',      'quad.thetaQuad')
        self.connect('nTubeQuad',      'quad.nTubeQuad')
        self.connect('lBiscuitQuad',   'quad.lBiscuitQuad')
        self.connect('flags.CFRPType', 'quad.CFRPType')
        self.connect('RQuad',          'quad.RQuad')
        self.connect('hQuad',          'quad.hQuad')

        self.add('mass', MassProperties())
        self.connect('flags',          'mass.flags')
        self.connect('b',              'mass.b')
        self.connect('spar.mSpar',     'mass.mSpar')
        self.connect('chord.mChord',   'mass.mChord')
        self.connect('chord.xCGChord', 'mass.xCGChord')
        self.connect('quad.mQuad',     'mass.mQuad')
        self.connect('xEA',            'mass.xEA')
        self.connect('ycmax',          'mass.ycmax')
        self.connect('zWire',          'mass.zWire')
        self.connect('yWire',          'mass.yWire')
        self.connect('tWire',          'mass.tWire')
        self.connect('mElseRotor',     'mass.mElseRotor')
        self.connect('mElseCentre',    'mass.mElseCentre')
        self.connect('mElseR',         'mass.mElseR')
        self.connect('R',              'mass.R')
        self.connect('mPilot',         'mass.mPilot')

        self.add('fem', FEM())
        self.connect('flags',        'fem.flags')
        self.connect('yN',           'fem.yN')
        self.connect('spar.EIx',     'fem.EIx')
        self.connect('spar.EIz',     'fem.EIz')
        self.connect('spar.EA',      'fem.EA')
        self.connect('spar.GJ',      'fem.GJ')
        self.connect('cE',           'fem.cE')
        self.connect('xEA',          'fem.xEA')
        self.connect('spar.mSpar',   'fem.mSpar')
        self.connect('chord.mChord', 'fem.mChord')
        self.connect('mass.xCG',     'fem.xCG')
        self.connect('zWire',        'fem.zWire')
        self.connect('yWire',        'fem.yWire')
        self.connect('TWire',        'fem.TWire')
        self.connect('Fblade',       'fem.Fblade')
        self.connect('presLoad',     'fem.presLoad')

        self.add('strains', Strains())
        self.connect('yN',    'strains.yN')
        self.connect('d',     'strains.d')
        self.connect('fem.k', 'strains.k')
        self.connect('fem.F', 'strains.F')
        self.connect('fem.q', 'strains.q')

        self.add('failure', Failures())
        self.connect('flags',             'failure.flags')
        self.connect('yN',                'failure.yN')
        self.connect('strains.Finternal', 'failure.Finternal')
        self.connect('strains.strain',    'failure.strain')
        self.connect('d',                 'failure.d')
        self.connect('theta',             'failure.theta')
        self.connect('nTube',             'failure.nTube')
        self.connect('nCap',              'failure.nCap')
        self.connect('yWire',             'failure.yWire')
        self.connect('zWire',             'failure.zWire')
        self.connect('joint.EIx',         'failure.EIxJ')
        self.connect('joint.EIz',         'failure.EIzJ')
        self.connect('lBiscuit',          'failure.lBiscuit')
        self.connect('dQuad',             'failure.dQuad')
        self.connect('thetaQuad',         'failure.thetaQuad')
        self.connect('nTubeQuad',         'failure.nTubeQuad')
        self.connect('lBiscuitQuad',      'failure.lBiscuitQuad')
        self.connect('RQuad',             'failure.RQuad')
        self.connect('hQuad',             'failure.hQuad')
        self.connect('quad.EIx',          'failure.EIQuad')  # confirm
        self.connect('quad.GJ',           'failure.GJQuad')  # confirm
        self.connect('tWire',             'failure.tWire')
        self.connect('TWire',             'failure.TWire')
        self.connect('TEtension',         'failure.TEtension')
        self.connect('b',                 'failure.b')
        self.connect('Fblade',            'failure.Fblade')
        self.connect('spar.mSpar',        'failure.mSpar')
        self.connect('chord.mChord',      'failure.mChord')
        self.connect('mElseRotor',        'failure.mElseRotor')

        # link up the outputs
        self.create_passthrough('mass.Mtot')
        self.create_passthrough('fem.q')
        self.create_passthrough('strains.Finternal')
        self.create_passthrough('strains.strain')
        self.create_passthrough('failure.fail')

        self.driver.workflow.add('chord')
        self.driver.workflow.add('fem')
        self.driver.workflow.add('joint')
        self.driver.workflow.add('mass')
        self.driver.workflow.add('quad')
        self.driver.workflow.add('spar')
        self.driver.workflow.add('strains')
        self.driver.workflow.add('failure')
