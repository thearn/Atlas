# pylint: disable=line-too-long, invalid-name, bad-whitespace, trailing-whitespace, too-many-locals, line-too-long
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int, Str, Enum
import numpy as np
from math import pi, sin, cos, floor, sqrt


# Material properties of the CFRP prepregs used in the HPH Project.
# Input Variables:
#   	type_flag - flag indicating prepreg to be used
#           1 - NCT301-1X HS40 G150 33 +/-2#RW
#           2 - HexPly 6376 HTS-12K 268gsm 35#RW
#           3 - MTM28-1/IM7-GP-145gsm 32 +/-2#RW
#           4 - HexPly 8552 IM7 160gsm 35#RW
#           5 - NCT304-1 HR40 G80 40 +/-2#RW
#           6 - MTM28-1B/M46J-140-37#RW
#
# Output Variables:
#
#       The properties contained in this function are:
#       RHO - Prepreg laminate density [Km/m^3]
#       T_PLY - Ply thickness [m]
#       E_11 - Elastic modulus in the material 11-axis [Pa]
#       E_22 - Elastic modulus in the material 22-axis [Pa]
#       G_12 - Shear modulus in the material 12-axis [Pa]
#       ULTIMATE_11_TENS - Failure strength in the material 11-axis, tension
#       ULTIMATE_11_COMP - Failure strength in the material 11-axis, compression
#       ULTIMATE_22_TENS - Failure strength in the material 22-axis, tension
#       ULTIMATE_22_COMP - Failure strength in the material 22-axis, compression
#       ULTIMATE_12 - Failure strength in shear in the material 12-axis
#       V_12 - Poisson's ratio of the material in the 12-axis [N/A]
prepregProperties =  {
    'NCT301-1X HS40 G150 33 +/-2%RW': {
        'RHO': 1.5806E+03,
        'T_PLY': 1.4173E-04,
        'E_11': 2.1066E+11,
        'E_22': 9.0101E+09,
        'G_12': 3.3426E+09,
        'V_12': 0.27,
        'ULTIMATE_11_TENS': 2.0469E+09,
        'ULTIMATE_11_COMP': 8.4593E+08,
        'ULTIMATE_22_TENS': 3.4903E+07,
        'ULTIMATE_22_COMP': 2.0806E+08,
        'ULTIMATE_12': 1.1974E+08,
    },
    'HexPly 6376 HTS-12K 268gsm 35%RW': {
       'RHO': 1.5763E+03,
       'T_PLY': 2.6157E-04,
       'E_11': 1.2329E+11,
       'E_22': 9.2393E+09,
       'G_12': 2.8958E+09,
       'V_12': 0.36,
       'ULTIMATE_11_TENS': 2.0342E+09,
       'ULTIMATE_11_COMP': 1.2689E+09,
       'ULTIMATE_22_TENS': 5.7517E+07,
       'ULTIMATE_22_COMP': 2.0884E+08,
       'ULTIMATE_12': 5.4882E+07,
    },
    'MTM28-1/IM7-GP-145gsm 32 +/-2%Rw': {
       'RHO': 1530.90379,
       'T_PLY': 1.3970E-04,
       'E_11': 1.3607E+11,
       'E_22': 7.4808E+09,
       'G_12': 2.8958E+09,
       'V_12': 0.36,
       'ULTIMATE_11_TENS': 2.2732E+09,
       'ULTIMATE_11_COMP': 1.0378E+09,
       'ULTIMATE_22_TENS': 3.7301E+07,
       'ULTIMATE_22_COMP': 1.8181E+08,
       'ULTIMATE_12': 5.4882E+07,
    },
    'HexPly 8552 IM7 160gsm 35%RW': {
       'RHO': 1575.7808,
       'T_PLY': 1.5669E-04,
       'E_11': 1.3410E+11,
       'E_22': 8.4806E+09,
       'G_12': 4.4816E+09,
       'V_12': 0.356,
       'ULTIMATE_11_TENS': 1.8805E+09,
       'ULTIMATE_11_COMP': 1.5406E+09,
       'ULTIMATE_22_TENS': 5.1021E+07,
       'ULTIMATE_22_COMP': 2.6745E+08,
       'ULTIMATE_12': 8.8598E+07,
    },
    'NCT304-1 HR40 G80 40 +/-2%RW': {
       'RHO': 1508.2873,
       'T_PLY': 9.5250E-05,
       'E_11': 1.9010E+11,
       'E_22': 9.0101E+09,
       'G_12': 3.3426E+09,
       'V_12': 0.27,
       'ULTIMATE_11_TENS': 2.1614E+09,
       'ULTIMATE_11_COMP': 1.0003E+09,
       'ULTIMATE_22_TENS': 3.4903E+07,
       'ULTIMATE_22_COMP': 2.0806E+08,
       'ULTIMATE_12': 1.1974E+08,
    },
    'MTM28-1B/M46J-140-37%RW': {
       'RHO': 1548.7788,
       'T_PLY': 1.4700E-04,
       'E_11': 2.1713E+11,
       'E_22': 6.7341E+09,
       'G_12': 3.3426E+09,
       'V_12': 0.27,
       'ULTIMATE_11_TENS': 1.9436E+09,
       'ULTIMATE_11_COMP': 6.0103E+08,
       'ULTIMATE_22_TENS': 1.6470E+07,
       'ULTIMATE_22_COMP': 1.4617E+08,
       'ULTIMATE_12': 1.0255E+08,
    }
}


class SparProperties(Component):
    """ Computes the structural properties of a CFRP spar given the diameter, d,
        wrap angle, theta, number of tube layers, nTube, and number of cap
        strips, nCap. Properties are computed for each element along the span,
        where the node locations are given by yN.
    """

    yN    = Array(iotype='in', desc='node locations for each element along the span')

    d     = Array(iotype='in', desc='diameter')
    theta = Array(iotype='in', desc='wrap angle')
    nTube = Array(iotype='in', desc='number of tube layers')
    nCap  = Array(iotype='in', desc='number of cap strips')

    lBiscuit = Array(iotype='in', desc='description')
    #CFRPType = Int(0, iotype='in', desc='description')
    CFRPType = Str('', iotype='in', desc='description')

    EIx = Array(iotype='out', desc='description')
    EIz = Array(iotype='out', desc='description')
    EA  = Array(iotype='out', desc='description')
    GJ  = Array(iotype='out', desc='description')
    mSpar = Array(iotype='out', desc='description')

    # yN = Array(np.zeros(2), iotype='in', desc='description')
    # d = Float(0, iotype='in', desc='description')
    # theta = Float(0, iotype='in', desc='description')
    # nTube = Int(0, iotype='in', desc='description')
    # nCap = Array(np.zeros(2), iotype='in', desc='description')
    # lBiscuit = Float(0, iotype='in', desc='description')
    # CFRPType = Int(0, iotype='in', desc='description')

    # EIx = Float(0, iotype='out', desc='description')
    # EIz = Float(0, iotype='out', desc='description')
    # EA = Float(0, iotype='out', desc='description')
    # GJ = Float(0, iotype='out', desc='description')
    # mSpar = Float(0, iotype='out', desc='description')

    def execute(self):
        # material properties for tube
        pp = prepregProperties[self.CFRPType]
        RHO_TUBE = pp['RHO']
        T_PLY_TUBE = pp['T_PLY']
        E_11_TUBE = pp['E_11']
        E_22_TUBE = pp['E_22']
        G_12_TUBE = pp['G_12']
        V_12_TUBE = pp['V_12']

        # material properties for cap
        RHO_CAP = pp['RHO']
        T_PLY_CAP = pp['T_PLY']
        E_11_CAP = pp['E_11']
        G_12_CAP = pp['G_12']

        # other material properties
        RHO_BALSA = 160.0           # Materials Database, Balsa
        RHO_STRUCTURAL_FOAM = 31.0  # Materials Database, Rohacell 31 IG-F

        # spar geometry constants
        ALPHA_BASE = 45 * (pi / 180)  # Half-angle of base cap ply width, based on 90-degree rule
        PLY_TAPER_RATIO = 0.05  # Based on HPO wing, 0.05 for 3" spar, ~0.01 for smaller diameters

        # sizing factors
        AF_biscuit = 1.133
        thickness_biscuit_core = (1.0 / 4) * (0.0254)
        thickness_biscuit_plate = 2 * 2 * (3.0 / 32) * (0.0254)
        biscuit_face_fraction = 0.7

        # determine length of each spar element
        Ns = max(self.yN.shape) - 1
        dy = np.zeros(Ns)
        for s in range(1, (Ns+1)):
            dy[(s-1)] = self.yN[(s)] - self.yN[(s-1)]

        # Pre-allocate EIx, EIz, EA, GJ, mSpar
        self.EIx = np.zeros(Ns)
        self.EIz = np.zeros(Ns)
        self.EA = np.zeros(Ns)
        self.GJ = np.zeros(Ns)
        self.mSpar = np.zeros(Ns)

        # Pre-populate Q matrix for composite angle transformation
        Q = np.zeros((3, 3))  # Preallocate Q-matrix

        # matrix of elastic constants
        Q[0, 0] = E_11_TUBE
        Q[1, 1] = E_22_TUBE
        Q[0, 1] = E_22_TUBE * V_12_TUBE
        Q[1, 0] = Q[0, 1]
        Q[2, 2] = G_12_TUBE

        # compute EIx, EIz, EA, GJ, mSpar for each element
        for s in range(1, (Ns+1)):
            # Determine E_xx_tube and G_xy_tube using Q_Transform Code
            x = self.theta[(s-1)]  # Composite angle (in radians)

            # Transformation matrix
            T = np.array([
                [cos(x)**2,       sin(x)**2,      2*sin(x)*cos(x)],
                [sin(x)**2,       cos(x)**2,     -2*sin(x)*cos(x)],
                [-sin(x)*cos(x),  sin(x)*cos(x), (cos(x)**2)-(sin(x)**2)]
            ])

            # Transform the elastic constants using the transformation matrix
            # to obtain the elastic constants at the composite angle.
            #
            # MATLAB version of this line is
            #    Qbar = (T\Q)/T';
            # backslash is left division
            # linalg.lstsq(x,y)
            # another place says
            #   linalg.solve(a,b)
            # b[np.newaxis].T
            #Qbar = (np.linalg.solve(T,Q)) / T.T # ???
            #
            # A/B is A*inv(B)
            # A\B is inv(A)*B
            #Qbar = (np.linalg.lstsq(T,Q)) / T[np.newaxis].T
            #Qbar = (np.linalg.inv(T) * Q ) *  np.linalg.inv(T[np.newaxis].T)
            Qbar = (np.linalg.inv(T) * Q) *  np.linalg.inv(T.transpose())
            #print T
            #print (np.linalg.inv(T) * Q )
            b = np.linalg.solve(T, Q)
            a = T.T
            QbarT = np.linalg.solve(a.T, b.T)
            #print QbarT.T
            #print T.transpose() * np.linalg.inv( np.linalg.solve(T,Q) )
            #print T.transpose()
            Qbar = QbarT.T

            # Breakout tube elastic constants at the transformed angle
            E_xx_tube = Qbar[0, 0]
            G_xy_tube = Qbar[2, 2]

            # I Tube
            I_tube = pi * (((self.d[(s-1)] / 2) + (1.0 / 2) * self.nTube[(s-1)] * T_PLY_TUBE) ** 3) * (self.nTube[(s-1)] * T_PLY_TUBE)

            # A Tube
            A_tube = pi * ((((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE) ** 2) - ((self.d[(s-1)] / 2) ** 2))

            # Initialize Ix_cap(s) (In-plane), Iz_cap(s) (Out-of-plane), A_cap
            Ix_cap_0 = 0
            Iz_cap_0 = 0
            A_cap_0  = 0
            Ix_cap_1 = 0
            Iz_cap_1 = 0
            A_cap_1  = 0

            # Linearly interpolate between discrete values of nCap
            if self.nCap[(s-1)] < 0:
                self.nCap[(s-1)] = 0
            nCap_0 = int(floor(self.nCap[(s-1)]))
            nCap_1 = nCap_0 + 1

            if nCap_0 !=  0:
                width_ply_0 = np.zeros(nCap_0)
            else:
                width_ply_0 = np.zeros(1)
            width_ply_1 = np.zeros(nCap_1)

            for i in range(1, (nCap_0+1)):
                # width_ply_0  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
                # width_ply_0  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths
                width_ply_0[(i-1)] = self.d[(s-1)] * ALPHA_BASE - (i - 1) * (PLY_TAPER_RATIO * self.d[(s-1)])  # Sets base cap width from 90 degree rule, tapers each subsequent layer based on tube diameter
                alpha_ply = width_ply_0[(i-1)] / (2 * ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE + (i - (1.0 / 2)) * T_PLY_CAP))
                Ix_cap_0 = Ix_cap_0 + 2 * (alpha_ply * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                Iz_cap_0 = Iz_cap_0 + 2 * ((alpha_ply + (1.0 / 2) * sin(2 * alpha_ply)) * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                A_cap_0 = A_cap_0 + 2 * (2 * alpha_ply * ((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) * T_PLY_CAP)

            for i in range(1, (nCap_1+1)):
                # width_ply_1  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
                # width_ply_1  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths
                width_ply_1[(i-1)] = self.d[(s-1)] * ALPHA_BASE - (i - 1) * (PLY_TAPER_RATIO * self.d[(s-1)])
                alpha_ply = width_ply_1[(i-1)] / (2 * ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE + (i - (1.0 / 2)) * T_PLY_CAP))
                Ix_cap_1 = Ix_cap_1 + 2 * (alpha_ply * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                Iz_cap_1 = Iz_cap_1 + 2 * ((alpha_ply + (1.0 / 2) * sin(2 * alpha_ply)) * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                A_cap_1 = A_cap_1 + 2 * (2 * alpha_ply * ((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) * T_PLY_CAP)

            Ix_cap = Ix_cap_0 + (Ix_cap_1 - Ix_cap_0) * (self.nCap[(s-1)] - nCap_0)
            Iz_cap = Iz_cap_0 + (Iz_cap_1 - Iz_cap_0) * (self.nCap[(s-1)] - nCap_0)
            A_cap = A_cap_0 + (A_cap_1 - A_cap_0) * (self.nCap[(s-1)] - nCap_0)

            # GJ spar
            r_tube_avg = (self.d[(s-1)] / 2) + (1.0 / 2) * self.nTube[(s-1)] * T_PLY_TUBE
            r_cap_avg = ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE) + (1.0 / 2) * self.nCap[(s-1)] * T_PLY_CAP
            width_cap_avg = (1.0 / 2) * (np.mean(width_ply_0) + np.mean(width_ply_1))
            t_cap_avg = self.nCap[(s-1)] * T_PLY_CAP  # Do not know if the average is actually necessary here, and this value is close regardless
            r_spar_avg = ((pi * r_tube_avg - width_cap_avg) / (pi * r_tube_avg)) * (r_tube_avg) + ((width_cap_avg) / (pi * r_tube_avg)) * (r_cap_avg)
            self.GJ[(s-1)] = (4 * (pi ** 2) * (r_spar_avg ** 4)) / (((2 * pi * r_tube_avg - 2 * width_cap_avg) / (G_xy_tube * self.nTube[(s-1)] * T_PLY_TUBE)) + ((2 * width_cap_avg) / (G_xy_tube * self.nTube[(s-1)] * T_PLY_TUBE + G_12_CAP * t_cap_avg)))

            # Biscuit mass
            mass_biscuit = (AF_biscuit) * (dy[(s-1)] / self.lBiscuit[(s-1)]) * ((pi * (self.d[(s-1)] / 2) ** 2) * (RHO_BALSA * thickness_biscuit_plate * biscuit_face_fraction + RHO_STRUCTURAL_FOAM * thickness_biscuit_core))
            #print AF_biscuit, dy[(s-1)], self.lBiscuit[(s-1)], self.d[(s-1)], RHO_BALSA, thickness_biscuit_plate, biscuit_face_fraction, RHO_STRUCTURAL_FOAM, thickness_biscuit_core

            # Determine Total Spar Properties
            self.EIx[(s-1)] = E_xx_tube * I_tube + E_11_CAP * Ix_cap
            self.EIz[(s-1)] = E_xx_tube * I_tube + E_11_CAP * Iz_cap
            self.EA[(s-1)] = E_xx_tube * A_tube + E_11_CAP * A_cap
            #print A_tube, RHO_TUBE , A_cap , RHO_CAP, dy[(s-1)] , mass_biscuit
            self.mSpar[(s-1)] = (A_tube * RHO_TUBE + A_cap * RHO_CAP) * dy[(s-1)] + mass_biscuit


class DiscretizedProperties(Component):
    """
    Discretize properties along rotor blade. Y defines the locations at which
    the properties are defined. Properties are linearly interpolated between
    Y locations.
    """

    #def DiscretizeProperties(Ns,ycmax,R,c_,Cl_,Cm_,t_,xtU_,xtL_,xEA_,yWire,d_,theta_,nTube_,nCap_,lBiscuit_):

    Ns = Int(0, iotype='in', desc='description')
    ycmax = Array(iotype='in', desc='description')
    R = Float(0, iotype='in', desc='description')
    c_ = Array(iotype='in', desc='description')
    Cl_ = Array(iotype='in', desc='description')
    Cm_ = Array(iotype='in', desc='description')
    t_ = Array(iotype='in', desc='description')
    xtU_ = Array(iotype='in', desc='description')
    xtL_ = Array(iotype='in', desc='description')
    xEA_ = Array(iotype='in', desc='description')
    yWire = Array(iotype='in', desc='description')
    d_ = Array(iotype='in', desc='description')
    theta_ = Array(iotype='in', desc='description')
    nTube_ = Array(iotype='in', desc='description')
    nCap_ = Array(iotype='in', desc='description')
    lBiscuit_ = Array(iotype='in', desc='description')

    #return cE,cN,c100,Cl,Cm,t,xtU,xtL,xEA,d,theta,nTube,nCap,lBiscuit,yN,yE
    cE = Array(iotype='out', desc='description')
    cN = Array(iotype='out', desc='description')
    c100 = Array(np.zeros(100), iotype='out', desc='description')
    Cl = Array(iotype='out', desc='description')
    Cm = Array(iotype='out', desc='description')
    t = Array(iotype='out', desc='description')
    xtU = Array(iotype='out', desc='description')
    xtL = Array(iotype='out', desc='description')
    xEA = Array(iotype='out', desc='description')
    d = Array(iotype='out', desc='description')
    theta = Array(iotype='out', desc='description')
    nTube = Array(iotype='out', desc='description')
    nCap = Array(iotype='out', desc='description')
    lBiscuit = Array(iotype='out', desc='description')
    yN = Array(iotype='out', desc='description')
    yE = Array(iotype='out', desc='description')

    def execute(self):

        self.yN = np.zeros(self.Ns + 1)
        self.cN = np.zeros(self.Ns + 1)
        sTrans = np.zeros(self.Ns + 1)
        self.yE = np.zeros(self.Ns)
        cR = np.zeros(self.Ns)
        cZ = np.zeros(self.Ns)
        self.cE = np.zeros(self.Ns)
        self.Cl = np.zeros(self.Ns)
        self.Cm = np.zeros(self.Ns)
        self.t = np.zeros(self.Ns)
        self.xtU = np.zeros(self.Ns)
        self.xtL = np.zeros(self.Ns)
        self.xEA = np.zeros(self.Ns)
        self.d = np.zeros(self.Ns)
        self.theta = np.zeros(self.Ns)
        self.nTube = np.zeros(self.Ns)
        self.nCap = np.zeros(self.Ns)
        self.lBiscuit = np.zeros(self.Ns)
        Ns100 = 100
        yN100 = np.zeros(Ns100 + 1)
        yE100 = np.zeros(Ns100)
        cR100 = np.zeros(100)
        cZ100 = np.zeros(100)
        self.c100 = np.zeros(100)

        # compute node locations
        for s in range(1, (self.Ns + 1+1)):
            self.yN[(s-1)] = self.R / self.Ns * (s - 1)
            if self.yN[(s-1)] < self.ycmax[0]:
                sTrans[0] = s
            if self.yN[(s-1)] < self.ycmax[1]:
                sTrans[1] = s

        # compute element locations
        for s in range(1, (self.Ns + 1)):
            self.yE[(s-1)] = 0.5 * (self.yN[(s-1)] + self.yN[(s + 1-1)])

        # compute chord lengths at elements
        for s in range(1, (self.Ns + 1)):
            if s < sTrans[0]:  # root section
                x = self.yE[(s-1)] / self.ycmax[0]
                self.cE[(s-1)] = self.c_[0] + x * (self.c_[1] - self.c_[0])
            else:  # chord section
                # compute curve component
                x = (self.yE[(s-1)] - self.ycmax[0]) / (self.R - self.ycmax[0])
                pStart = self.c_[2]
                pCurve = self.c_[3]
                pEnd = self.c_[4]
                xx = x * (1 - pCurve) + sin(x * pi / 2) * pCurve
                cZ[(s-1)] = pStart + (pEnd - pStart) * xx

                # compute 1/r component
                c3 = self.c_[2] / (self.c_[4] * self.R / self.ycmax[0])
                cR[(s-1)] = self.c_[4] * self.R / self.yE[(s-1)] * (c3 + (1 - c3) * x)

                # average based on c_(2)
                self.cE[(s-1)] = cR[(s-1)] + (cZ[(s-1)] - cR[(s-1)]) * self.c_[1]

            if self.cE[(s-1)] == 0:
                self.cE[(s-1)] = 0.001

        # compute chord for display purposes
        for s in range(1, (Ns100 + 1 + 1)):
            yN100[(s-1)] = self.R / Ns100 * (s - 1)
            if all(yN100[(s-1)] < self.ycmax):
                sTrans100 = s

        for s in range(1, (Ns100+1)):
            yE100[(s-1)] = 0.5 * (yN100[(s-1)] + yN100[(s + 1-1)])

        for s in range(1, (Ns100+1)):
            if s < sTrans100:  # root section
                x = yE100[(s-1)] / self.ycmax[0]
                self.c100[(s-1)] = self.c_[0] + x * (self.c_[1] - self.c_[0])
            else:  # chord section
                # compute chord component
                x = (yE100[(s-1)] - self.ycmax[0]) / (self.R - self.ycmax[0])
                pStart = self.c_[2]
                pCurve = self.c_[3]
                pEnd = self.c_[4]
                xx = x * (1 - pCurve) + sin(x * pi / 2) * pCurve
                cZ100[(s-1)] = pStart + (pEnd - pStart) * xx

                # compute 1/r component
                c3 = self.c_[2] / (self.c_[4] * self.R / self.ycmax[0])
                cR100[(s-1)] = self.c_[4] * self.R / yE100[(s-1)] * (c3 + (1 - c3) * x)

                # average based on c_(2)
                self.c100[(s-1)] = cR100[(s-1)] + (cZ100[(s-1)] - cR100[(s-1)]) * self.c_[1]

            if self.c100[(s-1)] == 0:
                self.c100[(s-1)] = 0.001

        # Compute aero properties for each element
        # Y = [ycmax(1) ycmax(2) R];
        Y = np.array([self.ycmax[0], self.ycmax[1], self.R])
        for s in range(1, (self.Ns+1)):
            if s < sTrans[0]:  # root section
                self.Cl[(s-1)]  = self.Cl_[0]
                self.Cm[(s-1)]  = self.Cm_[0]
                self.t[(s-1)]   = self.t_[0]
                self.xEA[(s-1)] = self.xEA_[0]
            else:
                # check which segment the element is on
                for j in range(1, (max(Y.shape)+1)):
                    if Y[(j-1)] > self.yE[(s-1)]:
                        break
                if s == sTrans[0]:
                    j = 2
                x = (self.yE[(s-1)] - Y[(j - 1-1)]) / (Y[(j-1)] - Y[(j - 1-1)])

                # linearly interpolate between Y(j) and Y(j-1)
                self.Cl[(s-1)]  = self.Cl_[(j - 1-1)] + x * (self.Cl_[(j-1)] - self.Cl_[(j - 1-1)])
                self.Cm[(s-1)]  = self.Cm_[(j - 1-1)] + x * (self.Cm_[(j-1)] - self.Cm_[(j - 1-1)])
                self.t[(s-1)]   = self.t_[(j - 1-1)] + x * (self.t_[(j-1)] - self.t_[(j - 1-1)])
                self.xEA[(s-1)] = self.xEA_[(j - 1-1)] + x * (self.xEA_[(j-1)] - self.xEA_[(j - 1-1)])

        self.Cl[(self.Ns-1)] = self.Cl[(self.Ns-1)] * 2 / 3

        # compute xtU and xtL for each element
        # changes instantly from xtU(1) to xtU(3) at point xtU(2)
        for s in range(1, (self.Ns + 1+1)):
            if self.yN[(s-1)] < self.ycmax[0]:
                sTrans[0] = s
            if self.yN[(s-1)] < self.xtU_[1]:
                sTrans[1] = s

        for s in range(1, (self.Ns + 1)):
            if s < sTrans[0]:  # root section
                self.xtU[(s-1)] = 0.05
                self.xtL[(s-1)] = 0.05
            else:
                if s < sTrans[1]:
                    self.xtU[(s-1)] = self.xtU_[0]
                    self.xtL[(s-1)] = self.xtL_[0]
                elif s == sTrans[1]:
                    x = (self.xtU_[1] - self.yN[(s-1)]) / (self.yN[(s-1)] - self.yN[(s - 1-1)])
                    self.xtU[(s-1)] = (1 - x) * self.xtU_[0] + x * self.xtU_[2]
                    self.xtL[(s-1)] = (1 - x) * self.xtL_[0] + x * self.xtL_[2]
                elif s > sTrans[1]:
                    self.xtU[(s-1)] = self.xtU_[2]
                    self.xtL[(s-1)] = self.xtL_[2]

        # compute str properties for each element
        Y = np.array([0, self.yWire[0], self.R])
        for s in range(1, (self.Ns+1)):
            # check which segment the lement is on
            for j in range(1, (max(Y.shape) + 1)):
                if Y[(j-1)] > self.yE[(s-1)]:
                    break
            x = (self.yE[(s-1)] - Y[(j - 1-1)]) / (Y[(j-1)] - Y[(j - 1-1)])

            # linearly interpolate between Y(j) and Y(j-1)
            self.d[(s-1)] = self.d_[(j - 1-1)] + x * (self.d_[(j-1)] - self.d_[(j - 1-1)])
            self.theta[(s-1)] = self.theta_[(j - 1-1)] + x * (self.theta_[(j-1)] - self.theta_[(j - 1-1)])
            self.nTube[(s-1)] = self.nTube_[(j - 1-1)] + x * (self.nTube_[(j-1)] - self.nTube_[(j - 1-1)])
            self.nCap[(s-1)] = self.nCap[(j - 1-1)] + x * (self.nCap[(j-1)] - self.nCap[(j - 1-1)])
            self.lBiscuit[(s-1)] = self.lBiscuit_[(j - 1-1)] + x * (self.lBiscuit_[(j-1)] - self.lBiscuit_[(j - 1-1)])


class WireProperties(Component):
    """
    Computes the material properties of wire options used in the HPH Project.
    """

    material = Enum('Pianowire', ('Pianowire', 'Vectran'), iotype='in',
        desc='Material to be used')

    RHO = Float(0, iotype='out',
        desc='Material density, calculated for woven material if necessary [Km/m^3]')

    E = Float(0, iotype='out', desc='Elastic modulus [Pa]')

    ULTIMATE = Float(0, iotype='out', desc='Failure strength [Pa]')

    D = Array(iotype='out', desc='Available diameters [m]')

    def execute(self):
        if self.material == 'Pianowire':  # High-Tensile Steel
            self.RHO = 7.85e3       # From ASTM228 Standard, Accessed Online at MatWeb
            self.E = 2.10e11        # From ASTM228 Standard, Accessed Online at MatWeb
            self.ULTIMATE = 2.62e9  # FortePiano.com
            self.D = np.array([0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.72]).reshape(1, -1)  # FortePiano.com
        elif self.material == 'Vectran':
            self.RHO = 1.1065e3
            self.E = 3.921e10        # Effective E of braided Vectran, averaged from data supplied by Cortland Cable
            self.ULTIMATE = 9.828e8  # Based on minimum tensile strength supplied by Cortland Cable
            self.D = np.array([1.5, 2]).reshape(1, -1)


class ChordProperties(Component):
    """
    Computes the mass of the ribs, trailing edge, LE-sheeting and covering
    given the chord. The mass is computed for each spanwise element. If
    flagTESpar is 1 then the mass of a trailing edge spar and in-plane truss
    is added as well.
    """

    yN = Array(iotype='in', desc='description')
    c = Array(iotype='in', desc='description')
    d = Array(iotype='in', desc='description')
    flagGWing = Int(iotype='in', desc='description')
    xtU = Array(iotype='in', desc='description')

    mChord = Array(iotype='out', desc='description')
    xCGChord = Array(iotype='out', desc='description')

    def execute(self):

        # Determine the length of each spar element
        Ns = max(self.yN.shape) - 1
        dy = np.zeros(Ns)
        for s in range(1, (Ns+1)):
            dy[(s-1)] = self.yN[(s)] - self.yN[(s-1)]  # length of each element
        self.mChord = np.zeros(Ns)
        self.xCGChord = np.zeros(Ns)

        # Material properties
        RHO_BALSA = 160.0
        RHO_BASSWOOD = 387.0
        RHO_CARBON = 1610.6
        T_PLY_CARBON = 0.000147
        RHO_EPS = 16.0
        RHO_KEVLAR = 1300.0
        T_PLY_KEVLAR = 0.000127
        RHO_MYLAR = 1400.0
        T_MYLAR = 1.2e-05
        RHO_STRUCTURAL_FOAM = 31.0
        RHO_XPS = 24.67
        RHO_STEEL_WIRE = 7850.0

        # Airfoil Geometry Properties (DT1068156, HPO root airfoil)
        AREA_AIRFOIL = 0.08233
        PERIMETER_AIRFOIL = 2.06814
        XCG_AIRFOIL = 0.39003

        # Wing Design Variables
        thickness_rib = 0.005
        rib_spacing = 12 * (0.0254)
        thickness_LE_sheeting = 0.003
        thickness_rib_caps = (1.0 / 32) * (0.0254)
        percent_rib_caps = 0.97
        thickness_spar_plate = (1.0 / 32) * (0.0254)
        percent_radius_spar_plate = 0.0632
        percent_LE_sheeting_top = 0.6
        percent_LE_sheeting_top_GWing = self.xtU
        percent_LE_sheeting_bottom = 0.1
        d_TE_spar = (3.0 / 4) * (0.0254)
        nTube_TE_spar = 4
        d_comp_member = (3.0 / 4) * (0.0254)
        nTube_comp_member = 4
        spacing_comp_member = 3
        percent_length_comp_member = 0.782 - 0.25
        d_cross_bracing = 0.001
        height_TE_foam = 0.007
        length_TE_foam = 0.0505
        percent_height_TE_plate = 0.05
        percent_length_TE_plate = 0.24
        thickness_TE_plate = (1.0 / 32) * (0.0254)
        diameter_piano_wire = (0.022) * (0.0254)
        spar_location = 0.25

        # Adjustment Factors
        AF_ribs = 1.569
        AF_TE = 1.602
        AF_leading_edge_sheeting = 1.302
        AF_leading_edge_sheeting_GWing = 1.21
        AF_covering = 0.972
        AF_TE_spar = 0.529
        AF_comp_member = 1.0
        AF_cross_bracing = 4.603

        # Compute mChord for each element
        for s in range(1, (Ns+1)):
            # Rib Mass & Xcg
            mass_rib_foam = RHO_EPS * (thickness_rib * ((self.c[(s-1)] ** 2) * AREA_AIRFOIL))
            mass_rib_caps = RHO_BASSWOOD * (thickness_rib * thickness_rib_caps * (percent_rib_caps * PERIMETER_AIRFOIL * self.c[(s-1)]))
            mass_rib_plate_spar = RHO_BALSA * (thickness_spar_plate * (pi * (((self.c[(s-1)] * percent_radius_spar_plate) ** 2) - ((self.d[(s-1)] / 2) ** 2))))
            if self.flagGWing == 0:
                mass_rib_plate_TE = RHO_BALSA * (thickness_TE_plate * ((1.0 / 2) * (self.c[(s-1)] ** 2) * percent_height_TE_plate * percent_length_TE_plate))
                mass_rib = AF_ribs * ((dy[(s-1)] / rib_spacing) * (mass_rib_foam + mass_rib_caps + 2 * mass_rib_plate_spar + 2 * mass_rib_plate_TE))
                Xcg_rib = XCG_AIRFOIL
            else:
                mass_rib = AF_ribs * (0.66) * ((dy[(s-1)] / rib_spacing) * (mass_rib_foam + mass_rib_caps + 2 * mass_rib_plate_spar))
                Xcg_rib = ((XCG_AIRFOIL + (2.0 / 3) * (spar_location)) / 2)

            # Trailing Edge Mass & Xcg
            if self.flagGWing == 0:
                mass_TE = AF_TE * (RHO_STRUCTURAL_FOAM * dy[(s-1)] * ((1.0 / 2) * (length_TE_foam * height_TE_foam)) + RHO_KEVLAR * dy[(s-1)] * T_PLY_KEVLAR * (length_TE_foam + height_TE_foam + sqrt(length_TE_foam ** 2 + height_TE_foam ** 2)))
                Xcg_TE = (self.c[(s-1)] - (2.0 / 3) * length_TE_foam) / self.c[(s-1)]
            else:
                mass_TE = RHO_STEEL_WIRE * dy[(s-1)] * pi * ((diameter_piano_wire / 2) ** 2)
                Xcg_TE = 1

            # Leading Edge Sheeting Mass
            if self.flagGWing == 0:
                mass_LE_sheeting = AF_leading_edge_sheeting * (RHO_XPS * dy[(s-1)] * self.c[(s-1)] * (percent_LE_sheeting_top + percent_LE_sheeting_bottom) * PERIMETER_AIRFOIL * thickness_LE_sheeting)
                Xcg_LE_sheeting = (1.0 / 2) * ((1.0 / 2) * (percent_LE_sheeting_top) + (1.0 / 2) * (percent_LE_sheeting_bottom))
            else:
                mass_LE_sheeting = AF_leading_edge_sheeting_GWing * (RHO_EPS * dy[(s-1)] * self.c[(s-1)] * (percent_LE_sheeting_top_GWing[(s-1)] + percent_LE_sheeting_bottom) * PERIMETER_AIRFOIL * thickness_LE_sheeting) + RHO_STEEL_WIRE * dy[(s-1)] * pi * ((diameter_piano_wire / 2) ** 2)
                Xcg_LE_sheeting = (1.0 / 2) * ((1.0 / 2) * (percent_LE_sheeting_top_GWing[(s-1)]) + (1.0 / 2) * (percent_LE_sheeting_bottom))

            # Covering Mass (Mylar, use perimeter estimates from HPO airfoils)
            mass_covering = AF_covering * (RHO_MYLAR * dy[(s-1)] * self.c[(s-1)] * PERIMETER_AIRFOIL * T_MYLAR)
            Xcg_covering = 0.5

            # Trailing Edge Spar and In-Plane Truss Mass & Xcg (TE Spar, Kevlar cross-bracing, compression members)
            if self.flagGWing == 0:
                mass_TE_spar = AF_TE_spar * (RHO_CARBON * dy[(s-1)] * (pi * (((d_TE_spar / 2) + nTube_TE_spar * T_PLY_CARBON) ** 2 - (d_TE_spar / 2) ** 2)))
                Xcg_TE_spar = 0.9
                mass_comp_member = AF_comp_member * (dy[(s-1)] / spacing_comp_member) * RHO_CARBON * (self.c[(s-1)] * percent_length_comp_member) * (pi * (((d_comp_member / 2) + nTube_comp_member * T_PLY_CARBON) ** 2 - (d_comp_member / 2) ** 2))
                Xcg_comp_member = 0.25 + (1.0 / 2) * (percent_length_comp_member)
                mass_cross_bracing = AF_cross_bracing * RHO_KEVLAR * (dy[(s-1)] / spacing_comp_member) * 2 * (pi * (d_cross_bracing / 2) ** 2) * ((((self.c[(s-1)] * percent_length_comp_member) ** 2) + ((spacing_comp_member) ** 2)) ** (1.0 / 2))
                Xcg_cross_bracing = 0.25 + (1.0 / 2) * (percent_length_comp_member)
            else:
                mass_TE_spar = 0
                Xcg_TE_spar = 0
                mass_comp_member = 0
                Xcg_comp_member = 0
                mass_cross_bracing = 0
                Xcg_cross_bracing = 0

            # Total Chord Mass
            self.mChord[(s-1)] = mass_rib + mass_TE + (mass_LE_sheeting) + (mass_covering) + (mass_TE_spar + mass_comp_member + mass_cross_bracing)
            self.xCGChord[(s-1)] = ((mass_rib * Xcg_rib) + (mass_TE * Xcg_TE) + (mass_LE_sheeting * Xcg_LE_sheeting) + (mass_covering * Xcg_covering) + (mass_TE_spar * Xcg_TE_spar) + (mass_comp_member * Xcg_comp_member) + (mass_cross_bracing * Xcg_cross_bracing)) / (mass_rib + mass_TE + mass_LE_sheeting + mass_covering + mass_TE_spar + mass_comp_member + mass_cross_bracing)
