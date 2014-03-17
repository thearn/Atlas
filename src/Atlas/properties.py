# pylint: disable=line-too-long, invalid-name, bad-whitespace, trailing-whitespace, too-many-locals, line-too-long
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int, Str
import numpy as np
import math

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
#       ULTIMATE_11_COMP - Failure strength in the material 11-axis,
#           compression
#       ULTIMATE_22_TENS - Failure strength in the material 22-axis, tension
#       ULTIMATE_22_COMP - Failure strength in the material 22-axis,
#           compression
#       ULTIMATE_12 - Failure strength in shear in the material
#           12-axis
#       V_12 - Poisson's ratio of the material in the 12-axis [N/A]
prepregProperties = \
                  {'NCT301-1X HS40 G150 33 +/-2%RW' :
                   {
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
                   'HexPly 6376 HTS-12K 268gsm 35%RW' :
                   {
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
                   'MTM28-1/IM7-GP-145gsm 32 +/-2%Rw' :
                   {
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
                   'HexPly 8552 IM7 160gsm 35%RW' :
                   {
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
                   'NCT304-1 HR40 G80 40 +/-2%RW' :
                   {
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
                   'MTM28-1B/M46J-140-37%RW' :
                   {
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









class sparProperties(Component):
    yN = Array( iotype='in', desc='description')
    d = Array( iotype='in', desc='description')
    theta = Array( iotype='in', desc='description')
    nTube = Array( iotype='in', desc='description')
    nCap = Array( iotype='in', desc='description')
    lBiscuit = Array(iotype='in', desc='description')
    #CFRPType = Int(0, iotype='in', desc='description')
    CFRPType = Str('', iotype='in', desc='description')
    
    EIx = Array(iotype='out', desc='description')
    EIz = Array(iotype='out', desc='description')
    EA = Array(iotype='out', desc='description')
    GJ = Array(iotype='out', desc='description')
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
        pp = prepregProperties[self.CFRPType]
        RHO_TUBE = pp['RHO']
        T_PLY_TUBE = pp['T_PLY']
        E_11_TUBE = pp['E_11']
        E_22_TUBE = pp['E_22']
        G_12_TUBE = pp['G_12']
        V_12_TUBE = pp['V_12']
        ULTIMATE_11_TENS = pp['ULTIMATE_11_TENS']
        ULTIMATE_11_COMP = pp['ULTIMATE_11_COMP']
        ULTIMATE_22_TENS = pp['ULTIMATE_22_TENS']
        ULTIMATE_22_COMP = pp['ULTIMATE_22_COMP']
        ULTIMATE_12 = pp['ULTIMATE_12']
        
        RHO_CAP = pp['RHO']
        T_PLY_CAP = pp['T_PLY']
        E_11_CAP = pp['E_11']
        E_22_CAP = pp['E_22']
        G_12_CAP = pp['G_12']
        V_12_CAP = pp['V_12']

        RHO_BALSA = 160.0 # Materials Database, Balsa
        RHO_STRUCTURAL_FOAM = 31.0 # Materials Database, Rohacell 31 IG-F
        ALPHA_BASE = 45 * (math.pi / 180) # Half-angle of base cap ply width, based on 90-degree rule
        PLY_TAPER_RATIO = 0.05
        AF_biscuit = 1.133
        thickness_biscuit_core = (1.0 / 4) * (0.0254)
        thickness_biscuit_plate = 2 * 2 * (3.0 / 32) * (0.0254)
        biscuit_face_fraction = 0.7
        Ns = max(self.yN.shape) - 1
        dy = np.zeros(Ns)
        for s in range(1,(Ns+1)):
            dy[(s-1)] = self.yN[(s)] - self.yN[(s-1)]

        # Pre-allocate EIx, EIz, EA, GJ, mSpar
        self.EIx = np.zeros(Ns)
        self.EIz = np.zeros(Ns)
        self.EA = np.zeros(Ns)
        self.GJ = np.zeros(Ns)
        self.mSpar = np.zeros(Ns)
        
        # Pre-populate Q matrix for composite angle transformation
        Q = np.zeros((3,3))  # Preallocate Q-matrix
    
        Q[0,0] = E_11_TUBE
        Q[1,1] = E_22_TUBE
        Q[0,1] = E_22_TUBE * V_12_TUBE
        Q[1,0] = Q[0,1]
        Q[2,2] = G_12_TUBE
        for s in range(1,(Ns+1)):
            # Determine E_xx_tube and G_xy_tube using Q_Transform Code
            x = self.theta[(s-1)] # Composite angle (in radians)

            # Transformation matrix
            T = np.array([
                            [
                              math.cos(x) ** 2,
                              math.sin(x) ** 2,
                              2 * math.sin(x) * math.cos(x),
                            ],
                            [
                              math.sin(x) ** 2,
                              math.cos(x) ** 2,
                              - 2 * math.sin(x) * math.cos(x),
                            ],
                            [
                              - math.sin(x) * math.cos(x),
                              math.sin(x) * math.cos(x),
                              (math.cos(x) ** 2) - (math.sin(x) ** 2)
                            ]
                          ]
                         )
            # Transform the elastic constants umath.sing the transformation matrix to obtain the
            #   elastic constants at the composite angle.
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
            Qbar = (np.linalg.inv(T) * Q ) *  np.linalg.inv(T.transpose())

            #print T
            #print (np.linalg.inv(T) * Q )

            b =  np.linalg.solve(T,Q)
            a = T.T
            QbarT = np.linalg.solve(a.T,b.T)
            print QbarT.T
            #print T.transpose() * np.linalg.inv( np.linalg.solve(T,Q) )
            #print T.transpose()
            Qbar = QbarT.T
            
            # Breakout tube elastic constants at the transformed angle
            E_xx_tube = Qbar[0,0]
            G_xy_tube = Qbar[2,2]

            # I Tube
            I_tube = math.pi * (((self.d[(s-1)] / 2) + (1.0 / 2) * self.nTube[(s-1)] * T_PLY_TUBE) ** 3) * (self.nTube[(s-1)] * T_PLY_TUBE)

            # A Tube
            A_tube = math.pi * ((((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE) ** 2) - ((self.d[(s-1)] / 2) ** 2))

            # Initialize Ix_cap(s) (In-plane), Iz_cap(s) (Out-of-plane), A_cap
            Ix_cap_0 = 0
            Iz_cap_0 = 0
            A_cap_0 = 0
            Ix_cap_1 = 0
            Iz_cap_1 = 0
            A_cap_1 = 0
            
            if self.nCap[(s-1)] < 0:
                self.nCap[(s-1)] = 0
            nCap_0 = int( math.floor(self.nCap[(s-1)]) )
            nCap_1 = nCap_0 + 1
            
            if nCap_0 !=  0:
                width_ply_0 = np.zeros(nCap_0)
            else:
                width_ply_0 = np.zeros(1)
            width_ply_1 = np.zeros(nCap_1)
            
            for i in range(1,(nCap_0+1)):
                # width_ply_0  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
                # width_ply_0  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths
                width_ply_0[(i-1)] = self.d[(s-1)] * ALPHA_BASE - (i - 1) * (PLY_TAPER_RATIO * self.d[(s-1)]) # Sets base cap width from 90 degree rule, tapers each subsequent layer based on tube diameter
                alpha_ply = width_ply_0[(i-1)] / (2 * ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE + (i - (1.0 / 2)) * T_PLY_CAP))
                Ix_cap_0 = Ix_cap_0 + 2 * (alpha_ply * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                Iz_cap_0 = Iz_cap_0 + 2 * ((alpha_ply + (1.0 / 2) * math.sin(2 * alpha_ply)) * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                A_cap_0 = A_cap_0 + 2 * (2 * alpha_ply * ((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) * T_PLY_CAP)
                
            for i in range(1,(nCap_1+1)):
                # width_ply_1  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
                # width_ply_1  =  [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths    
                width_ply_1[(i-1)] = self.d[(s-1)] * ALPHA_BASE - (i - 1) * (PLY_TAPER_RATIO * self.d[(s-1)])
                alpha_ply = width_ply_1[(i-1)] / (2 * ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE + (i - (1.0 / 2)) * T_PLY_CAP))
                Ix_cap_1 = Ix_cap_1 + 2 * (alpha_ply * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                Iz_cap_1 = Iz_cap_1 + 2 * ((alpha_ply + (1.0 / 2) * math.sin(2 * alpha_ply)) * (((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) ** 3) * T_PLY_CAP)
                A_cap_1 = A_cap_1 + 2 * (2 * alpha_ply * ((self.d[(s-1)] / 2) + (i - (1.0 / 2)) * (T_PLY_CAP)) * T_PLY_CAP)
                
            Ix_cap = Ix_cap_0 + (Ix_cap_1 - Ix_cap_0) * (self.nCap[(s-1)] - nCap_0)
            Iz_cap = Iz_cap_0 + (Iz_cap_1 - Iz_cap_0) * (self.nCap[(s-1)] - nCap_0)
            A_cap = A_cap_0 + (A_cap_1 - A_cap_0) * (self.nCap[(s-1)] - nCap_0)
            r_tube_avg = (self.d[(s-1)] / 2) + (1.0 / 2) * self.nTube[(s-1)] * T_PLY_TUBE
            r_cap_avg = ((self.d[(s-1)] / 2) + self.nTube[(s-1)] * T_PLY_TUBE) + (1.0 / 2) * self.nCap[(s-1)] * T_PLY_CAP
            width_cap_avg = (1.0 / 2) * (np.mean(width_ply_0) + np.mean(width_ply_1))
            t_cap_avg = self.nCap[(s-1)] * T_PLY_CAP # Do not know if the average is actually necessary here, and this value is close regardless
            r_spar_avg = ((math.pi * r_tube_avg - width_cap_avg) / (math.pi * r_tube_avg)) * (r_tube_avg) + ((width_cap_avg) / (math.pi * r_tube_avg)) * (r_cap_avg)
            self.GJ[(s-1)] = (4 * (math.pi ** 2) * (r_spar_avg ** 4)) / (((2 * math.pi * r_tube_avg - 2 * width_cap_avg) / (G_xy_tube * self.nTube[(s-1)] * T_PLY_TUBE)) + ((2 * width_cap_avg) / (G_xy_tube * self.nTube[(s-1)] * T_PLY_TUBE + G_12_CAP * t_cap_avg)))

            # Biscuit mass
            mass_biscuit = (AF_biscuit) * (dy[(s-1)] / self.lBiscuit[(s-1)]) * ((math.pi * (self.d[(s-1)] / 2) ** 2) * (RHO_BALSA * thickness_biscuit_plate * biscuit_face_fraction + RHO_STRUCTURAL_FOAM * thickness_biscuit_core))
            print AF_biscuit, dy[(s-1)], self.lBiscuit[(s-1)], self.d[(s-1)], RHO_BALSA, thickness_biscuit_plate, biscuit_face_fraction, RHO_STRUCTURAL_FOAM, thickness_biscuit_core
            # Determine Total Spar Properties
            
            self.EIx[(s-1)] = E_xx_tube * I_tube + E_11_CAP * Ix_cap
            self.EIz[(s-1)] = E_xx_tube * I_tube + E_11_CAP * Iz_cap
            self.EA[(s-1)] = E_xx_tube * A_tube + E_11_CAP * A_cap

            #print A_tube, RHO_TUBE , A_cap , RHO_CAP, dy[(s-1)] , mass_biscuit
            self.mSpar[(s-1)] = (A_tube * RHO_TUBE + A_cap * RHO_CAP) * dy[(s-1)] + mass_biscuit
            
# class sparProperties(Component):
#     yN = Array(np.zeros(2), iotype='in', desc='description')
#     d = Float(0, iotype='in', desc='description')
#     theta = Float(0, iotype='in', desc='description')
#     nTube = Int(0, iotype='in', desc='description')
#     nCap = Array(np.zeros(2), iotype='in', desc='description')
#     lBiscuit = Float(0, iotype='in', desc='description')
#     CFRPType = Int(0, iotype='in', desc='description')

#     EIx = Float(0, iotype='in', desc='description')
#     EIz = Float(0, iotype='in', desc='description')
#     EA = Float(0, iotype='in', desc='description')
#     GJ = Float(0, iotype='in', desc='description')
#     mSpar = Float(0, iotype='in', desc='description')

#     def execute(self):
#         pass

class discretizeProperties(Component):
    """
    Computes discretized Properties
    """
    Ns = Int(0, iotype='in', desc='description')
    ycmax = Array(np.zeros(2), iotype='in', desc='description')
    R = Float(0, iotype='in', desc='description')
    c_ = Array(np.zeros(2), iotype='in', desc='description')
    Cl_ = Array(np.zeros(2), iotype='in', desc='description')
    Cm_ = Array(np.zeros(2), iotype='in', desc='description')
    t_ = Array(np.zeros(2), iotype='in', desc='description')
    xtU_ = Array(np.zeros(2), iotype='in', desc='description')
    xtL_ = Array(np.zeros(2), iotype='in', desc='description')
    xEA_ = Array(np.zeros(2), iotype='in', desc='description')
    yWire = Float(0, iotype='in', desc='description')
    d_ = Array(np.zeros(2), iotype='in', desc='description')
    theta_ = Array(np.zeros(2), iotype='in', desc='description')
    nTube_ = Array(np.zeros(2), iotype='in', desc='description')
    nCap_ = Array(np.zeros(2), iotype='in', desc='description')
    lBiscuit_ = Array(np.zeros(2), iotype='in', desc='description')

    cE = Array(np.zeros(2), iotype='out', desc='description')
    cN = Array(np.zeros(2), iotype='out', desc='description')
    c100 = Array(np.zeros(100), iotype='out', desc='description')
    cl = Array(np.zeros(2), iotype='out', desc='description')
    Cm = Array(np.zeros(2), iotype='out', desc='description')
    t = Array(np.zeros(2), iotype='out', desc='description')
    xtU = Array(np.zeros(2), iotype='out', desc='description')
    xtL = Array(np.zeros(2), iotype='out', desc='description')
    xEA = Array(np.zeros(2), iotype='out', desc='description')
    d = Array(np.zeros(2), iotype='out', desc='description')
    theta = Array(np.zeros(2), iotype='out', desc='description')
    nTube = Array(np.zeros(2), iotype='out', desc='description')
    nCap = Array(np.zeros(2), iotype='out', desc='description')
    lBiscuit = Array(np.zeros(2), iotype='out', desc='description')


    def execute(self):
        pass

class wireProperties(Component):
    """
    Computes wire properties
    """
    type_flag = Int(0, iotype='in', desc='description')

    RHO = Float(0, iotype='out', desc='description')
    E = Float(0, iotype='out', desc='description')
    ULTIMATE = Float(0, iotype='out', desc='description')

    def execute(self):
        pass
