# pylint: disable=line-too-long, invalid-name, bad-whitespace, trailing-whitespace, too-many-locals, line-too-long
from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int, Str
import numpy as np
import math

#from __future__ import division

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
            

class discretizeProperties(Component):
    """
    Computes discretized Properties
    """

    #def DiscretizeProperties(Ns,ycmax,R,c_,Cl_,Cm_,t_,xtU_,xtL_,xEA_,yWire,d_,theta_,nTube_,nCap_,lBiscuit_):

    Ns = Int(0, iotype='in', desc='description')
    ycmax = Array( iotype='in', desc='description')
    R = Float(0, iotype='in', desc='description')
    c_ = Array( iotype='in', desc='description')
    Cl_ = Array( iotype='in', desc='description')
    Cm_ = Array( iotype='in', desc='description')
    t_ = Array( iotype='in', desc='description')
    xtU_ = Array( iotype='in', desc='description')
    xtL_ = Array( iotype='in', desc='description')
    xEA_ = Array( iotype='in', desc='description')
    yWire = Float(0, iotype='in', desc='description')
    d_ = Array( iotype='in', desc='description')
    theta_ = Array( iotype='in', desc='description')
    nTube_ = Array( iotype='in', desc='description')
    nCap_ = Array( iotype='in', desc='description')
    lBiscuit_ = Array( iotype='in', desc='description')

    #return cE,cN,c100,Cl,Cm,t,xtU,xtL,xEA,d,theta,nTube,nCap,lBiscuit,yN,yE
    cE = Array( iotype='out', desc='description')
    cN = Array( iotype='out', desc='description')
    c100 = Array(np.zeros(100), iotype='out', desc='description')
    cl = Array( iotype='out', desc='description')
    Cm = Array( iotype='out', desc='description')
    t = Array( iotype='out', desc='description')
    xtU = Array( iotype='out', desc='description')
    xtL = Array( iotype='out', desc='description')
    xEA = Array( iotype='out', desc='description')
    d = Array( iotype='out', desc='description')
    theta = Array( iotype='out', desc='description')
    nTube = Array( iotype='out', desc='description')
    nCap = Array( iotype='out', desc='description')
    lBiscuit = Array( iotype='out', desc='description')


    def execute(self):

        yN=np.zeros(Ns + 1,1)
        cN=np.zeros(Ns + 1,1)
        yE=np.zeros(Ns,1)
        cR=np.zeros(Ns,1)
        cZ=np.zeros(Ns,1)
        cE=np.zeros(Ns,1)
        Cl=np.zeros(Ns,1)
        Cm=np.zeros(Ns,1)
        t=np.zeros(Ns,1)
        xtU=np.zeros(Ns,1)
        xtL=np.zeros(Ns,1)
        xEA=np.zeros(Ns,1)
        d=np.zeros(Ns,1)
        theta=np.zeros(Ns,1)
        nTube=np.zeros(Ns,1)
        nCap=np.zeros(Ns,1)
        lBiscuit=np.zeros(Ns,1)
        Ns100=100
        yN100=np.zeros(Ns100 + 1,1)
        yE100=np.zeros(Ns100,1)
        cR100=np.zeros(100,1)
        cZ100=np.zeros(100,1)
        c100=np.zeros(100,1)
        for s in range(1,(Ns + 1+1)):
            yN[(s-1)]=R / Ns * (s - 1)
            if yN[(s-1)] < ycmax[0]:
                sTrans[0]=s
            if yN[(s-1)] < ycmax[1]:
                sTrans[1]=s
        for s in range(1,(Ns+1)):
            yE[(s-1)]=0.5 * (yN[(s-1)] + yN[(s + 1-1)])
        for s in range(1,(Ns+1)):
            if s < sTrans[0]:
                x=yE[(s-1)] / ycmax[0]
                cE[(s-1)]=c_[0] + x * (c_[1] - c_[0])
            else:
                x=(yE[(s-1)] - ycmax[0]) / (R - ycmax[0])
                pStart=c_[2]
                pCurve=c_[3]
                pEnd=c_[4]
                xx=x * (1 - pCurve) + sin(x * math.pi / 2) * pCurve
                cZ[(s-1)]=pStart + (pEnd - pStart) * xx
                c3=c_[2] / (c_[4] * R / ycmax[0])
                cR[(s-1)]=c_[4] * R / yE[(s-1)] * (c3 + (1 - c3) * x)
                cE[(s-1)]=cR[(s-1)] + (cZ[(s-1)] - cR[(s-1)]) * c_[1]
            if cE[(s-1)] == 0:
                cE[(s-1)]=0.001
        for s in range(1,(Ns100 + 1+1)):
            yN100[(s-1)]=R / Ns100 * (s - 1)
            if yN100[(s-1)] < ycmax:
                sTrans100=s
        for s in range(1,(Ns100+1)):
            yE100[(s-1)]=0.5 * (yN100[(s-1)] + yN100[(s + 1-1)])
        for s in range(1,(Ns100+1)):
            if s < sTrans100:
                x=yE100[(s-1)] / ycmax[0]
                c100[(s-1)]=c_[0] + x * (c_[1] - c_[0])
            else:
                x=(yE100[(s-1)] - ycmax[0]) / (R - ycmax[0])
                pStart=c_[2]
                pCurve=c_[3]
                pEnd=c_[4]
                xx=x * (1 - pCurve) + sin(x * math.pi / 2) * pCurve
                cZ100[(s-1)]=pStart + (pEnd - pStart) * xx
                c3=c_[2] / (c_[4] * R / ycmax[0])
                cR100[(s-1)]=c_[4] * R / yE100[(s-1)] * (c3 + (1 - c3) * x)
                c100[(s-1)]=cR100[(s-1)] + (cZ100[(s-1)] - cR100[(s-1)]) * c_[1]
            if c100[(s-1)] == 0:
                c100[(s-1)]=0.001
        Y=np.array([ycmax[0],ycmax[1],R]).reshape(1,-1)
        for s in range(1,(Ns+1)):
            if s < sTrans[0]:
                Cl[(s-1)]=Cl_[0]
                Cm[(s-1)]=Cm_[0]
                t[(s-1)]=t_[0]
                xEA[(s-1)]=xEA_[0]
            else:
                for j in range(1,(max(Y.shape)+1)):
                    if Y[(j-1)] > yE[(s-1)]:
                        break
                if s == sTrans[0]:
                    j=2
                x=(yE[(s-1)] - Y[(j - 1-1)]) / (Y[(j-1)] - Y[(j - 1-1)])
                Cl[(s-1)]=Cl_[(j - 1-1)] + x * (Cl_[(j-1)] - Cl_[(j - 1-1)])
                Cm[(s-1)]=Cm_[(j - 1-1)] + x * (Cm_[(j-1)] - Cm_[(j - 1-1)])
                t[(s-1)]=t_[(j - 1-1)] + x * (t_[(j-1)] - t_[(j - 1-1)])
                xEA[(s-1)]=xEA_[(j - 1-1)] + x * (xEA_[(j-1)] - xEA_[(j - 1-1)])
        Cl[(Ns-1)]=Cl[(Ns-1)] * 2 / 3
        for s in range(1,(Ns + 1+1)):
            if yN[(s-1)] < ycmax[0]:
                sTrans[0]=s
            if yN[(s-1)] < xtU_[1]:
                sTrans[1]=s
        for s in range(1,(Ns+1)):
            if s < sTrans[0]:
                xtU[(s-1)]=0.05
                xtL[(s-1)]=0.05
            else:
                if s < sTrans[1]:
                    xtU[(s-1)]=xtU_[0]
                    xtL[(s-1)]=xtL_[0]
                else:
                    if s == sTrans[1]:
                        x=(xtU_[1] - yN[(s-1)]) / (yN[(s-1)] - yN[(s - 1-1)])
                        xtU[(s-1)]=(1 - x) * xtU_[0] + x * xtU_[2]
                        xtL[(s-1)]=(1 - x) * xtL_[0] + x * xtL_[2]
                    else:
                        if s > sTrans[1]:
                            xtU[(s-1)]=xtU_[2]
                            xtL[(s-1)]=xtL_[2]
        Y=np.array([0,yWire[0],R]).reshape(1,-1)
        for s in range(1,(Ns+1)):
            for j in range(1,(max(Y.shape)+1)):
                if Y[(j-1)] > yE[(s-1)]:
                    break
            x=(yE[(s-1)] - Y[(j - 1-1)]) / (Y[(j-1)] - Y[(j - 1-1)])
            d[(s-1)]=d_[(j - 1-1)] + x * (d_[(j-1)] - d_[(j - 1-1)])
            theta[(s-1)]=theta_[(j - 1-1)] + x * (theta_[(j-1)] - theta_[(j - 1-1)])
            nTube[(s-1)]=nTube_[(j - 1-1)] + x * (nTube_[(j-1)] - nTube_[(j - 1-1)])
            nCap[(s-1)]=nCap_[(j - 1-1)] + x * (nCap_[(j-1)] - nCap_[(j - 1-1)])
            lBiscuit[(s-1)]=lBiscuit_[(j - 1-1)] + x * (lBiscuit_[(j-1)] - lBiscuit_[(j - 1-1)])
    










class wireProperties(Component):
    """
    Computes wire properties
    """
    type_flag = Int(0, iotype='in', desc='description')

    RHO = Float(0, iotype='out', desc='description')
    E = Float(0, iotype='out', desc='description')
    ULTIMATE = Float(0, iotype='out', desc='description')

    def execute(self):
        if self.type_flag == 1:
            self.RHO=7850.0
            self.E=2.1e+11
            self.ULTIMATE=2620000000.0
            #D=np.array([0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.72]).reshape(1,-1)
        elif self.type_flag == 2:
            self.RHO=1106.5
            self.E=39210000000.0
            self.ULTIMATE=982800000.0
            #D=np.array([1.5,2]).reshape(1,-1)







class chordProperties(Component):
    """
    Computes wire properties
    """

    yN = Array( iotype='in', desc='description')
    c = Array( iotype='in', desc='description')
    d = Array( iotype='in', desc='description')
    flagGWing = Int(iotype='in', desc='description')
    xtU = Array( iotype='in', desc='description')
    
    mChord = Array(iotype='out', desc='description')
    xCGChord = Array(iotype='out', desc='description')

    def execute(self):

        Ns=max(self.yN.shape) - 1
        dy=np.zeros(Ns)
        for s in range(1,(Ns+1)):
            dy[(s-1)]=self.yN[(s)] - self.yN[(s-1)]
        self.mChord=np.zeros(Ns)
        self.xCGChord=np.zeros(Ns)
        RHO_BALSA=160.0
        RHO_BASSWOOD=387.0
        RHO_CARBON=1610.6
        T_PLY_CARBON=0.000147
        RHO_EPS=16.0
        RHO_KEVLAR=1300.0
        T_PLY_KEVLAR=0.000127
        RHO_MYLAR=1400.0
        T_MYLAR=1.2e-05
        RHO_STRUCTURAL_FOAM=31.0
        RHO_XPS=24.67
        RHO_STEEL_WIRE=7850.0
        AREA_AIRFOIL=0.08233
        PERIMETER_AIRFOIL=2.06814
        XCG_AIRFOIL=0.39003
        thickness_rib=0.005
        rib_spacing=12 * (0.0254)
        thickness_LE_sheeting=0.003
        thickness_rib_caps=(1.0 / 32) * (0.0254)
        percent_rib_caps=0.97
        thickness_spar_plate=(1.0 / 32) * (0.0254)
        percent_radius_spar_plate=0.0632
        percent_LE_sheeting_top=0.6
        percent_LE_sheeting_top_GWing=self.xtU
        percent_LE_sheeting_bottom=0.1
        d_TE_spar=(3.0 / 4) * (0.0254)
        nTube_TE_spar=4
        d_comp_member=(3.0 / 4) * (0.0254)
        nTube_comp_member=4
        spacing_comp_member=3
        percent_length_comp_member=0.782 - 0.25
        d_cross_bracing=0.001
        height_TE_foam=0.007
        length_TE_foam=0.0505
        percent_height_TE_plate=0.05
        percent_length_TE_plate=0.24
        thickness_TE_plate=(1.0 / 32) * (0.0254)
        diameter_piano_wire=(0.022) * (0.0254)
        spar_location=0.25
        AF_ribs=1.569
        AF_TE=1.602
        AF_leading_edge_sheeting=1.302
        AF_leading_edge_sheeting_GWing=1.21
        AF_covering=0.972
        AF_TE_spar=0.529
        AF_comp_member=1.0
        AF_cross_bracing=4.603
        for s in range(1,(Ns+1)):
            mass_rib_foam=RHO_EPS * (thickness_rib * ((self.c[(s-1)] ** 2) * AREA_AIRFOIL))
            mass_rib_caps=RHO_BASSWOOD * (thickness_rib * thickness_rib_caps * (percent_rib_caps * PERIMETER_AIRFOIL * self.c[(s-1)]))
            mass_rib_plate_spar=RHO_BALSA * (thickness_spar_plate * (math.pi * (((self.c[(s-1)] * percent_radius_spar_plate) ** 2) - ((self.d[(s-1)] / 2) ** 2))))
            if self.flagGWing == 0:
                mass_rib_plate_TE=RHO_BALSA * (thickness_TE_plate * ((1.0 / 2) * (self.c[(s-1)] ** 2) * percent_height_TE_plate * percent_length_TE_plate))
                mass_rib=AF_ribs * ((dy[(s-1)] / rib_spacing) * (mass_rib_foam + mass_rib_caps + 2 * mass_rib_plate_spar + 2 * mass_rib_plate_TE))
                Xcg_rib=XCG_AIRFOIL
            else:
                mass_rib=AF_ribs * (0.66) * ((dy[(s-1)] / rib_spacing) * (mass_rib_foam + mass_rib_caps + 2 * mass_rib_plate_spar))
                Xcg_rib=((XCG_AIRFOIL + (2.0 / 3) * (spar_location)) / 2)
            if self.flagGWing == 0:
                mass_TE=AF_TE * (RHO_STRUCTURAL_FOAM * dy[(s-1)] * ((1.0 / 2) * (length_TE_foam * height_TE_foam)) + RHO_KEVLAR * dy[(s-1)] * T_PLY_KEVLAR * (length_TE_foam + height_TE_foam + sqrt(length_TE_foam ** 2 + height_TE_foam ** 2)))
                Xcg_TE=(self.c[(s-1)] - (2.0 / 3) * length_TE_foam) / self.c[(s-1)]
            else:
                mass_TE=RHO_STEEL_WIRE * dy[(s-1)] * math.pi * ((diameter_piano_wire / 2) ** 2)
                Xcg_TE=1
            if self.flagGWing == 0:
                mass_LE_sheeting=AF_leading_edge_sheeting * (RHO_XPS * dy[(s-1)] * self.c[(s-1)] * (percent_LE_sheeting_top + percent_LE_sheeting_bottom) * PERIMETER_AIRFOIL * thickness_LE_sheeting)
                Xcg_LE_sheeting=(1.0 / 2) * ((1.0 / 2) * (percent_LE_sheeting_top) + (1.0 / 2) * (percent_LE_sheeting_bottom))
            else:
                mass_LE_sheeting=AF_leading_edge_sheeting_GWing * (RHO_EPS * dy[(s-1)] * self.c[(s-1)] * (percent_LE_sheeting_top_GWing[(s-1)] + percent_LE_sheeting_bottom) * PERIMETER_AIRFOIL * thickness_LE_sheeting) + RHO_STEEL_WIRE * dy[(s-1)] * math.pi * ((diameter_piano_wire / 2) ** 2)
                Xcg_LE_sheeting=(1.0 / 2) * ((1.0 / 2) * (percent_LE_sheeting_top_GWing[(s-1)]) + (1.0 / 2) * (percent_LE_sheeting_bottom))
            mass_covering=AF_covering * (RHO_MYLAR * dy[(s-1)] * self.c[(s-1)] * PERIMETER_AIRFOIL * T_MYLAR)
            Xcg_covering=0.5
            if self.flagGWing == 0:
                mass_TE_spar=AF_TE_spar * (RHO_CARBON * dy[(s-1)] * (math.pi * (((d_TE_spar / 2) + nTube_TE_spar * T_PLY_CARBON) ** 2 - (d_TE_spar / 2) ** 2)))
                Xcg_TE_spar=0.9
                mass_comp_member=AF_comp_member * (dy[(s-1)] / spacing_comp_member) * RHO_CARBON * (self.c[(s-1)] * percent_length_comp_member) * (math.pi * (((d_comp_member / 2) + nTube_comp_member * T_PLY_CARBON) ** 2 - (d_comp_member / 2) ** 2))
                Xcg_comp_member=0.25 + (1.0 / 2) * (percent_length_comp_member)
                mass_cross_bracing=AF_cross_bracing * RHO_KEVLAR * (dy[(s-1)] / spacing_comp_member) * 2 * (math.pi * (d_cross_bracing / 2) ** 2) * ((((self.c[(s-1)] * percent_length_comp_member) ** 2) + ((spacing_comp_member) ** 2)) ** (1.0 / 2))
                Xcg_cross_bracing=0.25 + (1.0 / 2) * (percent_length_comp_member)
            else:
                mass_TE_spar=0
                Xcg_TE_spar=0
                mass_comp_member=0
                Xcg_comp_member=0
                mass_cross_bracing=0
                Xcg_cross_bracing=0
            self.mChord[(s-1)]=mass_rib + mass_TE + (mass_LE_sheeting) + (mass_covering) + (mass_TE_spar + mass_comp_member + mass_cross_bracing)
            self.xCGChord[(s-1)]=((mass_rib * Xcg_rib) + (mass_TE * Xcg_TE) + (mass_LE_sheeting * Xcg_LE_sheeting) + (mass_covering * Xcg_covering) + (mass_TE_spar * Xcg_TE_spar) + (mass_comp_member * Xcg_comp_member) + (mass_cross_bracing * Xcg_cross_bracing)) / (mass_rib + mass_TE + mass_LE_sheeting + mass_covering + mass_TE_spar + mass_comp_member + mass_cross_bracing)
