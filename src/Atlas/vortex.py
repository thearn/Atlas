from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np
from numpy import pi, cos, sin, mean, linspace, sqrt


class VortexRing(Component):
    """
    Vortex ring calculations
    Computes the induced velocity on the rotor blades given the
    thrust distribution on the rotor
    """
    # Inputs:
    yN     = Array(np.arange(11), iotype="in", desc='Node locations')
    rho    = Float(1.18, iotype="in", desc="Air density")
    vc     = Float(0, iotype="in", desc="Vertical velocity")
    Omega  = Float(1.0367, iotype="in", desc="Rotor angular velocity")
    b      = Int(2, iotype="in", desc="Number of blades")
    h      = Float(1.5, iotype="in", desc="Height of rotor")
    Nw     = Int(8, iotype="in")
    Ntt    = Int(3, iotype="in")
    Ntheta = Int(20, iotype="in")
    dT     = Array(np.array([0.0649,   1.5887,  6.8049, 11.2556, 15.7854,
                             19.8941, 23.2617, 25.7216, 27.2537, 18.6446]),
             iotype="in", desc="Thrust")
    qh = Array(np.array([0, -0.0179, -0.0427, -0.0727, -0.1049, -0.1331,
                            -0.1445, -0.1247, -0.0789, -0.0181, 0.0480]),
             iotype="in")

    # Outputs:
    dtheta = Float(1, iotype="out")
    cr     = Float(0, iotype="out")
    Ns     = Int(0, iotype="out")
    Gamma  = Array(np.zeros(10), iotype="out")
    z      = Array(np.zeros(10), iotype="out")
    r      = Array(np.zeros(10), iotype="out")
    vz     = Array(np.zeros(10), iotype="out")
    vr     = Array(np.zeros(10), iotype="out")
    thetaArray = Array(np.zeros(10), iotype="out")
    yE     = Array(np.zeros(10), iotype="out")

    def execute(self):
        self.Ns = max(self.yN.shape) - 1  # number of elements
        dy = np.zeros(self.Ns)
        self.yE = np.zeros(self.Ns)
        for s in range(self.Ns):
            dy[s] = self.yN[s+1] - self.yN[s]  # length of each element
            self.yE[s] = 0.5 * (self.yN[s] + self.yN[s+1])  # radial location of each element

        # R = self.yN[self.Ns]

        self.cr = 0.5 * mean(dy)
        self.dtheta = pi / self.Ntheta
        self.thetaArray = linspace(self.dtheta/2, pi - self.dtheta/2, self.Ntheta)

        # pre-allocate
        self.Gamma = np.zeros((self.Nw+1, self.Ns+1))
        self.z     = np.zeros((self.Nw+1, self.Ns+1))
        self.r     = np.zeros((self.Nw+1, self.Ns+1))
        self.vz    = np.zeros((self.Nw+1, self.Ns+1))
        self.vr    = np.zeros((self.Nw+1, self.Ns+1))
        self.vi    = np.zeros((self.Ns, 1))

        # create nacent vortex rings
        GammaBound = self.dT / (self.rho*(self.Omega*self.yE)*dy)

        self.Gamma[0, 0] = -GammaBound[0]
        for s in range(1, self.Ns):
            self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
        self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]

        self.r[0, :] = self.yN
        self.z[0, :] = self.qh[:]

        # free-wake time stepping
        for t in range(self.Nw+1):
            # proceed with substeps
            for tt in range(self.Ntt):
                # Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
                self.vz = np.zeros((self.Nw+1, self.Ns+1))
                self.vr = np.zeros((self.Nw+1, self.Ns+1))
                for i in range(t):                       # for each disk
                    for s in range(self.Ns+1):           # and for each ring on each disk
                        for ii in range(t):              # add the velocity induced from each disk
                            for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                                zr = self.z[ii, ss]
                                r  = self.r[ii, ss]
                                zp = self.z[i, s]
                                yp = self.r[i, s]
                                M  = self.Gamma[ii, ss] * r * self.dtheta / (2*pi)
                                X2 = (-r*sin(self.thetaArray))**2
                                Y2 = (yp - r*cos(self.thetaArray))**2
                                Z2 = (zp - zr)**2
                                Normal = sqrt(X2 + Y2 + Z2)
                                for iii in range(max(self.thetaArray.shape)):
                                    if Normal[iii] < self.cr:
                                        Normal[iii] = self.cr
                                Norm3 = Normal**3
                                self.vr[i, s] = self.vr[i, s] + np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[i, s] = self.vz[i, s] + np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M
                                zr = -2*self.h - self.z[ii, ss]
                                Z2 = (zp - zr)**2
                                Normal = sqrt(X2 + Y2 + Z2)
                                for iii in range(max(self.thetaArray.shape)):
                                    if Normal[iii] < self.cr:
                                        Normal[iii] = self.cr
                                Norm3 = Normal**3
                                self.vr[i, s] = self.vr[i, s] - np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[i, s] = self.vz[i, s] - np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M

                # Compute altitude and time power approximation
                # if tt == 0:
                #     PiApprox = 8 * np.sum(self.dT.dot(vi))
                # realtime = 2*pi / self.Omega / self.b * ((t-1) * self.Ntt + tt) / self.Ntt
                # altitude = self.vc * realtime

                # Convect rings downstream
                dt = 2*pi / self.Omega / self.b / self.Ntt
                for i in range(t):              # for each disk
                    for s in range(self.Ns+1):  # and for each ring on each disk
                        self.z[i, s] = self.z[i, s] + self.vz[i, s]*dt
                        self.r[i, s] = self.r[i, s] + self.vr[i, s]*dt
                        self.Gamma[i, s] = self.Gamma[i, s]

            # Shift elements in ring array
            for i in range(t)[::-1]:
                self.Gamma[i+1, :] = self.Gamma[i, :]
                self.r[i+1, :] = self.r[i, :]
                self.z[i+1, :] = self.z[i, :]

            GammaBound = self.dT / (self.rho*(self.Omega*self.yE)*dy)
            self.Gamma[0, 0] = -GammaBound[0]
            for s in range(1, self.Ns):
                self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
            self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]
            self.r[0, :] = self.yN
            self.z[0, :] = self.qh[:]


class InducedVelocity(Component):
    """
    Induced Velocity calculation
    """
    # Inputs:
    qh = Array(np.array([ 0, -0.0179, -0.0427, -0.0727, -0.1049, -0.1331,
                             -0.1445, -0.1247, -0.0789, -0.0181,  0.0480]),
               iotype="in")
    h      = Float(1.5, iotype="in", desc="Height of rotor")
    Gamma  = Array(np.zeros(10), iotype="in")
    z      = Array(np.zeros(10), iotype="in")
    r      = Array(np.zeros(10), iotype="in")
    thetaArray = Array(np.zeros(10), iotype="in")
    yE     = Array(np.zeros(10), iotype="in")
    cr     = Float(0, iotype="in")
    Ns     = Int(0, iotype="in")
    Nw     = Int(8, iotype="in")
    dtheta = Float(1, iotype="in")

    # Outputs:
    vi     = Array(np.zeros(10), iotype="out")

    def execute(self):
        for s in range(self.Ns):           # for each element
            self.vi[s] = 0
            for ii in range(self.Nw):      # add the velocity induced from each disk
                for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                    if (ii == 0):
                        ringFrac = 0.675
                    else:
                        ringFrac = 1

                    zr = self.z[ii, ss]
                    r  = self.r[ii, ss]
                    zp = (self.qh[s] + self.qh[s+1]) / 2
                    yp = self.yE[s]

                    M  = self.Gamma[ii, ss] * r * self.dtheta / (2*pi)
                    X2 = (-r*sin(self.thetaArray))**2
                    Y2 = (yp - r*cos(self.thetaArray))**2
                    Z2 = (zp - zr)**2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(max(self.thetaArray.shape)):
                        if Normal[iii] < self.cr:
                            Normal[iii] = self.cr
                    Norm3 = Normal**3

                    self.vi[s] = self.vi[s] + ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

                    # ground effect ring
                    zr = -2*self.h - self.z[ii, ss]
                    Z2 = (zp - zr)**2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(max(self.thetaArray.shape)):
                        if Normal[iii] < self.cr:
                            Normal[iii] = self.cr
                    Norm3 = Normal**3

                    self.vi[s] = self.vi[s] - ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

        # vi is positive downwards
        self.vi = -self.vi
