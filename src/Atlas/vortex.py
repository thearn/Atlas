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
    # inputs
    b        = Int(iotype='in', desc='number of blades')
    Ns       = Int(iotype='in', desc='number of elements')
    yN       = Array(iotype='in', desc='node locations')
    rho      = Float(iotype='in', desc='air density')
    vc       = Float(iotype='in', desc='vertical velocity')
    Omega    = Float(iotype='in', desc='rotor angular velocity')
    h        = Float(iotype='in', desc='height of rotor')
    dT       = Array(iotype='in', desc='thrust')
    q        = Array(iotype='in', desc='deformation')
    anhedral = Float(iotype='in')

    # outputs
    vi       = Array(iotype='out', desc='induced velocity')
    Gamma    = Array(iotype='out', desc='vortex strength')
    z        = Array(iotype='out', desc='')
    r        = Array(iotype='out', desc='')
    vz       = Array(iotype='out', desc='')
    vr       = Array(iotype='out', desc='')

    def execute(self):
        dy = np.zeros((self.Ns, 1))
        yE = np.zeros((self.Ns, 1))
        for s in range(self.Ns):
            dy[s] = self.yN[s+1] - self.yN[s]  # length of each element
            yE[s] = 0.5 * (self.yN[s] + self.yN[s+1])  # radial location of each element

        # set fidelity
        if self.Ns == 15:
            Nw = 15
            Ntt = 5
            Ntheta = 40
        elif self.Ns == 10:
            Nw = 8
            Ntt = 3
            Ntheta = 20
        else:
            Nw = 5
            Ntt = 1
            Ntheta = 15

        # Break out deformations
        qq = np.zeros((6, self.Ns+1))
        for s in range(1, self.Ns+1):
            qq[:, s] = self.q[s*6:s*6+6].T
        qh = (qq[2, :] - self.yN.T * self.anhedral).flatten()

        cr = 0.5 * mean(dy)
        dtheta = pi / Ntheta
        self.thetaArray = linspace(dtheta/2, pi - dtheta/2, Ntheta)

        # pre-allocate
        self.Gamma = np.zeros((Nw+1, self.Ns+1))
        self.z     = np.zeros((Nw+1, self.Ns+1))
        self.r     = np.zeros((Nw+1, self.Ns+1))
        self.vz    = np.zeros((Nw+1, self.Ns+1))
        self.vr    = np.zeros((Nw+1, self.Ns+1))
        self.vi    = np.zeros((self.Ns, 1))

        # create nacent vortex rings
        GammaBound = self.dT / (self.rho*(self.Omega*yE)*dy)

        self.Gamma[0, 0] = -GammaBound[0]
        for s in range(1, self.Ns):
            self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
        self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]

        self.r[0, :] = self.yN.T
        self.z[0, :] = qh[:]

        # free-wake time stepping
        for t in range(Nw+1):
            # proceed with substeps
            for tt in range(Ntt):
                # Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
                self.vz = np.zeros((Nw+1, self.Ns+1))
                self.vr = np.zeros((Nw+1, self.Ns+1))
                for i in range(t):                       # for each disk
                    for s in range(self.Ns+1):           # and for each ring on each disk
                        for ii in range(t):              # add the velocity induced from each disk
                            for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                                zr = self.z[ii, ss]
                                r  = self.r[ii, ss]
                                zp = self.z[i, s]
                                yp = self.r[i, s]
                                M  = self.Gamma[ii, ss] * r * dtheta / (2*pi)
                                X2 = (-r*sin(self.thetaArray))**2
                                Y2 = (yp - r*cos(self.thetaArray))**2
                                Z2 = (zp - zr)**2
                                Normal = sqrt(X2 + Y2 + Z2)
                                for iii in range(max(self.thetaArray.shape)):
                                    if Normal[iii] < cr:
                                        Normal[iii] = cr
                                Norm3 = Normal**3
                                self.vr[i, s] = self.vr[i, s] + np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[i, s] = self.vz[i, s] + np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M
                                zr = -2*self.h - self.z[ii, ss]
                                Z2 = (zp - zr)**2
                                Normal = sqrt(X2 + Y2 + Z2)
                                for iii in range(max(self.thetaArray.shape)):
                                    if Normal[iii] < cr:
                                        Normal[iii] = cr
                                Norm3 = Normal**3
                                self.vr[i, s] = self.vr[i, s] - np.sum(-cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[i, s] = self.vz[i, s] - np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M

                # Compute altitude and time power approximation
                # if tt == 0:
                #     PiApprox = 8 * np.sum(self.dT.dot(vi))
                # realtime = 2*pi / self.Omega / self.b * ((t-1) * Ntt + tt) / Ntt
                # altitude = self.vc * realtime

                # Convect rings downstream
                dt = 2*pi / self.Omega / self.b / Ntt
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

            # Create nacent vortex rings
            GammaBound = self.dT / (self.rho*(self.Omega*yE)*dy)
            self.Gamma[0, 0] = -GammaBound[0]
            for s in range(1, self.Ns):
                self.Gamma[0, s] = GammaBound[s-1] - GammaBound[s]
            self.Gamma[0, self.Ns] = GammaBound[self.Ns-1]
            self.r[0, :] = self.yN.T
            self.z[0, :] = qh[:]

        # Compute induced velocity on rotor (rp = [0 r(s) 0])
        for s in range(self.Ns):                # for each element
            self.vi[s] = 0
            for ii in range(Nw):                # add the velocity induced from each disk
                for ss in range(1, self.Ns+1):  # and each ring on each disk (inner ring cancels itself out)
                    if (ii == 0):
                        ringFrac = 0.675
                    else:
                        ringFrac = 1

                    zr = self.z[ii, ss]
                    r  = self.r[ii, ss]
                    zp = (qh[s] + qh[s+1]) / 2
                    yp = yE[s]

                    M  = self.Gamma[ii, ss] * r * dtheta / (2*pi)
                    X2 = (-r*sin(self.thetaArray))**2
                    Y2 = (yp - r*cos(self.thetaArray))**2
                    Z2 = (zp - zr)**2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(max(self.thetaArray.shape)):
                        if Normal[iii] < cr:
                            Normal[iii] = cr
                    Norm3 = Normal**3

                    self.vi[s] = self.vi[s] + ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

                    # ground effect ring
                    zr = -2*self.h - self.z[ii, ss]
                    Z2 = (zp - zr)**2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(max(self.thetaArray.shape)):
                        if Normal[iii] < cr:
                            Normal[iii] = cr
                    Norm3 = Normal**3

                    self.vi[s] = self.vi[s] - ringFrac * np.sum((cos(self.thetaArray)*yp - r) / Norm3) * M

        # vi is positive downwards
        self.vi = -self.vi
