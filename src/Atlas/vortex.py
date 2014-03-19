from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np
from numpy import pi, cos, sin, mean, linspace, sqrt

class vortexRing(Component):
    """
    Vortex ring calculations
    """
    # Inputs:
    yN = Array(np.arange(11), iotype="in", desc='Node locations')
    rho = Float(1.18, iotype="in", desc="Air density")
    vc = Float(0, iotype="in", desc="Vertical velocity")
    Omega = Float(1.0367, iotype="in", desc="Rotor angular velocity")
    b = Int(2, iotype="in", desc="Number of blades")
    h = Float(1.5, iotype="in", desc="Height of rotor")
    Nw = Int(8, iotype="in")
    Ntt = Int(3, iotype="in")
    Ntheta = Int(20, iotype="in")
    dT = Array(np.array([0.0649, 1.5887, 6.8049, 11.2556, 15.7854, 19.8941,
                         23.2617,25.7216, 27.2537, 18.6446]),
               iotype="in", desc="Thrust")
    qh = Array(np.array([0, -0.0179, -0.0427, -0.0727, -0.1049, -0.1331,
                         -0.1445, -0.1247, -0.0789, -0.0181, 0.0480]),
               iotype="in")

    # Outputs:
    dtheta = Float(1, iotype="out")
    cr = Float(0, iotype="out")
    Ns = Int(0, iotype="out")
    Gamma = Array(np.zeros(10), iotype="out")
    z = Array(np.zeros(10), iotype="out")
    r = Array(np.zeros(10), iotype="out")
    vz = Array(np.zeros(10), iotype="out")
    vr = Array(np.zeros(10), iotype="out")
    thetaArray = Array(np.zeros(10), iotype="out")
    yE = Array(np.zeros(10), iotype="out")


    def execute(self):
        self.Ns=max(self.yN.shape) - 1
        dy=np.zeros(self.Ns)
        self.yE=np.zeros(self.Ns)
        for s in range(1,(self.Ns+1)):
            dy[(s-1)]=self.yN[(s + 1-1)] - self.yN[(s-1)]
            self.yE[(s-1)]=0.5 * (self.yN[(s-1)] + self.yN[(s + 1-1)])
        R=self.yN[(self.Ns + 1-1)]
        self.cr=0.5 * mean(dy)
        self.dtheta=pi / self.Ntheta
        self.thetaArray=linspace(0 + self.dtheta / 2,pi - self.dtheta / 2,self.Ntheta)
        self.Gamma=np.zeros((self.Nw + 1,self.Ns + 1))
        self.z=np.zeros((self.Nw + 1,self.Ns + 1))
        self.r=np.zeros((self.Nw + 1,self.Ns + 1))
        self.vz=np.zeros((self.Nw + 1,self.Ns + 1))
        self.vr=np.zeros((self.Nw + 1,self.Ns + 1))
        vi=np.zeros((self.Ns,1))

        GammaBound = self.dT / (self.rho*(self.Omega*self.yE)*dy)


        self.Gamma[0,0]=- GammaBound[0]
        for s in range(2,(self.Ns+1)):
            self.Gamma[0,(s-1)]=GammaBound[(s - 1-1)] - GammaBound[(s-1)]
        self.Gamma[0,(self.Ns + 1-1)]=GammaBound[(self.Ns-1)]
        self.r[0,:]=self.yN
        self.z[0,:]=self.qh[:]

        for t in range(1,(self.Nw+1)):
            for tt in range(1,(self.Ntt+1)):
                self.vz=np.zeros((self.Nw + 1,self.Ns + 1))
                self.vr=np.zeros((self.Nw + 1,self.Ns + 1))
                for i in range(1,(t+1)):
                    for s in range(1,(self.Ns + 1+1)):
                        for ii in range(1,(t+1)):
                            for ss in range(2,(self.Ns + 1+1)):
                                zr=self.z[(ii-1),(ss-1)]
                                r=self.r[(ii-1),(ss-1)]
                                zp=self.z[(i-1),(s-1)]
                                yp=self.r[(i-1),(s-1)]
                                M=self.Gamma[(ii-1),(ss-1)] * r * self.dtheta / (2 * pi)
                                X2=(- r * sin(self.thetaArray)) ** 2
                                Y2=(yp - r * cos(self.thetaArray)) ** 2
                                Z2=(zp - zr) ** 2
                                Normal=sqrt(X2 + Y2 + Z2)
                                for iii in range(1,(max(self.thetaArray.shape)+1)):
                                    if Normal[(iii-1)] < self.cr:
                                        Normal[(iii-1)]=self.cr
                                Norm3=Normal ** 3
                                self.vr[(i-1),(s-1)]=self.vr[(i-1),(s-1)] + np.sum(- cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[(i-1),(s-1)]=self.vz[(i-1),(s-1)] + np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M
                                zr=- 2 * self.h - self.z[(ii-1),(ss-1)]
                                Z2=(zp - zr) ** 2
                                Normal=sqrt(X2 + Y2 + Z2)
                                for iii in range(1,(max(self.thetaArray.shape)+1)):
                                    if Normal[(iii-1)] < self.cr:
                                        Normal[(iii-1)]=self.cr
                                Norm3=Normal ** 3
                                self.vr[(i-1),(s-1)]=self.vr[(i-1),(s-1)] - np.sum(- cos(self.thetaArray) * (zp - zr) / Norm3) * M
                                self.vz[(i-1),(s-1)]=self.vz[(i-1),(s-1)] - np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M

                if tt == 1:
                    PiApprox=8 * np.sum(self.dT.dot(vi))
                realtime=2 * pi / self.Omega / self.b * ((t - 1) * self.Ntt + tt) / self.Ntt
                altitude=self.vc * realtime
                dt=2 * pi / self.Omega / self.b / self.Ntt

                for i in range(1,(t+1)):
                    for s in range(1,(self.Ns + 1+1)):
                        self.z[(i-1),(s-1)]=self.z[(i-1),(s-1)] + self.vz[(i-1),(s-1)] * dt
                        self.r[(i-1),(s-1)]=self.r[(i-1),(s-1)] + self.vr[(i-1),(s-1)] * dt
                        self.Gamma[(i-1),(s-1)]=self.Gamma[(i-1),(s-1)]
            for i in range(1, t+1)[::-1]:
                self.Gamma[(i + 1-1),:]=self.Gamma[(i-1),:]
                self.r[(i + 1-1),:]=self.r[(i-1),:]
                self.z[(i + 1-1),:]=self.z[(i-1),:]

            GammaBound = self.dT / (self.rho*(self.Omega*self.yE)*dy)
            self.Gamma[0,0]=- GammaBound[0]
            for s in range(2,(self.Ns+1)):
                self.Gamma[0,(s-1)]=GammaBound[(s - 1-1)] - GammaBound[(s-1)]
            self.Gamma[0,(self.Ns + 1-1)]=GammaBound[(self.Ns-1)]
            self.r[0,:]=self.yN
            self.z[0,:]=self.qh[:]


class inducedVelocity(Component):
    """
    Induced Velocity calculation
    """
    # Inputs:
    qh = Array(np.array([0, -0.0179, -0.0427, -0.0727, -0.1049, -0.1331,
                         -0.1445, -0.1247, -0.0789, -0.0181, 0.0480]),
               iotype="in")
    h = Float(1.5, iotype="in", desc="Height of rotor")
    Gamma = Array(np.zeros(10), iotype="in")
    z = Array(np.zeros(10), iotype="in")
    r = Array(np.zeros(10), iotype="in")
    thetaArray = Array(np.zeros(10), iotype="in")
    yE = Array(np.zeros(10), iotype="in")
    cr = Float(0, iotype="in")
    Ns = Int(0, iotype="in")
    Nw = Int(8, iotype="in")
    dtheta = Float(1, iotype="in")

    # Outputs:
    vi = Array(np.zeros(10), iotype="out")


    def execute(self):
        for s in range(1,(self.Ns+1)):
            self.vi[(s-1)]=0
            for ii in range(1,(self.Nw+1)):
                for ss in range(2,(self.Ns + 1+1)):
                    if (ii == 1):
                        ringFrac=0.675
                    else:
                        ringFrac=1
                    zr = self.z[(ii-1),(ss-1)]
                    r = self.r[(ii-1),(ss-1)]
                    zp = (self.qh[(s-1)] + self.qh[(s + 1-1)]) / 2
                    yp = self.yE[(s-1)]
                    M = self.Gamma[(ii-1),(ss-1)] * r * self.dtheta / (2 * pi)
                    X2 = (- r * sin(self.thetaArray)) ** 2
                    Y2 = (yp - r * cos(self.thetaArray)) ** 2
                    Z2 = (zp - zr) ** 2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(1,(max(self.thetaArray.shape)+1)):
                        if Normal[(iii-1)] < self.cr:
                            Normal[(iii-1)] = self.cr
                    Norm3 = Normal ** 3
                    self.vi[(s-1)] = self.vi[(s-1)] + ringFrac * np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M
                    zr = - 2 * self.h - self.z[(ii-1),(ss-1)]
                    Z2 = (zp - zr) ** 2
                    Normal = sqrt(X2 + Y2 + Z2)
                    for iii in range(1,(max(self.thetaArray.shape)+1)):
                        if Normal[(iii-1)] < self.cr:
                            Normal[(iii-1)] = self.cr
                    Norm3 = Normal ** 3
                    self.vi[(s-1)] = self.vi[(s-1)] - ringFrac * np.sum((cos(self.thetaArray) * yp - r) / Norm3) * M
        self.vi=- self.vi





