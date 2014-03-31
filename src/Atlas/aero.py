import numpy as np

from openmdao.main.api import Component, Assembly
from openmdao.main.datatypes.api import Int, Float, Array

from Atlas import Thrust, ActuatorDiskInducedVelocity, LiftDrag


class Configuration(Component):
    # inputs
    yN    = Array(iotype="in", desc='Node locations')

    # outputs
    Ns = Int(iotype="out",   desc="number of elements")
    dr = Array(iotype="out", desc="length of each element")
    r  = Array(iotype="out", desc="radial location of each element")
    R  = Float(iotype="out", desc="rotor radius")

    def execute(self):
        self.Ns = max(self.yN.shape) - 1  # number of elements

        self.self.dr = np.zeros(self.Ns)
        self.r = np.zeros(self.Ns)

        for s in range(self.Ns+1):
            self.dr[s-1] = self.yN[s] - self.yN[s-1]  # length of each element
            self.r[s-1] = 0.5*self.yN[s-1] + self.yN[s]

        self.R = self.yN[self.Ns+1]


class AeroCalc(Assembly):
    #yN,rho,visc,vw,vc,b,h,Omega,c,Cl,Cm,t,xtU,xtL,yWire,zWire,tWire,d,q,anhedral,ycmax,flags):
    def configure(self):
        self.add('thrust', Thrust())
        self.connect('yN', 'thrust.yN')
        self.connect('ycmax', 'thrust.ycmax')
        self.connect('Cl', 'thrust.Cl')
        self.connect('c', 'thrust.c')
        self.connect('rho', 'thrust.rho')
        self.connect('Omega', 'thrust.Omega')
        self.connect('dr', 'thrust.dr')
        self.connect('r', 'thrust.r')

        self.add('adiv', ActuatorDiskInducedVelocity())
        self.connect('vc',  'adiv.vc')
        self.connect('rho', 'adiv.rho')
        self.connect('b',   'adiv.b')
        self.connect('R',   'adiv.R')
        self.connect('h',   'adiv.h')
        self.connect('yN',  'adiv.yN')
        self.connect('r',   'adiv.r')
        self.connect('dr',  'adiv.dr')

        self.add('lift_drag', LiftDrag())
        self.connect('yN',    'lift_drag.yN')
        self.connect('rho',   'lift_drag.rho')
        self.connect('visc',  'lift_drag.visc')
        self.connect('vw',    'lift_drag.vw')
        self.connect('vc',    'lift_drag.vc')
        self.connect('Omega', 'lift_drag.Omega')
        self.connect('r',     'lift_drag.r')
        self.connect('t',     'lift_drag.t')
        self.connect('c',     'lift_drag.c')
        self.connect('cl',    'lift_drag.cl')
        self.connect('dr',    'lift_drag.dr')
        self.connect('d',     'lift_drag.d')
        self.connect('yWire', 'lift_drag.yWire')
        self.connect('zWire', 'lift_drag.zWire')
        self.connect('tWire', 'lift_drag.tWire')
        self.connect('Cm',    'lift_drag.Cm')
        self.connect('xtU',   'lift_drag.xtU')
        self.connect('xtL',   'lift_drag.xtL')

        self.connect('thrust.vi', 'lift_drag.vi')
        self.connect('thrust.chordFrac', 'lift_drag.chordFrac')

        self.connect('thrust.dT', 'dT')


    # Ns=max(yN.shape) - 1
    # dr=np.zeros(Ns,1)
    # r=np.zeros(Ns,1)
    # for s in range(1,(Ns+1)):
    #     dr[(s-1)]=yN[(s + 1-1)] - yN[(s-1)]
    #     r[(s-1)]=0.5 * (yN[(s-1)] + yN[(s + 1-1)])
    # R=yN[(Ns + 1-1)]

    # Fblade.Fx=np.zeros(Ns,1)
    # Fblade.Fz=np.zeros(Ns,1)
    # Fblade.My=np.zeros(Ns,1)
    # Fblade.Q=np.zeros(Ns,1)
    # Fblade.P=np.zeros(Ns,1)
    # Fblade.Pi=np.zeros(Ns,1)
    # Fblade.Pp=np.zeros(Ns,1)
    # phi=np.zeros(Ns,1)
    # vi=np.zeros(Ns,1)
    # dT=np.zeros(Ns,1)
    # Re=np.zeros(Ns,1)
    # Cd=np.zeros(Ns,1)
    # chordFrac=ones(Ns,1)
    # for s in range(1,(Ns + 1+1)):
    #     if yN[(s-1)] < ycmax:
    #         sTrans=s
    # chordFrac[(sTrans-1)]=(yN[(sTrans + 1-1)] - ycmax) / (yN[(sTrans + 1-1)] - yN[(sTrans-1)])
    # for s in range(1,(Ns+1)):
    #     dT[(s-1)]=chordFrac[(s-1)] * 0.5 * rho * (Omega * r[(s-1)]) ** 2 * Cl[(s-1)] * c[(s-1)] * dr[(s-1)]
    # if flags.FreeWake:
    #     if Ns == 15:
    #         Nw=15
    #         Ntt=5
    #         Ntheta=40
    #     else:
    #         if Ns == 10:
    #             Nw=8
    #             Ntt=3
    #             Ntheta=20
    #         else:
    #             Nw=5
    #             Ntt=1
    #             Ntheta=15
    #     qq(:,1)=np.array([0,0,0,0,0,0]).reshape(1,-1)
    #     for s in range(2,(Ns + 1+1)):
    #         qq[:,(s-1)]=q[((s - 1) * 6 + 1-1):(s - 1) * 6 + 6]
    #     qh=qq[2,:] - yN.T * anhedral
    #     vi,ring = VortexWakeCover(yN,rho,dT,vc,Omega,b,h,Nw,Ntt,Ntheta,qh,ycmax,flags.Cover,flags.PlotWake,flags.DynamicClimb)
    # else:
    #     for s in range(1,(Ns+1)):
    #         vi[(s-1)]=- 0.5 * vc + sqrt(0.25 * vc ** 2 + 0.25 * b * dT[(s-1)] / (pi * rho * r[(s-1)] * dr[(s-1)]))
    #     vi=vi / (1 + (R / h / 4) ** 2)
    #     ring=np.array([])
    # for s in range(1,(Ns+1)):
    #     U=sqrt((Omega * r[(s-1)] + vw) ** 2 + (vc + vi[(s-1)]) ** 2)
    #     if c[(s-1)] > 0.001:
    #         Re[(s-1)]=rho * U * c[(s-1)] / visc
    #         Cd[(s-1)]=dragCoefficientFit(Re[(s-1)],t[(s-1)],xtU[(s-1)],xtL[(s-1)])
    #         dL=0.5 * rho * U ** 2 * Cl[(s-1)] * c[(s-1)] * dr[(s-1)]
    #         dD=0.5 * rho * U ** 2 * Cd[(s-1)] * c[(s-1)] * dr[(s-1)]
    #     else:
    #         Re[(s-1)]=rho * U * d[(s-1)] / visc
    #         if Re < 3500:
    #             Cd[(s-1)]=- 1e-10 * Re[(s-1)] ** 3 + 7e-07 * Re[(s-1)] ** 2 - 0.0013 * Re[(s-1)] + 1.7397
    #         else:
    #             Cd[(s-1)]=1
    #         dL=0
    #         dD=0.5 * rho * U ** 2 * Cd[(s-1)] * d[(s-1)] * dr[(s-1)]
    #     for w in range(1,(max(yWire.shape)+1)):
    #         if yN[(s-1)] < yWire[(w-1)]:
    #             if yN[(s + 1-1)] < yWire[(w-1)]:
    #                 L=dr[(s-1)] * sqrt(zWire ** 2 + yWire[(w-1)] ** 2) / yWire[(w-1)]
    #             else:
    #                 L=(yWire[(w-1)] - yN[(s-1)]) * sqrt(zWire ** 2 + yWire[(w-1)] ** 2) / yWire[(w-1)]
    #             ReWire=rho * U * tWire / visc
    #             CdWire=- 1e-10 * ReWire ** 3 + 7e-07 * ReWire ** 2 - 0.0013 * ReWire + 1.7397
    #             dD=dD + 0.5 * rho * U ** 2 * CdWire * tWire * L
    #     phi[(s-1)]=atan2(vc + vi[(s-1)],vw + Omega * r[(s-1)])
    #     Fblade.Fz(s)=chordFrac[(s-1)] * (dL * cos(phi[(s-1)]) - dD * sin(phi[(s-1)]))
    #     Fblade.Fx(s)=chordFrac[(s-1)] * (dD * cos(phi[(s-1)]) + dL * sin(phi[(s-1)]))
    #     Fblade.My(s)=chordFrac[(s-1)] * (0.5 * rho * U ** 2 * Cm[(s-1)] * c[(s-1)] * c[(s-1)] * dr[(s-1)])
    #     Fblade.Q(s)=Fblade.Fx(s) * r[(s-1)]
    #     Fblade.P(s)=Fblade.Q(s) * Omega
    #     Fblade.Pi(s)=chordFrac[(s-1)] * (dL * sin(phi[(s-1)]) * r[(s-1)] * Omega)
    #     Fblade.Pp(s)=chordFrac[(s-1)] * (dD * cos(phi[(s-1)]) * r[(s-1)] * Omega)
    # return Fblade,vi,phi,Re,Cd,ring
