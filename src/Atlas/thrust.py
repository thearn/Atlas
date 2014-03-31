import numpy as np

from openmdao.main.api import Component
from openmdao.main.datatypes.api import Int, Float, Array

class Thrust(Component):

    #Aerocalc Inputs
    yN = Array(iotype="in")
    ycmax = Float(iotype="in")
    Cl = Array(iotype="in")
    c = Array(iotype="in")
    rho = Float(iotype="in")
    Omega = Float(iotype="in")

    #Aerocalc intermediate values
    Ns = Int(iotype="in")
    dr = Array(iotype="in")
    r = Array(iotype="in")
    dT = Array(iotype="out")
    chordFrac = Array(iotype="out")

    def execute(self):
        self.chordFrac = np.ones(self.Ns)
        self.dT = np.zeros(self.Ns)

        for index, element in enumerate(self.yN):
            if element < self.ycmax:
                sTrans = index

        self.chordFrac[sTrans]  = self.yN[sTrans + 1]
        self.chordFrac[sTrans] -= self.ycmax
        self.chordFrac[sTrans] /= (self.yN[sTrans + 1] - self.yN[sTrans])

        self.dT += self.chordFrac
        self.dT *= 0.5 
        self.dT *= self.rho
        self.dT *= (self.Omega * self.r) ** 2
        self.dT *= self.Cl 
        self.dT *= self.c
        self.dT *= self.dr

class InducedVelocity(Component):
    vc  = Float(iotype="in")
    rho = Float(iotype="in")
    b   = Float(iotype="in")
    R   = Float(iotype="in")
    h   = Float(iotype="in")

    Ns  = Int(iotype="in")

    dT = Array(iotype="in")
    r  = Array(iotype="in")
    dr = Array(iotype="in")

    vi = Array(iotype="out")

    def execute(self):
        self.vi = np.zeros(self.Ns)
        sqrt = np.zeros(self.Ns)

        sqrt = 0.25 * self.b * self.dT / (np.pi * self.rho * self.r * self.dr)
        sqrt += 0.25 * self.vc ** 2
        sqrt = np.sqrt(sqrt)

        self.vi += -0.5 * self.vc
        self.vi += sqrt

        self.vi /= (1. + (self.R / self.h / 4.) ** 2)



