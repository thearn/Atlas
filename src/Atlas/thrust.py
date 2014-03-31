import numpy as np

from openmdao.main.api import Component
from openmdao.main.datatypes.api import Int, Float, Array


class Thrust(Component):

    #Aerocalc Inputs
    yN    = Array(iotype="in", desc='Node locations')
    ycmax = Float(iotype="in")
    Cl    = Array(iotype="in", desc="lift coefficient distribution")
    c     = Array(iotype="in", desc="chord distribution")
    rho   = Float(iotype="in", desc="Air density")
    Omega = Float(iotype="in", desc="Rotor angular velocity")

    #Aerocalc intermediate values
    dr = Array(iotype="in", desc="length of each element")
    r  = Array(iotype="in", desc="radial location of each element")

    # outputs
    dT = Array(iotype="out", desc="Thrust")
    chordFrac = Array(iotype="out")

    def execute(self):
        Ns = max(self.yN.shape) - 1  # number of elements

        self.chordFrac = np.ones(Ns)
        self.dT = np.zeros(Ns)

        # Compute multiplyer for partial element
        for index, element in enumerate(self.yN):
            if element < self.ycmax:
                sTrans = index  # determine transitional partial element

        self.chordFrac[sTrans] = self.yN[sTrans+1] - self.ycmax  \
                               / (self.yN[sTrans+1] - self.yN[sTrans])

        # Compute thrust assuming small angles
        self.dT = self.chordFrac * 0.5 * self.rho \
                * (self.Omega * self.r)**2        \
                * self.Cl * self.c * self.dr


class ActuatorDiskInducedVelocity(Component):
    """
    Compute induced velocity using annual-ring actuator disk theory
    """

    # inputs
    vc  = Float(iotype="in", desc="vertical velocity")
    rho = Float(iotype="in", desc="Air density")
    b   = Float(iotype="in", desc="Number of blades")
    R   = Float(iotype="in", desc="rotor radius")
    h   = Float(iotype="in", desc="Height of rotor")

    yN  = Int(iotype="in")

    dT  = Array(iotype="in", desc="Thrust")
    r   = Array(iotype="in")
    dr  = Array(iotype="in")

    # outputs
    vi  = Array(iotype="out")

    def execute(self):
        Ns = max(self.yN.shape) - 1  # number of elements

        self.vi = np.zeros(Ns)
        sq = np.zeros(Ns)

        sq = 0.25 * self.vc**2 + \
             0.25 * self.b * self.dT / (np.pi * self.rho * self.r * self.dr)

        self.vi += -0.5 * self.vc + np.sqrt(sq)

        # Add ground effect Cheesemen & Benett's
        self.vi /= (1. + (self.R / self.h / 4.) ** 2)
