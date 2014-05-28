from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float
from numpy import sqrt


def frictionCoefficient(Re, xtc):
    """
    Computes Cf of a flat plate at Re, with xtc fraction of laminar flow
    """
    # Fully turbulent
    if xtc == 0:
        Cfflat = 0.072/(Re)**0.2
    # Fully laminar
    elif xtc == 1:
        Cfflat = (1.328/sqrt(Re))
    # Partially laminar
    else:
        Cflam = (1.328/sqrt(Re))*xtc**(-0.5)             # Cf of laminar part
        deltalamc = (5/sqrt(Re))*sqrt(xtc)               # boundary layer thickness, delta/c, of laminar part
        deltaturbc = (0.13/0.097)*deltalamc              # boundary layer thickness, delta/c, of turbulent part
        x0c = xtc - (Re**0.2*deltaturbc/0.375)**(1/0.8)  # imaginary start point of turbulent BL
        CfturbFull = 0.072/((1-x0c)*Re)**0.2             # Cf of flat plate of length c-x0
        CfturbStart = 0.072/((xtc-x0c)*Re)**0.2          # Cf of imaginary part of turb BL
        Cfturb = (CfturbFull*(1-x0c) - CfturbStart*(xtc-x0c))/(1-xtc)  # Cf of turbulent part
        Cfflat = Cflam*xtc + Cfturb*(1-xtc)

    return Cfflat


class DragCoefficient(Component):
    """
    Computes drag coefficient
    """
    Re = Float(0.0, iotype='in', desc='description')
    tc = Float(0.0, iotype='in', desc='description')
    xtcU = Float(0.0, iotype='in', desc='description')
    xtcL = Float(0.0, iotype='in', desc='description')

    Cd = Float(0.0, iotype='out', desc='description')

    def execute(self):
        CfU = frictionCoefficient(self.Re, self.xtcU)
        CfL = frictionCoefficient(self.Re, self.xtcL)
        Cfflat = (CfU + CfL)/2
        self.Cd = 2*Cfflat*(1 + 2*self.tc + 60*self.tc**4)
