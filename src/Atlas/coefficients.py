from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float


class dragCoefficient(Component):
    """
    Computes drag coefficient
    """
    Re = Float(0.0, iotype='in', desc='description')
    tc = Float(0.0, iotype='in', desc='description')
    xtcU = Float(0.0, iotype='in', desc='description')
    xtcL = Float(0.0, iotype='in', desc='description')

    Cd = Float(0.0, iotype='out', desc='description')

    def execute(self):
        pass
        
class dragCoefficientFit(Component):
    """
    Computes drag coefficient
    """
    Re = Float(0.0, iotype='in', desc='description')
    tc = Float(0.0, iotype='in', desc='description')
    xtcU = Float(0.0, iotype='in', desc='description')
    xtcL = Float(0.0, iotype='in', desc='description')

    Cd = Float(0.0, iotype='out', desc='description')

    def execute(self):
        pass

class frictionCoefficient(Component):
    """
    Computes friction coefficient
    """
    Re = Float(0.0, iotype='in', desc='description')
    xtc = Float(0.0, iotype='in', desc='description')

    Cfflat = Float(0.0, iotype='out', desc='description')

    def execute(self):
        pass