from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np

class sparProperties(Component):
    yN = Array(np.zeros(2), iotype='in', desc='description')
    d = Float(0, iotype='in', desc='description')
    theta = Float(0, iotype='in', desc='description')
    nTube = Int(0, iotype='in', desc='description')
    nCap = Array(np.zeros(2), iotype='in', desc='description')
    lBiscuit = Float(0, iotype='in', desc='description')
    CFRPType = Int(0, iotype='in', desc='description')

    EIx = Float(0, iotype='in', desc='description')
    EIz = Float(0, iotype='in', desc='description')
    EA = Float(0, iotype='in', desc='description')
    GJ = Float(0, iotype='in', desc='description')
    mSpar = Float(0, iotype='in', desc='description')

    def execute(self):
        pass

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