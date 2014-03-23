from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np


class torsionalBucklingFailure(Component):
    Ns = Int(0, iotype='in', desc='description')
    Finternal = Array(np.zeros(2), iotype='in', desc='description')
    d = Float(0, iotype='in', desc='description')
    theta = Float(0, iotype='in', desc='description')
    nTube = Int(0, iotype='in', desc='description')
    nCap = Int(0, iotype='in', desc='description')
    lBiscuit = Float(0, iotype='in', desc='description')
    CFRPType = Int(0, iotype='in', desc='description')

    failure = Float(0, iotype='out', desc='description')

    def execute(self):
        pass


class failureCalc(Component):
    """
    Computes failure
    """
    yN = Array(np.zeros(2), iotype='in', desc='description')
    Finternal = Array(np.zeros(2), iotype='in', desc='description')

    strain_top = Array(np.zeros(2), iotype='in', desc='description')
    strain_bottom = Array(np.zeros(2), iotype='in', desc='description')
    strain_back = Array(np.zeros(2), iotype='in', desc='description')
    strain_front = Array(np.zeros(2), iotype='in', desc='description')
    strain_bending_x = Array(np.zeros(2), iotype='in', desc='description')
    strain_bending_z = Array(np.zeros(2), iotype='in', desc='description')
    strain_axial_y = Array(np.zeros(2), iotype='in', desc='description')
    strain_torsion_y = Array(np.zeros(2), iotype='in', desc='description')

    d = Array(np.zeros(2), iotype='in', desc='description')
    theta = Array(np.zeros(2), iotype='in', desc='description')
    nTube = Array(np.zeros(2), iotype='in', desc='description')
    nCap = Array(np.zeros(2), iotype='in', desc='description')
    yWire = Float(0, iotype='in', desc='description')
    zWire = Float(0, iotype='in', desc='description')
    EIxJ = Float(0, iotype='in', desc='description')
    EIzJ = Float(0, iotype='in', desc='description')

    lBiscuit = Array(np.zeros(2), iotype='in', desc='description')
    dQuad = Float(0, iotype='in', desc='description')
    thetaQuad = Float(0, iotype='in', desc='description')
    nTubeQuad = Float(0, iotype='in', desc='description')
    lBiscuitQuad = Float(0, iotype='in', desc='description')
    RQuad = Float(0, iotype='in', desc='description')
    hQuad = Float(0, iotype='in', desc='description')
    TQuad = Float(0, iotype='in', desc='description')
    EIQuad = Float(0, iotype='in', desc='description')
    GJQuad = Float(0, iotype='in', desc='description')
    tWire = Float(0, iotype='in', desc='description')
    TWire = Float(0, iotype='in', desc='description')
    TEtension = Float(0, iotype='in', desc='description')

    def execute(self):
        pass
