from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Int
import numpy as np


class vortexWakeCover(Component):
    """
    Computes vortexWakeCover 
    """
    yN = Array(np.zeros(2), iotype='in', desc='description')
    rho = Float(0, iotype='in', desc='description')
    dT = Array(np.zeros(2), iotype='in', desc='description')
    vc = Float(0, iotype='in', desc='description')
    Omega = Float(0, iotype='in', desc='description')
    b = Float(0, iotype='in', desc='description')
    h = Float(0, iotype='in', desc='description')
    Nw = Int(0, iotype='in', desc='description')
    Ntt = Int(0, iotype='in', desc='description')
    Ntheta = Int(0, iotype='in', desc='description')
    qh = Array(np.zeros(2), iotype='in', desc='description')
    ycmax = Float(0, iotype='in', desc='description')
    flag_cover = Int(0, iotype='in', desc='description')
    flag_plot = Int(0, iotype='in', desc='description')
    flag_dynamic_climb = Int(0, iotype='in', desc='description')

    vi = Array(np.zeros(2), iotype='out', desc='description')
    ring_gamma = Array(np.zeros(2), iotype='out', desc='description')
    ring_z = Array(np.zeros(2), iotype='out', desc='description')
    ring_r = Array(np.zeros(2), iotype='out', desc='description')
    ring_vz = Array(np.zeros(2), iotype='out', desc='description')
    ring_vr = Array(np.zeros(2), iotype='out', desc='description')


    def execute(self):
        pass