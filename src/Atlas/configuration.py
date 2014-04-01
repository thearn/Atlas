import numpy as np

from openmdao.main.api import Component, VariableTree
from openmdao.main.datatypes.api import Int, Float, Array, Str, Enum, VarTree


class Flags(VariableTree):
    Opt          = Int(1, desc='0 - single run, 1 - optimization')
    ConFail      = Int(0, desc='1 - use structural failure as a constraint on optimization')
    ConWireCont  = Int(0, desc='1 - use wire length continuity as a constraint to set appropriate wire forces in multi-point optimizations')
    ConJigCont   = Int(0, desc='1 - use jig continuity')
    ConDef       = Int(0, desc='1 - constraints on maximum deformation of the rotor')
    MultiPoint   = Int(4, desc='0 - single point optimization, 1 - 4 point optimization (h=0.5, h=3, wind case, gravity load)')
    Quad         = Int(1, desc='0 - prop drive, 1 - quad rotor')
    FreeWake     = Int(1, desc='0 - momentum theory, 1 - free vortex ring wake')
    PlotWake     = Int(0, desc='0 - dont plot wake, 1 - plot wake ')
    DynamicClimb = Int(0, desc='0 - vc imposes downward velocity, 1 - vc represents climb (final altitude depends on Nw)')
    Cover        = Int(0, desc='0 - no cover over root rotor blades, 1 - cover')
    Load         = Int(0, desc='0 - normal run, 1 - gravity forces only, 2 - prescribed load from pLoad')
    Cdfit        = Int(1, desc='0 - analytic model for drag coefficient, 1 - curve fit on BE airfoils')
    GWing        = Int(1, desc='0 - Daedalus style wing, 1 - Gossamer style wing (changes amount of laminar flow)')
    AeroStr      = Int(1, desc='0 - Assume flat wing, 1 - take deformation into account')
    Movie        = Int(0, desc='0 - dont save animation, 1 - save animation')
    wingWarp     = Int(0, desc='0 - no twist constraint, >0 - twist constraint at wingWarp')

    CFRPType     = Str('NCT301-1X HS40 G150 33 +/-2%RW', desc='type of carbon fibre reinforced polymer')

    WireType     = Enum('Pianowire',
                        ('Pianowire', 'Vectran'),
                        desc='Material to be used for lift wire')


class AtlasConfiguration(Component):
    """ Atlas configuration
    """
    # inputs
    yN    = Array(iotype="in", desc='node locations')

    # outputs (manual configuration)
    flags = VarTree(Flags(), iotype='out')

    # outputs (calculated)
    Ns = Int(iotype="out",   desc="number of elements")
    dr = Array(iotype="out", desc="length of each element")
    r  = Array(iotype="out", desc="radial location of each element")
    R  = Float(iotype="out", desc="rotor radius")

    def execute(self):
        self.Ns = max(self.yN.shape) - 1  # number of elements

        self.dr = np.zeros(self.Ns)
        self.r = np.zeros(self.Ns)

        for s in range(self.Ns):
            self.dr[s] = self.yN[s+1] - self.yN[s]     # length of each element
            self.r[s] = 0.5*(self.yN[s] + self.yN[s+1])

        self.R = self.yN[self.Ns]
