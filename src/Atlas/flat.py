import numpy as np
from numpy import pi
from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.main.datatypes.api import Int, Float, Array, Str, Enum, VarTree
from openmdao.lib.drivers.api import SLSQPdriver
from openmdao.lib.drivers.api import FixedPointIterator

from Atlas import DragCoefficient, frictionCoefficient
from Atlas import prepregProperties, wireProperties, DiscretizeProperties, \
                       JointProperties, SparProperties, ChordProperties, \
                       JointSparProperties, QuadSparProperties
from Atlas import PrescribedLoad, Strain, \
                       MassProperties, FEM, Strains, Failures
from Atlas import LiftDrag, Fblade
from Atlas import VortexRing
from Atlas import VortexRingC
from Atlas import Thrust, ActuatorDiskInducedVelocity

from Atlas import AtlasConfiguration

class Results(Component):
    # inputs
    b          = Int(iotype='in', desc='number of blades')
    Ns         = Int(iotype='in', desc='number of elements')
    yN         = Array(iotype='in', desc='node locations')
    yE         = Array(iotype='in', desc='')
    cE         = Array(iotype='in', desc='chord of each element')
    Cl         = Array(iotype='in', desc='lift coefficient distribution')
    q          = Array(iotype='in', desc='deformation')
    phi        = Array(iotype='in', desc='')
    collective = Float(iotype='in', desc='collective angle in radians')
    fblade     = VarTree(Fblade(), iotype='in')
    Mtot       = Float(0.0, iotype='in', desc='total mass')

    # outputs
    di         = Array(iotype='out', desc='dihedral angle')
    alphaJig   = Array(iotype='out', desc='aerodynamic jig angle')
    Ttot       = Float(iotype='out', desc='')
    Qtot       = Float(iotype='out', desc='')
    MomRot     = Float(iotype='out', desc='')
    Ptot       = Float(iotype='out', desc='')

    def execute(self):
        # Compute aerodynamic jig angle
        self.alphaJig = np.zeros(self.cE.shape)

        qq = np.zeros((6, self.Ns+1))
        for s in range(1, self.Ns+1):
            qq[:, s] = self.q[s*6:s*6+6].T

        Clalpha = 2*pi
        for s in range(0, len(self.yN)-1):
            self.alphaJig[s] = self.Cl[s] / Clalpha            \
                             - (qq[4, s] + qq[4, s+1]) / 2     \
                             + self.phi[s] - self.collective

        # Compute dihedral angle
        self.di = np.zeros((self.Ns, 1))
        for s in range(0, self.Ns):
            self.di[s] = arctan2(qq[2, s+1] - qq[2, s], self.yN[s+1] - self.yN[s])

        # Compute totals
        # (Note: reshaping is due to numpy vs MATLAB 1D array shapes.. should scrub this)
        self.Ttot   = np.sum(self.fblade.Fz.reshape(-1, 1) * np.cos(self.di)) * self.b * 4
        self.MomRot = np.sum(self.fblade.Fz.reshape(-1, 1) * self.yE.reshape(-1, 1))
        self.Qtot   = np.sum(self.fblade.Q)  * self.b * 4
        Pitot       = np.sum(self.fblade.Pi) * self.b * 4
        Pptot       = np.sum(self.fblade.Pp) * self.b * 4
        self.Ptot   = Pptot + Pitot  # non-covered centre

class Switch(Component):
    """ select the appropriate source for blade force data """

    # inputs
    fblade_initial  = VarTree(Fblade(), iotype='in')
    fblade_updated  = VarTree(Fblade(), iotype='in')

    # outputs
    fblade          = VarTree(Fblade(), iotype='out')

    def __init__(self):
        super(Switch, self).__init__()
        self.initial = True
        self.force_execute = True

    def execute(self):
        if self.initial:
            self.fblade = self.fblade_initial
            self.initial = False
        else:
            self.fblade = self.fblade_updated


class HeliCalc(Assembly):

    def configure(self):


        #self.add('driver', SLSQPdriver())


        self.add('config', AtlasConfiguration())
        self.driver.workflow.add("config")
        self.add('discrete', DiscretizeProperties())
        self.driver.workflow.add("discrete")

        self.add('results', Results())
        self.driver.workflow.add("results")

        self.add('thrust', Thrust())
        self.driver.workflow.add("thrust")

        #self.add('induced', ActuatorDiskInducedVelocity()) # low
        self.add('induced', VortexRingC()) # high
        self.driver.workflow.add("induced")

        self.add('lift_drag', LiftDrag())
        self.driver.workflow.add("lift_drag")
        self.add('spar', SparProperties())
        self.driver.workflow.add("spar")
        self.add('joint', JointSparProperties())
        self.driver.workflow.add("joint")
        self.add('chord', ChordProperties())
        self.driver.workflow.add("chord")
        self.add('quad', QuadSparProperties())
        self.driver.workflow.add("quad")
        self.add('mass', MassProperties())
        self.driver.workflow.add("mass")
        self.add('fem', FEM())
        self.driver.workflow.add("fem")

        self.add('strains', Strains())
        self.driver.workflow.add("strains")
        self.add('failure', Failures())
        self.driver.workflow.add("failure")

        # self.connect('aero.Fblade',         'switch.fblade_initial')
        # self.connect('switch.fblade',       'struc.fblade')
        # # self.connect('struc.q',             'aero2.q')
        # self.connect('aero2.Fblade',        'switch.fblade_updated')

        self.connect('config.b',           'results.b')
        self.connect('config.Ns',  'results.Ns')
        self.connect('discrete.yN',            'results.yN')
        self.connect('discrete.yE',            'results.yE')
        self.connect('discrete.cE',           'results.cE')
        self.connect('config.Cl',           'results.Cl')
        self.connect('fem.q',            'results.q')
        self.connect('lift_drag.phi',          'results.phi')
        self.connect('config.collective',          'results.collective')
        self.connect('lift_drag.Fblade',          'results.fblade')
        self.connect('mass.Mtot',        'results.Mtot')


        self.connect('config.Ns',           'discrete.Ns')
        self.connect('config.ycmax_array',  'discrete.ycmax_array')
        self.connect('config.R',            'discrete.R')
        self.connect('config.c',            'discrete.c_in')
        self.connect('config.Cl',           'discrete.Cl_in')
        self.connect('config.Cm',           'discrete.Cm_in')
        self.connect('config.t',            'discrete.t_in')
        self.connect('config.xtU',          'discrete.xtU_in')
        self.connect('config.xtL',          'discrete.xtL_in')
        self.connect('config.xEA',          'discrete.xEA_in')
        self.connect('config.yWire',        'discrete.yWire')
        self.connect('config.d',            'discrete.d_in')
        self.connect('config.theta',        'discrete.theta_in')
        self.connect('config.nTube',        'discrete.nTube_in')
        self.connect('config.nCap',         'discrete.nCap_in')
        self.connect('config.lBiscuit',     'discrete.lBiscuit_in')


        self.connect('discrete.yN', 'chord.yN')
        self.connect('discrete.yN', 'strains.yN')
        self.connect('discrete.yN', 'failure.yN')
        self.connect('discrete.yN', 'lift_drag.yN')
        self.connect('discrete.yN', 'fem.yN')

        self.connect('discrete.yN', 'thrust.yN')
        self.connect('discrete.yN', 'spar.yN')
        #self.connect('discrete.yN', 'joint.yN')

        self.connect('discrete.lBiscuit', 'failure.lBiscuit')
        self.connect('discrete.lBiscuit', 'spar.lBiscuit')

        #self.connect('discrete.lBiscuit', 'joint.lBiscuit')

        self.connect('discrete.xtU', 'chord.xtU')
        self.connect('discrete.xtU', 'lift_drag.xtU')

        self.connect('discrete.xtL', 'lift_drag.xtL')

        self.connect('discrete.Cm', 'lift_drag.Cm')

        self.connect('discrete.Cl', 'lift_drag.Cl')
        self.connect('discrete.Cl', 'thrust.Cl')

        self.connect('discrete.nCap', 'failure.nCap')

        self.connect('discrete.nCap', 'spar.nCap')
        #self.connect('discrete.nCap', 'joint.nCap')

        self.connect('quad.EIx',          'failure.EIQuad')

        self.connect('joint.EIx',         'failure.EIxJ')
        self.connect('joint.EIz',         'failure.EIzJ')

        self.connect('spar.GJ',  'fem.GJ')
        self.connect('spar.GJ',  'failure.GJQuad')

        self.connect('discrete.nTube', 'failure.nTube')

        self.connect('discrete.nTube', 'spar.nTube')
        #self.connect('discrete.nTube', 'joint.nTube')

        self.connect('discrete.theta', 'failure.theta')

        self.connect('discrete.theta', 'spar.theta')
        #self.connect('discrete.theta', 'joint.theta')

        self.connect('discrete.xEA', 'mass.xEA')
        self.connect('discrete.xEA', 'fem.xEA')

        self.connect('discrete.d', 'strains.d')
        self.connect('discrete.d', 'chord.d')
        self.connect('discrete.d', 'failure.d')
        self.connect('discrete.d', 'lift_drag.d')
        self.connect('discrete.d', 'spar.d')

        #self.connect('discrete.d', 'joint.d')

        self.connect('spar.EIx',     'fem.EIx')
        self.connect('spar.EIz',     'fem.EIz')
        self.connect('spar.EA',      'fem.EA')
        self.connect('spar.mSpar',   'failure.mSpar')
        self.connect('spar.mSpar',   'fem.mSpar')
        self.connect('spar.mSpar',     'mass.mSpar')

        self.connect("mass.xCG" , "fem.xCG")
        self.connect("config.presLoad" , "fem.presLoad")
        self.connect("config.hQuad" , "failure.hQuad")
        self.connect("config.hQuad" , "quad.hQuad")
        self.connect("config.RQuad" , "failure.RQuad")
        self.connect("config.RQuad" , "quad.RQuad")
        self.connect("config.GWing" , "chord.GWing")
        self.connect("config.lBiscuitQuad" , "failure.lBiscuitQuad")
        self.connect("config.lBiscuitQuad" , "quad.lBiscuitQuad")
        self.connect("config.mElseR" , "mass.mElseR")

        self.connect("config.dr" , "induced.dr")
        self.connect("config.dr" , "lift_drag.dr")
        self.connect("config.dr" , "thrust.dr")

        self.connect("config.tWire" , "failure.tWire")
        self.connect("config.tWire" , "lift_drag.tWire")
        self.connect("config.tWire" , "mass.tWire")

        self.connect("config.mElseRotor" , "failure.mElseRotor")
        self.connect("config.mElseRotor" , "mass.mElseRotor")

        self.connect("config.zWire" , "failure.zWire")
        self.connect("config.zWire" , "lift_drag.zWire")
        self.connect("config.zWire" , "mass.zWire")
        self.connect("config.zWire" , "fem.zWire")

        self.connect("config.h" , "induced.h")
        self.connect("thrust.dT" , "induced.dT")
        self.connect("config.TWire" , "failure.TWire")
        self.connect("config.TWire" , "fem.TWire")
        self.connect("config.yWire" , "failure.yWire")
        self.connect("config.yWire" , "lift_drag.yWire")
        self.connect("config.yWire" , "mass.yWire")
        self.connect("config.yWire" , "fem.yWire")
        self.connect("config.mElseCentre" , "mass.mElseCentre")
        self.connect("config.thetaQuad" , "failure.thetaQuad")
        self.connect("config.thetaQuad" , "quad.thetaQuad")
        self.connect("chord.xCGChord" , "mass.xCGChord")
        self.connect("discrete.cE" , "lift_drag.c")
        self.connect("discrete.cE" , "thrust.c")
        self.connect("fem.k" , "strains.k")

        self.connect("lift_drag.Fblade" , "failure.Fblade")
        self.connect("lift_drag.Fblade" , "fem.Fblade")
        self.connect("config.TEtension" , "failure.TEtension")
        self.connect("config.nTubeQuad" , "failure.nTubeQuad")
        self.connect("config.nTubeQuad" , "quad.nTubeQuad")
        self.connect("quad.mQuad" , "mass.mQuad")
        self.connect("chord.mChord" , "failure.mChord")
        self.connect("chord.mChord" , "mass.mChord")
        self.connect("chord.mChord" , "fem.mChord")
        self.connect("config.Ns" , "induced.Ns")
        self.connect("config.Ns" , "lift_drag.Ns")
        self.connect("config.Ns" , "thrust.Ns")
        self.connect("discrete.cE" , "chord.cE")
        self.connect("discrete.cE" , "fem.cE")
        self.connect("induced.vi" , "lift_drag.vi")
        self.connect("config.ycmax" , "thrust.ycmax")
        self.connect("config.ycmax" , "mass.ycmax")
        self.connect("strains.Finternal" , "failure.Finternal")

        self.connect("config.R" , "induced.R")
        self.connect("config.R" , "mass.R")
        self.connect("config.b" , "failure.b")
        self.connect("config.b" , "induced.b")
        self.connect("config.b" , "mass.b")
        self.connect("config.r" , "induced.r")
        self.connect("config.r" , "lift_drag.r")
        self.connect("config.r" , "thrust.r")
        self.connect("config.vc" , "induced.vc")
        self.connect("config.vc" , "lift_drag.vc")

        self.connect("strains.strain" , "failure.strain")
        self.connect("config.vw" , "lift_drag.vw")
        self.connect("fem.F" , "strains.F")
        self.connect("config.mPilot" , "mass.mPilot")
        self.connect("config.Omega" , "lift_drag.Omega")
        self.connect("config.Omega" , "thrust.Omega")
        self.connect("config.Jprop" , "joint.Jprop")
        self.connect("thrust.chordFrac" , "lift_drag.chordFrac")
        self.connect("config.visc" , "lift_drag.visc")
        self.connect("config.dQuad" , "failure.dQuad")
        self.connect("config.dQuad" , "quad.dQuad")
        self.connect("config.rho" , "induced.rho")
        self.connect("config.rho" , "lift_drag.rho")
        self.connect("config.rho" , "thrust.rho")
        self.connect("config.CFRPType" , "spar.CFRPType")
        self.connect("config.CFRPType" , "quad.CFRPType")
        self.connect("config.CFRPType" , "joint.CFRPType")
        self.connect("fem.q" , "strains.q")

        self.connect("config.Load" , "fem.Load")
        self.connect("config.wingWarp" , "fem.wingWarp")
        self.connect("config.Cover" , "mass.Cover")
        self.connect("config.WireType" , "mass.WireType")

        self.connect("config.CFRPType" , "failure.CFRPType")
        self.connect("config.WireType" , "failure.WireType")
        self.connect("config.Quad" , "mass.Quad")

        self.connect('discrete.yN',        'induced.yN')
        self.connect('config.Omega',     'induced.Omega')
        self.connect('config.anhedral',  'induced.anhedral')

        #self.connect('fem.q',         'induced.q')

        fpi = self.add('fpi', FixedPointIterator())
        fpi.max_iteration = 10
        fpi.tolerance = 1e-10
        fpi.add_parameter('induced.q', low=-1e999, high=1e999)
        fpi.add_constraint('induced.q = fem.q')
        self.driver.workflow.add("fpi")


if __name__ == "__main__":

    top = HeliCalc()
    #from openmdao.util.dotgraph import plot_graphs
    #plot_graphs(top)

    top.run()
    print np.linalg.norm(top.fem.q - top.induced.q)
    print top.config.Omega
    print top.results.Ptot
