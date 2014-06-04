import numpy as np
from numpy import pi
from util import arctan2
import warnings
warnings.simplefilter("ignore", np.ComplexWarning)
from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.main.datatypes.api import Int, Float, Array, Str, Enum, VarTree

from openmdao.lib.drivers.api import SLSQPdriver
from openmdao.lib.drivers.api import FixedPointIterator
from pyopt_driver import pyopt_driver

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


class HeliCalc(Assembly):

    def configure(self):


        #optimizer = self.add('driver', SLSQPdriver())
        optimizer = self.add("driver", pyopt_driver.pyOptDriver())
        optimizer.optimizer = "SNOPT"

        fpi = self.add('fpi', FixedPointIterator())

        self.add('config', AtlasConfiguration())
        optimizer.workflow.add("config")

        self.add('discrete', DiscretizeProperties())
        optimizer.workflow.add("discrete")

        self.add('results', Results())
        optimizer.workflow.add("results")

        self.add('thrust', Thrust())
        optimizer.workflow.add("thrust")

        self.add('induced', VortexRingC())
        optimizer.workflow.add("induced")

        self.add('lift_drag', LiftDrag())
        optimizer.workflow.add("lift_drag")

        self.add('spar', SparProperties())
        optimizer.workflow.add("spar")

        self.add('joint', JointSparProperties())
        optimizer.workflow.add("joint")

        self.add('chord', ChordProperties())
        optimizer.workflow.add("chord")

        self.add('quad', QuadSparProperties())
        optimizer.workflow.add("quad")

        self.add('mass', MassProperties())
        optimizer.workflow.add("mass")

        self.add('fem', FEM())
        optimizer.workflow.add("fem")

        self.add('strains', Strains())
        optimizer.workflow.add("strains")

        self.add('failure', Failures())
        optimizer.workflow.add("failure")

        # Config
        self.connect("config.rho" , ["induced.rho", "lift_drag.rho",
                                     "thrust.rho"])
        self.connect("config.Ns" , ["induced.Ns", "lift_drag.Ns", "thrust.Ns",
                                    "results.Ns"])
        self.connect("config.lBiscuitQuad" , ["failure.lBiscuitQuad",
                                              "quad.lBiscuitQuad"])
        self.connect("config.tWire" , ["failure.tWire", "lift_drag.tWire",
                                       "mass.tWire"])
        self.connect("config.mElseRotor" , ["failure.mElseRotor",
                                            "mass.mElseRotor"])
        self.connect("config.zWire" , ["failure.zWire", "lift_drag.zWire",
                                       "mass.zWire", "fem.zWire"])
        self.connect("config.nTubeQuad" , ["failure.nTubeQuad",
                                           "quad.nTubeQuad"])
        self.connect("config.yWire" , ["failure.yWire", "lift_drag.yWire",
                                       "mass.yWire", "fem.yWire"])
        self.connect("config.thetaQuad" , ["failure.thetaQuad",
                                           "quad.thetaQuad"])
        self.connect('config.b', 'results.b')
        self.connect('config.collective', 'results.collective')
        self.connect('config.Ns', 'discrete.Ns')
        self.connect('config.ycmax_array',  'discrete.ycmax_array')
        self.connect('config.R', 'discrete.R')
        self.connect('config.c', 'discrete.c_in')
        self.connect('config.Cl', 'discrete.Cl_in')
        self.connect('config.Cm', 'discrete.Cm_in')
        self.connect('config.t', 'discrete.t_in')
        self.connect('config.xtU', 'discrete.xtU_in')
        self.connect('config.xtL', 'discrete.xtL_in')
        self.connect('config.xEA', 'discrete.xEA_in')
        self.connect('config.yWire', 'discrete.yWire')
        self.connect('config.d', 'discrete.d_in')
        self.connect('config.theta', 'discrete.theta_in')
        self.connect('config.nTube', 'discrete.nTube_in')
        self.connect('config.nCap', 'discrete.nCap_in')
        self.connect('config.lBiscuit', 'discrete.lBiscuit_in')
        self.connect("config.presLoad" , "fem.presLoad")
        self.connect("config.hQuad" , ["failure.hQuad", "quad.hQuad"])
        self.connect("config.RQuad" , ["failure.RQuad", "quad.RQuad"])
        self.connect("config.GWing" , "chord.GWing")
        self.connect("config.mElseR" , "mass.mElseR")
        self.connect("config.dr" , ["induced.dr", "lift_drag.dr", "thrust.dr"])
        self.connect("config.h" , "induced.h")
        self.connect("config.TWire" , ["failure.TWire", "fem.TWire"])
        self.connect("config.mElseCentre" , "mass.mElseCentre")
        self.connect("config.TEtension" , "failure.TEtension")
        self.connect("config.R" , ["induced.R", "mass.R"])
        self.connect("config.b" , ["failure.b", "induced.b", "mass.b"])
        self.connect("config.r" , ["induced.r", "lift_drag.r", "thrust.r"])
        self.connect("config.vc" , ["induced.vc", "lift_drag.vc"])
        self.connect("config.ycmax" , ["thrust.ycmax", "mass.ycmax"])
        self.connect("config.mPilot" , "mass.mPilot")
        self.connect("config.Omega" , ["lift_drag.Omega", "thrust.Omega"])
        self.connect("config.Jprop" , "joint.Jprop")
        self.connect("config.visc" , "lift_drag.visc")
        self.connect("config.dQuad" , ["failure.dQuad", "quad.dQuad"])
        self.connect("config.vw" , "lift_drag.vw")


        # Config - FLAGS
        self.connect("config.Load" , "fem.Load")
        self.connect("config.wingWarp" , "fem.wingWarp")
        self.connect("config.Cover" , "mass.Cover")
        self.connect("config.WireType" , "mass.WireType")
        self.connect("config.WireType" , "failure.WireType")
        self.connect("config.Quad" , "mass.Quad")
        self.connect('config.Omega',     'induced.Omega')
        self.connect('config.anhedral',  'induced.anhedral')
        self.connect("config.CFRPType" , ["spar.CFRPType", "quad.CFRPType",
                                          "joint.CFRPType", "failure.CFRPType"])

        # Discrete
        self.connect('discrete.yN', ['chord.yN', 'strains.yN', 'failure.yN',
                                     'lift_drag.yN', 'fem.yN', 'thrust.yN',
                                     'spar.yN', 'induced.yN', 'results.yN'])
        self.connect('discrete.d', ['strains.d', 'chord.d', 'failure.d',
                                    'lift_drag.d', 'spar.d'])
        self.connect("discrete.cE" , ["chord.cE", "fem.cE", "lift_drag.c",
                                      "thrust.c", "results.cE"])
        self.connect('discrete.lBiscuit', ['failure.lBiscuit', 'spar.lBiscuit'])
        self.connect('discrete.yE', 'results.yE')
        self.connect('discrete.xtU', ['chord.xtU', 'lift_drag.xtU'])
        self.connect('discrete.xtL', 'lift_drag.xtL')
        self.connect('discrete.Cm', 'lift_drag.Cm')
        self.connect('discrete.Cl', ['lift_drag.Cl', 'thrust.Cl', 'results.Cl'])
        self.connect('discrete.nCap', ['failure.nCap', 'spar.nCap'])
        self.connect('discrete.nTube', ['failure.nTube', 'spar.nTube'])
        self.connect('discrete.theta', ['failure.theta', 'spar.theta'])
        self.connect('discrete.xEA', ['mass.xEA', 'fem.xEA'])

        # Quad
        self.connect('quad.EIx', 'failure.EIQuad')
        self.connect("quad.mQuad" , "mass.mQuad")

        # Joint
        self.connect('joint.EIx', 'failure.EIxJ')
        self.connect('joint.EIz', 'failure.EIzJ')

        # Spar
        self.connect('spar.GJ',  ['fem.GJ', 'failure.GJQuad'])
        self.connect('spar.EIx',     'fem.EIx')
        self.connect('spar.EIz',     'fem.EIz')
        self.connect('spar.EA',      'fem.EA')
        self.connect('spar.mSpar',   ['failure.mSpar', 'fem.mSpar',
                                      'mass.mSpar'])

        # Mass
        self.connect("mass.xCG" , "fem.xCG")
        self.connect('mass.Mtot',        'results.Mtot')

        # Thrust
        self.connect("thrust.dT" , "induced.dT")
        self.connect("thrust.chordFrac" , "lift_drag.chordFrac")

        # Chord
        self.connect("chord.xCGChord" , "mass.xCGChord")
        self.connect("chord.mChord" , ["failure.mChord", "mass.mChord",
                                       "fem.mChord"])

        # FEM
        self.connect("fem.k" , "strains.k")
        self.connect("fem.F" , "strains.F")
        self.connect("fem.q" , ["strains.q", 'results.q'])

        # LiftDrag
        self.connect("lift_drag.Fblade" , ["failure.Fblade", "fem.Fblade",
                                           "results.fblade"])
        self.connect('lift_drag.phi',          'results.phi')

        # Induced
        self.connect("induced.vi" , "lift_drag.vi")

        # Strains
        self.connect("strains.Finternal" , "failure.Finternal")
        self.connect("strains.strain" , "failure.strain")


        fpi.max_iteration = 10
        fpi.tolerance = 1e-10
        fpi.add_parameter('induced.q', low=-1e999, high=1e999)
        fpi.add_constraint('induced.q = fem.q')
        optimizer.workflow.add("fpi")

        optimizer.add_parameter("config.Omega", low=0.15*2*pi, high=0.25*2*pi)
        optimizer.add_objective("results.Ptot")
        optimizer.add_constraint('results.Mtot * 9.8 - results.Ttot <= 0')


if __name__ == "__main__":
    from openmdao.util.dotgraph import plot_graphs
    top = HeliCalc()

    #top.driver.gradient_options.fd_form = 'complex_step'
    plot_graphs(top)
    print
    top.driver.run_iteration()
    print "BASELINE:"
    print "aero structural residual:", np.linalg.norm(top.fem.q - top.induced.q)
    print "lift weight residual:", top.results.Mtot*9.8-top.results.Ttot
    print "Omega:", top.config.Omega
    print "Ptot:", top.results.Ptot
    top.run()
    print
    print "OPTIMIZED:"
    print "aero structural residual:", np.linalg.norm(top.fem.q - top.induced.q)
    print "lift weight residual:", top.results.Mtot*9.8-top.results.Ttot
    print "Omega:", top.config.Omega
    print "Ptot:", top.results.Ptot
