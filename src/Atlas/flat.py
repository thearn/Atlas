import numpy as np
from numpy import pi
from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.main.datatypes.api import Int, Float, Array, Str, Enum, VarTree

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
        self.add('config', AtlasConfiguration())
        self.add('discrete', DiscretizeProperties())

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

        self.add('thrust', Thrust())

        self.add('induced', ActuatorDiskInducedVelocity()) # low
        #self.add('induced', VortexRing()) # high

        self.add('lift_drag', LiftDrag())
        self.add('spar', SparProperties())
        self.add('joint', JointSparProperties())
        self.add('chord', ChordProperties())
        self.add('quad', QuadSparProperties())
        self.add('mass', MassProperties())
        self.add('fem', FEM())

        self.add('strains', Strains())
        self.add('failure', Failures())

        self.connect('discrete.yN', 'chord.yN')
        self.connect('discrete.yN', 'strains.yN')
        self.connect('discrete.yN', 'failure.yN')
        self.connect('discrete.yN', 'lift_drag.yN')
        self.connect('discrete.yN', 'fem.yN')
        self.connect('discrete.yN', 'quad.yN')
        self.connect('discrete.yN', 'thrust.yN')
        self.connect('discrete.yN', 'spar.yN')
        self.connect('discrete.yN', 'joint.yN')

        self.connect('discrete.lBiscuit', 'failure.lBiscuit')
        self.connect('discrete.lBiscuit', 'spar.lBiscuit')
        self.connect('discrete.lBiscuit', 'quad.lBiscuit')
        self.connect('discrete.lBiscuit', 'joint.lBiscuit')

        self.connect('discrete.xtU', 'chord.xtU')
        self.connect('discrete.xtU', 'lift_drag.xtU')

        self.connect('discrete.xtL', 'lift_drag.xtL')

        self.connect('discrete.Cm', 'lift_drag.Cm')

        self.connect('discrete.Cl', 'lift_drag.Cl')
        self.connect('discrete.Cl', 'thrust.Cl')

        self.connect('discrete.nCap', 'failure.nCap')
        self.connect('discrete.nCap', 'quad.nCap')
        self.connect('discrete.nCap', 'spar.nCap')
        self.connect('discrete.nCap', 'joint.nCap')

        self.connect('quad.EIx',          'failure.EIQuad')

        self.connect('joint.EIx',         'failure.EIxJ')
        self.connect('joint.EIz',         'failure.EIzJ')

        self.connect('spar.GJ',  'fem.GJ')
        self.connect('spar.GJ',  'failure.GJQuad')

        self.connect('discrete.nTube', 'failure.nTube')
        self.connect('discrete.nTube', 'quad.nTube')
        self.connect('discrete.nTube', 'spar.nTube')
        self.connect('discrete.nTube', 'joint.nTube')

        self.connect('discrete.theta', 'failure.theta')
        self.connect('discrete.theta', 'quad.theta')
        self.connect('discrete.theta', 'spar.theta')
        self.connect('discrete.theta', 'joint.theta')

        self.connect('discrete.xEA', 'mass.xEA')
        self.connect('discrete.xEA', 'fem.xEA')

        self.connect('discrete.d', 'strains.d')
        self.connect('discrete.d', 'chord.d')
        self.connect('discrete.d', 'failure.d')
        self.connect('discrete.d', 'lift_drag.d')
        self.connect('discrete.d', 'spar.d')
        self.connect('discrete.d', 'quad.d')
        self.connect('discrete.d', 'joint.d')

        self.connect('spar.EIx',     'fem.EIx')
        self.connect('spar.EIz',     'fem.EIz')
        self.connect('spar.EA',      'fem.EA')
        self.connect('spar.mSpar',        'failure.mSpar')
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
        self.connect("config.c" , "lift_drag.c")
        self.connect("config.c" , "thrust.c")
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


    def make_connections(self):
        """
        Collects the names of all input and output variables for all
        components within the assembly (drivers excluded).
        Then establishes connections between
        any output variable and input variable that has the same name so
        long as the variable name does not exist as an output to more than
        a single component (so excludes default outputs).
        """

        inputs, outputs = {}, {}
        for compname in self.list_components():

            comp_inputs = self.get(compname).list_inputs()

            for input_name in comp_inputs:
                if input_name not in inputs:
                    inputs[input_name] = [compname]
                else:
                    inputs[input_name].append(compname)

            comp_outputs = self.get(compname).list_outputs()

            for output_name in comp_outputs:
                if output_name not in outputs:
                    outputs[output_name] = [compname]
                else:
                    outputs[output_name].append(compname)

        io = {}

        for compname in self.list_components():
            comp_io = self.get(compname).list_inputs() + \
                self.get(compname).list_outputs()

            for io_name in comp_io:
                if io_name not in ['derivative_exec_count',
                                   'directory',
                                   'exec_count',
                                   'external_vars',
                                   'fixed_external_vars',
                                   'force_execute',
                                   'init_state_var',
                                   'itername',
                                   'state_var',
                                   ]:

                    if io_name not in io:
                        io[io_name] = [compname]
                    else:
                        io[io_name].append(compname)

        for io_name in sorted( io.keys() ):
            for compname in io[ io_name ] :
                p = '.'.join( [compname, io_name] )
                metadata = self.get_metadata( p )

        assym_level = self.list_inputs()
        assym_level.remove('directory')
        assym_level.remove('force_execute')

        for var in assym_level:
            outputs[var] = ['']

        for varname in outputs.keys():
            comps = outputs[varname]
            if len(comps) > 1:
                continue
            if comps[0]:
                frompath = '.'.join([comps[0], varname])
            else:
                frompath = varname

            if varname in inputs:
                for compname in inputs[varname]:
                    topath = '.'.join([compname, varname])
                    try:
                        self.connect(frompath, topath)
                        print 'self.connect("'+frompath+'" , "'+topath+'")'
                    except:
                        pass

    def get_unconnected_inputs(self):
        unconnected_inputs = []
        connected_inputs = [i[1] for i in self.get_dataflow()['connections']]
        defaults = ['itername', 'force_execute', 'directory', 'exec_count',
                    'derivative_exec_count', 'fixed_external_vars']
        for compname in self.list_components():
            if compname == "driver":
                continue
            comp = self.get(compname)
            for var in comp.list_inputs():
                if var in defaults:
                    continue
                fullname = '.'.join([compname, var])
                if fullname not in connected_inputs:
                    unconnected_inputs.append(fullname)
        return unconnected_inputs

if __name__ == "__main__":
    # failure.flags
    # mass.flags
    # fem.flags


    top = HeliCalc()
    #top.make_connections()
    # from openmdao.util.dotgraph import plot_graphs
    # plot_graphs(top)

    # for i in top.get_unconnected_inputs():
    #     print i

    top.run()
