from openmdao.main.api import Assembly

from Atlas import AtlasConfiguration, DiscretizeProperties, Aero, Structures

# see: HeliCalc.m


class AeroStructural(Assembly):
    """
    Performs an aerodynamic and structural computation on a single
    configuration given the full set of design parameters. The aerodynamic
    computation returns the thrust, torque, power, and induced velocity. The
    structural computation first computes the mass of the helicopter based on
    the structural description of the spars and chord lengths. It then
    computes the deformation of the spars, the strains, and the resulting
    factor of safety for each of the failure modes.
    """

    def configure(self):
        self.add('config', AtlasConfiguration())

        self.add('discrete', DiscretizeProperties())
        self.connect('config.Ns',    'discrete.Ns')
        self.connect('config.ycmax', 'discrete.ycmax')
        self.connect('config.R',     'discrete.R')
        self.connect('config.c',     'discrete.c_in')
        self.connect('config.Cl',    'discrete.Cl_in')
        self.connect('config.Cm',    'discrete.Cm_in')
        self.connect('config.t',     'discrete.t_in')
        self.connect('config.xtU',   'discrete.xtU_in')
        self.connect('config.xtL',   'discrete.xtL_in')
        self.connect('config.xEA',   'discrete.xEA_in')
        self.connect('config.yWire', 'discrete.yWire')
        self.connect('config.d',     'discrete.d_in')
        self.connect('config.theta', 'discrete.theta_in')
        self.connect('config.nTube', 'discrete.nTube_in')
        self.connect('config.nCap',  'discrete.nCap_in')
        self.connect('config.lBiscuit', 'discrete.lBiscuit_in')

        self.add('aero', Aero())
        self.connect('config.b',        'aero.b')
        self.connect('config.R',        'aero.R')
        self.connect('config.Ns',       'aero.Ns')
        self.connect('discrete.yN',     'aero.yN')
        self.connect('config.dr',       'aero.dr')
        self.connect('config.r',        'aero.r')
        self.connect('config.h',        'aero.h')
        self.connect('config.ycmax[0]', 'aero.ycmax')
        self.connect('config.rho',      'aero.rho')
        self.connect('config.visc',     'aero.visc')
        self.connect('config.vw',       'aero.vw')
        self.connect('config.vc',       'aero.vc')
        self.connect('config.Omega',    'aero.Omega')
        self.connect('discrete.cE',     'aero.c')
        self.connect('discrete.Cl',     'aero.Cl')
        self.connect('discrete.d',      'aero.d')
        self.connect('config.yWire',    'aero.yWire')
        self.connect('config.zWire',    'aero.zWire')
        self.connect('config.tWire',    'aero.tWire')
        self.connect('discrete.Cm',     'aero.Cm')
        self.connect('discrete.xtU',    'aero.xtU')
        self.connect('discrete.xtL',    'aero.xtL')

        self.add('structures', Structures())
        self.connect('config.flags',        'structures.flags')
        self.connect('discrete.yN',         'structures.yN')
        self.connect('discrete.d',          'structures.d')
        self.connect('discrete.theta',      'structures.theta')
        self.connect('discrete.nTube',      'structures.nTube')
        self.connect('discrete.nCap',       'structures.nCap')
        self.connect('discrete.lBiscuit',   'structures.lBiscuit')
        self.connect('config.Jprop',        'structures.Jprop')
        self.connect('config.b',            'structures.b')
        self.connect('discrete.cE',         'structures.cE')
        self.connect('discrete.xEA',        'structures.xEA')
        self.connect('discrete.xtU',        'structures.xtU')
        self.connect('config.dQuad',        'structures.dQuad')
        self.connect('config.thetaQuad',    'structures.thetaQuad')
        self.connect('config.nTubeQuad',    'structures.nTubeQuad')
        self.connect('config.lBiscuitQuad', 'structures.lBiscuitQuad')
        self.connect('config.RQuad',        'structures.RQuad')
        self.connect('config.hQuad',        'structures.hQuad')
        self.connect('config.ycmax[0]',     'structures.ycmax')
        self.connect('config.yWire',        'structures.yWire')
        self.connect('config.zWire',        'structures.zWire')
        self.connect('config.tWire',        'structures.tWire')
        self.connect('config.TWire',        'structures.TWire')
        self.connect('config.TEtension',    'structures.TEtension')
        self.connect('config.mElseRotor',   'structures.mElseRotor')
        self.connect('config.mElseCentre',  'structures.mElseCentre')
        self.connect('config.mElseR',       'structures.mElseR')
        self.connect('config.R',            'structures.R')
        self.connect('config.mPilot',       'structures.mPilot')
        self.connect('aero.Fblade',         'structures.fblade')
        self.connect('config.presLoad',     'structures.presLoad')

        self.driver.workflow.add('config')
        self.driver.workflow.add('discrete')
        self.driver.workflow.add('aero')
        self.driver.workflow.add('structures')


if __name__ == "__main__":
    top = AeroStructural()
    top.run()
