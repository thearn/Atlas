
from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Float, Array

from Atlas import AtlasConfiguration, Thrust, ActuatorDiskInducedVelocity, LiftDrag


class Aero(Assembly):
    # inputs
    yN    = Array(iotype="in",   desc='node locations')

    ycmax = Float(iotype="in")
    b     = Float(iotype="in", desc="number of blades")
    h     = Float(iotype="in", desc="height of rotor")

    rho   = Float(0.0, iotype='in', desc='air density')
    visc  = Float(0.0, iotype='in', desc='air viscosity')
    vw    = Float(0.0, iotype='in', desc='wind velocity')
    vc    = Float(0.0, iotype='in', desc='vertical velocity')
    Omega = Float(0.0, iotype='in', desc='rotor angular velocity')

    c     = Array(iotype='in', desc='chord distribution')
    Cl    = Array(iotype='in', desc='lift coefficient distribution')
    d     = Array(iotype='in', desc='spar diameter distribution')

    yWire = Array(iotype='in', desc='location of wire attachment along span')
    zWire = Float(iotype='in', desc='depth of wire attachement')
    tWire = Float(iotype='in', desc='thickness of wire')

    Cm    = Array(iotype='in', desc='')
    xtU   = Array(iotype='in', desc='fraction of laminar flow on the upper surface')
    xtL   = Array(iotype='in', desc='fraction of laminar flow on the lower surface')

    def configure(self):
        self.add('config', AtlasConfiguration())
        self.connect('yN', 'config.yN')

        self.add('thrust', Thrust())
        self.connect('yN',        'thrust.yN')
        self.connect('config.Ns', 'thrust.Ns')
        self.connect('config.dr', 'thrust.dr')
        self.connect('config.r',  'thrust.r')
        self.connect('ycmax',     'thrust.ycmax')
        self.connect('Cl',        'thrust.Cl')
        self.connect('c',         'thrust.c')
        self.connect('rho',       'thrust.rho')
        self.connect('Omega',     'thrust.Omega')

        self.add('adiv', ActuatorDiskInducedVelocity())
        self.connect('yN',        'adiv.yN')
        self.connect('config.Ns', 'adiv.Ns')
        self.connect('config.r',  'adiv.r')
        self.connect('config.dr', 'adiv.dr')
        self.connect('config.R',  'adiv.R')
        self.connect('b',         'adiv.b')
        self.connect('h',         'adiv.h')
        self.connect('vc',        'adiv.vc')
        self.connect('rho',       'adiv.rho')
        self.connect('thrust.dT', 'adiv.dT')

        self.add('lift_drag', LiftDrag())
        self.connect('yN',        'lift_drag.yN')
        self.connect('config.Ns', 'lift_drag.Ns')
        self.connect('config.dr', 'lift_drag.dr')
        self.connect('config.r',  'lift_drag.r')
        self.connect('rho',       'lift_drag.rho')
        self.connect('visc',      'lift_drag.visc')
        self.connect('vw',        'lift_drag.vw')
        self.connect('vc',        'lift_drag.vc')
        self.connect('Omega',     'lift_drag.Omega')
        self.connect('c',         'lift_drag.c')
        self.connect('Cl',        'lift_drag.Cl')
        self.connect('d',         'lift_drag.d')
        self.connect('yWire',     'lift_drag.yWire')
        self.connect('zWire',     'lift_drag.zWire')
        self.connect('tWire',     'lift_drag.tWire')
        self.connect('Cm',        'lift_drag.Cm')
        self.connect('xtU',       'lift_drag.xtU')
        self.connect('xtL',       'lift_drag.xtL')
        self.connect('adiv.vi',   'lift_drag.vi')
        self.connect('thrust.chordFrac', 'lift_drag.chordFrac')

        self.create_passthrough('adiv.vi')
        self.create_passthrough('lift_drag.phi')
        self.create_passthrough('lift_drag.Re')
        self.create_passthrough('lift_drag.Cd')
        self.create_passthrough('lift_drag.Fblade')

        self.driver.workflow.add('config')
        self.driver.workflow.add('thrust')
        self.driver.workflow.add('adiv')
        self.driver.workflow.add('lift_drag')


if __name__ == '__main__':
    comp = Aero()

    # set inputs
    comp.yN    = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    comp.ycmax = 1.4656
    comp.b     = 2
    comp.h     = 1.5000

    comp.rho   = 1.1800
    comp.visc  = 1.7800e-005
    comp.vw    = 0
    comp.vc    = 0
    comp.Omega = 1.0367

    comp.c     = [0.2729, 1.3903, 1.1757, 1.0176, 0.8818, 0.7602, 0.6507, 0.5528, 0.4666, 0.3925]
    comp.Cl    = [1.5000, 1.4987, 1.4604, 1.4239, 1.3940, 1.3642, 1.3344, 1.3046, 1.2747, 0.8299]
    comp.d     = [0.0843, 0.0780, 0.0718, 0.0655, 0.0592, 0.0530, 0.0477, 0.0431, 0.0384, 0.0338]

    comp.yWire = [5.8852]
    comp.zWire = 1
    comp.tWire = 0.0016

    comp.Cm    = [-0.1500, -0.1494, -0.1330, -0.1200, -0.1200, -0.1200, -0.1200, -0.1200, -0.1200, -0.1200]
    comp.xtU   = [0.0500, 0.1500, 0.1500, 0.1500, 0.1500, 0.1500, 0.1500, 0.1500, 0.1500, 0.1500]
    comp.xtL   = [0.0500, 0.3000, 0.3000, 0.3000, 0.3000, 0.3000, 0.3000, 0.3000, 0.3000, 0.300]

    # run
    comp.run()

    # check outputs
    tol = 5e-4

    print 'comp.vi:\n', comp.vi
