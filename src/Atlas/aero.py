
from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Int, Float, Array

from Atlas import Thrust, ActuatorDiskInducedVelocity, LiftDrag, VortexRing


class Aero(Assembly):
    # inputs
    b     = Int(iotype="in",   desc="number of blades")
    R     = Float(iotype='in', desc='rotor radius')
    Ns    = Int(iotype='in',   desc='number of elements')
    yN    = Array(iotype="in", desc='node locations')
    dr    = Array(iotype="in", desc="length of each element")
    r     = Array(iotype="in", desc="radial location of each element")
    h     = Float(iotype="in", desc="height of rotor")

    ycmax = Float(iotype="in")

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
        self.add('thrust', Thrust())
        self.connect('Ns',        'thrust.Ns')
        self.connect('yN',        'thrust.yN')
        self.connect('dr',        'thrust.dr')
        self.connect('r',         'thrust.r')
        self.connect('ycmax',     'thrust.ycmax')
        self.connect('Cl',        'thrust.Cl')
        self.connect('c',         'thrust.c')
        self.connect('rho',       'thrust.rho')
        self.connect('Omega',     'thrust.Omega')

        self.add('induced', ActuatorDiskInducedVelocity())
        self.connect('Ns',        'induced.Ns')
        self.connect('r',         'induced.r')
        self.connect('dr',        'induced.dr')
        self.connect('R',         'induced.R')
        self.connect('b',         'induced.b')
        self.connect('h',         'induced.h')
        self.connect('vc',        'induced.vc')
        self.connect('rho',       'induced.rho')
        self.connect('thrust.dT', 'induced.dT')

        self.add('lift_drag', LiftDrag())
        self.connect('yN',         'lift_drag.yN')
        self.connect('Ns',         'lift_drag.Ns')
        self.connect('dr',         'lift_drag.dr')
        self.connect('r',          'lift_drag.r')
        self.connect('rho',        'lift_drag.rho')
        self.connect('visc',       'lift_drag.visc')
        self.connect('vw',         'lift_drag.vw')
        self.connect('vc',         'lift_drag.vc')
        self.connect('Omega',      'lift_drag.Omega')
        self.connect('c',          'lift_drag.c')
        self.connect('Cl',         'lift_drag.Cl')
        self.connect('d',          'lift_drag.d')
        self.connect('yWire',      'lift_drag.yWire')
        self.connect('zWire',      'lift_drag.zWire')
        self.connect('tWire',      'lift_drag.tWire')
        self.connect('Cm',         'lift_drag.Cm')
        self.connect('xtU',        'lift_drag.xtU')
        self.connect('xtL',        'lift_drag.xtL')
        self.connect('induced.vi', 'lift_drag.vi')
        self.connect('thrust.chordFrac', 'lift_drag.chordFrac')

        self.create_passthrough('induced.vi')
        self.create_passthrough('lift_drag.phi')
        self.create_passthrough('lift_drag.Re')
        self.create_passthrough('lift_drag.Cd')
        self.create_passthrough('lift_drag.Fblade')

        self.driver.workflow.add('thrust')
        self.driver.workflow.add('induced')
        self.driver.workflow.add('lift_drag')


class Aero2(Assembly):
    # inputs
    b     = Int(iotype="in",   desc="number of blades")
    R     = Float(iotype='in', desc='rotor radius')
    Ns    = Int(iotype='in',   desc='number of elements')
    yN    = Array(iotype="in", desc='node locations')
    dr    = Array(iotype="in", desc="length of each element")
    r     = Array(iotype="in", desc="radial location of each element")
    h     = Float(iotype="in", desc="height of rotor")

    ycmax = Float(iotype="in")

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

    q     = Array(iotype='in', desc='deformation')

    anhedral = Float(iotype='in')

    def configure(self):
        self.add('thrust', Thrust())
        self.connect('Ns',        'thrust.Ns')
        self.connect('yN',        'thrust.yN')
        self.connect('dr',        'thrust.dr')
        self.connect('r',         'thrust.r')
        self.connect('ycmax',     'thrust.ycmax')
        self.connect('Cl',        'thrust.Cl')
        self.connect('c',         'thrust.c')
        self.connect('rho',       'thrust.rho')
        self.connect('Omega',     'thrust.Omega')

        self.add('induced', VortexRing())
        self.connect('yN',        'induced.yN')
        self.connect('Ns',        'induced.Ns')
        self.connect('b',         'induced.b')
        self.connect('h',         'induced.h')
        self.connect('vc',        'induced.vc')
        self.connect('rho',       'induced.rho')
        self.connect('thrust.dT', 'induced.dT')
        self.connect('Omega',     'induced.Omega')
        self.connect('q',         'induced.q')
        self.connect('anhedral',  'induced.anhedral')

        self.add('lift_drag', LiftDrag())
        self.connect('yN',         'lift_drag.yN')
        self.connect('Ns',         'lift_drag.Ns')
        self.connect('dr',         'lift_drag.dr')
        self.connect('r',          'lift_drag.r')
        self.connect('rho',        'lift_drag.rho')
        self.connect('visc',       'lift_drag.visc')
        self.connect('vw',         'lift_drag.vw')
        self.connect('vc',         'lift_drag.vc')
        self.connect('Omega',      'lift_drag.Omega')
        self.connect('c',          'lift_drag.c')
        self.connect('Cl',         'lift_drag.Cl')
        self.connect('d',          'lift_drag.d')
        self.connect('yWire',      'lift_drag.yWire')
        self.connect('zWire',      'lift_drag.zWire')
        self.connect('tWire',      'lift_drag.tWire')
        self.connect('Cm',         'lift_drag.Cm')
        self.connect('xtU',        'lift_drag.xtU')
        self.connect('xtL',        'lift_drag.xtL')
        self.connect('induced.vi', 'lift_drag.vi')
        self.connect('thrust.chordFrac', 'lift_drag.chordFrac')

        self.create_passthrough('induced.vi')
        self.create_passthrough('lift_drag.phi')
        self.create_passthrough('lift_drag.Re')
        self.create_passthrough('lift_drag.Cd')
        self.create_passthrough('lift_drag.Fblade')

        self.driver.workflow.add('thrust')
        self.driver.workflow.add('induced')
        self.driver.workflow.add('lift_drag')
