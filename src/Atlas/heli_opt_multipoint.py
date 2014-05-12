
from openmdao.main.api import Assembly, set_as_top
from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Float, Array
from openmdao.lib.drivers.api import SLSQPdriver

try:
    from pyopt_driver import pyopt_driver
except ImportError:
    pyopt_driver = None

from openmdao.util.log import enable_trace, disable_trace

from numpy import pi

from Atlas import AtlasConfiguration, AeroStructural


class ConfigLow(AtlasConfiguration):
    """ Atlas configuration for low altitude """

    Omega_opt = Float(iotype='in', desc='rotor angular velocity')

    H_opt     = Float(iotype='in', desc='height of aircraft')

    def execute(self):
        super(ConfigLow, self).execute()

        # use optimizer provided values
        self.Omega = self.Omega_opt

        self.H     = self.H_opt + self.zWire

        print self.parent.name, '\t', self.H, '\t',  self.Omega, '\t', self.vw, '\t', self.TWire
        print '\t', self.Cl


class AeroStructuralLow(AeroStructural):
    """ AeroStructural assembly for multipoint optimization """

    def configure(self):
        super(AeroStructuralLow, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigLow())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')

        self.create_passthrough('config.H_opt')

        self.create_passthrough('struc.Mtot')
        self.create_passthrough('results.Ttot')
        self.create_passthrough('results.Ptot')


class ConfigHigh(AtlasConfiguration):
    """ Atlas configuration for high altitude """

    Omega_opt = Float(iotype='in', desc='rotor angular velocity')

    Cl0_opt   = Float(iotype='in')
    Cl1_opt   = Float(iotype='in')

    H_opt     = Float(iotype='in', desc='height of aircraft')

    TWire_opt = Float(iotype='in', desc='')

    def execute(self):
        super(ConfigHigh, self).execute()

        # use optimizer provided values
        self.Omega = self.Omega_opt

        self.Cl[0] = self.Cl0_opt
        self.Cl[1] = self.Cl1_opt

        self.H = self.H_opt + self.zWire

        self.TWire = [self.TWire_opt]

        print self.parent.name, '\t', self.H, '\t',  self.Omega, '\t', self.vw, '\t', self.TWire
        print '\t', self.Cl


class AeroStructuralHigh(AeroStructural):
    """ AeroStructural assembly for multipoint optimization """

    def configure(self):
        super(AeroStructuralHigh, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigHigh())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')
        self.create_passthrough('config.Cl0_opt')
        self.create_passthrough('config.Cl1_opt')

        self.create_passthrough('config.H_opt')
        self.create_passthrough('config.TWire_opt')

        self.create_passthrough('struc.Mtot')
        self.create_passthrough('results.Ttot')
        self.create_passthrough('results.Ptot')


class ConfigWind(AtlasConfiguration):
    """ Atlas configuration for wind case """

    Omega_opt = Float(iotype='in', desc='rotor angular velocity')
    OmegaRatio = Float(iotype='in')

    Cl_opt    = Array(iotype='in')

    H_opt     = Float(iotype='in', desc='height of aircraft')

    TWire_opt = Float(iotype='in', desc='')

    vw_opt    = Float(iotype='in', desc='wind velocity')

    def execute(self):
        super(ConfigWind, self).execute()

        # use optimizer provided values
        self.Omega = (self.Omega_opt**3 * self.OmegaRatio)**(1./3.)

        self.Cl = self.Cl_opt

        self.H = self.H_opt + self.zWire

        self.TWire = [self.TWire_opt]

        self.vw = self.vw_opt

        print self.parent.name, '\t', self.H, '\t',  self.Omega, '\t', self.vw, '\t', self.TWire
        print '\t', self.Cl


class AeroStructuralWind(AeroStructural):
    """ AeroStructural assembly for multipoint optimization """

    def configure(self):
        super(AeroStructuralWind, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigWind())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')
        self.create_passthrough('config.OmegaRatio')

        self.create_passthrough('config.Cl_opt')
        self.create_passthrough('config.H_opt')
        self.create_passthrough('config.TWire_opt')
        self.create_passthrough('config.vw_opt')

        self.create_passthrough('struc.Mtot')
        self.create_passthrough('results.Ttot')
        self.create_passthrough('results.Ptot')


class ConfigGravity(AtlasConfiguration):
    """ Atlas configuration for gravity case """

    Omega_opt = Float(iotype='in', desc='rotor angular velocity')
    OmegaRatio = Float(iotype='in')

    Cl_opt    = Array(iotype='in')

    H_opt     = Float(iotype='in', desc='height of aircraft')

    TWire_opt = Float(iotype='in', desc='')

    def execute(self):
        super(ConfigGravity, self).execute()

        # use optimizer provided values
        self.Omega = (self.Omega_opt**3 * self.OmegaRatio)**(1./3.)

        self.Cl = self.Cl_opt

        self.H = self.H_opt + self.zWire

        self.TWire = [self.TWire_opt]

        print self.parent.name, '\t', self.H, '\t',  self.Omega, '\t', self.vw, '\t', self.TWire
        print '\t', self.Cl


class AeroStructuralGravity(AeroStructural):
    """ AeroStructural assembly for multipoint optimization """

    def configure(self):
        super(AeroStructuralGravity, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigGravity())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')
        self.create_passthrough('config.OmegaRatio')

        self.create_passthrough('config.Cl_opt')
        self.create_passthrough('config.H_opt')
        self.create_passthrough('config.TWire_opt')

        self.create_passthrough('struc.Mtot')
        self.create_passthrough('results.Ttot')
        self.create_passthrough('results.Ptot')


class Multipoint(Assembly):
    """ Assembly for multipoint optimization.

        Evaluates AeroStructural for four cases:
            low altitude
            high altitude
            wind
            gravity only
    """

    # configuration inputs
    alt_low   = Float(iotype='in', desc='low altitude')
    alt_high  = Float(iotype='in', desc='high altitude')
    alt_ratio = Float(iotype='in', desc='proportion of time near ground')

    TWire_high = Float(iotype='in')
    TWire_wind = Float(iotype='in')
    TWire_grav = Float(iotype='in')

    OmegaRatio = Float(iotype='in')

    vw         = Float(iotype='in', desc='wind velocity')

    Cl_max     = Array(iotype='in')

    # optimizer parameters
    Omega_low  = Float(iotype='in', desc='rotor angular velocity, low altitude')
    Omega_high = Float(iotype='in', desc='rotor angular velocity, high altitude')
    Cl0_high   = Float(iotype='in')
    Cl1_high   = Float(iotype='in')

    # outputs
    P       = Float(iotype='out', desc='')

    def configure(self):
        # low altitude
        self.add('low', AeroStructuralLow())

        self.connect('Omega_low', 'low.Omega_opt')
        self.connect('alt_low',   'low.H_opt')

        self.create_passthrough('low.Mtot', 'Mtot_low')
        self.create_passthrough('low.Ttot', 'Ttot_low')

        # high altitude
        # need a different rotor speed and lift distribution at altitude
        self.add('high', AeroStructuralHigh())

        self.connect('Omega_high', 'high.Omega_opt')
        self.connect('Cl0_high',   'high.Cl0_opt')
        self.connect('Cl1_high',   'high.Cl1_opt')
        self.connect('alt_high',   'high.H_opt')
        self.connect('TWire_high', 'high.TWire_opt')

        self.create_passthrough('high.Mtot', 'Mtot_high')
        self.create_passthrough('high.Ttot', 'Ttot_high')

        # wind case
        self.add('wind', AeroStructuralWind())

        self.connect('Omega_high', 'wind.Omega_opt')
        self.connect('OmegaRatio', 'wind.OmegaRatio')
        self.connect('Cl_max',     'wind.Cl_opt')
        self.connect('alt_high',   'wind.H_opt')
        self.connect('TWire_wind', 'wind.TWire_opt')
        self.connect('vw',         'wind.vw_opt')
        # TODO: verify that the following flags are respected
        self.high.config.flags.FreeWake = 0  # momentum theory
        self.high.config.flags.AeroStr  = 0  # assume flat wing (no deformation)

        # gravity case
        self.add('grav', AeroStructuralGravity())

        self.connect('Omega_high', 'grav.Omega_opt')
        self.connect('OmegaRatio', 'grav.OmegaRatio')
        self.connect('Cl_max',     'grav.Cl_opt')
        self.connect('alt_high',   'grav.H_opt')
        self.connect('TWire_grav', 'grav.TWire_opt')
        # TODO: verify that the following flags are respected
        self.grav.config.flags.Load     = 1  # gravity and wire forces only
        self.grav.config.flags.FreeWake = 0  # momentum theory
        self.grav.config.flags.AeroStr  = 0  # assume flat wing (no deformation)

        self.connect('alt_ratio*low.Ptot + (1 - alt_ratio)*high.Ptot', 'P')

        self.driver.workflow.add(['low', 'high', 'wind', 'grav'])


class HeliOptM(Assembly):
    """ Multipoint aero-structural optimization """

    def configure(self):
        # add an optimizer and a multi-point AeroStructural assembly
        if pyopt_driver and 'SNOPT' in pyopt_driver._check_imports():
            self.add("driver", pyopt_driver.pyOptDriver())
            self.driver.optimizer = "SNOPT"
            self.driver.options = {
                # any changes to default SNOPT options?
            }
        else:
            print 'SNOPT not available, using SLSQP'
            self.add('driver', SLSQPdriver())

        self.add('mp', Multipoint())

        self.mp.alt_low = 0.5         # low altitude
        self.mp.alt_high = 3.5        # high altitude
        self.mp.alt_ratio = 35./60.   # proportion of time near ground

        self.mp.TWire_high = 900
        self.mp.TWire_wind = 2100
        self.mp.TWire_grav = 110

        self.mp.OmegaRatio  = 2

        self.mp.vw = 0/3.6   # zero

        self.mp.Cl_max = [1.4, 1.35, 1.55]    # max control

        # objective: minimize total power
        self.driver.add_objective('mp.P')

        # parameter: rotor speed
        self.driver.add_parameter('mp.Omega_low',
                                  low=0.15*2*pi, high=0.25*2*pi)
        self.mp.Omega_low = 0.20*2*pi  # initial value

        self.driver.add_parameter('mp.Omega_high',
                                  low=0.15*2*pi, high=0.19*2*pi)
        self.mp.Omega_high = 0.17*2*pi  # initial value

        # parameter: lift distribution at high altitude
        self.driver.add_parameter('mp.Cl0_high',
                                  low=0.8, high=1.4)
        self.mp.Cl0_high = 1.

        self.driver.add_parameter('mp.Cl1_high',
                                  low=0.8, high=1.3)
        self.mp.Cl1_high = 1.

        # constraint: lift >= weight
        self.driver.add_constraint('mp.Mtot_low*9.8-mp.Ttot_low<=0')
        self.driver.add_constraint('mp.Mtot_high*9.8-mp.Ttot_high<=0')

        # TODO: optional constraints
        #    if flags.ConFail:
        #       Structural Failure in Rotor Spar (ConFail)
        #       Buckling failure of spar (ConFailBuck)
        #       Tensile failure in wire (ConFailWire)
        #
        #    if flags.ConDef:
        #       Constraints on Maximum Deformation (ConDelta)
        #
        #    if flags.MultiPoint && flags.ConJigCont:
        #       Consistent jig twist (ConAlphaJig)
        #
        #    if flags.MultiPoint && flags.ConWireCont
        #       Wire stretch consistency (conWire)

        # Optimization Constraints  (not used... yet)
        vrCon = VariableTree()
        vrCon.MaxDelta    = -0.1
        vrCon.MinDelta    = 0.1
        vrCon.FOSmat      = 0.55    # 1.3
        vrCon.FOSbuck     = 0.5     # 1.3
        vrCon.FOSquadbuck = 5.
        vrCon.FOStorbuck  = 0.5     # 1.5
        vrCon.FOSwire     = 0.5     # 2


if __name__ == '__main__':
    # TODO: create units tests for the following

    if True:
        print '====== Multipoint ======'
        mp = set_as_top(Multipoint())

        mp.alt_low = 0.5         # low altitude
        mp.alt_high = 3.5        # high altitude
        mp.alt_ratio = 35./60.   # proportion of time near ground

        mp.TWire_high = 900
        mp.TWire_wind = 2100
        mp.TWire_grav = 110

        mp.OmegaRatio = 2

        mp.vw = 0/3.6

        mp.Cl_max = [1.4, 1.35, 1.55]  # max control

        mp.Omega_low  = 1.0512
        mp.Omega_high = 1.0771
        mp.Cl0_high   = 1.4000
        mp.Cl1_high   = 1.3000

        mp.run()

        print 'low Ptot  =', mp.low.Ptot,  ' (reference:  421.3185)'
        print 'low Mtot  =', mp.Mtot_low,  ' (reference:  126.1670)'
        print 'low Ttot  =', mp.Ttot_low,  ' (reference: 1236.4369)'
        print
        print 'high Ptot =', mp.high.Ptot, ' (reference:  846.6429)'
        print 'high Mtot =', mp.Mtot_high, ' (reference:  126.1670)'
        print 'high Ttot =', mp.Ttot_high, ' (reference: 1236.4368)'
        print
        print 'wind Ptot =', mp.wind.Ptot, ' (reference: 1815.7617)'
        print 'wind Mtot =', mp.wind.Mtot, ' (reference:  126.1670)'
        print 'wind Ttot =', mp.wind.Ttot, ' (reference: 2239.7556)'
        print
        print 'grav Ptot =', mp.grav.Ptot, ' (reference:    0.00000'
        print 'grav Mtot =', mp.grav.Mtot, ' (reference:  126.1670)'
        print 'grav Ttot =', mp.grav.Ttot, ' (reference:    0.0000)'
        print
        print 'P =', mp.P,                 ' (reference: 598.537)'

    if False:
        print '====== HeliOptM ======'
        opt = set_as_top(HeliOptM())

        # enable_trace()
        opt.run()

        print 'Objective:  P =', opt.mp.P

        print 'Constraint: Low Weight-Lift =', opt.mp.Mtot_low*9.8-opt.mp.Ttot_low
        print 'Constraint: High Weight-Lift =', opt.mp.Mtot_high*9.8-opt.mp.Ttot_high

        print 'Parameter:  Omega (Low) =', opt.mp.Omega_low
        print 'Parameter:  Omega (High) =', opt.mp.Omega_high
