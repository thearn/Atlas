
from openmdao.main.api import Assembly, set_as_top
from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Float, Array
from openmdao.lib.drivers.api import SLSQPdriver

from openmdao.util.log import enable_trace  # , disable_trace

from numpy import pi

from Atlas import AtlasConfiguration, AeroStructural


class ConfigOptM(AtlasConfiguration):
    """ Atlas configuration for multipoint optimization """

    # inputs for optimizer
    H_opt     = Float(iotype='in', desc='height of aircraft')
    Omega_opt = Float(iotype='in', desc='rotor angular velocity')
    Cl_opt    = Float(iotype='in', desc='lift coefficient distribution')

    TWire_opt = Float(iotype='in', desc='')

    vw_opt    = Float(0.0,  iotype='in', desc='wind velocity')

    def execute(self):
        super(ConfigOptM, self).execute()

        # use optimizer provided values
        self.Omega = self.Omega_opt
        self.H     = self.H_opt + self.zWire
        self.Cl    = self.Cl_opt

        if self.TWire:
            self.TWire = [self.TWire_opt]

        if self.vw_opt:
            self.vw = [self.vw_opt]


class AeroStructuralOptM(AeroStructural):
    """ AeroStructural assembly for multipoint optimization """

    def configure(self):
        super(AeroStructuralOptM, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigOptM())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')
        self.create_passthrough('config.H_opt')
        self.create_passthrough('config.Cl_opt')
        self.create_passthrough('config.TWire_opt')
        self.create_passthrough('config.vw_opt')

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

    # inputs
    alt_low   = Float(iotype='in', desc='low altitude')
    alt_high  = Float(iotype='in', desc='high altitude')
    alt_ratio = Float(iotype='in', desc='proportion of time near ground')

    TWire_high = Float(iotype='in')
    TWire_wind = Float(iotype='in')
    TWire_grav = Float(iotype='in')

    OmegaRatio = Float(iotype='in')

    vw         = Float(iotype='in', desc='wind velocity')

    Cl_max     = Float(iotype='in')

    # optimizer parameters
    Omega_opt  = Float(iotype='in', desc='rotor angular velocity')
    Cl_opt     = Array(iotype='in', desc='lift coefficient distribution')

    # outputs
    P       = Float(iotype='out', desc='')

    def configure(self):
        # low altitude
        self.add('low', AeroStructuralOptM())
        self.low.replace('config', ConfigOptM())

        self.connect('alt_low',   'low.H_opt')
        self.connect('Omega_opt', 'low.Omega_opt')
        self.connect('Cl_opt',    'low.Cl_opt')

        self.create_passthrough('low.Mtot', 'Mtot_low')
        self.create_passthrough('low.Ttot', 'Ttot_low')

        # high altitude
        self.add('high', AeroStructuralOptM())
        self.high.replace('config', ConfigOptM())

        self.connect('alt_high',   'high.H_opt')
        self.connect('Omega_opt',  'high.Omega_opt')
        self.connect('Cl_opt',     'high.Cl_opt')
        self.connect('TWire_high', 'high.TWire_opt')

        self.create_passthrough('high.Mtot', 'Mtot_high')
        self.create_passthrough('high.Ttot', 'Ttot_high')

        # wind case
        self.add('wind', AeroStructuralOptM())
        self.high.replace('config', ConfigOptM())

        self.connect('alt_high',   'wind.H_opt')
        omega_expr = 'math.pow(Omega_opt**3 * OmegaRatio, float(1)/3)'
        self.connect(omega_expr,   'wind.Omega_opt')  # FIXME
        self.connect('Cl_max',     'wind.Cl_opt')
        self.connect('TWire_wind', 'wind.TWire_opt')
        self.connect('vw',         'wind.vw_opt')
        # TODO: verify that the following flags are respected
        self.high.config.flags.FreeWake = 0  # momentum theory
        self.high.config.flags.AeroStr  = 0  # assume flat wing (no deformation)

        # gravity case
        self.add('grav', AeroStructuralOptM())
        self.high.replace('config', ConfigOptM())

        self.connect('alt_high',   'grav.H_opt')
        self.connect('Omega_opt',  'grav.Omega_opt')
        self.connect('Cl_opt',     'grav.Cl_opt')
        self.connect('TWire_grav', 'grav.TWire_opt')
        # TODO: verify that the following flags are respected
        self.high.config.flags.Load     = 1  # gravity and wire forces only
        self.high.config.flags.FreeWake = 0  # momentum theory
        self.high.config.flags.AeroStr  = 0  # assume flat wing (no deformation)

        self.connect('alt_ratio*low.Ptot + (1 - alt_ratio)*high.Ptot', 'P')


class HeliOptM(Assembly):
    """ Multipoint aero-structural optimization """

    def configure(self):
        # add an optimizer and an AeroStructural assembly
        self.add('driver', SLSQPdriver())
        self.add('mp', Multipoint())

        self.mp.alt_low = 0.5         # low altitude
        self.mp.alt_high = 3.5        # high altitude
        self.mp.alt_ratio = 35./60.   # proportion of time near ground

        self.mp.TWire_high = 900
        self.mp.TWire_wind = 2100
        self.mp.TWire_grav = 110

        self.mp.OmegaRatio  = 2

        self.mp.vw = 0/3.6

        # objective: minimize total power
        self.driver.add_objective('mp.P')

        # parameter: rotor speed
        self.driver.add_parameter('mp.Omega_opt',
                                  low=0.15*2*pi, high=0.19*2*pi)
        self.mp.Omega_opt = 0.17*2*pi  # initial value

        # parameter: lift coefficient distribution
        self.driver.add_parameter('mp.Cl_opt',
                                  low=[0.8, 0.8], high=[1.4, 1.3])
        self.mp.Cl_opt = [1.0, 1.0]  # initial value

        # constraint: lift >= weight
        self.driver.add_constraint('mp.M_low*9.8-mp.Ttot_low<=0')
        self.driver.add_constraint('mp.M_high*9.8-mp.Ttot_high<=0')

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
        vrCon.ClMax       = [1.4, 1.25, 0.85]  # min control
        vrCon.ClMax       = [1.4, 1.35, 1.55]  # max control
        vrCon.FOSmat      = 0.55    # 1.3
        vrCon.FOSbuck     = 0.5     # 1.3
        vrCon.FOSquadbuck = 5.
        vrCon.FOStorbuck  = 0.5     # 1.5
        vrCon.FOSwire     = 0.5     # 2


if __name__ == '__main__':
    opt = set_as_top(HeliOptM())

    # enable_trace()
    opt.run()

    print 'Objective:  Ptot =', opt.mp.Ptot

    print 'Constraint: Weight-Lift =', (opt.mp.Mtot*9.8-opt.mp.Ttot)

    print 'Parameter:  Omega =', opt.mp.config.Omega
