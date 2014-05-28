from openmdao.main.api import Assembly, set_as_top
from openmdao.lib.datatypes.api import Float
from openmdao.lib.drivers.api import SLSQPdriver

try:
    from pyopt_driver import pyopt_driver
except ImportError:
    pyopt_driver = None

from openmdao.util.log import enable_trace  # , disable_trace

from numpy import pi

from Atlas import AtlasConfiguration, AeroStructural


class ConfigOpt(AtlasConfiguration):
    """ Atlas configuration for single point optimization """

    # inputs for optimizer
    Omega_opt = Float(iotype='in', desc='rotor angular velocity')

    def execute(self):
        super(ConfigOpt, self).execute()

        # use optimizer provided value for Omega
        self.Omega = self.Omega_opt


class AeroStructuralOpt(AeroStructural):
    """ AeroStructural assembly for single point optimization """

    def configure(self):
        super(AeroStructuralOpt, self).configure()

        # replace config with optimizer driven config
        self.replace('config', ConfigOpt())

        # create passthroughs for variables used by the optimizer
        self.create_passthrough('config.Omega_opt')
        self.create_passthrough('struc.Mtot')
        self.create_passthrough('results.Ttot')
        self.create_passthrough('results.Ptot')


class HeliOpt(Assembly):
    """ Single point aero-structural optimization """

    def configure(self):
        # add an optimizer and an AeroStructural assembly
        if pyopt_driver and 'SNOPT' in pyopt_driver._check_imports():
            self.add("driver", pyopt_driver.pyOptDriver())
            self.driver.optimizer = "SNOPT"
            self.driver.options = {
                # any changes to default SNOPT options?
            }
        else:
            print 'SNOPT not available, using SLSQP'
            self.add('driver', SLSQPdriver())

        self.add('aso', AeroStructuralOpt())

        # objective: minimize total power
        self.driver.add_objective('aso.Ptot')

        # parameter: rotor speed
        self.driver.add_parameter('aso.Omega_opt',
                                  low=0.15*2*pi, high=0.25*2*pi)
        self.aso.Omega_opt = 0.2*2*pi  # initial value

        # constraint: lift >= weight
        self.driver.add_constraint('aso.Mtot*9.8-aso.Ttot<=0')

        # TODO: optional constraints
        #
        #    if flags.ConFail:
        #       Structural Failure in Rotor Spar (ConFail)
        #       Buckling failure of spar (ConFailBuck)
        #       Tensile failure in wire (ConFailWire)
        #
        #    if flags.ConDef:
        #       Constraints on Maximum Deformation (ConDelta)


if __name__ == '__main__':
    opt = set_as_top(HeliOpt())
    opt.aso.Omega_opt = 1.0512

    opt.driver.run_iteration()

    print 'Parameter:  Omega =', opt.aso.config.Omega

    print 'Constraint: Weight-Lift =', (opt.aso.Mtot*9.8-opt.aso.Ttot)

    print 'Objective:  Ptot =', opt.aso.Ptot

    # enable_trace()
    exit()
    opt.run()


    print 'Parameter:  Omega =', opt.aso.config.Omega

    print 'Constraint: Weight-Lift =', (opt.aso.Mtot*9.8-opt.aso.Ttot)

    print 'Objective:  Ptot =', opt.aso.Ptot

    # for reference, MATLAB solution:
    #    Omega: 1.0512
    #    Ptot: 421.3185
