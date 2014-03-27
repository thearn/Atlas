from Atlas import dragCoefficient, frictionCoefficient
import unittest
from openmdao.util.testutil import assert_rel_error


class AtlasTestCase(unittest.TestCase):

    def setUp(self):
        self.tol = 0.0001

    def test_dragCoefficient(self):
        comp = dragCoefficient()
        comp.Re = 4.7632e5
        comp.tc = 0.15
        comp.xtcU = 0.15
        comp.xtcL = 0.15

        comp.run()

        assert_rel_error(self, comp.Cd, 0.013329, self.tol)

        comp.Re = 4.7632e5
        comp.tc = 0.15
        comp.xtcU = 0.5
        comp.xtcL = 1

        comp.run()
        assert_rel_error(self, comp.Cd, 0.0077468, self.tol)

    # def test_dragCoefficientFit(self):
    #     comp = dragCoefficientFit()
    #     comp.Re = 2.9078e5
    #     comp.tc = 0.14
    #     comp.xtcU = 0.15
    #     comp.xtcL = 0.3

    #     comp.run()
    #     assert_rel_error(self, comp.Cd, 0.018375, self.tol)

    #     comp.Re = 2.8504e5
    #     comp.tc = 0.14
    #     comp.xtcU = 0.15
    #     comp.xtcL = 0.3

    #     comp.run()
    #     assert_rel_error(self, comp.Cd, 0.018487, self.tol)

    def test_frictionCoefficient(self):
        Re = 500000
        xtc = 0.5
        Cfflat = frictionCoefficient(Re, xtc)

        assert_rel_error(self, Cfflat, 0.003846, self.tol)

        Re = 4.7632e5
        xtc = 1.
        Cfflat = frictionCoefficient(Re, xtc)

        assert_rel_error(self, Cfflat, 0.0019241, self.tol)


if __name__ == "__main__":
    unittest.main()
