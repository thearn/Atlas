from Atlas import dragCoefficient, dragCoefficientFit, frictionCoefficient
import unittest


class AtlasTestCase(unittest.TestCase):

    
    def test_dragCoefficient(self):
        comp = dragCoefficient()
        comp.Re = 4.7632e5
        comp.tc = 0.15
        comp.xtcU = 0.15
        comp.xtcL = 0.15

        comp.run()
        assert comp.Cd == 0.0133

        comp.Re = 4.7632e5
        comp.tc = 0.15
        comp.xtcU = 0.5
        comp.xtcL = 1

        comp.run()
        assert comp.Cd == 0.0077

    def test_dragCoefficientFit(self):
        comp = dragCoefficientFit()
        comp.Re = 2.9078
        comp.tc = 0.14
        comp.xtcU = 0.15
        comp.xtcL = 0.3

        comp.run()
        assert comp.Cd == 0.0184

        comp.Re = 2.8504e5
        comp.tc = 0.14
        comp.xtcU = 0.15
        comp.xtcL = 0.3

        comp.run()
        assert comp.Cd == 0.0185

    def test_frictionCoefficient(self):
        comp = frictionCoefficient()
        comp.Re = 500000
        comp.xtc = 0.5

        comp.run()
        assert comp.Cfflat == 0.0038

        comp.Re = 4.7632e5
        comp.xtc = 1.

        comp.run()
        assert comp.Cfflat == 0.0019

        
if __name__ == "__main__":
    unittest.main()
    