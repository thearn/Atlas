import unittest

import numpy as np

from Atlas.configuration import AtlasConfiguration


def relative_err(x, y):
    return (np.abs(x-y)/np.linalg.norm(x)).max()


def absolute_err(x, y):
    return np.abs(x-y).max()


class AtlasTestCase(unittest.TestCase):
    def setup(self):
        pass

    def test_configuration(self):
        comp = AtlasConfiguration()

        # run with default values
        comp.run()

        # check calculated outputs
        tol = 1e-10

        self.assertLess(relative_err(comp.yN,
            [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]), tol)

        self.assertLess(relative_err(comp.dr,
            [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]), tol)

        self.assertLess(relative_err(comp.r,
            [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]), tol)


if __name__ == "__main__":
    unittest.main()
