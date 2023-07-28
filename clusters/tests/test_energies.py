import unittest
import warnings
import numpy as np
import clusters.energies as energies

class TestEnergies(unittest.TestCase):
    def test_pe(self):
        pos = np.array([
            [0,0,0],
            [0,0,10],
            [0,0,5]
        ])
        masses = np.array([
            5,
            8,
            15
        ])
        with warnings.catch_warnings():
            # Ignore numpy's division warnings
            warnings.filterwarnings('ignore', r'invalid value encountered')
            pe = energies.calculate_pe(pos, masses)
            self.assertAlmostEqual(pe[0], (-1.055642595e51 + -3.958659731e51)/(3.085678e16), -35+6)

if __name__ == '__main__':
    unittest.main()