import unittest
import numpy as np
from IMF.maschberger_imf import MaschbergerIMF

class TestMaschbergerIMF(unittest.TestCase):
    def test_cdf(self):
        # Test that the CDF is 0 at lower limit, and 1 at higher limit
        imf = MaschbergerIMF(m_lower=0.001, m_upper=300)
        self.assertEqual(imf.cdf(0.001), 0.0)
        self.assertEqual(imf.cdf(300), 1.0)

    def test_pdf(self):
        # Test that the (approximate) integral over all space under the PDF is (roughly) equal to 1
        imf = MaschbergerIMF(m_lower=0.001, m_upper=300)
        x = np.logspace(np.log10(0.001), np.log10(300), 1000)
        y = imf.pdf(x)
        self.assertAlmostEqual(np.trapz(y, x), 1.0, 3)

        imf = MaschbergerIMF(m_lower=1, m_upper=150)
        x = np.logspace(np.log10(1), np.log10(150), 1000)
        y = imf.pdf(x)
        self.assertAlmostEqual(np.trapz(y, x), 1.0, 3)

    def test_quantile_function(self):
        # Test that the quantile function acts like the inverse of the CDF
        imf = MaschbergerIMF(m_lower=0.001, m_upper=300)
        self.assertAlmostEqual(imf.quantile_function(0), 0.001)
        self.assertAlmostEqual(imf.quantile_function(1), 300)
        self.assertAlmostEqual(imf.quantile_function(imf.cdf(100)), 100)

if __name__ == '__main__':
    unittest.main()