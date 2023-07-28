import numpy as np
from scipy.stats import chisquare, ks_1samp, ks_2samp, cramervonmises

class MaschbergerIMF():
    """
    Recreating the IMF from Maschberger (2011)
    """
    # Parameters
    m_upper = None  # Upper mass limit (in solar masses)
    m_lower = None  # Lower mass limit (in solar masses)
    alpha   = None  # High-mass exponent
    beta    = None  # Low-mass exponent
    mu      = None  # Scale parameter

    def __init__(self, m_lower=0.01, m_upper=150, alpha=2.3, beta=1.4, mu=0.2):
        self.m_upper = m_upper
        self.m_lower = m_lower
        self.alpha = alpha
        self.beta = beta
        self.mu = mu

    # Auxiliary function
    def _G(self, mass: float) -> float:
        return np.power((1 + np.power(mass/self.mu, 1-self.alpha)), 1-self.beta)
    
    # Constant for PDF
    def _A(self) -> float:
        a = ( (1-self.alpha)*(1-self.beta) )/self.mu
        b = 1 / ( self._G(self.m_upper)-self._G(self.m_lower) )
        return a*b

    def cdf(self, mass: float) -> float:
        """
        Cumulative distribution function (CDF)
        
        Parameters
        ----------
        mass : float, list or np.ndarray

        Returns
        -------
        float or np.ndarray
        """
        return (self._G(mass) - self._G(self.m_lower))/(self._G(self.m_upper)-self._G(self.m_lower))

    def quantile_function(self, u: float) -> float:
        """
        Quantile function - inverse of the CDF
        
        Parameters
        ----------
        u : float, list or np.ndarray

        Returns
        -------
        float or np.ndarray
        """
        inner = u * (self._G(self.m_upper) - self._G(self.m_lower)) + self._G(self.m_lower)
        middle = np.power(inner, 1/(1-self.beta)) - 1
        outer = self.mu * np.power(middle, 1/(1-self.alpha))
        return outer
    
    def pdf(self, mass: float) -> float:
        """
        Probability density function (PDF)
        
        Parameters
        ----------
        mass : float, list or np.ndarray

        Returns
        -------
        float or np.ndarray
        """
        return self._A() * np.power(mass/self.mu, -self.alpha) * np.power(
            1 + np.power(mass/self.mu, 1-self.alpha), -self.beta
        )
    
    def histogram(self, N: int, bins: int, start: float = None, end: float = None) -> tuple:
        """
        Generate a histogram of the no. of stars across the IMF
        
        Parameters
        ----------
        N : int
            Total number of stars.
        bins : int
            Number of bins.
        start : float, default: `self.m_upper`
            Start value
        end : float, default: `self.m_lower`
            End value
        
        Returns
        -------
        no_of_stars : list
            Number of stars in each bin.
        bin_edges : list
            Edges of bins.
        """
        hist_start = self.m_upper if (start == None) else start
        hist_end = self.m_upper if (start == None) else end

        # Calculate the width of each bin
        width = (hist_end-hist_start)/bins

        no_of_stars = []    # No. of stars in each bin
        bin_edges = []      # Centre points of bins

        bin_edges.append(hist_start)

        for i in range(0,bins):
            # Calculate bin start & end
            bin_start = hist_start + i*width
            bin_end = hist_start + (i+1)*width

            # Calculate probability
            x = np.linspace(bin_start, bin_end, 1000)
            P = np.trapz(self.pdf(x),x)
            
            # Append values to lists
            no_of_stars.append(N*P)
            bin_edges.append(bin_end)
        
        return (no_of_stars, bin_edges)
    
    def ks_test(self, mass_data: list):
        """
        Perform a 1 way ks-test of the data against the IMF

        Parameters
        ----------
        mass_data : list
            Data to test

        Returns
        -------
        KstestResult
            See docs for scipy.stats.ks_1samp.
        """
        return ks_1samp(mass_data, self.cdf)

    def chisquare(self, mass_data: list, n: int = 20):
        """
        Perform a 1 way ks-test of the data against the IMF

        Parameters
        ----------
        mass_data : list
            Data to test
        n : int, default: 20
            Number of bins to split the data into.

        Returns
        -------
        chisq : float
            The chi-squared test statistic.
        p : float
            The p-value of the test.
        """
        cluster_hist, cluster_edges = np.histogram(mass_data, n)
        standard_hist, standard_edges = self.histogram(len(mass_data), n, start=min(mass_data), end=max(mass_data))

        if (min(standard_hist) < 5):
            pass
            # print("Some counts are less than 5. Result may be invalid.")
        
        # This is a hack to stop scipy from screaming about the sums of the data not matching up!
        difference = sum(cluster_hist)-sum(standard_hist)
        if (np.abs(difference) > 1e-8):
            standard_hist[0] = standard_hist[0] + difference

        return chisquare(cluster_hist, standard_hist)

    def cvm_test(self, mass_data: list):
        """
        Perform a one-sample Cramér-von Mises test of the data against the IMF

        Parameters
        ----------
        mass_data : list
            Data to test

        Returns
        -------
        statistic : float
            The cvm test statistic.
        p : float
            The p-value of the test.
        """
        return cramervonmises(mass_data, self.cdf)