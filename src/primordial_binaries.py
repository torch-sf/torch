"""
Binary generation algorithm, making use of statistics by Moe & Di Stefano (2017), 
Winters et al. (2019) and Offner et al. (2022)
Claude Cournoyer-Cloutier, McMaster University, 2020-2025
If used, please cite Cournoyer-Cloutier et al. (2024), ApJ 977, Issue 2, id. 203
https://ui.adsabs.harvard.edu/abs/2024ApJ...977..203C/abstract
"""

import numpy as np
from amuse.lab import units
from amuse.ext.orbital_elements import *

def get_multiplicity(m, mult_frac='field'):
    print("importing correct function...")
    """
    Function to select primaries and singles
    args:
        m: array of masses in MSun
        mult_frac: mass-dependant multiplicity fraction,
                   default is 'field'
    returns:
        primaries_IDs: indices for primaries
        singles_IDs: indices for singel stars
    """
    
    def companion_frequency(m):
        """
        Use a piecewise function to return the companion frequency 
        as a function of primary mass. Use the average between 
        mass ranges.
        args:
            m: array of masses in MSun
        returns:
            step function for sampling
        """
        mass_ranges = [m < 0.15, 
                       (m >= 0.15) & (m < 0.30), 
                       (m >= 0.30) & (m < 0.60), 
                       (m >= 0.60) & (m < 0.80),
                       (m >= 0.80) & (m < 1.2), 
                       (m >= 1.2) & (m < 2), 
                       (m >= 2) & (m < 5), 
                       (m >= 5) & (m < 9), 
                       (m >= 9) & (m < 16), 
                       m >= 16]
        frequencies = [0.19, 0.23, 0.30, 0.35, 0.40, 0.50, 0.59, 0.76, 0.84, 0.94]
        
        return np.piecewise(m, mass_ranges, frequencies)
    
    frq = companion_frequency(m)
    rnd = np.random.uniform(size=len(m))
    
    primaries_IDs = np.zeros(len(frq), dtype=bool)
    singles_IDs   = np.zeros(len(frq), dtype=bool)
    primaries_IDs[np.where(rnd < frq)] = np.ones(len(np.where(rnd < frq)), dtype=bool)
    singles_IDs[np.where(rnd >= frq)]  = np.ones(len(np.where(rnd >= frq)), dtype=bool)

    return primaries_IDs, singles_IDs



def get_periods(m, pdist='inner'):
    """
    Sample the period based on the mass-dependant period distributions 
    reported by Moe & di Stefano (2017). For M-dwarfs, the period distributions
    are derived from the semi-major distributions compiled by Winters et al. 2019, 
    using the mass ratio distributions reported in Offner et al. 2023.
    args:
        m: array of masses in MSun
        pdist: period distribution, possible values are 'field' and 'inner'
    returns:
        p: array of orbital periods in days
    """
        
    def pdf(log_p, m_range=-1, pdist=pdist):
        """
        Generate the period probability distribution function for a given
        primary mass range.
        args:
            log_p: values of log(period/days) over which the pdf is evaluated
            m_range: mass range between -5 (solar-type) and -1 (O stars)
            p_dist: period distribution, possible values are 'field' and 'inner'
        returns:
            probabilities: value of the pdf for each log_p value
        """
        
        log_p_values = [1., 3., 5., 7.]
        
        # Companion frequencies at log_p = 1, 3, 5, 7 from Moe & di Stefano 2017,
        # for each primary mass bin
        frequencies = np.array([[0.027, 0.057, 0.095, 0.075], #solar
                                [0.07, 0.12, 0.13, 0.09],     # AB
                                [0.14, 0.22, 0.20, 0.11],     # mid B
                                [0.19, 0.26, 0.23, 0.13],     # early B
                                [0.29, 0.32, 0.30, 0.18]])    # O

        # Generate the pdf
        probabilities = np.interp(log_p, log_p_values, frequencies[m_range]) # Flat beyond bounds
        
        if pdist == 'inner':
            # Correct the pdf for the inner companion
            # 1./2 leaves the probabilities unchanged
            frac_close = np.array([1./2, 1./2, 63./76, 80./84, 1.]) # From close binary frequency
            _inner = np.where(log_p <= 3.7)[0]
            _outer = np.where(log_p > 3.7)[0]
            probabilities[_inner] *= frac_close[m_range]
            probabilities[_outer] *= (1-frac_close)[m_range]
            
        return probabilities
    
    def pdf_Mdwarf(log_p, m_range=-1):
        """
        Generate the period probability distribution function for a given
        primary mass range for M-dwarfs.
        args:
            log_p: values of log(period/days) over which the pdf is evaluated
            m_range: mass range between 0 (< 0.15 MSun) and 2 (0.3-0.6 MSun)
        returns:
            probabilities: value of the pdf for each log_p value
        """
        
        log_p_values = [1., 2., 3., 4., 5., 6., 7.]
        
        # Companion frequencies for log_p = 1, 2, 3, 4, 5, 6, 7, derived from the semi-major
        # axis distributions of Winters et al. 2019 using 1e6 stars and a random seed of 0
        frequencies = np.array([[0.020,  0.117,  0.308,  0.346,  0.170,  0.035,  0.003],  # < 0.15 MSun
                                [0.053,  0.125,  0.211,  0.245,  0.203,  0.116,  0.047],  # 0.15-0.30 MSun
                                [0.049,  0.106,  0.171,  0.217,  0.212,  0.154,  0.091]]) # 0.30-0.60 MSun
        
        probabilities = np.interp(log_p, log_p_values, frequencies[m_range]) # Flat beyond bounds
            
        return probabilities
        
    def inv_transform_sampling(m, m_range=-1, pdist=pdist):
        """
        Use inverse transform sampling to sample the period pdf for a given mass range.
        args:
            m: array of masses in MSun
            m_range: mass range between -5 and 2
            p_dist: period distribution, possible values are 'field' and 'inner'
        returns:
            return_log_p: array of log(period/days) of the same length as m
        """
            
        log_p = np.arange(0.2, 7.501, 0.001)
        
        if m_range < 0:
            _pdf = pdf(log_p, m_range, pdist=pdist)
        else:
            _pdf = pdf_Mdwarf(log_p, m_range)
            
        cdf = np.cumsum(_pdf)
        cdf = cdf/cdf[-1] # Normalize the cdf
            
        percentiles = np.random.uniform(size=len(m))
        args = abs(percentiles[:, None] - cdf[None, :]).argmin(axis=-1)
        return_log_p = log_p[args]
            
        return return_log_p

    p = np.zeros(len(m))
    # Use negative ranges for M17 and positive for M dwarfs
    # O stars
    _mask = np.where(m >= 16)[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=-1, pdist=pdist)
    # Early B stars
    _mask = np.where((m >= 9) & (m < 16))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=-2, pdist=pdist)
    # Mid B stars
    _mask = np.where((m >= 5) & (m < 9))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=-3, pdist=pdist)
    # AB stars
    _mask = np.where((m >= 1.6) & (m < 5))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=-4, pdist=pdist)
    # Solar-type stars and below
    _mask = np.where((m >= 0.6) & (m < 1.6))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=-5, pdist=pdist)
    # M-dwarfs
    _mask = np.where(m < 0.15)[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=0, pdist=pdist)
    _mask = np.where((m >= 0.15) & (m < 0.3))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=1, pdist=pdist)
    _mask = np.where((m >= 0.3) & (m < 0.6))[0]
    p[_mask] = 10**inv_transform_sampling(m[_mask], m_range=2, pdist=pdist)

    # Verify that all period values are within bounds
    assert np.min(p) >= 10**0.2
    assert np.max(p) <= 10**7.5
    
    return p


def get_mass_ratios(m, p, qdist='field', mmin=0.08):
    """
    Sample the mass ratio based on the distributions reported by Moe & di Stefano (2017),
    for the specified mass range, for a given period. For M dwarfs, slopes are from
    Offner et al. 2023, with the 0.15-0.30 MSun also used below 0.15 MSun.
    args:
        m: array of masses in MSun
        p: array of orbital periods in days
        q_dist: distribution of mass ratios, only possible choice is 'field'
        mmin: minimum companion mass
    returns: 
        q: array of mass ratios of the same length as m and p
    """

    def pdf(q, p, m_range=-1):
        """
        Generate the mass ratio probability distribution function for a given
        period and mass range.
        args:
            q: values of the mass ratio q over which the pdf is evaluated
            p: array of periods in days
            m_range: mass range between -5 and 2
        returns:
            pdfs: probability distribution functions for each of the four period bins
        """
        
        q_ranges = [q < 0.1, (q >= 0.1) & (q < 0.3), (q >= 0.3) & (q < 0.95), (q >= 0.95) & (q <= 1), q > 1]

        # Power-law exponent gamma for q < 0.3
        gamma_S = np.array([[0.7, 0.7, 0.7, 0.7],
                            [0.1, 0.1, 0.1, 0.1],
                            [0.3, 0.3, 0.3, 0.3],
                            [0.2, 0.1, -0.5, -1.0],
                            [0.1, -0.2, -1.2, -1.5],
                            [0.1, -0.2, -1.2, -1.5],
                            [0.1, -0.2, -1.2, -1.5]])

        # Power-law exponent gamma for q > 0.3
        gamma_L = np.array([[0.7, 0.7, 0.7, 0.7],
                            [0.1, 0.1, 0.1, 0.1],
                            [-0.5, -0.5, -0.5, -1.1],
                            [-0.5, -0.9, -1.4, -2.0],
                            [-0.5, -1.7, -2.0, -2.0],
                            [-0.5, -1.7, -2.0, -2.0],
                            [-0.5, -1.7, -2.0, -2.0]])

        # Excess twin fraction
        excess_twin = np.array([[0., 0., 0., 0.,],
                                [0., 0., 0., 0.,],
                                [0.30, 0.2, 0.1, 0.],
                                [0.22, 0.10, 0., 0.],
                                [0.17, 0., 0., 0.],
                                [0.14, 0., 0., 0.],
                                [0.08, 0., 0., 0.]])
            
        pdfs = np.zeros((4, len(q))) # Create an array for the different pdfs, with four period bins
        
        for i in range(4): # 4 period bins here
            frequencies = [lambda q: 0, 
                           lambda q: q**gamma_S[m_range][i] * 0.3**(gamma_L[m_range][i]-gamma_S[m_range][i]), 
                           lambda q: q**gamma_L[m_range][i], 
                           lambda q: q**gamma_L[m_range][i] + excess_twin[m_range][i], 
                           lambda q: 0]
            
            probabilities = np.piecewise(q, q_ranges, frequencies) #m, mass_ranges, frequencies
            pdfs[i] = probabilities
        
        return pdfs
    
    def inv_transform_sampling(p, m):
        """
        Use inverse transform sampling to sample the mass ratio pdf for a given mass and period range.
        args:
            p: array of orbital periods in days
            m: array of masses in MSun
        returns:
            return_q: array of mass ratios of the same length as m and p
        """
            
        q = np.arange(0.1, 1.001, 0.001)
        
        return_q = np.zeros(len(p))
        
        def get_q(q, _p, m_range=-1):
            
            _return_q = np.zeros(len(_p))

            pdfs = pdf(q, _p, m_range) # One per period range
        
            log_p_centres = np.array([1, 3, 5, 7])
            args = abs(np.log10(_p)[:, None] - log_p_centres[None, :]).argmin(axis=-1)
        
            # Per period range
            for i, arg in enumerate(np.unique(args)):
                _pdf = pdfs[i]
                p_mask = np.where(args == arg)[0] 
            
                cdf = np.cumsum(_pdf)
                cdf = cdf/cdf[-1] # Normalize the cdf
            
                percentiles = np.random.uniform(size=len(p[p_mask]))
                q_args = abs(percentiles[:, None] - cdf[None, :]).argmin(axis=-1)
            
                _return_q[p_mask] = q[q_args]
            
            return _return_q
                
        # O stars
        _mask = np.where(m >= 16)[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=-1)            
        # Early B stars
        _mask = np.where((m >= 9) & (m < 16))[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=-2) 
        # Mid B stars
        _mask = np.where((m >= 5) & (m < 9))[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=-3) 
        # AB stars
        _mask = np.where((m >= 1.6) & (m < 5))[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=-4) 
        # Solar-type stars
        _mask = np.where((m >= 0.6) & (m < 1.6))[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=-5) 
        # 0.3-0.6 MSun M dwarfs
        _mask = np.where((m >= 0.3) & (m < 0.6))[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=1) 
        # < 0.3 MSun M dwarfs
        _mask = np.where(m < 0.3)[0]
        return_q[_mask] = get_q(q, p[_mask], m_range=0) 

        # We can have companions below the hydrogen burning limit
        # or below the minimum stellar mass in the simulation.
        # For those stars, set their mass to the minimum mass.
        _mask = np.where(m*return_q < mmin)[0]
        return_m = np.ones(len(_mask)) * mmin
        return_q[_mask] = return_m / m[_mask]

        return return_q
    
    q = inv_transform_sampling(p, m)

    # Verify that all mass ratios and companion masses are within bounds
    assert np.min(q) >=	0.1
    assert np.max(q) <=	1
    assert np.allclose(np.min(m*q), mmin)

    return q


def get_eccentricities(m, p, edist='field'):
    """
    Sample the eccentricity based on the distributions reported by Moe & di Stefano (2017),
    for the specified mass range, for a given period. 
    args:
        m: array of masses in MSun
        p: array of orbital periods in days
        e_dist: distribution of eccentricities, only possible choice is 'field'
    returns: 
        ecc: array of eccentricities of the same length as m and p
    """
    
    def ecc_max(p_range=0):
        """
        Calculate the maximum eccentricity for the minimum period in a period bin
        Equation evaluated at logP=0.55, 0.6, 0.7, 0.85, 1, 1.3, 1.6, 2, 2.5, 3.5, 4.5, 5.5
        """
        e_max = np.array([0.318, 0.368, 0.458, 0.569, 0.658, 0.748, 0.815, 0.864, 0.926, 0.966, 0.993, 0.998, 1])
    
        return e_max[p_range]
    
    def pdf(ecc_values, p_range=0, m_range=-1):
        """
        Calculate the power law exponent eta of the e^eta distribution for a period bin
        Equation evaluated at logP=0.55, 0.6, 0.7, 0.85, 1, 1.3, 1.6, 2, 2.5, 3.5, 4.5, 5.5
        """
        # Mass ranges: <= 5 MSun and > 5 MSun
        eta = np.array([[-13.4, -6.4, -2.9, -1.4, -0.8, -0.4, -0.178, -0.036, 0.133, 0.25, 0.367, 0.425, 0.46],
                        [-3.1, -1.1, -0.1, 0.329, 0.5, 0.614, 0.678, 0.718, 0.767, 0.8, 0.833, 0.85, 0.86]])
        
        probabilities = ecc_values**eta[m_range][p_range]
            
        return probabilities
        
    def inv_transform_sampling(p, p_range=0, m_range=-1):
        """
        Use inverse transform sampling to sample the eccentricity pdf for a given mass and period range.
        args:
            p: array of orbital periods in days
            p_range: period range between 0 and 11
            m_range: mass range between 0 (< 5 MSun) and 1
        returns:
            ecc: array of eccentricities of the same length as m and p
        """
            
        ecc = np.arange(0.001, 1., 0.001)
        emax = ecc_max(p_range)
        
        _pdf = pdf(ecc[ecc <= emax], p_range, m_range)
            
        cdf = np.cumsum(_pdf)
        cdf = cdf/cdf[-1] # Normalize the cdf
            
        percentiles = np.random.uniform(size=len(p))
        args = abs(percentiles[:, None] - cdf[None, :]).argmin(axis=-1)
        return_ecc = ecc[args]
            
        return return_ecc

    ecc = np.zeros(len(m))
    
    # Late-type stars
    # Leave periods below logP = 0.55 circularized
    # logP < 0.6
    _mask = np.where((m <= 5) & (np.log10(p) < 0.6) & (np.log10(p) >= 0.55))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=0, m_range=0)
    # 0.6 < logP < 0.7
    _mask = np.where((m <= 5) & (0.6 <= np.log10(p)) & (np.log10(p) < 0.7))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=1, m_range=0)
    # 0.7 < logP < 0.85
    _mask = np.where((m <= 5) & (0.7 <= np.log10(p)) & (np.log10(p) < 0.85))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=2, m_range=0)
    # 0.85 < logP < 1
    _mask = np.where((m <= 5) & (0.85 <= np.log10(p)) & (np.log10(p) < 1))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=3, m_range=0)
    # 1 < logP < 1.2
    _mask = np.where((m <= 5) & (1 <= np.log10(p)) & (np.log10(p) < 1.2))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=4, m_range=0)
    # 1.2 < logP < 1.4
    _mask = np.where((m <= 5) & (1.2 <= np.log10(p)) & (np.log10(p) < 1.4))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=5, m_range=0)
    # 1.4 < logP < 1.6
    _mask = np.where((m <= 5) & (1.4 <= np.log10(p)) & (np.log10(p) < 1.6))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=6, m_range=0)
    # 1.6 < logP < 2
    _mask = np.where((m <= 5) & (1.6 <= np.log10(p)) & (np.log10(p) < 2))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=7, m_range=0)
    # 2 < logP < 2.5
    _mask = np.where((m <= 5) & (2 <= np.log10(p)) & (np.log10(p) < 2.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=8, m_range=0)
    # 2.5 < logP < 3.5
    _mask = np.where((m <= 5) & (2.5 <= np.log10(p)) & (np.log10(p) < 3.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=9, m_range=0)
    # 3.5 < logP < 4.5
    _mask = np.where((m <= 5) & (3.5 <= np.log10(p)) & (np.log10(p) < 4.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=10, m_range=0)
    # logP > 4.5
    _mask = np.where((m <= 5) & (4.5 <= np.log10(p)))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=-1, m_range=0)
    
    # Early-type stars
    # Leave periods below logP = 0.55 circularized
    # logP < 0.6
    _mask = np.where((m > 5) & (np.log10(p) < 0.6) & (np.log10(p) >= 0.55))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=0, m_range=1)
    # 0.6 < logP < 0.7
    _mask = np.where((m > 5) & (0.6 <= np.log10(p)) & (np.log10(p) < 0.7))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=1, m_range=1)
    # 0.7 < logP < 0.85
    _mask = np.where((m > 5) & (0.7 <= np.log10(p)) & (np.log10(p) < 0.85))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=2, m_range=1)
    # 0.85 < logP < 1
    _mask = np.where((m > 5) & (0.85 <= np.log10(p)) & (np.log10(p) < 1))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=3, m_range=1)
    # 1 < logP < 1.2
    _mask = np.where((m > 5) & (1 <= np.log10(p)) & (np.log10(p) < 1.2))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=4, m_range=1)
    # 1.2 < logP < 1.4
    _mask = np.where((m > 5) & (1.2 <= np.log10(p)) & (np.log10(p) < 1.4))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=5, m_range=1)
    # 1.4 < logP < 1.6
    _mask = np.where((m > 5) & (1.4 <= np.log10(p)) & (np.log10(p) < 1.6))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=6, m_range=1)
    # 1.6 < logP < 2
    _mask = np.where((m > 5) & (1.6 <= np.log10(p)) & (np.log10(p) < 2))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=7, m_range=1)
    # 2 < logP < 2.5
    _mask = np.where((m > 5) & (2 <= np.log10(p)) & (np.log10(p) < 2.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=8, m_range=1)
    # 2.5 < logP < 3.5
    _mask = np.where((m > 5) & (2.5 <= np.log10(p)) & (np.log10(p) < 3.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=9, m_range=1)
    # 3.5 < logP < 4.5
    _mask = np.where((m > 5) & (3.5 <= np.log10(p)) & (np.log10(p) < 4.5))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=10, m_range=1)
    # logP > 4.5
    _mask = np.where((m > 5) & (4.5 <= np.log10(p)))[0]
    ecc[_mask] = inv_transform_sampling(p[_mask], p_range=-1, m_range=1)

    # Verify that all eccentricity values are within bounds
    assert np.min(ecc) >= 0
    assert np.max(ecc) < 1 # Note that we require e < 1 rather than e <= 1
    
    return ecc



def orbits(mass_array, mult_frac='field', pdist='inner', qdist='field', edist='field', min_mass=0.08):
    """
    Generate binaries from array of stellar masses.
    args:
        mass_array: stellar masses in MSun
        mult_frac: choice of multiplicity fraction 
        pdist: choice of period distribution
        qdist: choice of semi-major axis distribution
        edist: choice of eccentricity distribution
        min_mass: minimum companion mass in MSun
    returns:
        masses: array of primary masses in MSun
        system_masses: array of system masses in MSun
        positions: array of positions relative to the COM, in cm
        velocities: array of velocities relative to the COM, in cm/s
    """

    p_IDs, s_IDs  = get_multiplicity(mass_array, mult_frac=mult_frac)

    # Insert companions with binaries
    which_s = np.ones(len(mass_array))
    # Leave primaries as 1's
    # Set singles to 0's
    which_s[s_IDs] = np.zeros(len(which_s[s_IDs]))
    # Set companions to -1's
    which_s = np.insert(which_s, np.arange(len(p_IDs))[p_IDs] + 1, -1) # Insert companions after primaries
    _mask_s = np.where(which_s == 0)
    _mask_p = np.where(which_s == 1)
    _mask_c = np.where(which_s == -1)
    
    masses        = np.zeros(len(which_s))
    system_masses = np.zeros(len(which_s))
    positions     = np.zeros((len(which_s), 3))
    velocities    = np.zeros((len(which_s), 3))
    
    primaries      = mass_array[p_IDs]
    singles        = mass_array[s_IDs]
    periods        = get_periods(primaries, pdist=pdist)
    mass_ratios    = get_mass_ratios(primaries, periods, qdist=qdist, mmin=min_mass)
    companions     = primaries * mass_ratios
    semimajor_axes = orbital_period_to_semimajor_axis(periods | units.day, primaries | units.MSun, companions | units.MSun)
    eccentricities = get_eccentricities(primaries, periods, edist=edist)

    # Place the binary at the random time in its orbit, and choose a random orientation
    E = np.random.uniform(-1 * np.pi, np.pi, size=len(primaries))
    true_anomalies                   = true_anomaly_from_eccentric_anomaly(E, eccentricities) | units.rad
    inclinations                     = np.random.uniform(-np.pi / 2, np.pi / 2, size=len(primaries)) | units.rad
    longitudes_of_the_ascending_node = np.random.vonmises(np.pi, 0, size=len(primaries)) | units.rad
    arguments_of_periapsis           = np.random.vonmises(np.pi, 0, size=len(primaries)) | units.rad
    
    rel_pos, rel_vel = rel_posvel_arrays_from_orbital_elements(primaries | units.MSun, companions | units.MSun, 
                                                               semimajor_axes, eccentricities, true_anomalies, 
                                                               inclinations, longitudes_of_the_ascending_node, 
                                                               arguments_of_periapsis, G = units.constants.G)
    
    # Offset by COM
    COM_pos = center_of_mass_array(rel_pos, primaries | units.MSun, companions | units.MSun)
    COM_vel = center_of_mass_array(rel_vel, primaries | units.MSun, companions | units.MSun)

    # Set the masses
    masses[_mask_s] = singles
    masses[_mask_p] = primaries
    masses[_mask_c] = companions
    # Set the system masses to Mtot for primaries and 0 for companions
    system_masses[_mask_s] = singles
    system_masses[_mask_p] = primaries + companions
    system_masses[_mask_c] = np.zeros(len(_mask_c))
    # Set the positions and velocities for binaries only
    positions[_mask_p]  = (-1*COM_pos).value_in(units.cm)
    positions[_mask_c]  = (rel_pos - COM_pos).value_in(units.cm)
    velocities[_mask_p] = (-1*COM_vel).value_in(units.cm/units.s)
    velocities[_mask_c] = (rel_vel - COM_vel).value_in(units.cm/units.s)
    
    return masses, system_masses, positions, velocities


if __name__ == '__main__':
    pass
