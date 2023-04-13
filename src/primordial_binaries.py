"""
Binary generation algorithm, making use of statistics by Moe & Di Stefano (2017), Winters et al. (2019) and Offner et al. (2022)
Claude Cournoyer-Cloutier, McMaster University, 2020, 2021, 2023
Version used in CCC+21 took the log twice in the period distribution, resulting in an absence of close massive binaries --> Fixed, 11/2021, CCC
"""

import numpy as np
import random
from amuse.lab import units
from amuse.ext.orbital_elements import generate_binaries, true_anomaly_from_eccentric_anomaly


def get_multiplicity(m_arr, binaries=True, mult_frac='field'):

    def interpolate(m_low, m_high, CF_low, CF_high, m):
        a = (CF_high - CF_low) / (m_high - m_low)
        b = CF_high - a * m_high
        return a * m + b

    def companion_frequency(m):

        if mult_frac == 'field':
    
            """
            For masses below 0.6 MSun, we use the primary mass dependent binary fraction from Winters et al. (2019)
            as reported and corrected in Offner et al. (2022). Above 0.8 MSun, we use the multiple star fraction
            (binary + triple/quad fraction). Between mass bins, we interpolate.
            """
            
            if m < 0.15:
                CF = 0.19
            elif 0.15 <= m < 0.3:
                CF = 0.23
            elif 0.3 <= m < 0.6:
                CF = 0.30
            elif 0.6 <= m < 0.8:
                CF = interpolate(0.6, 0.8, 0.282, 0.4, m)                        
            elif 0.8 <= m < 1.2:
                CF = 0.4                                
            elif 1.2 <= m < 2:
                CF = interpolate(1.2, 2, 0.4, 0.59, m)                                                
            elif 2 <= m < 5:
                CF = 0.59                                                
            elif 5 <= m < 9:
                CF = 0.76                                                        
            elif 9 <= m < 16:
                CF = 0.84                                                                
            elif m >= 16:
                CF = 0.94
                                                                                        
            return CF

        else:
            print('Please select a valid argument for the multiplicity. Options are \'field\' and TBD.')

    multiplicity = []
    singles      = []
    primaries    = []

    if (binaries):
        for m in m_arr:
            mult_prob = random.random()
            if mult_prob <= companion_frequency(m):
                multiplicity.append(1)
                primaries.append(m)
            else:
                multiplicity.append(0)
                singles.append(m)
                
    else:
        for m in m_arr:
            multiplicity.append(0)
            singles.append(m)

    return multiplicity, singles, primaries



def get_period(mass, pdist='field'):
    
    """
    For masses below 0.6 MSun, we use the lognormal period distributions from Winters et al. (2019).
    The means are those reported in Offner et al. (2022) and the standard deviations are obtained
    from the FWHM = 2.355 sigma, from the figures shown in Winters et al. (2019).
    For masses above 0.6 MSun, we use the period distributions from Moe & Di Stefano. We extend the
    period distrbutions up and down to 1.6 MSun from 1.2 and 2 MSun.
    """
    
    def interpolate(p_high, p_low, prob_high, prob_low, p):
        a = (prob_high - prob_low) / (10 ** p_high - 10 ** p_low)
        b = prob_high - a * p_high
        return a * p + b
    
    
    def m_dwarfs(m):
        
        def probability(m):
            
            if m < 0.15:
                log_period = np.log10(365.25) + (3./2) * np.random.normal(0.845, 1.04)
            elif 0.15 <= m < 0.30:
                log_period = np.log10(365.25) + (3./2) * np.random.normal(1.04, 1.21)
            elif 0.30 <= m < 0.60:
                log_period = np.log10(365.25) + (3./2) * np.random.normal(1.69, 1.48)
            
            return log_period
        
        def period(m):
            
            p = 0
            while (p < 0.5) or (p > 7.5):
                p = probability(m)
            
            return p
        
        period_m_dwarf = period(m)
        return period_m_dwarf
    
    
    def solar_and_above(m):
        
        def probability(m, x):
            
            if 0.6 <= m < 1.6:

                if 0.5 <= x < 1.5:
                    prob = 0.027
                elif 1.5 <= x < 2.5:
                    prob = interpolate(2.5, 1.5, 0.057, 0.027, x)
                elif 2.5 <= x < 3.5:
                    prob = 0.057
                elif 3.5 <= x < 4.5:
                    prob = interpolate(4.5, 3.5, 0.095, 0.057, x)
                elif 4.5 <= x < 5.5:
                    prob = 0.095
                elif 5.5 <= x < 6.5:
                    prob = interpolate(6.5, 5.5, 0.075, 0.095, x)
                elif 6.5 <= x < 7.5:
                    prob = 0.075
                else:
                    prob = 0
        
        
            elif 1.6 <= m < 5:
                if 0.5 <= x < 1.5:
                    prob = 0.07
                elif 1.5 <= x < 2.5:
                    prob = interpolate(2.5, 1.5, 0.12, 0.07, x)
                elif 2.5 <= x < 3.5:
                    prob = 0.12
                elif 3.5 <= x < 4.5:
                    prob = interpolate(4.5, 3.5, 0.13, 0.12, x)
                elif 4.5 <= x < 5.5:
                    prob = 0.13
                elif 5.5 <= x < 6.5:
                    prob = interpolate(6.5, 5.5, 0.09, 0.13, x)
                elif 6.5 <= x < 7.5:
                    prob = 0.09
                else:
                    prob = 0
    
            elif 5 <= m < 9:
                if 0.5 <= x < 1.5:
                    prob = 0.14
                elif 1.5 <= x < 2.5:
                    prob = interpolate(2.5, 1.5, 0.22, 0.14, x)
                elif 2.5 <= x < 3.5:
                    prob = 0.22
                elif 3.5 <= x < 4.5:
                    prob = interpolate(4.5, 3.5, 0.20, 0.22, x)
                elif 4.5 <= x < 5.5:
                    prob = 0.20
                elif 5.5 <= x < 6.5:
                    prob = interpolate(6.5, 5.5, 0.11, 0.20, x)
                elif 6.5 <= x < 7.5:
                    prob = 0.11
                else:
                    prob = 0
                        
            elif 9 <= m < 16:
                if 0.5 <= x < 1.5:
                    prob = 0.19
                elif 1.5 <= x < 2.5:
                    prob = interpolate(2.5, 1.5, 0.26, 0.19, x)
                elif 2.5 <= x < 3.5:
                    prob = 0.26
                elif 3.5 <= x < 4.5:
                    prob = interpolate(4.5, 3.5, 0.23, 0.26, x)
                elif 4.5 <= x < 5.5:
                    prob = 0.23
                elif 5.5 <= x < 6.5:
                    prob = interpolate(6.5, 5.5, 0.13, 0.23, x)
                elif 6.5 <= x < 7.5:
                    prob = 0.13
                else:
                    prob = 0
                                                                                            
            elif m >= 16:
                if 0.5 <= x < 1.5:
                    prob = 0.29
                elif 1.5 <= x < 2.5:
                    prob = interpolate(2.5, 1.5, 0.32, 0.29, x)
                elif 2.5 <= x < 3.5:
                    prob = 0.32
                elif 3.5 <= x < 4.5:
                    prob = interpolate(4.5, 3.5, 0.30, 0.32, x)
                elif 4.5 <= x < 5.5:
                    prob = 0.30
                elif 5.5 <= x < 6.5:
                    prob = interpolate(6.5, 5.5, 0.18, 0.30, x)
                elif 6.5 <= x < 7.5:
                    prob = 0.18
                else:
                    prob = 0
                                                                                                                                                                
            return prob
                                                                                                                                                                    
                                                                                                                                                                    
        def period(m):
                                                                                                                                                                        
            p = 1
            h = 1
            while probability(m, p) < h:
                 p = random.uniform(0.5, 7.5)
                 h = random.uniform(0, 0.32)
                                                                                                                                                                                        
            return p

        period_above_solar = period(m)
        return period_above_solar

    if pdist == 'field':
        if mass < 0.6:
            p = m_dwarfs(mass)
        else:
            p = solar_and_above(mass)
    else:
        print('Please select a valid argument for the period distribution. Options are \'field\' and TBD.')

    return p


def get_companion_mass(mass, period, qdist='field'):
    
    """
    From Moe & di Stefano (2017), except for M-dwarfs from Winters et al. (2019)
    M-dwarfs have a uniform distribution of mass ratios; the mass ratio is already
    selected in the period step for M-dwarfs. The mass bins are extended, there is
    no interpolation.
    """
    
    # Between 0.6 and 1.6 MSun
    def prob_solar(P, x):
        
        # Period ranges
        if 0.5 <= P < 6:
            g_small = 0.3
            g_large = -0.5
        else:
            g_small = 0.3
            g_large = -1.1
        
        # Match the values at q=0.3
        corr_fac = 0.3**(g_small-g_large)
        
        # Excess twin fraction
        if x >= 0.95:
            if P <= 2:
                twin = 0.3
            elif 2 < P <= 4:
                twin = 0.2
            elif 4 < P <= 6:
                twin = 0.1
            else:
                twin = 0
        else:
            twin = 0
        
        # Probabilities
        if x < 0.3:
            prob = x**g_small + twin
        else:
            prob = corr_fac * x**g_large + twin
        
        return prob
    
    # Between 1.6 and 5 MSun
    def prob_AB(P, x):
        
        # Period ranges
        if 0.5 <= P < 2:
            g_small = 0.2
            g_large = -0.5
        elif 2 <= P < 4:
            g_small = 0.1
            g_large = -0.9
        elif 4 <= P < 6:
            g_small = -0.5
            g_large = -1.4
        else:
            g_small = -1.0
            g_large = -2.0
        
        # Match the values at q=0.3
        corr_fac = 0.3**(g_small-g_large)
        
        # Excess twin fraction
        if x >= 0.95:
            if P <= 2:
                twin = 0.22
            elif 2 < P <= 4:
                twin = 0.10
            else:
                twin = 0
        else:
            twin = 0
        
        # Probabilities
        if x < 0.3:
            prob = x**g_small + twin
        else:
            prob = corr_fac * x**g_large + twin
        
        return prob
    
    # Between 5 and 9 MSun
    def prob_midB(P, x):
        
        # Period ranges
        if 0.5 <= P < 2:
            g_small = 0.1
            g_large = -0.5
        elif 2 <= P < 4:
            g_small = -0.2
            g_large = -1.7
        elif 4 <= P < 6:
            g_small = -1.2
            g_large = -2.0
        else:
            g_small = -1.5
            g_large = -2.0
        
        # Match the values at q=0.3
        corr_fac = 0.3**(g_small-g_large)
        
        # Excess twin fraction
        if x >= 0.95:
            if P <= 2:
                twin = 0.17
            else:
                twin = 0
        else:
            twin = 0
        
        # Probabilities
        if x < 0.3:
            prob = x**g_small + twin
        else:
            prob = corr_fac * x**g_large + twin
        
        return prob
    
    # Between 9 and 16 MSun
    def prob_earlyB(P, x):
        
        # Period ranges
        if 0.5 <= P < 2:
            g_small = 0.1
            g_large = -0.5
        elif 2 <= P < 4:
            g_small = -0.2
            g_large = -1.7
        elif 4 <= P < 6:
            g_small = -1.2
            g_large = -2.0
        else:
            g_small = -1.5
            g_large = -2.0
        
        # Match the values at q=0.3
        corr_fac = 0.3**(g_small-g_large)
        
        # Excess twin fraction
        if x >= 0.95:
            if P <= 2:
                twin = 0.14
            else:
                twin = 0
        else:
            twin = 0
        
        # Probabilities
        if x < 0.3:
            prob = x**g_small + twin
        else:
            prob = corr_fac * x**g_large + twin
        
        return prob
    
    # Above 16 MSun
    def prob_O(P, x):
        
        # Period ranges
        if 0.5 <= P < 2:
            g_small = 0.1
            g_large = -0.5
        elif 2 <= P < 4:
            g_small = -0.2
            g_large = -1.7
        elif 4 <= P < 6:
            g_small = -1.2
            g_large = -2.0
        else:
            g_small = -1.5
            g_large = -2.0
        
        # Match the values at q=0.3
        corr_fac = 0.3**(g_small-g_large)
        
        # Excess twin fraction
        if x >= 0.95:
            if P <= 2:
                twin = 0.08
            else:
                twin = 0
        else:
            twin = 0
        
        # Probabilities
        if x < 0.3:
            prob = x**g_small + twin
        else:
            prob = corr_fac * x**g_large + twin
        
        return prob
    
    def mass_ratio(m, P, q_temp):
        
        q = 0
        h = 10
        prob = 0
        
        while prob < h:
            low = np.max([0.1, 0.08 / m])
            q = rng.uniform(low, 1)
            if m < 0.6:
                q = q_temp
                # Exit loop
                prob = 1
                h = 0
            elif 0.6 < m <= 1.6:
                prob = prob_solar(P, q)
                h = rng.uniform(0, np.max([prob_solar(P, 0.1), prob_solar(P, 0.3), prob_solar(P, 1)]))
            elif 1.6 < m <= 5:
                prob = prob_AB(P, q)
                h = rng.uniform(0, np.max([prob_AB(P, 0.1), prob_AB(P, 0.3), prob_AB(P, 1)]))
            elif 5 < m < 9:
                prob = prob_midB(P, q)
                h = rng.uniform(0, np.max([prob_midB(P, 0.1), prob_midB(P, 0.3), prob_midB(P, 1)]))
            elif 9 < m < 16:
                prob = prob_earlyB(P, q)
                h = rng.uniform(0, np.max([prob_earlyB(P, 0.1), prob_earlyB(P, 0.3), prob_earlyB(P, 1)]))
            else:
                prob = prob_O(P, q)
                h = rng.uniform(0, np.max([prob_O(P, 0.1), prob_O(P, 0.3), prob_O(P, 1)]))
            
            mass_ratio = q
        
        return mass_ratio

    
    if qdist == 'field':
        mr = mass_ratio(mass, period)
        cm = mass * mr
    else:
        print('Please select a valid argument for the mass ratio distribution. Options are \'field\' and TBD.')
    
    return cm



def get_eccentricity(mass, period, edist = 'field'):
    
    def ecc_max(p):
        e_max = 1 - (10 ** p / 2) ** (-2 / 3)
        return e_max
    
    def prob_eta(m, p):
        if m <= 5:
            eta = 0.6 - 0.7 / (p - 0.5)
        else:
            eta = 0.9 - 0.2 / (p - 0.5)
        return eta
    
    def get_ecc(m, P):
        if P <= 0.55:
            ecc = random.uniform(0, ecc_max(P))
        else:
            t = 0
            e = 1
            h = 10
            n = 0
            while (e ** n) < h:
                emax = ecc_max(P)
                ecc  = random.uniform(0, emax)
                n    = prob_eta(m, P)
                h    = random.uniform(0, emax ** n)
        return ecc
    
    if edist == 'field':
        eccentricity = get_ecc(mass, period)
    else:
        print('Please select a valid argument for the eccentricity distribution. Options are \'field\' and TBD.')
        
    return eccentricity



def orbits(mass_array, binaries=True, mult_frac='field', pdist='field', qdist='field', edist='field'):

    def semi_major_axis_from_period(primary_mass, companion_mass, log_period):
        """
        returns the semi-major axes from the primary and companion masses and the period
        uses eq. 2.45 from van de Kamp (1964) chap. 8.4
        """
        primary_mass = primary_mass
        companion_mass = companion_mass
        period = 10 ** log_period | units.day
        semi_major_axis = (units.constants.G * (primary_mass + companion_mass)) ** (1./3) * period ** (2./3)

        return semi_major_axis.as_quantity_in(units.AU)


    def generate_binaries_with_orientation(mass_array = mass_array, binaries = binaries, mult_frac = mult_frac, pdist = pdist, qdist = qdist, edist = edist):
        """
        Generates binaries from arrays of masses
        The input array has no units
        """

        multiplicity  = get_multiplicity(mass_array, binaries, mult_frac)[0]
        masses        = []
        system_masses = []
        positions     = []
        velocities    = []

        for i in range(len(multiplicity)):

            if multiplicity[i] == 1:
                primary_mass    = mass_array[i]
                log_period      = get_period(primary_mass, pdist)
                companion_mass  = get_companion_mass(primary_mass, log_period, qdist) | units.MSun
                eccentricity    = get_eccentricity(primary_mass, log_period, edist)
                primary_mass    = primary_mass | units.MSun
                semi_major_axis = semi_major_axis_from_period(primary_mass, companion_mass, log_period)

                E = random.uniform(-1 * np.pi, np.pi)
                true_anomaly                    = true_anomaly_from_eccentric_anomaly(E, eccentricity) | units.rad
                inclination                     = random.uniform(-np.pi / 2, np.pi / 2) | units.rad
                longitude_of_the_ascending_node = random.vonmisesvariate(np.pi, 0) | units.rad
                argument_of_periapsis           = random.vonmisesvariate(np.pi, 0) | units.rad

                binary = generate_binaries(primary_mass, companion_mass, semi_major_axis, eccentricity, true_anomaly, inclination, longitude_of_the_ascending_node, argument_of_periapsis, G = units.constants.G)

                primary_set = binary[0]
                primary = primary_set[0]
                companion_set = binary[1]
                companion = companion_set[0]

                masses.append(primary.mass.value_in(units.MSun))
                masses.append(companion.mass.value_in(units.MSun))
                system_masses.append(primary.mass.value_in(units.MSun) + companion.mass.value_in(units.MSun))
                system_masses.append(0)
                positions.append(primary.position.value_in(units.cm))
                positions.append(companion.position.value_in(units.cm))
                velocities.append(primary.velocity.value_in(units.cm / units.s))
                velocities.append(companion.velocity.value_in(units.cm / units.s))

            else:
                mass = mass_array[i]
                masses.append(mass)
                system_masses.append(mass)
                positions.append([0, 0, 0])
                velocities.append([0, 0 ,0])

        return masses, system_masses, positions, velocities
    
    binaries = generate_binaries_with_orientation()
    
    return binaries


if __name__ == '__main__':
    pass
