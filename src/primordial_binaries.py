"""
Binary generation algorithm, making use of statistics by Moe & Di Stefano (2017) and Winters et al. (2019)
Claude Cournoyer-Cloutier, McMaster University, 2020
The functions get_multiplicity, get_period and get_eccentricity have been reworked
There is still some work to be done on get_companion_masses, which still defaults the gamma (June 3, 2020)
"""

import numpy as np
import random
from amuse.lab import units
from amuse.ext.orbital_elements import generate_binaries, true_anomaly_from_eccentric_anomaly


def get_multiplicity(m_arr, binaries='field'):

    def interpolate(m_low, m_high, CF_low, CF_high, m):
        a = (CF_high - CF_low) / (m_high - m_low)
        b = CF_high - a * m_high
        return a * m + b

    def companion_frequency(m):
    
        """
        For masses below 0.6 M_sun, we use the primary mass dependent binary fraction from Winters et al. (2020)
        Above 0.8 M_sun, we use these from Moe & Di Stefano (2017) -- we use binary fraction + triple/quad frac
        Between mass bins, we interpolate
        """
            
            
        if m < 0.15:
            CF = 0.16
                
        if 0.15 <= m < 0.3:
            CF = 0.214
                        
        if 0.3 <= m < 0.6:
            CF = 0.282
                                
        if 0.6 <= m < 0.8:
            CF = interpolate(0.6, 0.8, 0.282, 0.4, m)
                                        
        if 0.8 <= m < 1.2:
            CF = 0.4
                                                
        if 1.2 <= m < 2:
            CF = interpolate(1.2, 2, 0.4, 0.59, m)
                                                        
        if 2 <= m < 5:
            CF = 0.59
                                                                
        if 5 <= m < 9:
            CF = 0.76
                                                                        
        if 9 <= m < 16:
            CF = 0.84
                                                                                
        if m >= 16:
            CF = 0.94
                                                                                        
        return CF


    multiplicity = []
    singles      = []
    primaries    = []

    if binaries in ['field', 'Field', 'FIELD']:
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



def get_period(mass):
    
    """
    For masses below 0.6 M_sun, we use the lognormal period distributions from Winters et al.
    For masses above 0.6 M_sun, we use the period distributions from Moe & Di Stefano
    We extend the period distrbutions up and down to 1.6 M_sun from 1.2 and 2 M_sun
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
        
        # Out of all the probabilities defined below, the maximum is 0.32
        
        def probability(m, x):
            
            if 0.6 <= m < 1.6:
                if 0.5 <= np.log10(x) < 1.5:
                    prob = 0.027
                elif 1.5 <= np.log10(x) < 2.5:
                    prob = interpolate(2.5, 1.5, 0.057, 0.027, x)
                elif 2.5 <= np.log10(x) < 3.5:
                    prob = 0.057
                elif 3.5 <= np.log10(x) < 4.5:
                    prob = interpolate(4.5, 3.5, 0.095, 0.057, x)
                elif 4.5 <= np.log10(x) < 5.5:
                    prob = 0.095
                elif 5.5 <= np.log10(x) < 6.5:
                    prob = interpolate(6.5, 5.5, 0.075, 0.095, x)
                elif 6.5 <= np.log10(x) < 7.5:
                    prob = 0.075
                else:
                    prob = 0
        
        
            elif 1.6 <= m < 5:
                if 0.5 <= np.log10(x) < 1.5:
                    prob = 0.07
                elif 1.5 <= np.log10(x) < 2.5:
                    prob = interpolate(2.5, 1.5, 0.12, 0.07, x)
                elif 2.5 <= np.log10(x) < 3.5:
                    prob = 0.12
                elif 3.5 <= np.log10(x) < 4.5:
                    prob = interpolate(4.5, 3.5, 0.13, 0.12, x)
                elif 4.5 <= np.log10(x) < 5.5:
                    prob = 0.13
                elif 5.5 <= np.log10(x) < 6.5:
                    prob = interpolate(6.5, 5.5, 0.09, 0.13, x)
                elif 6.5 <= np.log10(x) < 7.5:
                    prob = 0.09
                else:
                    prob = 0
    
            elif 5 <= m < 9:
                if 0.5 <= np.log10(x) < 1.5:
                    prob = 0.14
                elif 1.5 <= np.log10(x) < 2.5:
                    prob = interpolate(2.5, 1.5, 0.22, 0.14, x)
                elif 2.5 <= np.log10(x) < 3.5:
                    prob = 0.22
                elif 3.5 <= np.log10(x) < 4.5:
                    prob = interpolate(4.5, 3.5, 0.20, 0.22, x)
                elif 4.5 <= np.log10(x) < 5.5:
                    prob = 0.20
                elif 5.5 <= np.log10(x) < 6.5:
                    prob = interpolate(6.5, 5.5, 0.11, 0.20, x)
                elif 6.5 <= np.log10(x) < 7.5:
                    prob = 0.11
                else:
                    prob = 0
                        
            elif 9 <= m < 16:
                if 0.5 <= np.log10(x) < 1.5:
                    prob = 0.19
                elif 1.5 <= np.log10(x) < 2.5:
                    prob = interpolate(2.5, 1.5, 0.26, 0.19, x)
                elif 2.5 <= np.log10(x) < 3.5:
                    prob = 0.26
                elif 3.5 <= np.log10(x) < 4.5:
                    prob = interpolate(4.5, 3.5, 0.23, 0.26, x)
                elif 4.5 <= np.log10(x) < 5.5:
                    prob = 0.23
                elif 5.5 <= np.log10(x) < 6.5:
                    prob = interpolate(6.5, 5.5, 0.13, 0.23, x)
                elif 6.5 <= np.log10(x) < 7.5:
                    prob = 0.13
                else:
                    prob = 0
                                                                                            
            elif m >= 16:
                if 0.5 <= np.log10(x) < 1.5:
                    prob = 0.29
                elif 1.5 <= np.log10(x) < 2.5:
                    prob = interpolate(2.5, 1.5, 0.32, 0.29, x)
                elif 2.5 <= np.log10(x) < 3.5:
                    prob = 0.32
                elif 3.5 <= np.log10(x) < 4.5:
                    prob = interpolate(4.5, 3.5, 0.30, 0.32, x)
                elif 4.5 <= np.log10(x) < 5.5:
                    prob = 0.30
                elif 5.5 <= np.log10(x) < 6.5:
                    prob = interpolate(6.5, 5.5, 0.18, 0.30, x)
                elif 6.5 <= np.log10(x) < 7.5:
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

    if mass < 0.6:
        p = m_dwarfs(mass)
    else:
        p = solar_and_above(mass)
            
    return p


def get_companion_mass(mass, period):
    
    """
    All from Moe & Di Stefano, due to incompleteness of mass ratio distribution for M-dwarfs
    The solar-type distribution is extended down to 0.08 M_sun
    """
    
    def interpolate(p_low, p_high, g_low, g_high, p):
        a = (g_high - g_low) / (10 ** p_high - 10 ** p_low)
        b = g_high - a * p_high
        return a * 10**p + b
    
    
    def gamma_solar(P, x):
        """
        Returns the power law exponent for stellar masses between 0.08 MSun and 1.6 MSun
        """        
        if 0.3 <= x:
            if P <= 5:
                gamma = -0.5
            elif 5 < P:
                gamma = interpolate(5, 7, -0.5, -1.1, P)

        elif x < 0.3:
            gamma = 0.3
        
        else:
            print('Err. in gamma solar, default to gamma=0 !')
            print('Period =', P, 'and q =', x)
            gamma = 0

        return gamma
    
    
    def gamma_AB(P, x):
        """
        Returns the power law exponent for stellar masses between 1.6 MSun and 5 MSun
        """        
        if 0.3 <= x < 1:
            if P < 1.5:
                gamma = -0.5
            elif 1.5 <= P < 2.5:
                gamma = interpolate(1.5, 2.5, -0.5, -0.9, P)
            elif 2.5 <= P < 3.5:
                gamma = -0.9
            elif 3.5 <= P < 4.5:
                gamma = interpolate(3.5, 4.5, -0.9, -1.4, P)
            elif 4.5 <= P < 5.5:
                gamma = -1.4
            elif 5.5 <= P < 6.5:
                gamma = interpolate(5.5, 6.5, -1.4, -2.0, P)
            elif 6 <= P <= 7.5:
                gamma = -2.0
    
        elif 0.1 <= x < 0.3:
            if P < 1.5:
                gamma = 0.2
            elif 1.5 <= P < 2.5:
                gamma = interpolate(1.5, 2.5, 0.2, 0.1, P)
            elif 2.5 <= P < 3.5:
                gamma = 0.1
            elif 3.5 <= P < 4.5:
                gamma = interpolate(3.5, 4.5, 0.1, -0.5, P)
            elif 4.5 <= P < 5.5:
                gamma = -0.5
            elif 5.5 <= P < 6.5:
                gamma = interpolate(5.5, 6.5, -0.5, -1.0, P)
            elif 6 <= P <= 7.5:
                gamma = -1.0

        else:
            print('Err. in gamma AB, default to gamma=0 !')
            print('Period =', P, 'and q =', x)
            gamma = 0

        return gamma
                                
                                
    def gamma_early(P, x):
        """
        Returns the power law exponent for stellar masses above 5 MSun
        """                                    
        if 0.3 <= x < 1:
            if 0.5 <= P < 1.5:
                gamma = -0.5
            elif 1.5 <= P < 2.5:
                gamma = interpolate(1.5, 2.5, -0.5, -1.7, P)
            elif 2.5 <= P < 3.5:
                gamma = -1.7
            elif 3.5 <= P < 4.5:
                gamma = interpolate(3.5, 4.5, -1.7, -2.0, P)
            elif 4.5 <= P <= 7.5:
                gamma = -2.0
                                                                                
        elif 0.1 <= x < 0.3:
            if 0.5 <= P < 1.5:
                gamma = 0.1
            elif 1.5 <= P < 2.5:
                gamma = interpolate(1.5, 2.5, 0.1, -0.2, P)
            elif 2.5 <= P < 3.5:
                gamma = -0.2
            elif 3.5 <= P < 4.5:
                gamma = interpolate(3.5, 4.5, -0.2, -1.2, P)
            elif 4.5 <= P < 5.5:
                gamma = -1.2
            elif 5.5 <= P < 6.5:
                gamma = interpolate(5.5, 6.5, -1.2, -1.5, P)
            elif 6.5 <= P <= 7.5:
                gamma = -1.5
                                                                                                                                            
        else:
            print('Err. in gamma early, default to gamma=0 !')
            print('Period =', P, 'and q =', x)
            gamma = 0
                                                                                                                                                    
        return gamma


    def twin_solar(m, P):
        """
        Returns the excess twin fraction for stellar masses between 0.08 MSun and 1.6 MSun
        """
        if P <= 1.5:
            twin = 0.3
        elif 1.5 < P <= 2.5:
            twin = interpolate(2.5, 1.5, 0.2, 0.3, P)
        elif 2.5 < P <= 3.5:
            twin = 0.2
        elif 3.5 < P <= 4.5:
            twin = interpolate(4.5, 3.5, 0.1, 0.2, P)
        elif 4.5 < P <= 5.5:
            twin = 0.1
        elif 5.5 < P <= 6.5:
            twin = interpolate(6.5, 5.5, 0, 0.1, P)
        else:
            twin = 0
        return twin
    
    
    def twin_AB(m, P):
        """
        Returns the excess twin fraction for stellar masses between 1.6 MSun and 5 MSun
        """
        if P <= 1.5:
            twin = 0.22
        elif 1.5 < P <= 2.5:
            twin = interpolate(2.5, 1.5, 0.1, 0.22, P)
        elif 2.5 < P <= 3.5:
            twin = 0.1
        elif 3.5 < P <= 4.5:
            twin = interpolate(4.5, 3.5, 0, 0.1, P)
        else:
            twin = 0
        return twin
    
    
    def twin_midB(m, P):
        """
        Returns the excess twin fraction for stellar masses between 5 MSun and 9 MSun
        """
        if P <= 1.5:
            twin = 0.17
        elif 1.5 < P <= 2.5:
            twin = interpolate(2.5, 1.5, 0, 0.17, P)
        else:
            twin = 0
        return twin
    
    
    def twin_earlyB(m, P):
        """
        Returns the excess twin fraction for stellar masses between 9 MSun and 16 MSun
        """
        if P <= 1.5:
            twin = 0.14
        elif 1.5 < P <= 2.5:
            twin = interpolate(2.5, 1.5, 0, 0.14, P)
        else:
            twin = 0
        return twin
    
    
    def twin_O(m, P):
        """
        Returns the excess twin fraction for stellar masses above 16 MSun
        """
        if P <= 1.5:
            twin = 0.08
        elif 1.5 < P <= 2.5:
            twin = interpolate(2.5, 1.5, 0, 0.08, P)
        else:
            twin = 0
        return twin


    def mass_ratio(m, P):
        q = 1
        h = 10
        g = 0
        t = 0

        def prob(q, g, t):
            if q < 0.95:
                prob = q ** g
            else:
                prob = q ** g + t
            return prob

        while prob(q, g, t) < h:
            low = np.max([0.1, 0.08 / m])
            q = random.uniform(low, 1)
            if m <= 1.6:
                g = gamma_solar(P, q)
                t = twin_solar(m, P)
                h = random.uniform(0, low ** g)
            if 1.6 < m <= 5:
                g = gamma_AB(P, q)
                t = twin_AB(m, P)
                h = random.uniform(0, low ** g)
            if 5 < m <= 9:
                g = gamma_early(P, q)
                t = twin_midB(m, P)
                h = random.uniform(0, low ** g)
            if 9 < m <= 16:
                g = gamma_early(P, q)
                t = twin_earlyB(m, P)
                h = random.uniform(0, low ** g)
            if 16 < m:
                g = gamma_early(P, q)
                t = twin_O(m, P)
                h = random.uniform(0, low ** g)
            mass_ratio = q
        
        return mass_ratio

    
    mr = mass_ratio(mass, period)
    cm = mass * mr
    
    return cm



def get_eccentricity(mass, period):
    
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
    
    eccentricity = get_ecc(mass, period)

    return eccentricity



def orbits(mass_array, binaries='field'):

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


    def generate_binaries_with_orientation(mass_array = mass_array, binaries = binaries):
        """
        Generates binaries from arrays of masses
        The input array has no units
        """

        multiplicity  = get_multiplicity(mass_array, binaries)[0]
        masses        = []
        system_masses = []
        positions     = []
        velocities    = []

        for i in range(len(multiplicity)):

            if multiplicity[i] == 1:
                primary_mass    = mass_array[i]
                log_period      = get_period(primary_mass)
                companion_mass  = get_companion_mass(primary_mass, log_period) | units.MSun
                eccentricity    = get_eccentricity(primary_mass, log_period)
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
