"""
Binary generation algorithm, making use of statistics by Moe & Di Stefano (2017), Duchene & Kraus (2013) and Winters et al. (2019)
Claude Cournoyer-Cloutier, McMaster University, 2020
"""

import numpy as np
import random
from amuse.lab import units
from amuse.ext.orbital_elements import generate_binaries, true_anomaly_from_eccentric_anomaly


def get_multiplicity(m_arr):

    def companion_frequency(m):

        if 0.01 <= m < 0.1:
            CF = 0.22  # DK

        if 0.1 <= m < 0.5:
            CF = 0.33  # DK

        if 0.5 <= m < 0.8:
            a = (0.5 - 0.33) / (0.8 - 0.5)
            b = 0.5 - a * 0.8
            CF = a * m + b

        if 0.8 <= m < 1.2:
            CF = 0.5  # MDS

        if 1.2 <= m < 2:
            a = (0.84 - 0.5) / (2 - 1.2)
            b = 0.84 - a * 2
            CF = a * m + b

        if 2 <= m < 5:
            CF = 0.84  # MDS

        if 5 <= m < 9:
            CF = 1.3  # MDS

        if 9 <= m < 16:
            CF = 1.6  # MDS

        if m >= 16:
            CF = 2.1  # MDS

        return round(CF, 2)

    multiplicity = []
    for i in m_arr:
        N    = 100
        K    = int(companion_frequency(i) * 100)  # K ones, N-K zeros
        rd   = random.randint(0, 99)
        bin  = np.array([0] * (N - K) + [1] * (K))[rd]
        multiplicity.append(bin)

    indices_b = np.nonzero(multiplicity)
    binaries = []
    for j in indices_b:
        binaries.extend(m_arr[j])

    minus_mult = []
    for k in range(len(multiplicity)):
        minus_mult.append(multiplicity[k - 1] - 1)

    indices_s = np.nonzero(minus_mult)
    singles = []
    for j in indices_s:
        singles.extend(m_arr[j])

    return multiplicity, singles, binaries



def get_period(mass):

    def period(m):

        if 0.01 <= m < 0.15:

            P_years = (10 ** (np.random.normal(0.845, 1.04))) ** (3. / 2)
            P = np.log10(P_years * 365.25)

            if 0.5 <= P <= 7.5:
                period = P

            elif P < 0.5:
                period = 0.5

            elif P > 7.5:
                P_years = (10 ** (np.random.normal(0.845, 1.04))) ** (3. / 2)
                P = np.log10(P_years * 365.25)

                if 0.5 <= P <= 7.5:
                    period = P

                elif P < 0.5:
                    period = 0.5

                elif P > 7.5:
                    P_years = (10 ** (np.random.normal(0.845, 1.04))) ** (3. / 2)
                    P = np.log10(P_years * 365.25)

                    if 0.5 <= P <= 7.5:
                        period = P

                    elif P < 0.5:
                        period = 0.5

                    elif P > 7.5:
                        # random within 2 std
                        P_years = (10 ** (random.randint(-124, 293) / 100)) ** (3. / 2)
                        P = np.log10(P_years * 365.25)
                        period = P


        elif 0.15 <= m < 0.30:

            P_years = (10 ** (np.random.normal(1.04, 1.21))) ** (3. / 2)

            P = np.log10(P_years * 365.25)

            if 0.5 <= P <= 7.5:

                period = P


            elif P < 0.5:

                period = 0.5


            elif P > 7.5:

                P_years = (10 ** (np.random.normal(1.04, 1.21))) ** (3. / 2)

                P = np.log10(P_years * 365.25)

                if 0.5 <= P <= 7.5:

                    period = P


                elif P < 0.5:

                    period = 0.5


                elif P > 7.5:

                    P_years = (10 ** (np.random.normal(1.04, 1.21))) ** (3. / 2)

                    P = np.log10(P_years * 365.25)

                    if 0.5 <= P <= 7.5:

                        period = P


                    elif P < 0.5:

                        period = 0.5


                    elif P > 7.5:

                        # random within 2 std

                        P_years = (10 ** (random.randint(-138, 346) / 100)) ** (3. / 2)

                        P = np.log10(P_years * 365.25)

                        period = P


        elif 0.30 <= m < 0.60:

            P_years = (10 ** (np.random.normal(1.69, 1.48))) ** (3. / 2)
            P = np.log10(P_years * 365.25)

            if 0.5 <= P <= 7.5:
                period = P

            elif P < 0.5:
                period = 0.5

            elif P > 7.5:
                P_years = (10 ** (np.random.normal(1.69, 1.48))) ** (3. / 2)
                P = np.log10(P_years * 365.25)

                if 0.5 <= P <= 7.5:
                    period = P

                elif P < 0.5:
                    period = 0.5

                elif P > 7.5:
                    P_years = (10 ** (np.random.normal(1.69, 1.48))) ** (3. / 2)
                    P = np.log10(P_years * 365.25)

                    if 0.5 <= P <= 7.5:
                        period = P

                    elif P < 0.5:
                        period = 0.5

                    elif P > 7.5:
                        # random within 2 std
                        P_years = (10 ** (random.randint(-135, 474) / 100)) ** (3. / 2)
                        P = np.log10(P_years * 365.25)
                        period = P

        else:
            prob = period_prob(m)
            rd = random.randint(0, len(prob) - 1)
            period = prob[rd]

        return period

    def period_func(m, x):

        # if 0.6 <= m < 0.8:

        if 0.6 <= m < 2:

            if 0.6 <= np.log10(x) < 1.5:
                prob = 0.027

            elif 1.5 <= np.log10(x) < 2.5:
                a = (0.057 - 0.027) / (10 ** (2.5) - 10 ** (1.5))
                b = 0.057 - a * 10 ** (2.5)
                prob = a * x + b

            elif 2.5 <= np.log10(x) < 3.5:
                prob = 0.057

            elif 3.5 <= np.log10(x) < 4.5:
                a = (0.095 - 0.057) / (10 ** (4.5) - 10 ** (3.5))
                b = 0.095 - a * 10 ** (4.5)
                prob = a * x + b

            elif 4.5 <= np.log10(x) < 5.5:
                prob = 0.095

            elif 5.5 <= np.log10(x) < 6.5:
                a = (0.075 - 0.095) / (10 ** (6.5) - 10 ** (5.5))
                b = 0.075 - a * 10 ** (6.5)
                prob = a * x + b

            elif 6.5 <= np.log10(x) < 7.5:
                prob = 0.075

            else:
                prob = 0

            return round(prob, 3) * 1000


        # elif 1.2 <= m < 2:

        elif 2 <= m < 5:

            if 0.5 <= np.log10(x) < 1.5:
                prob = 0.07

            elif 1.5 <= np.log10(x) < 2.5:
                a = (0.12 - 0.07) / (10 ** (2.5) - 10 ** (1.5))
                b = 0.12 - a * 10 ** (2.5)
                prob = a * x + b

            elif 2.5 <= np.log10(x) < 3.5:
                prob = 0.12

            elif 3.5 <= np.log10(x) < 4.5:
                a = (0.13 - 0.12) / (10 ** (4.5) - 10 ** (3.5))
                b = 0.13 - a * 10 ** (4.5)
                prob = a * x + b

            elif 4.5 <= np.log10(x) < 5.5:
                prob = 0.13

            elif 5.5 <= np.log10(x) < 6.5:
                a = (0.09 - 0.13) / (10 ** (6.5) - 10 ** (5.5))
                b = 0.09 - a * 10 ** (6.5)
                prob = a * x + b

            elif 6.5 <= np.log10(x) < 7.5:
                prob = 0.09

            else:
                prob = 0

            return round(prob, 2) * 100



        elif 5 <= m < 9:

            if 0.5 <= np.log10(x) < 1.5:
                prob = 0.14

            elif 1.5 <= np.log10(x) < 2.5:
                a = (0.22 - 0.14) / (10 ** (2.5) - 10 ** (1.5))
                b = 0.22 - a * 10 ** (2.5)
                prob = a * x + b

            elif 2.5 <= np.log10(x) < 3.5:
                prob = 0.22

            elif 3.5 <= np.log10(x) < 4.5:
                a = (0.20 - 0.22) / (10 ** (4.5) - 10 ** (3.5))
                b = 0.20 - a * 10 ** (4.5)
                prob = a * x + b

            elif 4.5 <= np.log10(x) < 5.5:
                prob = 0.20

            elif 5.5 <= np.log10(x) < 6.5:
                a = (0.11 - 0.20) / (10 ** (6.5) - 10 ** (5.5))
                b = 0.11 - a * 10 ** (6.5)
                prob = a * x + b

            elif 6.5 <= np.log10(x) < 7.5:
                prob = 0.11

            else:
                prob = 0

            return round(prob, 2) * 100



        elif 9 <= m < 16:

            if 0.5 <= np.log10(x) < 1.5:
                prob = 0.19

            elif 1.5 <= np.log10(x) < 2.5:
                a = (0.26 - 0.19) / (10 ** (2.5) - 10 ** (1.5))
                b = 0.26 - a * 10 ** (2.5)
                prob = a * x + b

            elif 2.5 <= np.log10(x) < 3.5:
                prob = 0.26

            elif 3.5 <= np.log10(x) < 4.5:
                a = (0.23 - 0.26) / (10 ** (4.5) - 10 ** (3.5))
                b = 0.23 - a * 10 ** (4.5)
                prob = a * x + b

            elif 4.5 <= np.log10(x) < 5.5:
                prob = 0.23

            elif 5.5 <= np.log10(x) < 6.5:
                a = (0.13 - 0.23) / (10 ** (6.5) - 10 ** (5.5))
                b = 0.13 - a * 10 ** (6.5)
                prob = a * x + b

            elif 6.5 <= np.log10(x) < 7.5:
                prob = 0.13

            else:
                prob = 0

            return round(prob, 2) * 100



        elif m >= 16:

            if 0.5 <= np.log10(x) < 1.5:
                prob = 0.29

            elif 1.5 <= np.log10(x) < 2.5:
                a = (0.32 - 0.29) / (10 ** (2.5) - 10 ** (1.5))
                b = 0.29 - a * 10 ** (2.5)
                prob = a * x + b

            elif 2.5 <= np.log10(x) < 3.5:
                prob = 0.32

            elif 3.5 <= np.log10(x) < 4.5:
                a = (0.30 - 0.32) / (10 ** (4.5) - 10 ** (3.5))
                b = 0.30 - a * 10 ** (4.5)
                prob = a * x + b

            elif 4.5 <= np.log10(x) < 5.5:
                prob = 0.30

            elif 5.5 <= np.log10(x) < 6.5:
                a = (0.18 - 0.30) / (10 ** (6.5) - 10 ** (5.5))
                b = 0.18 - a * 10 ** (6.5)
                prob = a * x + b

            elif 6.5 <= np.log10(x) < 7.5:
                prob = 0.18

            else:
                prob = 0

            return round(prob, 2) * 100



        else:
            prob = 0
            return prob

    def period_prob(m):

        if 0.6 <= m:

            # Getting weird decimals here
            periods_err = np.arange(0.5, 7.6, 0.1)
            # Doing weird stuff to get rid of the weird decimals
            periods = []
            for i in periods_err:
                value = int(i * 10) / 10.
                periods.append(value)

            prob = []

            for i in periods:
                value = period_func(m, 10 ** i)
                val = int(value)
                norm = np.array([i] * val)
                prob.extend(norm)

            return prob

        else:
            prob = 0
            return prob

    #periods = []
    #for i in range(len(m_arr)):
        #periods.append(period(m_arr[i]))

    return period(mass)



def get_companion_mass(mass, period):

    q = np.arange(0.1, 1.01, 0.01)


    # ---------------------------------------------------------
    # Solar-type, extended from 0.5 M_sun to 2 M_sun
    # ---------------------------------------------------------

    # Power-law slope for large and small mass ratios
    def gamma_solar(P, x):
    # We extend the 3-1 slope down to 0.5, and the 7-5 up to 7.5
        if 0.3 <= x < 1:
            # Large q
            #if 0.5 <= P <= 5:
            if P <= 5:
                gamma = -0.5
                return gamma

            #elif 5 < P <= 7.5:
            elif 5 < P:
                a = (-1.1+0.5)/(10**7 - 10**5)
                b = -1.1 - a * 10**7
                gamma = a * x + b
                return gamma

        if 0.1 <= x < 0.3:
            # Small q
            gamma = 0.3
            return gamma

    # Non-normalized probability, as a function of period and mass ratio
    def prob_func_solar(P,x):

        prob = []
        if x > 0.94:
            if P <=2 :
                prob.append(np.round(130 * x ** (gamma_solar(P,x))))
            elif 2 < P <= 4:
                prob.append(np.round(120 * x ** (gamma_solar(P, x))))
            elif 4 < P <= 6:
                prob.append(np.round(110 * x ** (gamma_solar(P, x))))
            else:
                prob.append(np.round(100 * x ** (gamma_solar(P, x))))
        else:
            prob.append(np.round(100 * x ** (gamma_solar(P, x))))
        return prob

    # Non-normalized probability, as a function of period and for the complete array of mass ratios
    def prob_solar(P):
        prob = []
        for i in range(len(q)):
            prob.append(prob_func_solar(P, q[i]))
        return prob



    # ---------------------------------------------------------
    # A/late B, from 2 M_sun to 5 M_sun
    # ---------------------------------------------------------

    # Power-law slope for large and small mass ratios
    def gamma_AB(P, x):
    # (Maybe) temporarily, we use step functions
    # We extend P=1 from 0.5 to 2, P=3 from 2 to 4...
        if 0.3 <= x < 1:
            # Large q
            #if 0.5 <= P < 2:
            if P < 2:
                gamma = -0.5
                return gamma

            elif 2 <= P < 4:
                gamma = -0.9
                return gamma

            elif 4 <= P < 6:
                gamma = -1.4
                return gamma

            elif 6 <= P < 7.5:
                gamma = -2.0
                return gamma

        if 0.1 <= x < 0.3:
            # Small q
            #if 0.5 <= P < 2:
            if P < 2:
                gamma = 0.2
                return gamma

            elif 2 <= P < 4:
                gamma = 0.1
                return gamma

            elif 4 <= P < 6:
                gamma = -0.5
                return gamma

            elif 6 <= P < 7.5:
                gamma = -1.0
                return gamma

    # Non-normalized probability, as a function of period and mass ratio
    def prob_func_AB(P,x):

        prob = []
        if x > 0.94:
            if P <= 2:
                prob.append(np.round(122 * x ** (gamma_AB(P,x))))
            elif 2 < P <= 4:
                prob.append(np.round(110 * x ** (gamma_AB(P,x))))
            else:
                prob.append(np.round(100 * x ** (gamma_AB(P, x))))
        else:
            prob.append(np.round(100 * x ** (gamma_AB(P, x))))
        return prob

    # Non-normalized probability, as a function of period and for the complete array of mass ratios
    def prob_AB(P):
        prob = []
        for j in range(len(q)):
            prob.append(prob_func_AB(P, q[j]))
        return prob



    # ---------------------------------------------------------
    # Mid-B, Early-B and O, above 5 M_sun
    # ---------------------------------------------------------

    # Power-law slope for large and small mass ratios
    def gamma_early(P, x):
    # (Maybe) temporarily, we use step functions
    # We extend P=1 from 0.5 to 2, P=3 from 2 to 4...
        if 0.3 <= x < 1:
            # Large q
            if 0.5 <= P < 2:
                gamma = -0.5
                return gamma

            elif 2 <= P < 4:
                gamma = -1.7
                return gamma

            elif 4 <= P < 7.5:
                gamma = -2.0
                return gamma

        if 0.1 <= x < 0.3:
            # Small q
            if 0.5 <= P < 2:
                gamma = 0.1
                return gamma

            elif 2 <= P < 4:
                gamma = -0.2
                return gamma

            elif 4 <= P < 6:
                gamma = -1.2
                return gamma

            elif 6 <= P < 7.5:
                gamma = -1.5
                return gamma

    # Non-normalized probability, as a function of period and mass ratio
    def prob_func_early(P,x):

        prob = []
        if x > 0.94:
            if P <= 2:
                prob.append(np.round(113 * x ** (gamma_early(P,x))))
            else:
                prob.append(np.round(100 * x ** (gamma_early(P,x))))
        else:
            prob.append(np.round(100 * x ** (gamma_early(P, x))))
        return prob

    # Non-normalized probability, as a function of period and for the complete array of mass ratios
    def prob_early(P):
        prob = []
        for j in range(len(q)):
            prob.append(prob_func_early(P,q[j]))
        return prob



    # ---------------------------------------------------------
    # Express as function of mass, and as a probability function
    # ---------------------------------------------------------
    def non_norm_prob(P,m):

        if m < 2:
            return np.array(prob_solar(P))[:,0]

        elif 2 <= m <= 5:
            return np.array(prob_AB(P))[:,0]

        elif m > 5:
            return np.array(prob_early(P))[:,0]

    def prob(P,m):

        prob = []

        for k in range(len(non_norm_prob(P,m))):
            np.array(prob.append([k] * int(non_norm_prob(P,m)[k])))

        flat = []
        for sublist in prob:
            for item in sublist:
                flat.append(item)

        rd  = random.randint(0, len(flat)-1)
        val = flat[rd]

        mass_ratio = q[val]

        return mass_ratio

    return prob(period, mass) * mass


def get_eccentricity(mass, period):

    ecc_x = np.arange(0.005, 1.00, 0.005)

    def ecc_max(p):

        e_max = 1 - (10 ** p / 2.0) ** (-2 / 3.0)
        return e_max

    def prob_ecc(e, m, p):

        if m <= 3:
            eta = 0.6 - 0.7 / (p - 0.5)

        elif m > 7:
            eta = 0.9 - 0.2 / (p - 0.5)

        elif 3 < m <= 7:
            eta = -0.225 - 0.075 * m + (-0.125 + 0.475) / (p - 0.5)

        prob = int(100 * e ** eta)
        return prob

    # Make sure all arrays have a reasonable size
    # We want to avoid memory overflow errors

    if period > 0.6:
        prob_arr = []
        for d in range(len(ecc_x)):
            prob_arr.append(prob_ecc(ecc_x[d], mass, period))
        digits = int(np.log10(np.max(prob_arr))) + 1
        if digits <= 2:
            all_prob = prob_arr
        else:
            all_prob = []
            for c in prob_arr:
                all_prob.append(c / (10 ** (digits - 2)))

        all_prob = np.array(all_prob)
        prob_arr = []
        for dd in range(len(ecc_x)):
            arr = np.array([ecc_x[dd]] * all_prob)
            prob_arr.extend(arr)
        rd = random.randint(0, len(prob_arr) - 1)
        ecc = prob_arr[rd]
        if ecc <= ecc_max(period):
            eccentricity = ecc
        else:
            eccentricity = ecc_max(period)

    elif period <= 0.6:
        dd = np.random.randint(0, 100 * ecc_max(period))
        eccentricity = dd / 100

    return eccentricity


##################################
# Fifth function to call, orbits #
##################################
def orbits(mass_array):

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


    def generate_binaries_with_orientation(mass_array = mass_array):
        """
        Generates binaries from arrays of masses and orbital elements
        The input arrays have no units
        The function returns sets of two particles, containing the primary and companion
        """

        multiplicity  = get_multiplicity(mass_array)[0]
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
    #np.savetxt('masses.txt', binaries[0])
    #np.savetxt('system_masses.txt', binaries[1])
    #np.savetxt('positions.txt', binaries[2])
    #np.savetxt('velocities.txt', binaries[3])
    
    return binaries


if __name__ == '__main__':
    pass
