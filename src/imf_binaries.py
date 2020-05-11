"""
Generate samples from initial mass function (IMF).

Joshua Wall, Drexel University
"""

from __future__ import division

import numpy as np
from scipy.integrate import quad


def sample_stars(sample_imf_mass, num_bins=10, min_samp_mass=1.0,
                              max_samp_mass=150.0, sum_small=False):

    [n_stars, bins, lam, norm] = sample_stars_poisson(sample_imf_mass, min_samp_mass, max_samp_mass, num_bins)

    # Now use that to sample the IMF.
    masses     = np.zeros(n_stars.sum())
    positions  = np.zeros((n_stars.sum(), 3))
    velocities = np.zeros((n_stars.sum(), 3)) #Make sure positions and velocities have 3 components
    k = 0
    for i, n in enumerate(n_stars):
        #print "Pulling ", n, "stars from ranges ", bins[i], "to ", bins[i+1]
        for j in range(n):

            while (masses[k] == 0):

                m = np.random.uniform(low=bins[i], high=bins[i+1])
                r = np.random.uniform()
                p = mkroupa(m, norm)

                if (p/r > 1.0):
                    masses[k] = m

            k+=1

    # Sum all stars < 1 MSun into stars > 1 MSun.
    if (sum_small):
        masses = collect_small_stars_mass(masses)

    np.random.shuffle(masses)

    system_masses = masses

    return masses, system_masses, positions, velocities


def sample_stars_poisson(sink_mass, M_min, M_max, num_bins):
    """
    Return a poisson random sampling from the Kroupa IMF of sink total
    mass from M_min to M_max separated into num_bins in logspace.

    Returns:
        n_stars: Number of stars in each logarithmic bin
        binsL:   The bin edges, including the right most bin edge
        lam:     The average number of stars in each bin that the
                 Poisson sample is centered around.
        norm:    Norm to be used to sample the Kroupa IMF
                 using the n_stars array.
    """

    norm_inv = quad(kroupa,M_min,M_max,args=(1))[0]

    norm = 1/norm_inv

    binsL = np.logspace(np.log10(M_min),np.log10(M_max),num_bins+1)
    mass_per_bin = []
    frac_per_bin = []

    avg_mass = quad(mkroupa, M_min, M_max, args=(norm))[0]

    for i in range(num_bins):
        mass_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0]
                             / quad(kroupa, binsL[i], binsL[i+1], args=(norm))[0] )
        frac_per_bin.append( quad(mkroupa, binsL[i], binsL[i+1], args=(norm))[0]
                             / avg_mass )

    mass_per_bin = np.array(mass_per_bin)
    frac_per_bin = np.array(frac_per_bin)

    lam = sink_mass*frac_per_bin/mass_per_bin

    n_stars = np.random.poisson(lam=lam)

    return n_stars, binsL, lam, norm



def collect_small_stars_mass(masses):

    # Here we move all the stars smaller than 1 MSun into particles
    # that are at least 1 MSun. To do this we do a bit of fancy
    # footwork with the arrays.

    small_masses = masses[np.where(masses < 1.0)] # Smaller than 1.0 MSun.
    masses = masses[np.where(masses >= 1.0)]  # Everyone else.

    b = 0
    # If there are any left smaller than 1.0 MSun, sum with others
    # that are smaller than 1.0 MSun until there are none left.

    if (len(small_masses) > 1):
        while(small_masses[-1] < 1.0 and len(small_masses[b:])>1):

            small_masses[b] = small_masses[b]+small_masses[-1]
            small_masses = np.delete(small_masses, -1)
            if (len(small_masses[b:]) > 1):
                if(small_masses[b] >= 1.0):
                    b += 1

        # If the last one is smaller than 1.0 MSun, lump that bit into
        # the last star.
        if (small_masses[-1] < 1.0):
            small_masses[-2] = small_masses[-2] + small_masses[-1]
            small_masses = np.delete(small_masses, -1)

    masses = np.append(masses, small_masses)

    return masses


def m_max_star(m_max_clust):
    # The max stellar mass for sampling the "normal" IMF
    # calculated from Weidner et. al. 2013 eqn 1,
    # based on the integrated galatic IMF of Weidner and Kroupa 2004.

    # m_max_clust is the maximum cluster mass
    # and should figure in losses due to jets and other
    # feedback. Generally, I just assume a SFE of 0.5.

    a0 = -0.66
    a1 =  1.08
    a2 = -0.15
    a3 = 0.0084

    Lmclust = np.log10(m_max_clust)

    if (m_max_clust <= 2.5E5):
        m_max = a0 + a1*Lmclust + a2*Lmclust**2. + a3*Lmclust**3.0
    else:
        m_max = np.log10(150.0)

    return 10**m_max


def kroupa(m,a):

    if (0.001 <= m < 0.08):
        k = a*m**(-0.3)
    elif (0.08 <= m < 0.5):
        k = a*(0.08)*m**(-1.3)
    elif (0.5 <= m):
        k = a*(0.08*0.5)*m**(-2.3)
    else:
        print "Invalid mass range!"
        k=0
    return k


def mkroupa(m,a):

    if (0.001 <= m < 0.08):
        k = m*a*m**(-0.3)
    elif (0.08 <= m < 0.5):
        k = m*a*(0.08)*m**(-1.3)
    elif (0.5 <= m):
        k = m*a*(0.08*0.5)*m**(-2.3)
    else:
        print "Invalid mass range!"
        k=0
    return k


if __name__ == '__main__':
    pass
