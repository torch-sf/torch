"""
physics/constants
==========

Holds constants in various unit bases.
"""

class cgs:
    """
    Constants in centimeter, grams, and seconds.
    """

    # Mass [grams]
    Msun = Msol = 1.989e33  # Mass of Sun
    mH = 1.6736e-24         # Mass of proton

    # Distance [cenitmeter]
    km = 1e5
    pc = km*3.086e13
    A = 1e-8                # Angstrom
    au = 1.496e+13          # Astronomical unit

    # Time [second]
    hour = 3600
    day = 24*hour
    year = 365.25*day
    Myr = 1e6*year

    # Other
    G = 6.674e-8            # [cm**3 g**-1 s**-2] - Gravitational constant
    kB = 1.380e-16          # [erg/K] -  Boltzmann constant
    h = 6.626176e-27        # [erg s] - Plack constant
    eV = 1.602e-12          # [erg] - Electron volt
    c = 29979245800         # [cm/s] - Speed of light
    Lsol = Lsun = 3.839e33  # [erg/s] - Solar luminosity
    Zsol = Zsun = 0.02      # Solar metallicity

