a_dust_mass_weighted = 3.5355339059327366e-02 #in micron

def tau_sp(rhog, Tg):
    # rhog: gas density in g cm-3
    # Tg: gas temperature in K
    # returns sputtering timescale in Gyr
    #
    # Reference: Tsai & Mathews 1995, ApJ, 448, 84, Equation 14
    #            Based on Tielens et al. 1994 and Draine & Salpeter 1979
    return 0.17*(a_dust_mass_weighted/0.1) * (1e-27/rhog) * ((2e6/Tg)**2.5 + 1)



def sublimation_radius(L_star, T_sub):
    # L_star: stellar luminosity in L_sun
    # T_sub: sublimation temperature in K
    # returns sublimation radius in AU
    #
    # Reference: Monnier & Millan-Gabet 2002, ApJ, 579, 694, Equation 1
    return 1.1 * (L_star/1e3)**0.5 * (1500/T_sub)**2