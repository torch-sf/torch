a_dust_mass_weighted = 3.5355339059327366e-02 #in micron

def tau_sp(rhog, Tg):
    # rhog in g cm-3
    # Tg in K
    # returns sputtering timescale in Gyr
    return 0.17*(a_dust_mass_weighted/0.1) * (1e-27/rhog) * ((2e6/Tg)**2.5 + 1)



def sublimation_radius(L_star, T_sub):
    # L_star in L_sun
    # T_sub in K
    # returns sublimation radius in AU
    return 1.1 * (L_star/1e3)**0.5 * (1500/T_sub)**2