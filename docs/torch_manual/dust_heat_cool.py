#!/usr/bin/env python
"""
Compute dust heat/cool, find equilibrium dust temperature
and compute dust cooling time
"""
# Try Goldsmith 2001 values



from __future__ import division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

kB = 1.38e-16  # boltzmann constant, CGS
#G = 6.6725985e-8  # Newton's grav constant, CGS
#pi = 3.1415926535897932384
mP = 1.67e-24  # proton mass, g
#s2yr = 1/3.155e7  # convert seconds to years


for Tg in [10, 1e4]:

    n_vec = np.logspace(0, 6, 400)
    T_vec = np.logspace(1, 2, 200)

    n, T = np.meshgrid(n_vec, T_vec)

    fig, axes = plt.subplots(2, 2, figsize=(14,10))

    # -----------------------------------------------------------------------
    # Goldsmith 2001, eqns (7), (13), (15)
    # let n be H nucleon density, so n(H2) = n/2 for pure molec gas
    Gamma_ext = 3.9e-24 * (n/2)  # assume no attenuation; i.e., \chi = 1
    Gamma_gd = 2e-33 * (n/2)**2 * (Tg - T) * (Tg/10)**0.5
    Lambda = 6.8e-33 * T**6 * (n/2)
    du_dt = Gamma_ext + Gamma_gd - Lambda
    # dust-to-gas mass ratio is 0.01
    # assume dust grain mass ~ 1e-14 g from https://ned.ipac.caltech.edu/level5/Sept05/Li/Li2.html#note24
    # then, n_dust = 0.01 * m_H n_H / m_D
    t_unity = (3./2)*(0.01*mP*n/1e-14)*kB*T / np.abs(du_dt)

    plt.sca(axes[0,0])
    plt.title('goldsmith 2001, Tg={}, de/dt (erg cm**-3 s**-1)'.format(Tg))
    plt.imshow(du_dt, origin='lower',
               extent=(np.log10(n_vec[0]), np.log10(n_vec[-1]),
                       np.log10(T_vec[0]), np.log10(T_vec[-1])),
               aspect='auto',
               norm=mpl.colors.SymLogNorm(1e-23, vmin=-1e-15,vmax=1e-15),
               cmap='RdBu_r'
    )
    plt.colorbar()

    cs = plt.contour(np.log10(n_vec), np.log10(T_vec), du_dt,
            colors='k',
            levels=[-1e-15, -1e-16, -1e-17, -1e-18, 0, 1e-18, 1e-17, 1e-16, 1e-15]
    )

    plt.xlabel("log10(n(H))")
    plt.ylabel("log10(T)")

    plt.sca(axes[0,1])
    plt.title('goldsmith 2001, dust cooling time, s')
    plt.imshow(t_unity, origin='lower',
               extent=(np.log10(n_vec[0]), np.log10(n_vec[-1]),
                       np.log10(T_vec[0]), np.log10(T_vec[-1])),
               aspect='auto',
               norm=mpl.colors.LogNorm(),#vmin=1e0, vmax=1e4),
    )
    plt.colorbar()
    plt.xlabel("log10(n(H))")
    plt.ylabel("log10(T)")
    # -----------------------------------------------------------------------
    # Glover and Clark 2012
    Gamma_ext = 5.6e-24 * n
    Gamma_gd = 3.8e-33 * Tg**0.5 * (1 - 0.8*np.exp(-75./Tg)) * (Tg - T) * n**2
    Lambda = 4.68e-31 * T**6 * n
    du_dt = Gamma_ext + Gamma_gd - Lambda
    # dust-to-gas mass ratio is 0.01
    # assume dust grain mass ~ 1e-14 g from https://ned.ipac.caltech.edu/level5/Sept05/Li/Li2.html#note24
    # then, n_dust = 0.01 * m_H n_H / m_D
    t_unity = (3./2)*(0.01*mP*n/1e-14)*kB*T / np.abs(du_dt)

    plt.sca(axes[1,0])
    plt.title('glover clark 2012, Tg={}'.format(Tg))
    plt.imshow(du_dt, origin='lower',
               extent=(np.log10(n_vec[0]), np.log10(n_vec[-1]),
                       np.log10(T_vec[0]), np.log10(T_vec[-1])),
               aspect='auto',
               norm=mpl.colors.SymLogNorm(1e-23, vmin=-1e-15,vmax=1e-15),
               cmap='RdBu_r'
    )
    plt.colorbar()
    cs = plt.contour(np.log10(n_vec), np.log10(T_vec), du_dt,
            colors='k',
            levels=[-1e-15, -1e-16, -1e-17, -1e-18, 0, 1e-18, 1e-17, 1e-16, 1e-15]
    )
    plt.xlabel("log10(n(H))")
    plt.ylabel("log10(T)")

    plt.sca(axes[1,1])
    plt.title('glover clark 2012, dust cooling time, s')
    plt.imshow(t_unity, origin='lower',
               extent=(np.log10(n_vec[0]), np.log10(n_vec[-1]),
                       np.log10(T_vec[0]), np.log10(T_vec[-1])),
               aspect='auto',
               norm=mpl.colors.LogNorm(),#vmin=1e0, vmax=1e4),
    )
    plt.colorbar()
    plt.xlabel("log10(n(H))")
    plt.ylabel("log10(T)")
    # -----------------------------------------------------------------------

    plt.tight_layout()

    #plt.show()

    fname = 'fig/dust_heat_cool_Tgas{:g}.pdf'.format(Tg)
    plt.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()

    print("Saved", fname)
