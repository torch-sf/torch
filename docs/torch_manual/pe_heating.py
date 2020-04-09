#!/usr/bin/env python
"""
Decipher why Josh uses 1.0 - eff in equation for PE heating of dust
Spring 2019
Modified June 2019
"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt


def epsilon_wd01(Gfactor, T):
    """
    Weingartner Draine 2001 equation (44) and Table 2, R_V=3.1 b_c=6e-5
    Rescaled by factor 100x.
    They state that this provides the highest heating; the dust has the most
    small grains.
    """
    c0 = 7.64 / 100  # rescaling to match Bakes and Tielens (1994)
    c1 = 4.52 / 100
    c2 = 0.04371
    c3 = 0.00557
    c4 = 0.132
    c5 = 0.452
    c6 = 0.675
    return (c0 + c1*T**c4) /(1.0 + (c2*Gfactor**c5)*(1.0 + c3*Gfactor**c6))

def epsilon_bt94(Gfactor, T):
    """
    Bakes/Tielens 1994, equation (43)
    """
    eps0 = 4.87e-2 / (1 + 4e-3 * Gfactor**0.73)
    eps1 = 3.65e-2 * (T/1e4)**0.7 / (1 + 2e-4 * Gfactor)
    return eps0 + eps1

# ----------------------------------------------------------------------
# Check Gfactor expressions for varying T
# Then also plot heating rate...
# ... per H atom per G
# ... per H atom
#
# To compare with Fig 12 of Bakes/Tielens 1994, Fig 16 of Weingartner/Draine
# 2001.

G = 1.69
ne = np.logspace(-4,2,1000)
# http://colorbrewer2.org/#type=diverging&scheme=Spectral&n=4
T_range = [10, 100, 1000, 10000]
color_range = ['#d7191c', '#fdae61', '#abdda4', '#2b83ba']
assert len(T_range) == len(color_range)

fig, axes = plt.subplots(1, 3, figsize=(15,5), sharex=True)

for T, color in zip(T_range, color_range):
    Gfactor = G * np.sqrt(T) / ne

    plt.sca(axes[0])
    plt.plot(Gfactor, epsilon_bt94(Gfactor, T), '--', color=color,
             label='BT94, T={}'.format(T))
    plt.plot(Gfactor, epsilon_wd01(Gfactor, T), '-', color=color,
             label='WD01, T={}'.format(T))
    plt.xscale('log')
    plt.xlabel('G factor')
    plt.yscale('log')
    plt.ylabel('epsilon')
    plt.legend()

    plt.sca(axes[1])
    plt.plot(Gfactor, 1e-24 * epsilon_bt94(Gfactor, T), '--', color=color,
             label='BT94, T={}'.format(T))
    plt.plot(Gfactor, 1e-24 * epsilon_wd01(Gfactor, T), '-', color=color,
             label='WD01, T={}'.format(T))
    plt.xscale('log')
    plt.xlabel('G factor')
    plt.yscale('log')
    plt.ylabel(r'$\Gamma_\mathrm{pe} / G / n_\mathrm{H}$ [erg/s/Hatom/G]')

    plt.sca(axes[2])
    plt.plot(Gfactor, 1e-24 * G * epsilon_bt94(Gfactor, T), '--', color=color,
             label='BT94, T={}'.format(T))
    plt.plot(Gfactor, 1e-24 * G * epsilon_wd01(Gfactor, T), '-', color=color,
             label='WD01, T={}'.format(T))
    plt.xscale('log')
    plt.xlabel('G factor')
    plt.yscale('log')
    plt.ylabel(r'$\Gamma_\mathrm{pe} / n_\mathrm{H}$ [erg/s/Hatom]')

plt.tight_layout()
plt.savefig('fig/pe_heating.pdf', bbox_inches='tight')
#plt.show()
plt.clf()
plt.close()
