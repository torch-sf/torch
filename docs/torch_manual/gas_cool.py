#!/usr/bin/env python
"""
Plot gas cooling functions and stuff
2019 June 20, A Tran
"""

from __future__ import division, print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

kB = 1.38e-16  # boltzmann constant, CGS
mP = 1.67e-24  # proton mass, g


def lin_interp(x1, y1, x2, y2, x):
    """Linearly interpolate between points (x1,y1) and (x2,y2)"""
    assert x >= x1 and x <= x2
    ms = (y2-y1)/(x2-x1)
    return ms*(x-x1) + y1


def find(x, x0):
    """
    Given a monotonically increasing table x and a test value
    x0, return the index i of the largest table value less than
    or equal to x0 (or 0 if x0 < x(1)).  Use binary search.
    """
    return np.searchsorted(x, x0)


def radloss(T):
    """
    Joung/M-MML 2006 cooling curve
    which uses Dalgarno/McCray 1972 for T < 2e4, xHp = 1e-2
    and Sutherland/Dopita 1993 for T > 2e4, xHp = equilibrium ionization

    Units are erg cm^3 s^-1
    """
    # Possible to do linear if/elif chain,
    # but this binary-search-like structure is probably faster to evaluate.
    if T >= 1.70e4:
        if T >= 5.62e5:
            if T >= 2.75e6:
                if T >= 3.16e7:
                    return 3.090e-27*T**0.5  # >= 3.16e7
                else:
                    return 5.188e-21/T**0.33  # 2.75e6 to 3.16e7

            else:
                if T >= 1.78e6:
                    return 3.890e-4/T**2.95  # 1.78e6 to 2.75e6
                else:
                    return 1.3e-22*(T/1.5e5)**0.01  # 5.62e5 to 1.78e6
        else:
            if T >= 7.94e4:
                if T >= 2.51e5:
                    return 3.98e-11/T**2  # 2.51e5 to 5.62e5
                else:
                    return 6.31e-22*(T/1.e6)**0.01  # 7.94e4 to 2.51e5

            else:
                if T >= 3.98e4:
                    return 1.e-31*T**2  # 3.98e4 to 7.94e4
                else:
                    return 1.479e-21/(T**0.216)  # 1.70e4 to 3.98e4

    else:
        if T >= 1.e3:
            if T >= 6.31e3:
                if T >= 1.0e4:
                    return 7.63e-81*(T**13.8)  # 1.0e4 to 1.70e4, cliff
                else:
                    return 3.13e-27*(T**0.396)  # 6.31e3 to 1.0e4

            else:
                if T >= 3.16e3:
                    return 2.64e-26*(T**0.152)  # 3.16e3 to 6.31e3
                else:
                    return 5.28e-27*(T**0.352)  # 1e3 to 3.16e3
        else:
            if T >= 3.98e1:
                if T >= 2.00e2:
                    return 3.06e-27*(T**0.431)  # 200 to 1e3
                else:
                    return 1.52e-28*(T**0.997)  # 39.8 to 200

            else:
                if T >= 2.51e1:
                    return 2.39e-29*(T**1.50)  # 25.1 to 39.8
                else:
                    return 1.095e-32*(T**3.885)  # <= 25.1


def dustloss(T, Tdust):
    """Compute collisional cooling due to gas colliding with dust,
    taken from Torch code,
    units are erg cm^3 s^-1
    """
    if T <= Tdust:
        # Dust can heat gas, but we don't have correct physics/chemistry at low
        # temperature.  Therefore, just turn off cooling.
        return 0.
    return 3.8e-33 * (1.0 - 0.8*np.exp(-75.0/T)) * T**0.5 * (T - Tdust)


def get_cooling_data():
    """Ported from get_cooling_data.F90"""
    cool_dat = np.zeros((DENS_PTS, TEMP_PTS, 3))

    dat = np.loadtxt('cool.dat', skiprows=1)
    assert dat.shape[0] == TEMP_PTS * DENS_PTS

    for j in range(TEMP_PTS):
        for i in range(DENS_PTS):
            temp = dat[i+j*DENS_PTS, 0]  # flatten index
            dens = dat[i+j*DENS_PTS, 1]
            cool_pwr = dat[i+j*DENS_PTS, 5]

            cool_dat[i, j, 0] = 10**temp
            cool_dat[i, j, 1] = 10**dens
            cool_dat[i, j, 2] = 10**cool_pwr

    return cool_dat


def molecularloss(T, ndens):
    """Molecular cooling method, re-implemented
    units are erg cm^3 s^-1
    """
    ixdens = find(COOL_DAT[:,0,1], ndens)  # uses the first "block" to search
    ixtemp = find(COOL_DAT[0,:,0], T)
    ixdens = min(ixdens, DENS_PTS - 1)  # adjust bound... we are doing nearest neighbor instead of interpolating
    ixtemp = min(ixtemp, TEMP_PTS - 1)
    # heatCool.F90 returns COOL_DAT / mu is erg s^-1 g^-1
    # here we return COOL_DAT / ndens is erg cm^3 s^-1
    # to match other cooling functions in this python code
    return COOL_DAT[ixdens, ixtemp, 2] / ndens


DENS_PTS = 190  # Cool_data.F90
TEMP_PTS = 320  # Cool_data.F90
COOL_DAT = get_cooling_data()  # global . . .
print(COOL_DAT.shape)

radloss = np.vectorize(radloss)
dustloss = np.vectorize(dustloss)
molecularloss = np.vectorize(molecularloss)

# -----------------------------------------------
# Plot the "DGionfix" from Josh.

##def radloss_dgionfix(T, xHp):
##    """Interpolate for varying ionization"""
##    emin_out = radloss(T)
##    if T < 8e3:
##        if xHp <= 1e-1:  # Roughly 2-10x diff b/t xHp=1e-1 and xHp=1e-2 curve.
##            emin_out = emin_out / lin_interp(1e-2, 5.0, 1e-1, 1.0, max(1e-2,xHp))
##            if xHp <= 1e-2:  # Roughly 2-5x diff b/t xHp=1e-2 and xHp=1e-4 curve.
##                emin_out = emin_out / lin_interp(1e-4, 2.0, 1e-2, 1.0, max(1e-4,xHp))
##    return emin_out
##
##radloss_dgionfix = np.vectorize(radloss_dgionfix)
##
##T_vec = np.logspace(1, 8, 10000)
##atomic_cool = radloss(T_vec)
##plt.plot(T_vec, atomic_cool, 'k-', label='Joung/M-MML', linewidth=4, alpha=0.5)
##
##for xHp in [1e-1, 1e-2, 1e-3, 1e-4]:
##    T_vec = np.logspace(1, 8, 10000)
##    atomic_cool = radloss_dgionfix(T_vec, xHp)
##    plt.plot(T_vec, atomic_cool, '-', label='xHp={:e}'.format(xHp))
##
##plt.xscale('log')
##plt.yscale('log')
##plt.xlabel('T [K]')
##plt.ylabel(r'$\Lambda$ [erg cm${}^3$ s${}^{-1}$]')
##plt.xlim(1e1, 1e8)
##plt.ylim(1e-27, 1e-21)
##plt.legend()
##plt.show()

# ----------------------------------------------------------
# Plot proposed correction to Josh's implementation:
# (1) fix bug in linear interpolation;
# (2) adjust temperature/xHp ranges to ensure monotonic, piecewise-continuous
# behavior near T ~ 1e4 K

def radloss_dgionfix_bugfix(T, xHp):
    """Interpolate for varying ionization"""
    emin_out = radloss(T)
    if xHp > 1e-2 and T < 1.2e4:
        # Roughly 2-10x diff b/t xHp=1e-1 and xHp=1e-2 curve.
        factor = lin_interp(1e-2, 1.0, 1e-1, 5.0, min(1e-1,xHp))
        emin_out = emin_out * factor
        # piecewise-continuous, monotonic transition to CIE curve
        if T > 1e4:
            emin_out = emin_out / lin_interp(1e4, 1.0, 1.2e4, factor**0.5, T)**2
    elif xHp < 1e-2 and T < 1e4:
        # Roughly 2-5x diff b/t xHp=1e-2 and xHp=1e-4 curve.
        factor = lin_interp(1e-4, 0.5, 1e-2, 1.0, max(1e-4,xHp))
        emin_out = emin_out * factor
        # piecewise-continuous, monotonic transition to CIE curve
        if T > 9e3:
            emin_out = emin_out / lin_interp(9e3, 1.0, 1e4, factor, T)
    return emin_out

radloss_dgionfix_bugfix = np.vectorize(radloss_dgionfix_bugfix)

T_vec = np.logspace(1, 8, 10000)
atomic_cool = radloss(T_vec)
plt.plot(T_vec, atomic_cool, '-k', label='Joung/M-MML', linewidth=5, alpha=0.3)

for xHp in [1e-1, 1e-2, 1e-3, 1e-4]:
    T_vec = np.logspace(1, 8, 10000)
    atomic_cool = radloss_dgionfix_bugfix(T_vec, xHp)
    plt.plot(T_vec, atomic_cool, '-', label='xHp={:.0e}'.format(xHp))

# also show dust cooling...
T_vec = np.logspace(1, 8, 10000)
dust_cool = dustloss(T_vec, 10)  # Tdust=10K assumed
plt.plot(T_vec, dust_cool, ':k', label='dust')

# also show molec cooling...
for ndens in [1e1, 1e3, 1e6]:
    T_vec = np.logspace(1, 3.4, 10000)  # special bound for cool.dat
    mol_cool = molecularloss(T_vec, ndens)
    plt.plot(T_vec, mol_cool, '--', label='mol,n={:.0e}'.format(ndens))

plt.xscale('log')
plt.yscale('log')
plt.xlabel('T [K]')
plt.ylabel(r'$\Lambda$ [erg cm${}^3$ s${}^{-1}$]')
plt.xlim(1e1, 1e8)
#plt.ylim(1e-27, 1e-21)  # dalgarno/mccray limits
plt.ylim(1e-33, 1e-21)
plt.legend()

#plt.show()

plt.tight_layout()
plt.savefig('fig/gas_cool.pdf', bbox_inches='tight')
plt.clf()
plt.close()
