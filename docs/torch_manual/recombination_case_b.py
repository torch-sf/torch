#!/usr/bin/env python
"""
Quickly plot expression from Draine chapter 14
Aaron Tran - 2019 April 23
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def k_recb(T):
    return 2.54e-13 * (T/10000.0)**(-0.8163-0.0208* np.log(T/10000.0))

def k_recb_crude(T):
    return 2.54e-13 * (T/10000.0)**(-0.8163)

def k_recb_christian(T):
    """brech0 in Christian's calc_ionization"""
    return 2.59e-13 * (T/10000.0)**(-0.7)

t = np.logspace(np.log10(30), np.log10(3e4),100)

plt.plot(t, k_recb(t), label='Draine case B recombination')
plt.plot(t, k_recb_crude(t), label='Draine eqn, cruder')
plt.plot(t, k_recb_christian(t), label=r'C. Baczynski old $\alpha_\mathrm{B}$')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel(r'$\alpha_\mathrm{B}$ (cm${}^3$ s${}^{-1}$)')
plt.legend()

#plt.show()

plt.tight_layout()
plt.savefig('fig/recombination_case_b.pdf', bbox_inches='tight')
plt.clf()
plt.close()
