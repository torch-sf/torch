#!/usr/bin/env python
"""
Compare dust opacities, in attempt to track down difference in dust radiative
cooling rates.

* Goldsmith 2001 (quoting Goldsmith et al. 1997)

* Glover/Clark 2012 appendix (quoting Ossenkopf/Henning 1994)

Aaron Tran
2019 May 17
"""

from __future__ import division, print_function

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

c = 3e10  # cm s^-1, light speed
h = 6.626e-27  # erg Hz^-1, planck constant
mP = 1.67e-24  # g, proton mass


# goldsmith 2001 opacities
lamb = np.logspace(-6, -3, 1000) * 1e2  # 1 to 1000 micron, stated in cm
nu = c/lamb
nu0 = c / (790e-6 * 1e2)  # reference freq is 790 micron, stated here in cm
kappa = 3.3e-26 / (2*mP) * (nu/nu0)**2  # converted to cm^2 g^-1

plt.plot(lamb/1e2*1e6, kappa, '-k.', label='goldsmith')

del lamb, nu, nu0, kappa

# ossenkopf/henning 1994 opacity table 1, column 7
# MRN model non-coagulated, thick ice mantle grains
lamb = np.array([
    1, 2.15, 3.04, 5.01, 10, 20, 46.4,
    100, 226, 500, 1000,

]) * 1e-6 * 1e2  # micron to cm
nu = c/lamb
kappa = np.array([
    1.11e4, 4.42e3, 3.31e4, 1.27e3, 2.90e3, 2.46e3, 1.41e3,
    1.27e2, 1.75e1, 3.79e0, 1.01e0,
])  # cm^2 g^-1

plt.plot(lamb/1e2*1e6, kappa, '-o', label='ossenkopf/henning')

del lamb, nu, kappa

# display all together

plt.xlabel('wavelength (micron)')
plt.ylabel('kappa (cm2 g-1)')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')

plt.tight_layout()

#plt.show()

plt.savefig('fig/dust_opacity.pdf', bbox_inches='tight')
plt.clf()
plt.close()
