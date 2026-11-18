from amuse.units import units
from amuse.datamodel import Particles, Particle
from amuse.community.seba import Seba
import numpy as np
from matplotlib import pyplot as plt
from amuse.io import read_set_from_file

star_10Msun = read_set_from_file('data/seba_amuse_10Msun.amuse')
star_20Msun = read_set_from_file('data/seba_amuse_20Msun.amuse')
star_40Msun = read_set_from_file('data/seba_amuse_40Msun.amuse')
star_80Msun = read_set_from_file('data/seba_amuse_80Msun.amuse')

plt.figure()
plt.plot(star_10Msun.age.value_in(units.Myr), star_10Msun.mass.value_in(units.MSun), label=r'10 M$_\odot$')
plt.plot(star_20Msun.age.value_in(units.Myr), star_20Msun.mass.value_in(units.MSun), label=r'20 M$_\odot$')
plt.plot(star_40Msun.age.value_in(units.Myr), star_40Msun.mass.value_in(units.MSun), label=r'40 M$_\odot$')
plt.plot(star_80Msun.age.value_in(units.Myr), star_80Msun.mass.value_in(units.MSun), label=r'80 M$_\odot$')
plt.xlabel('Age (Myr)')
plt.ylabel('Mass (MSun)')
plt.legend()
plt.savefig('./seba_amuse_mass.png')

# Masks to only plot the non-remnant values for the HR diagram
_mask_10_Msun = np.where(star_10Msun.stellar_type < 13 | units.stellar_type)[0]
_mask_20_Msun = np.where(star_20Msun.stellar_type < 13 | units.stellar_type)[0]
_mask_40_Msun = np.where(star_40Msun.stellar_type < 13 | units.stellar_type)[0]
_mask_80_Msun = np.where(star_80Msun.stellar_type < 13 | units.stellar_type)[0]

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e6, 1e3)
plt.ylim(1e3, 1e7)
plt.plot(star_10Msun.temperature[_mask_10_Msun].value_in(units.K), star_10Msun.luminosity[_mask_10_Msun].value_in(units.LSun), label=r'10 M$_\odot$')
plt.plot(star_20Msun.temperature[_mask_20_Msun].value_in(units.K), star_20Msun.luminosity[_mask_20_Msun].value_in(units.LSun), label=r'20 M$_\odot$')
plt.plot(star_40Msun.temperature[_mask_40_Msun].value_in(units.K), star_40Msun.luminosity[_mask_40_Msun].value_in(units.LSun), label=r'40 M$_\odot$')
plt.plot(star_80Msun.temperature[_mask_80_Msun].value_in(units.K), star_80Msun.luminosity[_mask_80_Msun].value_in(units.LSun), label=r'80 M$_\odot$')
plt.xlabel('Temperature (K)')
plt.ylabel('Luminosity (LSun)')
plt.legend()
plt.savefig('./seba_amuse_HR.png')