from amuse.units import units
from amuse.datamodel import Particles, Particle
from amuse.community.seba import Seba
import numpy as np
from matplotlib import pyplot as plt
from amuse.io import read_set_from_file
import yt

# Avoid unnecessary yt output
yt.set_log_level(50)

def find_snapshots(data_dir, sim_name):
    from glob import glob 
    import os
    pattern = os.path.join(data_dir, f"stars*.amuse")
    amuse_files = sorted(glob(pattern))
    pattern = os.path.join(data_dir, f"{sim_name}_hdf5_plt_cnt_*")
    flash_files = sorted(glob(pattern))
    return amuse_files, flash_files

def find_time(filename):
    ds = yt.load(filename)
    t = ds.current_time.in_units("Myr").value | units.Myr
    return t

if __name__=="__main__":

    star_20Msun = read_set_from_file('data/seba_amuse_20Msun.amuse')
    star_40Msun = read_set_from_file('data/seba_amuse_40Msun.amuse')
    star_80Msun = read_set_from_file('data/seba_amuse_80Msun.amuse')
    amuse_files, flash_files = find_snapshots('data', 'seba')

    plt.figure()
    plt.plot(star_20Msun.age.value_in(units.Myr), 
             star_20Msun.mass.value_in(units.MSun), 
             label=r'20 M$_\odot$')
    plt.plot(star_40Msun.age.value_in(units.Myr), 
             star_40Msun.mass.value_in(units.MSun), 
             label=r'40 M$_\odot$')
    plt.plot(star_80Msun.age.value_in(units.Myr), 
             star_80Msun.mass.value_in(units.MSun), 
             label=r'80 M$_\odot$')
    for i in range(len(amuse_files)):
        amuse_file = amuse_files[i]
        flash_file = flash_files[i]
        if amuse_file.endswith('stars0000.amuse'):
            continue #skip the empty particle set at the start
        else: 
            time = find_time(flash_file)
            stars = read_set_from_file(amuse_file)
            plt.scatter(time.value_in(units.Myr)*np.ones(len(stars)), 
                        stars.mass.value_in(units.Msun), 
                        color='k', s=10, alpha=0.5)
    plt.xlabel('Age (Myr)')
    plt.ylabel('Mass (MSun)')
    plt.legend()
    plt.savefig('./seba_amuse_mass.png')

    plt.figure()
    plt.plot(star_20Msun.age.value_in(units.Myr), 
             star_20Msun.stellar_type.value_in(units.stellar_type), 
             label=r'20 M$_\odot$')
    plt.plot(star_40Msun.age.value_in(units.Myr), 
             star_40Msun.stellar_type.value_in(units.stellar_type), 
             label=r'40 M$_\odot$')
    plt.plot(star_80Msun.age.value_in(units.Myr), 
             star_80Msun.stellar_type.value_in(units.stellar_type), 
             label=r'80 M$_\odot$')
    for i in range(len(amuse_files)):
            amuse_file = amuse_files[i]
            flash_file = flash_files[i]
            if amuse_file.endswith('stars0000.amuse'):
                continue #skip the empty particle set at the start
            else: 
                time = find_time(flash_file)
                stars = read_set_from_file(amuse_file)
                plt.scatter(time.value_in(units.Myr)*np.ones(len(stars)), 
                            stars.stellar_type.value_in(units.stellar_type), 
                            color='k', s=10, alpha=0.5)
    plt.xlabel('Age (Myr)')
    plt.ylabel('Stellar type')
    plt.legend()
    plt.savefig('./seba_amuse_type.png')

    # Masks to only plot the non-remnant values for the HR diagram
    _mask_20_Msun = np.where(star_20Msun.stellar_type < 13 | units.stellar_type)[0]
    _mask_40_Msun = np.where(star_40Msun.stellar_type < 13 | units.stellar_type)[0]
    _mask_80_Msun = np.where(star_80Msun.stellar_type < 13 | units.stellar_type)[0]

    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e6, 1e3)
    plt.ylim(1e3, 1e7)
    plt.plot(star_20Msun.temperature[_mask_20_Msun].value_in(units.K), 
             star_20Msun.luminosity[_mask_20_Msun].value_in(units.LSun), 
             label=r'20 M$_\odot$')
    plt.plot(star_40Msun.temperature[_mask_40_Msun].value_in(units.K), 
             star_40Msun.luminosity[_mask_40_Msun].value_in(units.LSun), 
             label=r'40 M$_\odot$')
    plt.plot(star_80Msun.temperature[_mask_80_Msun].value_in(units.K), 
             star_80Msun.luminosity[_mask_80_Msun].value_in(units.LSun), 
             label=r'80 M$_\odot$')
    for i in range(len(amuse_files)):
            amuse_file = amuse_files[i]
            flash_file = flash_files[i]
            if amuse_file.endswith('stars0000.amuse'):
                continue #skip the empty particle set at the start
            else: 
                stars = read_set_from_file(amuse_file)
                plt.scatter(stars.temperature.value_in(units.K), 
                            stars.luminosity.value_in(units.LSun),
                            color='k', s=10, alpha=0.5)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Luminosity (LSun)')
    plt.legend()
    plt.savefig('./seba_amuse_HR.png')