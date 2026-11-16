from amuse.units import units
from amuse.io import read_set_from_file
from amuse.ext.orbital_elements import get_orbital_elements_from_binary
import numpy as np
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

def find_binary(filename):
    binary = read_set_from_file(filename)
    _, _, a, e, _, _, _, _ = get_orbital_elements_from_binary(binary)
    return a, e

def find_time(filename):
    ds = yt.load(filename)
    t = ds.current_time.in_units("Myr").value | units.Myr
    return t

if __name__=="__main__":

    import matplotlib.pyplot as plt
    # plot analytic solution
    t = np.linspace(0, 0.1, 101) | units.Myr 
    a = np.ones(len(t)) | units.au
    e = np.zeros(len(t))

    fig, ax1 = plt.subplots()
    ax1.plot(t.value_in(units.Myr), a.value_in(units.au), color='k', alpha=0.5, label='Kepler orbit')

    amuse_files, flash_files = find_snapshots('data', 'nbody')
    for i in range(len(amuse_files)):
        amuse_file = amuse_files[i]
        flash_file = flash_files[i]
        if amuse_file.endswith('stars0000.amuse'):
            continue #skip the empty particle set at the start
        elif amuse_file.endswith('stars0001.amuse'): 
            a_nbody, e_nbody = find_binary(amuse_file)
            t_nbody = find_time(flash_file)
            ax1.scatter(t_nbody.value_in(units.Myr), a_nbody.value_in(units.au), c='k', marker='o', label='N-body')
        else: 
            a_nbody, e_nbody = find_binary(amuse_file)
            t_nbody = find_time(flash_file)
            ax1.scatter(t_nbody.value_in(units.Myr), a_nbody.value_in(units.au), c='k', marker='o')
    
    ax1.set_ylabel('Semi-major axis (au)')
    ax1.set_xlabel('Time (Myr)')
    #ax1.set_xlim(0, 0.1)
    #ax1.set_ylim(-1, 2)

    fig.tight_layout()
    plt.legend(frameon=False)
    plt.savefig('r_bin_test.png')
