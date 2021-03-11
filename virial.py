"""
Multiplies the velocities in x, y and z to change the virial parameter of a checkpoint file
Multiples v^2 by virial_new/virial_old
Overwrites the file - CAREFUL, MAKE SURE TO COPY FILE BEFORE USE
"""

import h5py
import numpy as np
filename = "M4V10-C/turbsph_hdf5_chk_0000"

# Old and new virial parameters
v_i = 0.2
v_f = 1.0
vir_r = v_f/v_i
vel_r = vir_r ** (1./2)
print("Multiplying every velocity component by %.2f" % vel_r)

with h5py.File(filename, "r+") as f:
    # List all groups
    #print("Keys: %s" % f.keys())

    for key in f:
        if key in ['velx', 'vely', 'velz']:
            vel = f[key][...]
            print(np.max(vel))
            f[key][...] = vel * vel_r
            print(np.max(f[key][...]))
        else:
            pass
    
    
