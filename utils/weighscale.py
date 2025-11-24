#!/usr/bin/env python
"""
Inspect cloud properties
"""
from __future__ import division, print_function

from datetime import datetime
import numpy as np
from os import path
import sys

np.set_printoptions(precision=3)

# hardcode to save time.  yt import is slow.
G = 6.67384e-08  # Grav constant in cgs units; yt value. This is 2010 NIST recommendation
Msun = 1.98841586e+33  # solar mass in g; yt value
cmpc = 3.0856775809623245e+18  # 1 parsec in cm; yt value
mH = 1.6737352238051868e-24  # hydrogen mass in g; yt value
year = 3.155e7  # year in seconds

# --------------------
# argparse / file read
# --------------------

fname = sys.argv[1]
if not fname:
    fname = 'cube128'

#if path.exists(fname+".dat"):
try:
    meta = np.loadtxt(fname+".dat")
    Msph, Rsph, box = meta[0,:3]
    print("from metadata: Msph = {:e} Msun, Rsph = {:g} pc, box = {:g} pc".format(
            Msph, Rsph, box))
    r0_vals = np.array([Rsph, 5, 10, 50, 500]) * cmpc  # values pc, store as cm
    box_sizes = np.array([box, 7, 12.5, 55, 600]) * cmpc  # values pc, store as cm
except Exception as e:
    print("Couldn't get target mass, size from .dat file, exception:", e)
    r0_vals = np.array([5, 10, 50, 500]) * cmpc  # values pc, store as cm
    box_sizes = np.array([7, 12.5, 55, 600]) * cmpc  # values pc, store as cm

print("loading", fname)
started = datetime.now()
cube = np.loadtxt(fname)
rho = cube[:,3]
print("done loading! elapsed:", datetime.now()-started)

# ---------------------
# report useful numbers
# ---------------------

nblock = 128
assert nblock**3 == cube.shape[0]

total = np.sum(rho) * (box_sizes*2/nblock)**3  # g

# assumes uniform density ambient medium and no weird behavior like cloud being
# bigger than cube
amb = np.amin(rho) * (box_sizes*2)**3 \
      * (1 - 4*np.pi/3*r0_vals**3/(box_sizes*2)**3)

print("Rsph (pc)", r0_vals/cmpc)
print("box halfwidth (pc)", box_sizes/cmpc)
print("Total mass (Msun):".format(fname), total / Msun)
print("Ambient mass (Msun):", amb / Msun)
print("Total - ambient mass (Msun):", (total - amb) / Msun)

print("t_freefall (Myr):", (3*np.pi/32/G * r0_vals**3/(total-amb))**0.5 / (1e6*year))
