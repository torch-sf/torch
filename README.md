README
======

Quick-start
-----------

Set the environment variables:

    export AMUSE_DIR=/path/to/amuse
    export FLASH_DIR=/path/to/FLASH4.5
    export TORCH_DIR=/path/to/torch
    export PYTHONPATH=$PYTHONPATH:$TORCH_DIR/src

Then, run the script:

    ./install.sh

You will be prompted to confirm your environment variables, and continue with
installation.  The installer will then automatically copy interface code from
`$TORCH_DIR/src` to their correct locations in your AMUSE and FLASH repos.
It will also ensure that the FLASH4.5-a patch for HDF5 >1.10.x compatibility is
applied in your FLASH repo.


Contributing
------------

Torch is an open-source code!  Please see [CONTRIBUTING.md](CONTRIBUTING.md) to
learn how to get involved.


What's in this repository?
--------------------------

The top-level directory provides FLASH/AMUSE bridge code to perform
simulations:

    bridge_multiples_for_rad_testing.py
    bridge_multiples.py
    cool.dat
    cube128
    flash.par.radtest
    flash.par.turbsph_standard
    torch_user.py   # beta-release update of bridge_multiples.py

And, some code to create turbulent initial conditions (dens, pres, vel, etc.)
for your simulations:

    hAc_b_2.0E-17_e_0.021_FUV_1.69.dat
    turb-sphere.py
    turb-velbox.py
    weighscale.py

The path `src/` holds Python modules used by the top-level Torch code, which
perfoms a coupled FLASH and N-body simulation using the AMUSE framework.
Some stand-alone utility modules used by Torch are also included.

    src/
        torch_mainloop.py
        torch_param.py
        torch_se.py
        ...

The path `src/amuse/` holds AMUSE interface files, which are installed to
`$AMUSE_DIR/src/amuse/community/flash/`.

    src/amuse/
        base_grid_interface.F90
        interface.F90
        interface.py
        Makefile.prototype
        src/*

The path `src/flash/` holds FLASH interface files and add-on units, which are
installed to `$FLASH_DIR/`.  The directory structure mimics the FLASH4.5
source tree.

    src/flash/
        bin/
            setup_shortcuts.txt
        sites/
            Aliases
            cartesius.surfsara.nl/
            ...
        source/
            Driver/
            Grid/
            ...

The path `src/voramr/` holds Python-side utilities for the VorAMR
grid converter tool which are called by `torch_mainloop.py`.

    src/voramr/
        hdf5_convert.py
        kdtree.py
        ...
        

The path `utils/` provides user-contributed tools (AS-IS, support not
guaranteed) that may be of use for setting up, running, and/or analyzing Torch
simulations.


Simulation init condition setup
-------------------------------

`turb-sphere.py` creates a turbulent sphere in a 128^3 array for use with the
FLASH simulation `setup-cube-USM`.  Originally by Richard Wuensch (Astron.
Inst., Czech Acad. Sci.), with updates by J. Wall (Drexel).

`turb-velbox.py` creates a turbulent box (just velocity data) in a 128^3 array
for use with the FLASH simulation `Cube_AT` -- this is mostly for testing
purposes, and may be removed in the future.

`weighscale.py` reports useful information about the 128^3 plaintext arrays
output by `turb-sphere.py`.

__Example:__ create a 1e3 Msun sphere with radius 5 pc, virial ratio 0.2, T =
10 K, and ambient density 0.1 cm^-3 in a box of size 7^3 pc.

    ./turb-sphere.py --mass 1e3 --radius 5 --box_side 7 -v 0.2 -Ts 10 -n 0.1 --filename cube.m1e3.r5

This will output several files:

    cube.m1e3.r5        // Data array, ~200 MB
    cube.m1e3.r5.dat    // Metadata w/ input parameters, sphere properties
    cube.m1e3.r5velx.png    // x-y slice plot of x-velocity
    cube.m1e3.r5pres.png    // x-y slice plot of pressure
    cube.m1e3.r5temp.png    // x-y slice plot of temperature
    cube.m1e3.r5dens.png    // x-y slice plot of mass density

__Example:__ create a 128^3 cube of turbulent velocities.

    ./velbox.py -f vel128 --clobber

This will output the data array `vel128`, overwriting any existing file.

__Example:__ report some information about the cube array file.

    ./weighscale.py cube.m1e3.r5

This reports to STDOUT something like:

    from metadata: Msph = 1.000000e+03 Msun, Rsph = 5 pc, box = 7 pc
    loading cube.m1e3.r5
    done loading! elapsed: 0:02:43.201788
    Rsph (pc) [  5.   5.  10.  50. 500.]
    box halfwidth (pc) [  7.    7.   12.5  55.  600. ]
    Total mass (Msun): [1.005e+03 1.005e+03 5.723e+03 4.875e+05 6.329e+08]
    Ambient mass (Msun): [7.138e+00 7.138e+00 3.677e+01 2.596e+03 3.872e+06]
    Total - ambient mass (Msun): [9.978e+02 9.978e+02 5.686e+03 4.849e+05 6.290e+08]
    t_freefall (Myr): [2.865 2.865 3.394 4.109 3.608]

VorAMR
======
VorAMR is a utility developed within the Torch framework which allows for the conversion 
of output data from one hydrodynamical software suite to another. Specifically, VorAMR
takes data from the Vornoi mesh code AREPO and interpolates the data into FLASH for use
in Torch.

In principle, VorAMR can be extended to interpolate data from ANY type of hydrodynamical
code into FLASH/Torch including other AMR grid codes, or SPH codes.

How VorAMR works
----------------
1. Voronoi cell position data sent to FLASH. Positions and field value (density, internal energy, velocity, etc.) data sent to AMUSE.
2. FLASH views the Voronoi positional data as particles and builds an AMR grid that satisfies the refine-on-particle criteria. 
The particles used for refinement are erased from memory.
3. AMUSE constructs a KDtree with field values assigned to leaf nodes.
4. For each block data-structure, FLASH sends AMUSE an empty 16x16x16 matrix.
5. AMUSE passes the entire matrix into a N-dimensional nearest neighbor interpolation routine which maps field data 
corresponding to the closest KDtree leaf for each matrix cell, creating a 4 dimensional matrix of matrices.
6. AMUSE sends the populated block matrix back to FLASH via the interface, mapping each field value to the cells within the 
FLASH AMR grid structure.
7. FLASH outputs a refined AMR grid with all cells populated with field data ready for use in Torch, Torch then proceeds as normal.

See the [Quickstart Guide](https://torch-sf.bitbucket.io/quickstart.pdf) for how to use VorAMR.

Credits
=======

The Torch code includes contributions by:

* Eric Andersson
* Sabrina Appel
* Claude Cournoyer-Cloutier
* Will Farner
* Joseph Glaser
* Ralf Klessen
* Sean Lewis
* Steve McMillan
* Mordecai-Mark Mac Low
* Andrew Pellegrino
* Brooke Polak
* Simon Portegies Zwart
* Steven Rieder
* Aaron Tran
* Joshua Wall
* Maite Wilhelm

In addition to the FLASH and AMUSE codes, Torch also builds upon software by:

* Christian Bacyznski
* Robi Banerjee
* Moo Kwang Ryan Joung
* Juan Camilo Ibáñez-Mejía
* Daniel Seifried
* Long Wang
* Richard Wünsch

Torch also acknowledges the contributions of everyone who helped build FLASH and AMUSE.
