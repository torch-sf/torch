README
======

Quick-start
-----------

Set the environment variables:

    export AMUSE_DIR=/path/to/amuse
    export FLASH_DIR=/path/to/FLASH4.5
    export TORCH_DIR=/path/to/torch

Then, run the script:

    ./install.sh

You will be prompted to confirm your environment variables, and continue with
installation.  The installer will then automatically copy interface code from
`$TORCH_DIR/src` to their correct locations in your AMUSE and FLASH repos.
It will also ensure that the FLASH4.5-a patch for HDF5 >1.10.x compatibility is
applied in your FLASH repo.

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
    ionizingflux.py

And, some code to create turbulent initial conditions (dens, pres, vel, etc.)
for your simulations:

    hAc_b_2.0E-17_e_0.021_FUV_1.69.dat
    turb-sphere.py
    turb-velbox.py
    weighscale.py

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

Contributing
------------

Torch is an open-source code!  Please see [CONTRIBUTING.md](CONTRIBUTING.md) to
learn how to get involved.

The Torch code includes contributions by:

* Sabrina Appel
* Raúl Domínguez
* Ralf Klessen
* Sean Lewis
* Steve McMillan
* Mordecai-Mark Mac Low
* Andrew Pellegrino
* Simon Portegies Zwart
* Aaron Tran
* Joshua Wall

In addition to the FLASH and AMUSE codes, Torch also builds upon software by:

* Christian Bacyznski
* Robi Banerjee
* Moo Kwang Ryan Joung
* Juan Camilo Ibáñez-Mejía
* Daniel Seifried
* Richard Wünsch
