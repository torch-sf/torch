# Detailed Installation Guide

## VETTAM

VETTAM requires PETSC version 3.19.x. You can use a module if available, or download and install it with the following steps:
```
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.6.tar.gz
cd petsc-3.19.6
./configure --with-cc=mpicc --with-cxx=mpic++ --with-fc=mpif90
```
Then use the make command that petsc tells you to use after configure. It should look something like this:

```
make PETSC_DIR=/path/to/petsc-3.19.6 PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=/path/to/petsc-3.19.6 PETSC_ARCH=arch-linux-c-debug check
```

Once you have petsc installed/loaded, define the following environment variables. These must be defined in your run script as well.
```
export PETSC_DIR=/path/to/petsc-3.19.x
export PETSC_ARCH=arch-linux-c-debug
export LD_LIBRARY_PATH=$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

Compiling VETTAM:

In the FLASH Makefile.h:
- Add ``` -I${PETSC_ARCH}/include -I${PETSC_DIR}/include -L${PETSC_ARCH}/lib -lpetsc ``` to ``` FFLAGS_OPT```.
- Add LIB_PETSC  = -L${PETSC_ARCH}/lib -lpetsc -lz

In amuse/src/community/flash/Makefile.h:
- Add ``` -I${PETSC_ARCH}/include -I${PETSC_DIR}/include -L${PETSC_ARCH}/lib -lpetsc ``` to ``` FCFLAGS```.

Finally, in the setup step of the FLASH compilation, replace ```+rayPE``` with ```+vettam```. For example:

```
./setup Cube -auto -3d +amuseusm +amusemg +vettam +handcnew +amusesinksandstars +amusewinds +enerinj +supportppmupwind +pm4dev +cube16 -maxblocks=100
```

Running VETTAM: 

The default parameters for VETTAM can be found in ```flash.par.vettam```. The file ```opacity.inp``` needs to be linked to your run directory for VETTAM to work.  Below is a table of the relevant vettam parameters that are useful to know and change, and some helpful notes about their uses.

| Parameter        | Notes              |
| -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| rt_rtol                                                                    | Tolerance of implicit radiation update solve. The maximum value to set this if petsc isn't converging is 1e-2.                                                                                                                                                                                                                                                                    |
| rt_maxits                                                                  | Maximum iterations before attempting subcycling. This can also be raised to help with convergence, but less helpful than raising the tolerance.                                                                                                                                                                                                                                   |
| rt_healpix_nSide                                                           | This sets the number of healpix tessellations for raytracing from each cell. Must be multiple of 2. 2 is recommended; higher values are more accurate but much slower and memory intensive.                                                                                                                                                                                       |
| rt_nrOfAngleGroups                                                         | This controls how many angle groups are used for the raytracing step. If you have many blocks and VETTAM runs out of memory, this splits the raytrace steps into multiple, dividing the raytracing step over the specified number of angle groups. raytrace*.log file will inform the memory consumption of raytracing per cpu and will help in setting this parameter correctly. |
| sigma_star                                                                 | This is the Gaussian source size in units of dx for radiation injection from each star. Keeping this at 0.1 ensures that VETTAM still treats the stars as point sources of radiation.                                                                                                                                                                                             |
| noSink_VETTAM                                                              | This flag sets whether to use VETTAM if there are no sinks. This should be kept as true, because you can have stars with no sinks.                                                                                                                                                                                                                                                |
| noStar_VETTAM                                                              | This flag sets whether to use VETTAM if there are no stars. This can be set to false, since we can skip the expensive raytracing step if there are no sources on the grid.                                                                                                                                                                                                        |
| rt_etens                                                                   | This sets the closure method for VETTAM. For production runs, this should always be "vet".                                                                                                                                                                                                                                                                                        |
| rt_linearsolver                                                            | This sets the linear solver for the PETSC solve step. This can be changed if there are convergence issues, but I've found this rarely helps.                                                                                                                                                                                                                                      |
| preconditioner                                                             | This sets the preconditioner for the PETSC solve step. This can be changed if there are convergence issues, but I've found this rarely helps.                                                                                                                                                                                                                                     |
| rt_freqbands = 2<br>rt_band1 = "EUV"<br>rt_band2 = "FUV"                   | This sets the number and type of frequency bands to solve.                                                                                                                                                                                                                                                                                                                        |
| dust_euvopac_value                                                         | This sets the opacity of dust to EUV photons.                                                                                                                                                                                                                                                                                                                                     |
| dust_fuvopac_value                                                         | This sets the opacity of dust to FUV photons.                                                                                                                                                                                                                                                                                                                                     |
| dusttoGasRatio                                                             | This sets the dust to gas ratio is units of the solar value.                                                                                                                                                                                                                                                                                                                      |
| ion_sigmaH                                                                 | This sets the cross section of H atoms to EUV Photons.                                                                                                                                                                                                                                                                                                                            |
| xl_rad_bc<br>xr_rad_bc<br>yl_rad_bc<br>yr_rad_bc<br>zl_rad_bc<br>zr_rad_bc | This sets the radiation boundary conditions. The default value of 6 is Marshak boundary conditions. This must stay at 6 to prevent convergence issues in the FUV band.                                                                                                                                                                                                            |


## Debugging
