# FAQ

## General

### What resolution should I run at? 

Torch is considered converged at ~0.3 pc resolution. This is the maximum recommended cell size for production simulations. Runs at lower resolution than this produce **unphysical results**.
 
### How can I run a vanilla FLASH problem with Torch installed?

Torch requires a dependency on `rndMT`. 
You will need to add `rndMT.o` somewhere in your `FLASH` Makefile. 

### Why is my compilation failing with an HDF5 related error?

The automatic HDF5 library detection fails on some HPC systems. This can be due to unique system path naming conventions. For example, on TACC systems, the HDF5 path is set to `TACC_HDF5_DIR`, but autoconf expects it to be at `HDF5_DIR` or `HDF5_ROOT`. Set one of these environment variables to your installation, and recompile. 

## `FERVENT`

### Why is my run stalling in "Entering RadRay"?

You likely exceeded the max number of particles and rays. 
Raise the `pt_maxPerProc` parameter in your `flash.par` by a factor of 10-100, and restart.

## `VETTAM`

### Why am I running out of memory (OOM Error) when runnning with VETTAM?

When running at high resolution with a lot of blocks, the ray tracing step can become expensive. 
This is because ray tracing in VETTAM is done from every cell. The OOM will likely occur during the 
3DRT ray tracing step. You can either use more processors, or increase `rt_nrOfAngleGroups` in your `flash.par`.
This basically splits the raytracing step into the desired number of steps. For example, increasing `rt_nrOfAngleGroups`
from 1 to 2 halves the memory load, but can increase runtime due to reallocation steps.

### Why does my run keeps aborting with message on "[VETTAM]: ... subcycling limit reached."?

It is likely that extreme gradients in opacity have made it difficult for `VETTAM` to converge, likely caused by floored opacities at dust sputtering temperature. Try increasing `rt_maxits` if its lower than 2000, or decreasing your cfl. 

If none of the above methods work, try raising the opacity floors set in `torch/src/flash/source/physics/RadTrans/RadTransMain/VETTAM/rt_setOpacity.F90`:
```fortran
#ifdef TEMP_VAR
          !Destroy dust opacities due to thermal sputtering if gas temperatures > 10^6 K
          !TODO: Implement non-thermal sputtering
          if(solnData(TEMP_VAR,i,j,k) .gt. he_dust_sputter_temp) then
            opac_planck = 1.0e-27 ! set floor to prevent NANs in boundary conditions
            opac_rosseland = 1.0e-27
          endif
#endif
```
The maximum value you should use is ~$1.0e-26$. Setting the floor too high will result in too much radiation being absorbed by dust that should be sputtered. 