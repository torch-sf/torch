# FAQ

## What resolution should I run at? 

Torch is considered converged at ~0.3 pc resolution. This is the maximum recommended cell size for production simulations. Runs at lower resolution than this produce **unphysical results**.
 
## Why is my run stalling in "Entering RadRay"?

You likely exceeded the max number of particles and rays. 
Raise the `pt_maxPerProc` parameter in your `flash.par` by a factor of 10-100, and restart.

## How can I run a vanilla FLASH problem with Torch installed?

Torch requires a dependency on `rndMT`. 
You will need to add `rndMT.o` somewhere in your `FLASH` Makefile. 

