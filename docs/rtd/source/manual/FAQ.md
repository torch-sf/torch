# FAQ

## What resolution should I run at? 

Torch is considered converged at ~0.3 pc resolution. This is the maximum recommended cell size for production simulations. Runs at lower resolution than this produce **unphysical results**.

## How do I restart a simulation run?

<!--- To do: include some common pitfalls when restarting after a crash. -->

## Why did my run crash?

<!--- To do: list some common reasons for a run crashing and how to identify them.
Then point towards other FAQs specific to those problems.  -->

## How do I analyze my run?

<!--- To do: list some brief options for analyzing Torch runs. -->

##  Why is my simulation so slow?

<!--- To do: list typical run speeds for different categories of runs and 
some common pitfalls that may be slowing down your run. -->

## I am running on an HPC that is not listed in mpi_setup.sh. What do I do?

<!--- To do: describe process for updating mpi_setup.sh for a new HPC.  -->
 
## Why is my run stalling in "Entering RadRay"?

You likely exceeded the max number of particles and rays. 
Raise the `pt_maxPerProc` parameter in your `flash.par` by a factor of 10-100, and restart.

## How can I run a vanilla FLASH problem with Torch installed?

Torch requires a dependency on `rndMT`. 
You will need to add `rndMT.o` somewhere in your `FLASH` Makefile. 

## Can I expect feedback stars to always be in maximally refined regions?

No, you cannot. While the wind injection routine will refine the injection region to the maximum refinement level, other sources of feedback (e.g., radiation) are not guaranteed to do so. This is by design.
