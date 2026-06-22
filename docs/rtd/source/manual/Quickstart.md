# Quickstart
Torch has many moving parts. You’ll likely have to adjust the procedure below 
for your personal setup. In any case, don’t despair! Torch has been successfully 
built and run on a variety of systems.
We assume that you’re using (1) the bash shell, and (2) conda or a virtualenv 
for Python package management. But, the setup procedure should work for any 
shell and Python package manager with some adjustments.

The deprecated pdf version of the Quickstart can be found [here](https://torch-sf.github.io/quickstart.pdf). 

## Overview
Torch consists of a collection of Fortran files that extend the FLASH MHD code, 
bindings for AMUSE for that extended FLASH, and a set of Python scripts that use
AMUSE to run a simulation combining the extended FLASH with the other codes.
Installing Torch therefore entails 1) creating an environment and installing the
dependencies, including AMUSE, 2) installing Torch, with its extended FLASH and
Python parts, and 3) installing the FLASH AMUSE bindings to tie it all together.

## Environment
First, start a new shell session. You may want to check your .bashrc for any old
MPI or HDF5 configuration that could conflict with Torch setup – module load
commands, environment variable exports, et cetera – and comment them out.
Next, it is probably a good idea to make a project directory in which you can
download and install all the pieces:

```
mkdir torch_project
```

You can call it something different or arrange the files differently if you
want. In particular if you’re on an HPC cluster with a limited amount of space
available in your home directory, then you’ll want to use another location
instead, for example a scratch directory. The documentation for the machine
should tell you what is available. You’ll need about 2GB for the base
installation, plus space for model output which can also easily grow to multiple
gigabytes.

Next, we’ll need an environment to install the Python bits into, and possibly
the native tools and dependencies as well. There are two options here. If you’re
on an HPC machine, then you’ll probably want to use native tools and
dependencies, and a Python virtual environment for the Python bits. On a local
machine, you can do that as well, or use conda instead to get everything
contained in a single environment. Both options are described here.

## Local, with conda
We’ll assume here that you already have conda installed and available. If not,
we recommend installing the miniforge distribution from 
[https://github.com/conda-forge/miniforge][1].
First, we’ll create and activate a new conda environment, and set it to use
conda-forge:

[1]: https://github.com/conda-forge/miniforge

```
conda create -n Torch
conda activate Torch
```

Note that AMUSE is designed to work with the conda-forge package channel.
Conda-forge contains a very large variety of scientific and other packages, so
you’ll find everything you need for Torch as well as other tools you may want to
use.

Mixing conda channels is not recommended, as it leads to package
incompatibilities, long waits for the solver to try to deal with incongruent
sets of packages when installing, and hard to diagnose issues. So we strongly
recommend sticking with conda-forge only. Miniforge will do that by default.

If you open a new shell, then the environment needs to be reactivated using the
`conda activate Torch` command before you do anything else.

Once we have a conda environment, we can install the prerequisites inside of it,
once again making sure our packages come from conda-forge (some versions of
conda can be a bit stubborn on this):

```
conda install -c conda-forge --override-channels pip wheel c-compiler cxx-compiler fortran-compiler 'gfortran<14' python coreutils patch curl make openmpi hdf5 zlib 'docutils>=0.6' 'mpi4py>=1.1.0' 'numpy>=1.2.2' 'h5py>=1.1.0' scipy astropy jupyter pandas seaborn matplotlib yt
```

AMUSE is not yet available directly from conda-forge, so we’ll have to install
it from source. First, we need to download and unpack it (into our project
directory):

```
git clone https://github.com/amusecode/amuse.git
```

AMUSE has had a major upgrade to its installer recently, so be sure to get
v2025.9.0 or a later version, older ones won’t work.

Next, we can move into the newly created AMUSE source directory, and run the
installer:

```
cd amuse
./setup
```

This will produce quite a bit of output, including an overview of installable
codes, disabled codes and the reasons why, and follow-up instructions. If all
went well then you should now be able to install the AMUSE packages you need
into your conda environment using

```
./setup install framework kepler petar ph4 seba smalln
```

Note that Torch doesn’t necessarily use all of these, so if you know which ones
you want, then you can install only those. At any rate you can return to the
AMUSE directory later and use `./setup` to install more codes.

That concludes setting up the environment and the dependencies, so you can
scroll past the next section (one environment is enough, really!) and continue
with installing Torch.

## HPC, with virtualenv and modules
As an alternative to Conda, you can also use a Python virtual environment to
install the Python dependencies into. Dependencies for the Fortran and C code
in Torch will then have to be installed in some other way. Linux computers
typically have a built-in package manager like `apt` or `dnf`, while on a Mac
Homebrew and MacPorts are commonly used. On HPC machines, you typically load
the necessary modules.

If you’re on a local machine and have one of the above package managers
available, then you should continue with the instructions below, but skip the
loading of modules and just make the virtual environment and install the Python
dependencies into it. If you then continue with installing AMUSE, the AMUSE
installer will guide you through installing the dependencies using whichever
package manager you have available. (Note that you won’t need everything it
suggests, because Torch doesn’t use all of AMUSE. You can save some disk space
by limiting your selection.)

On an HPC cluster or similar environment, the main tools and dependencies are
likely already installed and available as modules. To use them, we’ll have to
activate them using the module command.

First, let’s see what’s installed:

```
module avail
```

This will show a long (potentially very long) list of installed packages, or
modules.

On some systems, this will only show a few years, e.g. `2023` `2024` `2025`. In that
case, use `module load 2025` (the latest is usually best) and then try 
`module avail` again to get the list.

Have a look through the list and see if you can find modules named `GCC`,
`HDF5`, `OpenMPI` (or some other MPI implementation, like `MPICH`), `make`, and
`Python`. The exact names will vary from machine to machine. If the list is
long, it’s better to let the computer search for you using e.g.
`module keyword HDF5`.

Often, there are multiple versions of a package available. The codes in Torch
(and AMUSE) are fairly old, so probably any version will work. If you see
different packages with names like `gompi` or `iimpi`, use the `gompi` ones, as
those match GCC and OpenMPI.

You may not need make if it’s already available. If the command
```
make --version
```
prints some text starting with GNU Make, then you’re good to go.

Once you have found the needed module names, you can load the modules, for
example for Snellius using the 2025 software stack (you’ll need to modify the
exact module names to suit your machine):

```
module load 2025
module load GCC/14.2.0
module load OpenMPI/5.0.7-GCC-14.2.0
module load HDF5/1.14.6-gompi-2025a
module load make/4.4.1-GCCcore-14.2.0
module load Python/3.13.1-GCCcore-14.2.0
```

Next, we can make a Python virtual environment. This is a directory, which can
go into your project directory or wherever you want it.

```
python3 -m venv /path/to/torch_project/Torch-env
```

We can then activate the environment, and install the Python dependencies:

```
. /path/to/torch_project/Torch-env/bin/activate
pip3 install -U pip wheel scipy astropy jupyter pandas seaborn matplotlib yt
```

At this point, it’s probably convenient to create a file you can load that will
activate the modules and the virtual environment whenever you want to work with
Torch. To do that, create a file named Torch.env in your project directory
containing the module load commands and the activation of the environment, like
so:

```
module load 2025
module load GCC/14.2.0
module load OpenMPI/5.0.7-GCC-14.2.0
module load HDF5/1.14.6-gompi-2025a
module load make/4.4.1-GCCcore-14.2.0
module load Python/3.13.1-GCCcore-14.2.0
. /path/to/torch_project/Torch-env/bin/activate
```

with appropriate modifications for your site. Now you can activate everything in
one command using `. Torch.env` in your project directory (note the period and
the space at the beginning, they’re required).

With the environment set up and activated, we can install AMUSE and most of the
codes that Torch uses into it. First, we need to download and unpack it (into
our project directory):

```
git clone https://github.com/amusecode/amuse.git
```

AMUSE has had a major upgrade to its installer recently, so be sure to get
v2025.9.0 or a later version, older ones won’t work.

Next, we can move into the newly created AMUSE source directory, and run the
installer:

```
cd amuse
./setup
```

This will produce quite a bit of output, including an overview of installable
codes, disabled codes and the reasons why, and follow-up instructions.

If you’re on a local machine and don’t already have the dependencies installed,
then ./setup will give you the required command to install them using your local
package manager. Note that you probably don’t need all the dependencies it
suggests installing, as that is the complete set and Torch only uses a few of
AMUSE’s many codes. Compilers for C, C++ and Fortran, MPI, HDF5, make and Python
should go a long way.

If all went well then you should now be able to install the AMUSE packages you
need into your virtual environment using

```
./setup install framework kepler petar ph4 seba smalln
```


Note that Torch doesn’t necessarily use all of these, so if you know which ones
you want, then you can install only those. At any rate you can return to the
AMUSE directory later and use ./setup to install more codes.

With that, we have a virtual environment that is ready to install Torch into.

## Installing Torch

The first step to installing Torch is to get a copy of FLASH. To do that, you
need to register, which may take a few days, and then download FLASH 4.6.2. The
FLASH homepage [(http://flash.uchicago.edu/site/flashcode)][2] will show you how.

[2]: http://flash.uchicago.edu/site/flashcode

Once you have a FLASH4.6.2.tar.gz file in your project directory, you need to
unpack it, and then set the FLASH DIR environment variable so that the Torch
installer can find it:

```
cd torch_project
tar xzf FLASH4.6.2.tar.gz
export FLASH_DIR=${PWD}/FLASH4.6.2
```

Next, we can download Torch, and set TORCH DIR:

```
git clone https://github.com/torch-sf/torch.git export TORCH_DIR=${PWD}/torch
```

Then we can move into the Torch directory, and run its installer:

```
cd ${TORCH_DIR} ./install.sh
```

This will add the Fortran parts of Torch to FLASH. Next, we need to make the
Python parts of Torch available:

```
export PYTHONPATH=$PYTHONPATH:$TORCH_DIR/src
```

You will want to add the exports of TORCH DIR and FLASH DIR to your Torch.env,
as well as the PYTHONPATH, like so:

```
module load 2025
module load GCC/14.2.0
module load OpenMPI/5.0.7-GCC-14.2.0
module load HDF5/1.14.6-gompi-2025a
module load make/4.4.1-GCCcore-14.2.0
module load Python/3.13.1-GCCcore-14.2.0
FLASH_DIR=/path/to/FLASH4.6.2
TORCH_DIR=/path/to/torch
. /path/to/torch_project/Torch-env/bin/activate export PYTHONPATH=$PYTHONPATH:$TORCH_DIR/src
```

With Torch installed into FLASH, we can now go to the FLASH directory and
configure it:

```
cd $FLASH_DIR
./setup Cube -auto -3d +amuseUSM +amuseMG +rayPE +HandCnew +amuseSinksAndStars +amuseWind +enerinj +supportPPMupwind +pm4dev -maxblocks=100 +cube16 --site=torch_auto
```

The final piece of the puzzle are Torch’s AMUSE FLASH bindings. These make it
possible for the Python parts of Torch to talk to FLASH, including its Torch
additions. Check that you still have your virtual environment active, and then
you can compile and install FLASH and the bindings like this:

```
cd $TORCH_DIR
make -C src/amuse/torch_amuse_flash develop-torch-amuse-flash
```

With that, your environment should contain a working Torch installation. Time to
give it a try!

## Running the turbsph test
You are now ready to run a Torch simulation! We will simulate a 
$3 \times 10^3\,{\rm M}_{\odot}$, $5\,{\rm pc}$ radius gas sphere in 
$\pm 7\,{\rm pc}$ cubic domain with outflow boundary conditions. In about a 
free-fall time ($\approx 2\,{\rm Myr}$), the cloud will collapse under its own 
gravity, create a sink particle, and begin forming stars.

Create a simulation directory and use to the `setup_simulation.sh` script to
setup the necessary files.

```
mkdir /path/to/test_simulation
cd /path/to/test_simulation
cp $TORCH_DIR/ic/turbsph/setup_simulation.sh
bash setup_simulation.sh
```

Modify the `run.sh` according to you HPC system and submit the job
```
sbatch run.sh
```

The simulation should take 1–10 minutes to reach 2 Myr. It will use 6 threads: 
3 for FLASH, 1 for n-body integration (petar), 1 for stellar evolution, and 1 
for AMUSE. At around 1 Myr, a sink particle will form and begin creating stars. 
One or a few $>7\,{\rm M}_{\odot}$ stars may begin blowing a wind, so the time
step may drop; you may wish to end the simulation early.

If you increase the number of tasks requested, torch user.py will automatically 
assign more threads (workers) to FLASH. Log files will be dumped in your current 
working directory, and simulation outputs will be written to the data directory. 
If you re-start a simulation, Torch (FLASH) will overwrite existing data, but 
append to (some) existing logs. You may want to rename or delete your old logs 
to help keep track of your runs.

## Analyse the output

To analyze the outputs, we recomend using yt to load the data and create plots.
There are also scripts available to make basic analysis. First up we can create
movie frames of all the outputs. For this we need the python package `yt`. See
their [website][4] for documenation and examples.

[4]:https://yt-project.org/

```
conda install yt
conda clean --all
```

To make the movie frames use `movie.py` available in the `utils`. The default 
parameters are set for the turbsph test, but can be changed using arguments.

```
python $TORCH_DIR/utils/movie.py . -q rho -ap -f
```

Enjoy creating your own Torch setups and simulations!
