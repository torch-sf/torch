#!/bin/sh
#
#===============================================================================
# Example job script for Torch
#
# Unless changed, this example runs a low-resolution turbulent sphere test on 
# 6 cpus for around 10 minutes. 
#
# Copy this file, set the system flag to the relevant machine and submit. The
# script will copy the files needed for the run, create a copy of the source
# code, and execute the simulation.
#
# Slurm variables are set using the #SBATCH flag. The resources and runtime
# variables have been preset for this script but additional variables are
# included but commented out below.
#
# To learn more about slurm, check out https://slurm.schedmd.com/sbatch.html
#===============================================================================

## Name of job
#SBATCH --job-name=turbsph_test

#---------------------------
# Resources
#---------------------------

## Number of nodes
#SBATCH --nodes=1

## Number of tasks
#SBATCH --ntasks=6

## Number of tasks per node
#SBATCH --ntasks-per-node=6

## Number of cpus per task
#SBATCH --cpus-per-task=1

#---------------------------
# Queue / runtime
#---------------------------

## Partition to run on
#SBATCH --partition=genoa

## Job run time in HH:MM:SS
#SBATCH --time=00:10:00

#===============================================================================
# Other settings that might be useful 
#===============================================================================

## Split log and error messages
##SBATCH --output=slurm-%j.out
##SBATCH --error=slurm-%j.err

## Reserve the entire node
##SBATCH --exclusive

## Request a specific amount of memory
##SBATCH --mem=8G

## Memory per CPU
##SBATCH --mem-per-cpu=2G

## Send email notifications
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.com

## Run on multiple nodes
##SBATCH --nodes=2
##SBATCH --ntasks-per-node=32

## Use OpenMP threads
##SBATCH --cpus-per-task=8

## Run on a different partition
##SBATCH --partition=debug

## Use a specific account/project
##SBATCH --account=my_project

## Specify a QoS
##SBATCH --qos=normal

#===============================================================================
# Job setup, update these to match you system.
#===============================================================================
## Available systems include: default, snellius
SYSTEM="default"

## Set the path to the torch environment you created.
TORCH_ENV="path/to/torch.env"

#===============================================================================
# System configurations 
#===============================================================================
set -e

if [[ ! -f "$TORCH_ENV" ]]; then
    echo "Could not find $TORCH_ENV"
    exit 1
fi

. "$TORCH_ENV"
. $TORCH_DIR/utils/mpi_setups.sh

#===============================================================================
# Run the application
#===============================================================================
echo "System        : $SYSTEM"
echo "Torch path    : $TORCH_DIR"
echo "Job ID        : $SLURM_JOB_ID"
echo "Job name      : $SLURM_JOB_NAME"
echo "Nodes         : $SLURM_JOB_NODELIST"
echo "MPI tasks     : $SLURM_NTASKS"
echo "CPUs per task : $SLURM_CPUS_PER_TASK"
echo "Started at    : $(date)"
echo

# Create turbulent sphere test
cd "$SLURM_SUBMIT_DIR"
bash "$TORCH_DIR/utils/setup_simulation.sh"
python3 "$TORCH_DIR/utils/ic-generator/turb-sphere.py" -np -s 42 -b 10 -m 1e4 -r 7 -f cube128

# Launch simulation
$MPI_LAUNCHER $MPI_ARGS -n 1 python torch_user.py

echo
echo "Finished at $(date)"
