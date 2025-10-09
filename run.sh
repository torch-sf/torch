#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=0-0:05:00
#SBATCH --mail-user=cournoyc@mpa-garching.mpg.de
#SBATCH --mail-type=NONE
#SBATCH --output=torch.out
#SBATCH --error=torch.err
#SBATCH --ntasks=113
#SBATCH --cpus-per-task=1
#SBATCH --partition=p.exclusive

. /u/cournoyc/code/bubble/env.sh

export UCX_LOG_LEVEL=FATAL

ulimit -s unlimited
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=8

export UCX_TLS=rc,sm,self #Needed to prevent MPI_Comm_Spawn errors
mpiexec -x UCX_TLS -n 1 python torch_user.py
#mpiexec -n 1 python torch_user.py
