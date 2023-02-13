#!/bin/sh
#SBATCH --job-name=torch-test
#SBATCH -n 6
#SBATCH --time=0-01:00:00

mpiexec --mca orte_base_help_aggregate 0 -n 1 python bridge_multiples.py

# Some system/hardware-specific call options

# On Rome nodes (on e.g. Snellius) an error can occur in wireup.c when starting a worker,
# this can be fixed with these options:
#mpiexec --mca orte_base_help_aggregate 0 -x UCX_ATOMIC_MODE=cpu -x UCX_NET_DEVICES=mlx5_0:1 -x UCX_RC_MLX5_MAX_NUM_EPS=inf -n 1 python bridge_multiples.py

# Some MPI call options for experimenting

# OpenMPI 4.x or newer defaults to UCX rather than infiniband ports,
# you may need to override that policy
#mpiexec  --mca btl_openib_allow_ib 1 -n 1 python bridge_multiples.py

# OpenMPI
#export OMPI_MCA_mpi_warn_on_fork=0
#mpirun --mca orte_base_help_aggregate 0 --mca btl_openib_warn_no_device_params_found 0 -n 1 python bridge_multiples.py

# Intel MPI
#mpiexec.hydra -n 1 python bridge_multiples.py
#-hostfile ./local_host.txt 
#srun --mpi=pmi2 -n 1 python bridge_multiples.py

# MVAPICH2
#mpirun_rsh -np 1 -hostfile ./local_host.txt MV2_SUPPORT_DPM=1 MV2_ON_DEMAND_THRESHOLD=128 MV2_VBUF_TOTAL_SIZE=128 MV2_IBA_EAGER_THRESHOLD=128 python ./bridge_multiples.py
##-hostfile local_host.txt 
