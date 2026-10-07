#!/bin/bash
#SBATCH --partition=PARTITION_NAME
#SBATCH --nodes=1
#SBATCH --tasks-per-node=12
#SBATCH --time=01-00:00:00
#SBATCH --mem=5120MB

#SBATCH --output=slurm_output/polaris_run_array%a.out
#SBATCH --error=slurm_output/polaris_run_array%a.err

#SBATCH --array=0-1

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

echo "Job started on $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Number of MPI tasks: $SLURM_NTASKS"

## ACTIVATE CONDA OR PYTHON ENVIORNMENT
echo "Conda environment loaded"

## No trailing slash for the following two paths
export POLARIS_DIR="/PATH/TO/POLARIS/REPO/NO/TRAILING/SLASH"
export TORCH_DIR="/PATH/TO/TORCH/REPO/NO/TRAILING/SLASH"

python3 "${TORCH_DIR}/utils/polaris-torch/run_polaris.py" -r "T_dust" -p "/PATH/TO/TORCH/DATA" -s 189 -o "/PATH/TO/OUTPUT/FOLDER/" -ep 1e-2 -el "n2Z_n2D"

python3 "${TORCH_DIR}/utils/polaris-torch/run_polaris.py" -r "T_dust" -p "/PATH/TO/TORCH/DATA" -f "${TORCH_DIR}/utils/polaris-torch/example_run/example_list_snapshot.txt" -l $SLURM_ARRAY_TASK_ID -o "/PATH/TO/OUTPUT/FOLDER/" -ep 1e-2 -el "n2Z_n2D"

python3 "${TORCH_DIR}/utils/polaris-torch/run_polaris.py" -r "detectors_spherical"  -ns 7 -sl "sph7_Mmin20" -p "/PATH/TO/TORCH/DATA" -f "${TORCH_DIR}/utils/polaris-torch/example_run/example_list_snapshot.txt" -l $SLURM_ARRAY_TASK_ID -o "/PATH/TO/OUTPUT/FOLDER/" -ep 1e-2 -el "n2Z_n2D"

