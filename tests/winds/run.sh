#!/bin/sh
#SBATCH --job-name=wind_test
#SBATCH --nodes=1
#SBATCH --ntasks=83
#SBATCH --ntasks-per-node=83
#SBATCH --cpus-per-task=1
#SBATCH --partition=genoa
#SBATCH --time=48:00:00

SYSTEM="snellius"
TORCH_ENV="/home/eanderson/code/2026/torch-sf/torch.env"

. "$TORCH_ENV"
. "$TORCH_DIR/utils/mpi_setups.sh"

#===============================================================================
# Run the wind test
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

# Create run directory setup
cd "$SLURM_SUBMIT_DIR"

ln -sf ${TORCH_DIR}/cool.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/opacity.inp
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kP_h2001.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kR_h2001.dat

mkdir data

cp $TORCH_DIR/tests/winds/flash.par.wind flash.par
cp $TORCH_DIR/torch_user.py .

# Make initial conditions file for wind test
python $TORCH_DIR/utils/ic-generator/turb-uniform.py -rho 10 -temp=100 -o wind_ic

# Launch simulation
$MPI_LAUNCHER $MPI_ARGS -n 1 python $TORCH_DIR/tests/winds/run_test.py

python3 $TORCH_DIR/tests/winds/test_winds.py ./ --plot

echo
echo "Finished at $(date)"
