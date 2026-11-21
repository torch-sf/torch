#!/bin/bash
#SBATCH --job-name=nbody-tests
#SBATCH --time=0-1:00:00
#SBATCH --array=1-1%1
#SBATCH --account=rrg-asills-ac
#SBATCH --mail-user=cournoyc@mpa-garching.mpg.de
#SBATCH --mail-type=ALL
#SBATCH --output=nbody.out
#SBATCH --error=nbody.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --ntasks=192
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug

. /home/cournoyc/new-repo/torch.env

# link files to avoid hdf5 errors                                                                                                                                                                                                   
ln -sf ${TORCH_DIR}/cool.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/opacity.inp
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kP_h2001.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kR_h2001.dat

export UCX_LOG_LEVEL=FATAL

ulimit -s unlimited
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=1

mkdir data_binary1au
mkdir data_binary55au
mkdir data_binary600au

# Create uniform n=0.1 cm^-3, T=100K ICs for N-body tests
python $TORCH_DIR/utils/ic-generator/turb-uniform.py -rho 0.1 -temp=100 -o nbody_ic

# Run test for binaries in regularization, direct N-body, and hybrid regimes of PeTar
cp $TORCH_DIR/tests/nbody/flash.par.binary1au flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/nbody/run_binary1au.py
cp $TORCH_DIR/tests/nbody/flash.par.binary55au flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/nbody/run_binary55au.py
cp $TORCH_DIR/tests/nbody/flash.par.binary600au flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/nbody/run_binary600au.py

python $TORCH_DIR/tests/nbody/plot_nbody.py
