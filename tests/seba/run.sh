#!/bin/bash
#SBATCH --job-name=seba-test
#SBATCH --time=0-1:00:00
#SBATCH --array=1-1%1
#SBATCH --account=rrg-asills-ac
#SBATCH --mail-user=cournoyc@mpa-garching.mpg.de
#SBATCH --mail-type=ALL
#SBATCH --output=seba.out
#SBATCH --error=seba.err
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

mkdir data

# Run test for SeBa in AMUSE (without Torch)
python $TORCH_DIR/tests/seba/run_seba_amuse.py

# Create uniform n=0.1 cm^-3, T=100K ICs for stellar evolution tests
#python $TORCH_DIR/utils/ic-generator/turb-uniform.py -rho 0.1 -temp=100 -o seba_ic

# Run test for SeBa in Torch, without restart and with restarts
#cp $TORCH_DIR/tests/seba/flash.par.seba flash.par
#mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/seba/run_seba_torch.py

# Plot the results
python $TORCH_DIR/tests/seba/plot_seba.py