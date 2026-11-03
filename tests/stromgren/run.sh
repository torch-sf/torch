#!/bin/sh
#SBATCH --job-name=rad_tests
#SBATCH -n 10
#SBATCH --ntasks-per-node=10       # adjust to fill your CPUs
#SBATCH --cpus-per-task=1
#SBATCH -p genoa
#SBATCH --time=5:00:00

#source /full/path/to/torch.env

# link files
ln -sf ${TORCH_DIR}/cool.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/opacity.inp
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kP_h2001.dat
ln -sf ${TORCH_DIR}/src/flash/source/physics/materialProperties/Opacity/RadTrans/Semenov/kR_h2001.dat

ulimit -s unlimited
export OMP_STACKSIZE=128M

mkdir data

# D-type expansion
# create uniform n=100 cm^-3, T=100K ICs for dtype test
python $TORCH_DIR/utils/ic-generator/turb-uniform.py -rho 100 -temp=100 -o dtype_ic
cp $TORCH_DIR/tests/stromgren/flash.par.dtype flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/stromgren/run_test.py

# R-type expansion
# create uniform n=100 cm^-3, T=1e4K isothermal ICs for rtype test
python $TORCH_DIR/utils/ic-generator/turb-uniform.py -rho 100 -temp=10000 -o rtype_ic
cp $TORCH_DIR/tests/stromgren/flash.par.rtype flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python $TORCH_DIR/tests/stromgren/run_test.py

python $TORCH_DIR/tests/stromgren/plot_stromgren.py
