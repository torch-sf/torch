#!/bin/sh
#SBATCH --job-name=rad_tests
#SBATCH -n 50
#SBATCH --ntasks-per-node=50       # adjust to fill your CPUs
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH -p genoa
#SBATCH --time=5:00:00

source ~/torch-new-build/torch.env

ulimit -s unlimited
export OMP_STACKSIZE=128M

# D-type expansion
# create uniform n=100 cm^-3, T=100K ICs for dtype test
python turb-uniform.py -rho 100 -temp=100 -o dtype_ic
cp flash.par.dtype flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python torch_user.py

# R-type expansion
# create uniform n=100 cm^-3, T=1e4K isothermal ICs for rtype test
python turb-uniform.py -rho 100 -temp=10000 -o rtype_ic
cp flash.par.rtype flash.par
mpirun --mca orte_base_help_aggregate 0 -n 1 python torch_user.py

python plot_stromgren.py
