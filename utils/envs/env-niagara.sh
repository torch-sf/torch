module --force purge
module load CCEnv arch/avx512 StdEnv/2020
module load gcc/9.3.0 #Reloads correct version of opnempi 4.0.3 automatically
module load python/3.8.10 #Default python 3.8
module load hdf5-mpi/1.12.1 #Loads h5py automatically
module load gsl/2.6
module load gmp/6.2.0
module load mpi4py/3.1.3
module load scipy-stack/2023a

source /home/a/asills/cournoyc/binaries/torch_env/bin/activate

export AMUSE_DIR=/home/a/asills/cournoyc/binaries/amuse
export FLASH_DIR=/home/a/asills/cournoyc/binaries/FLASH4.6.2
export TORCH_DIR=/home/a/asills/cournoyc/binaries/torch

export PYTHONPATH=$PYTHONPATH:$TORCH_DIR
export PYTHONPATH=$PYTHONPATH:$TORCH_DIR/src
export PYTHONPATH=$PYTHONPATH:$AMUSE_DIR/test
export PYTHONPATH=$PYTHONPATH:$AMUSE_DIR/src
