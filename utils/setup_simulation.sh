#!/usr/bin/env bash

# Example bash script to setup a Torch simulation.
# Requires AMUSE_DIR, FLASH_DIR, and TORCH_DIR to be declared.
# Execute this from the directory in which you want to run Torch.
# User should manually specify a FLASH object directory, and customize this
# script according to their needs.
# -- Aaron Tran, 2020 Feb 19

objdir="object"
#objdir="objectcube"
#objdir="objecteinj"

if [[ -e "amr_runtime_parameters.dump" ]]; then
  echo "ERROR: FLASH already run, exiting!"
  exit 255
fi

# -------------------------
# Run versioning / metadata
# -------------------------

# terse, report commit hash only
#git -C ${AMUSE_DIR} log --format='%H' -n 1 > commit-amuse
#git -C ${TORCH_DIR} log --format='%H' -n 1 > commit-torch

# commit hash, date, message
git -C ${AMUSE_DIR} log --format='medium' -n 1 > commit-amuse
git -C ${TORCH_DIR} log --format='medium' -n 1 > commit-torch

# FLASH setup call is also stored in chk/plt HDF5 files, so formally redundant
rsync -a "${FLASH_DIR}/${objdir}/setup_call" .

tar -cz -C ${FLASH_DIR} --checkpoint=100 --totals \
  --dereference --exclude=flash4 --exclude=*.o --exclude=*.mod \
  --file=code.tar.gz "${objdir}"

# ----------------
# torch code setup
# ----------------

ln -sf ${TORCH_DIR}/cool.dat
cp -nv ${TORCH_DIR}/bridge_multiples.py .
cp -nv ${TORCH_DIR}/flash.par.turbsph_standard flash.par
cp -nv ${TORCH_DIR}/run.sh .
