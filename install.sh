#!/usr/bin/env bash

# Install: copy AMUSE coupling code from Torch into separate FLASH and AMUSE
# repositories.

errstr="Stopping, Torch install failed!"

echo "Got AMUSE_DIR: $AMUSE_DIR"
echo "Got FLASH_DIR: $FLASH_DIR"
echo "Got TORCH_DIR: $TORCH_DIR"

# https://stackoverflow.com/a/226724
while true; do
    read -p "Install Torch? [y/n] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer y (yes) or n (no).";;
    esac
done

#### ----------------------------
#### Prepare the FLASH repository

cd ${FLASH_DIR} || { echo $errstr; exit 255; }

# for FLASH4.5 and earlier only, fixed in FLASH4.6
# -O prevents wget from writing fname.1, fname.2, etc on successive calls
wget -nv http://flash.uchicago.edu/site/flashcode/user_support/FLASH4.5-a.diff -O FLASH4.5-a.diff || { echo $errstr; exit 255; }
# this will return exit code >0 when patch is skipped over, so cannot exit on
# error code here, and cannot use "set -e" in this script.
patch -p0 -r - --forward < FLASH4.5-a.diff

# Patch for parallel HDF5 1.10.x, fixed in FLASH4.6
# http://flash.uchicago.edu/pipermail/flash-users/2018-May/002626.html
wget -nv http://flash.uchicago.edu/pipermail/flash-users/attachments/20180519/fdf3cd9b/attachment.obj -O FLASH4.5_parallelHDF5.diff || { echo $errstr; exit 255; }
patch -p0 -r - --forward < FLASH4.5_parallelHDF5.diff

rsync -avh "${TORCH_DIR}/src/flash/" "${FLASH_DIR}/." || { echo $errstr; exit 255; }

#### ----------------------------
#### Prepare the AMUSE repository

asrc="${TORCH_DIR}/src/amuse" || { echo $errstr; exit 255; }
adest="${AMUSE_DIR}/src/amuse/community/flash" || { echo $errstr; exit 255; }

# to allow python imports
touch ${adest}/__init__.py || { echo $errstr; exit 255; }

cp -v ${asrc}/base_grid_interface.F90   ${adest}/ || { echo $errstr; exit 255; }
cp -v ${asrc}/interface.F90             ${adest}/ || { echo $errstr; exit 255; }
cp -v ${asrc}/interface.py              ${adest}/ || { echo $errstr; exit 255; }
cp -v ${asrc}/Makefile.prototype        ${adest}/Makefile || { echo $errstr; exit 255; }

mkdir -p ${adest}/src
cp -v ${asrc}/src/*                     ${adest}/src/ || { echo $errstr; exit 255; }

# Point AMUSE Makefile to FLASH directory via symlink
# TODO will not work? if user placed FLASH4.5/ in ${adest}/src
echo -n "Linking: "
ln -sfTv ${FLASH_DIR}                    ${adest}/src/FLASH4.5 || { echo $errstr; exit 255; }
echo "Setting drive_loc in ${adest}/Makefile"
# "sed -i" doesn't work for BSD sed (e.g. on OS X), use a workaround from
# https://stackoverflow.com/a/44877280
sed "s/^drive_loc\s*=.*/drive_loc = src\/FLASH4.5\/object/" ${adest}/Makefile > ${adest}/Makefile.$$ \
  && mv "${adest}/Makefile.$$" "${adest}/Makefile" \
  || { echo $errstr; exit 255; }

#### -----------
#### All done!!!

echo "Torch install complete!  Ready to configure and make FLASH and AMUSE."
