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

(
    cd ${FLASH_DIR} || { echo $errstr; exit 255; }

    # Patch obsolete Python and C code to work with modern tools
    patch -p1 --forward <"${TORCH_DIR}/support/flash-4.6.2.patch"
)


# Symlink Torch files into FLASH
./relink.sh || { echo $errstr; exit 255; }


#### ----------------------------------------------
#### Detect dependencies and set up torch_auto site

echo
echo "Detecting dependencies and creating FLASH Makefile:"
mkdir -p "${FLASH_DIR}"/sites/torch_auto
(cd "${TORCH_DIR}"/support && ./configure) || { echo "Error detecting dependencies"; }
cp "${TORCH_DIR}"/support/Makefile.h "${FLASH_DIR}"/sites/torch_auto/Makefile.h


#### --------------------------
#### Prepare the AMUSE bindings

adest="${TORCH_DIR}/src/amuse"

# Point AMUSE Makefile to FLASH directory via symlink
# TODO will not work? if user placed FLASH4.6.2/ in ${adest}/src

echo
echo "Preparing torch-amuse-flash package for installation:"
ln -sfTv ${FLASH_DIR} ${adest}/torch_amuse_flash/src/FLASH4.6.2 || { echo $errstr; exit 255; }


#########################
## VorAMR-Lite Install ##
#########################
## Dev: Sean C. Lewis ##
#echo "Installing VorAMR-Lite interface"
#asrc="${TORCH_DIR}/src/amuse/voramrLite" || { echo $errstr; exit 255; }
#adest="${AMUSE_DIR}/src/amuse/community/voramr" || { echo $errstr; exit 255; }

#mkdir -v ${adest}

#touch ${adest}/__init__.py || { echo $errstr; exit 255; }
#cp -v ${asrc}/base_grid_interface.F90   ${adest}/ || { echo $errstr; exit 255; }
#cp -v ${asrc}/interface.F90             ${adest}/ || { echo $errstr; exit 255; }
#cp -v ${asrc}/interface.py              ${adest}/ || { echo $errstr; exit 255; }
#cp -v ${asrc}/Makefile.prototype        ${adest}/Makefile || { echo $errstr; exit 255; }

#mkdir -p ${adest}/src
#cp -v ${asrc}/src/*                     ${adest}/src/ || { echo $errstr; exit 255; }

# Point AMUSE Makefile to FLASH directory via symlink
# TODO will not work? if user placed FLASH4.6.2/ in ${adest}/src
#echo -n "Linking: "
#ln -sfTv ${FLASH_DIR}                    ${adest}/src/FLASH4.6.2 || { echo $errstr; exit 255; }
#echo "Setting drive_loc in ${adest}/Makefile"
# "sed -i" doesn't work for BSD sed (e.g. on OS X), use a workaround from
# https://stackoverflow.com/a/44877280
#sed "s/^drive_loc\s*=.*/drive_loc = src\/FLASH4.6.2\/object/" ${adest}/Makefile > ${adest}/Makefile.$$ \
#  && mv "${adest}/Makefile.$$" "${adest}/Makefile" \
#  || { echo $errstr; exit 255; }


#### -----------
#### All done!!!

echo
echo "Torch install complete!  Ready to configure and make FLASH."
