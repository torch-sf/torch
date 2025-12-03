#!/bin/bash


# Get paths from command line if not set in the environment
if [ "x$TORCH_DIR" == "x" ] || [ "x$FLASH_DIR" == "x" ] ; then
    TORCH_DIR="$1"
    FLASH_DIR="$2"
fi

# Check that we have paths
if [ "x$TORCH_DIR" == "x" ] || [ "x$FLASH_DIR" == "x" ] ; then
    echo "Usage: $0 /path/to/torch /path/to/FLASH4.6.2"
    echo "Or set TORCH_DIR and FLASH_DIR in your environment and try again."
    exit 1
fi

# Try to catch incorrect paths
if [ ! -d "${TORCH_DIR}/src/flash" ] ; then
    echo "TORCH_DIR is set to $TORCH_DIR, which does not seem to be a torch repository"
    exit 1
fi

if [ ! -d "${FLASH_DIR}/source/flashUtilities" ] ; then
    echo "FLASH_DIR is set to $FLASH_DIR, which does not seem to be a FLASH directory"
    exit 1
fi

echo
echo "Relinking Torch Fortran files into FLASH..."

# Link Torch files into the FLASH directory
# Existing files will be renamed and replaced by a link
while IFS= read -r -d $'\0' file <&3; do
    source="$(realpath "${file}")"
    target="${FLASH_DIR}/${file#${TORCH_DIR}/src/flash/}"
    mkdir -p "$(dirname "${target}")"
    if [ ! -L "${target}" ] ; then
        if [ -e "${target}" ] ; then
            mv "${target}" "${target}.replaced-by-torch"
        fi
        ln -s "${source}" "${target}"
    fi
done 3< <(find "${TORCH_DIR}"/src/flash -type f -print0)


# Restore original FLASH file if the Torch version no longer exists
while IFS= read -r -d $'\0' file <&3; do
    backup="$(realpath "${file}")"
    link="${backup%.replaced-by-torch}"
    if [ ! -e "${link}" ] ; then
        rm -f "${link}"
        mv "${backup}" "${link}"
    fi
done 3< <(find "${FLASH_DIR}" -name '*.replaced-by-torch' -print0)

