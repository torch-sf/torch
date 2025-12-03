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
# Existing links will be replaced unless they point into TORCH_DIR already
while IFS= read -r -d $'\0' file <&3; do
    source="$(realpath "${file}")"
    target="${FLASH_DIR}/${file#${TORCH_DIR}/src/flash/}"
    mkdir -p "$(dirname "${target}")"
    if [ -L "${target}" ] ; then
        link_target=$(readlink -n "${target}")
        if [ ${link_target#${TORCH_DIR}/src/flash/} = ${link_target} ] ; then
            # link does not point into TORCH_DIR
            rm "${target}"
        fi
    fi

    if [ -f "${target}" ] ; then
        mv "${target}" "${target}.replaced-by-torch"
    fi

    if [ ! -e "${target}" ] ; then
        ln -s "${source}" "${target}"
    fi
done 3< <(find "${TORCH_DIR}"/src/flash -type f -print0)


# Restore original FLASH file if the Torch version no longer exists
# Also removes the symlink if there was no original file
while IFS= read -r -d $'\0' link <&3; do
    link_target=$(readlink -n "${link}")
    if printf '%s' "${link_target}" | grep '/' 2>&1 >/dev/null ; then
        if [ "${link_target#${TORCH_DIR}/src/flash/}" = "${link_target}" ] ; then
            # This link points outside the FLASH tree, so was added by us
            # FLASH has some internal symlinks, which we don't touch
            rm "${link}"
            backup_file="${link}.replaced-by-torch"
            if [ -f "${backup_file}" ] ; then
                mv "${backup_file}" "${link}"
            fi
        fi
    fi
done 3< <(find "${FLASH_DIR}" -type l -print0)

