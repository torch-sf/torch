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

# macOS fix: printf does not work for macOS
if find --version >/dev/null 2>&1 ; then
    FIND_CMD='find "${TORCH_DIR}/src/flash" -type f -printf "%P\0"'
    BACKUP_FIND='find "${FLASH_DIR}" -name "*.replaced-by-torch" -printf "%P\0"'
    IS_GNU=1
else
    echo
    echo "System is not GNU, will try alternative solution for find"
    FIND_CMD='find "${TORCH_DIR}/src/flash" -type f -print0'
    BACKUP_FIND='find "${FLASH_DIR}" -name "*.replaced-by-torch" -print0'
    IS_GNU=0
fi

echo
echo "Relinking Torch Fortran files into FLASH..."

# Link Torch files into the FLASH directory
# Existing files will be renamed and replaced by a link
while IFS= read -r -d $'\0' file <&3; do
    
    # macOS fix: strip path prefix
    if [ "$IS_GNU" -eq 0 ]; then
        file="${file#${TORCH_DIR}/src/flash/}"
    fi

    source="$(realpath "${TORCH_DIR}/src/flash/${file}")"
    target="${FLASH_DIR}/${file}"
    mkdir -p "$(dirname "${target}")"
    if [ ! -L "${target}" ] ; then
        if [ -e "${target}" ] ; then
            mv "${target}" "${target}.replaced-by-torch"
        fi
        ln -s "${source}" "${target}"
    fi
done 3< <(eval "$FIND_CMD")


# Restore original FLASH file if the Torch version no longer exists
while IFS= read -r -d $'\0' file <&3; do

    # macOS fix: strip path prefix
    if [ "$IS_GNU" -eq 0 ]; then
        file="${file#${FLASH_DIR}/}"
    fi

    backup="$(realpath "${FLASH_DIR}/${file}")"
    link="${backup%.replaced-by-torch}"
    if [ ! -e "${link}" ] ; then
        rm -f "${link}"
        mv "${backup}" "${link}"
    fi
done 3< <(eval "$BACKUP_FIND")

