#!/bin/sh
# This script includes MPI setups for different systems. 
# Used by run.sh

case "$SYSTEM" in
    default)
        MPI_LAUNCHER="mpirun"
        MPI_ARGS=""
        ;;
    snellius)
        MPI_LAUNCHER="mpirun"
        MPI_ARGS="--mca orte_base_help_aggregate 0"
        ;;
    *)
        echo "Unknown system: $SYSTEM"
        exit 1
        ;;
esac