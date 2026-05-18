#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --output=configpack_to_h5_%A.out
#SBATCH --error=configpack_to_h5_%A.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

set -e

# --------------------------------
# Activate spack and load libraries
# --------------------------------
source ~/spack/share/spack/setup-env.sh
spack load boost gsl hdf5

echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# --------------------------------
# Run
# --------------------------------
./configpack_to_h5

echo "Done."
