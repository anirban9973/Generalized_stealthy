#!/bin/bash
# Simple local runner for ground_search.py on a laptop (no SLURM, no module load).
# Usage:  ./run_local.sh [seed]        (default seed = 1)
# Needs:  python3 with numpy, scipy, h5py, and ground_search.py in this directory.
set -e

seed=${1:-1}
threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
beg_idx=0
Verbosity=3
fname=GS_local                 # output -> GS_local.h5

# --------------------------------
# Parameters (edit here) -- same meaning as the slurm script
# --------------------------------
d=2
M=1                            # exactly ONE stealthy annulus (auto-Rmax)
K1s=(1.0)                      # annulus [1.0, 4.59], Keff=3.59 -> Rmax~1.31 (binds)
deltas=(3.59)
S0s=(0.0)
vareps0=0.0
phi_fake=0.15                  # (do not touch) fictitious packing fraction
sigma=0.20

timelimit=4                    # hours
N=500                          # smaller than the cluster N for laptop speed
Nc=1

eps0="tolerance 1e-16"
max_eval="maxsteps 100000"
algorithm="algorithm LBFGS"

lambda_h=1.0                   # hole penalty weight
epsilon_grid=0.1               # coarser grid than the cluster -> faster on a laptop

# --------------------------------
# a-scaling (identical to the slurm): script multiplies k by a; the code divides back
# --------------------------------
pi=$(echo "4*a(1)" | bc -l)
if   [ "$d" = 1 ]; then v=2.0
elif [ "$d" = 2 ]; then v=${pi}
else                    v=$(echo "4.*${pi}/3." | bc -l); fi
a=$(echo "e(l(${phi_fake}/${v})/${d})" | bc -l)

shell_str="${M}"
for (( n=0; n<M; n++ )); do
    K1a=$(echo   "${K1s[$n]}*${a}"   | bc -l)
    deltaa=$(echo "${deltas[$n]}*${a}" | bc -l)
    shell_str="${shell_str} ${K1a} ${deltaa} ${S0s[$n]}"
done

# numpy/BLAS threads
export OMP_NUM_THREADS=${threads}
export OPENBLAS_NUM_THREADS=${threads}
export MKL_NUM_THREADS=${threads}

echo "seed=${seed}, threads=${threads}, N=${N}  ->  ${fname}.h5"

python3 ground_search.py ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< \
"${d} ${shell_str} ${vareps0} ${sigma} ${phi_fake} ${threads} ${N} ${Nc} \
random ${fname} ground ${eps0} ${max_eval} ${algorithm} run ${lambda_h} ${epsilon_grid}"
