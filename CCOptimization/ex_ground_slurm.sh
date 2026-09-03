#!/bin/bash
# IMPORTANT: run  `mkdir -p logs`  BEFORE `sbatch`, because SLURM opens the
# --output/--error paths below before this script executes; the in-script mkdir
# would be too late if logs/ does not already exist.
#SBATCH --job-name=CCO_ground
#SBATCH --array=1-10               # array indices; each index = one independent run
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=4          # OpenMP threads per task; also sets $SLURM_CPUS_PER_TASK
#SBATCH --time=04:00:00            # wall-clock limit (hh:mm:ss); MUST match timelimit below
#SBATCH --mem=4G                   # memory per task
#SBATCH --output=logs/%x_%A_%a.out # stdout: jobname_arrayjobid_taskid.out
#SBATCH --error=logs/%x_%A_%a.err  # stderr

mkdir -p logs   # kept for reruns; the pre-sbatch mkdir above is the one that matters

# --------------------------------
# 1. Derived from SLURM environment
# --------------------------------
# seed = seed_base * array task ID  (each task gets a unique, reproducible seed)
seed_base=1000
seed=$(( seed_base * SLURM_ARRAY_TASK_ID ))

# threads = number of CPUs allocated to this task by SLURM
threads=${SLURM_CPUS_PER_TASK}

# each array task writes to its own save file to avoid conflicts
beg_idx=0
Verbosity=3
fname=GS_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Parameters of pair potential.
# --------------------------------
d=2                            # space dimension

# Generalized stealthy: M annular shells, each defined by (K1, delta) in unit number density.
# S(k) = 0 for k in union of [K1^(n), K1^(n) + delta^(n)], n = 1..M.
# To use a single standard stealthy shell: set M=1, K1s=(0.0), deltas=(your_delta).
M=1                            # number of shells
K1s=(0.5)                      # lower bound k1 (unit number density); K1>0 => stealthy
deltas=(0.2)                   # width delta => annulus [0.5, 0.7]        NONHYPERUNIFORM
S0s=(0.0)                      # S0=0 stealthy. (K1>0 leaves k->0 free, so NOT hyperuniform,
                               #  which is the regime where the hole penalty matters.)
vareps0=0.0                    # relative strength of the soft-core repulsion. (0 means no soft-core repulsions.)
phi_fake=0.15                  # (Do not touch) fictitious packing fraction
sigma=0.20                     # exclusion radius of soft-core repulsion (unit number density)

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=4                    # simulation time limit in hours (MUST match --time=04:00:00)
N=4000                          # number of particles
Nc=1                          # number of configurations per array task

# Minimization parameters
eps0="tolerance 1e-16"
max_eval="maxsteps 100000"
algorithm="algorithm LBFGS"

# Hole-constraint parameters (used by ground_search.py)
Rmax=1.2                       # max allowed hole radius (unit number density)
lambda_h=1.0                   # hole penalty weight
Mgrid="64 128 256"             # probe-grid refinement schedule

# ------------- do not touch ---------
pi=`echo "4*a(1)" | bc -l`
if [ "$d" == 1 ]; then
v=2.0
elif [ "$d" == 2 ]; then
v=${pi}
elif [ "$d" == 3 ]; then
v=`echo "4.*${pi}/3."| bc -l`
fi
a=`echo "e(l(${phi_fake}/${v})/${d})" | bc -l`

# Build shell string: "M  K1a_0 deltaa_0 S0_0  K1a_1 deltaa_1 S0_1 ..."
shell_str="${M}"
for (( n=0; n<M; n++ )); do
    K1a=`echo "${K1s[$n]}*${a}"       | bc -l`
    deltaa=`echo "${deltas[$n]}*${a}" | bc -l`
    shell_str="${shell_str} ${K1a} ${deltaa} ${S0s[$n]}"
done
# ------------------------------------

exe="python3 ground_search.py"        # was ./CCO.out

# numpy/BLAS honor these for the vectorized k-space term (matches old -fopenmp)
export OMP_NUM_THREADS=${threads}
export OPENBLAS_NUM_THREADS=${threads}
export MKL_NUM_THREADS=${threads}

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Seed          : ${seed}"
echo "Threads       : ${threads}"
echo "Save prefix   : ${fname}_GS"

time_start=$(date +%s)

# Same CLI + stdin as CCO.out, with the hole params (Rmax lambda_h grid) appended.
${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_GS \
ground $eps0 $max_eval $algorithm run ${Rmax} ${lambda_h} ${Mgrid}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
