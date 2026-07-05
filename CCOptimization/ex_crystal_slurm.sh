#!/bin/bash
#SBATCH --job-name=CCO_crystal
#SBATCH --array=1-10               # array indices; each index = one independent seed
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=20         # OpenMP threads per task; also sets $SLURM_CPUS_PER_TASK
#SBATCH --time=20:00:00            # wall-clock limit (hh:mm:ss); should match timelimit below
#SBATCH --mem=4G                   # memory per task
#SBATCH --output=logs/%x_%A_%a.out # stdout: jobname_arrayjobid_taskid.out
#SBATCH --error=logs/%x_%A_%a.err  # stderr

mkdir -p logs

# --------------------------------
# 1. Derived from SLURM environment
# --------------------------------
seed_base=1000
seed=$(( seed_base * SLURM_ARRAY_TASK_ID ))
threads=${SLURM_CPUS_PER_TASK}
Verbosity=3
fname=crystal_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Parameters read from params.dat (written by setup_crystal_runs.py)
# --------------------------------
# line 1: N   = points per side (total points = N*N)
# line 2: chi = target stealthiness; the code derives K from the discrete k-count
N=`awk 'NR==1' params.dat`
chi=`awk 'NR==2' params.dat`

# --------------------------------
# 3. Fixed crystal parameters
# --------------------------------
Nc=10000                       # number of perturbed-lattice attempts
sigma_pert=0.03                # perturbation std (in units of the lattice constant a)
max_steps=1000000              # L-BFGS eval ceiling per attempt (optimizer stops early at a
                               # local minimum; then the next attempt is a fresh perturbed lattice)
timelimit=20                   # simulation time limit in hours (matches --time above)

# --------------------------------
# 4. Run
# --------------------------------
exe=./CCO.out

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Seed          : ${seed}"
echo "Threads       : ${threads}"
echo "N (per side)  : ${N}   (total = $((N*N)))"
echo "chi (target)  : ${chi}"
echo "Nc            : ${Nc}"
echo "Save prefix   : ${fname}"

time_start=$(date +%s)

# stdin order matches GetCrystalCCO: N Nc sigma_pert chi threads savename max_steps
${exe} crystal ${timelimit} ${seed} ${Verbosity} <<< "${N} ${Nc} ${sigma_pert} ${chi} ${threads} ${fname} ${max_steps}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
