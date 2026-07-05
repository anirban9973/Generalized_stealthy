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
# 2. Crystal potential parameters
# --------------------------------
d=2                            # space dimension (GetCrystalCCO supports d=2 only)

# Commensurate rhombic (triangular) cell: Nx x Ny particles, N = Nx*Ny.
# Constraint radius K is read from params.dat (line 1: K, line 2: chi_actual),
# precomputed by compute_K_crystal.py / setup_crystal_runs.py so that chi lands
# on the desired plateau and K stays below the first Bragg peak.
K=`awk 'NR==1' params.dat`
chi=`awk 'NR==2' params.dat`

Nx=50                          # particles along a1
Ny=50                          # particles along a2  (N = Nx*Ny = 2500)
val=0.0                        # soft-core repulsion strength (0 = off, pure stealthy)
sigma=0.20                     # soft-core exclusion radius (unused when val=0)

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=20                   # simulation time limit in hours
Nc=1                           # target number of accepted (Phi<1e-16) configs per task
max_steps=10000000             # L-BFGS steps per relaxation attempt
sigma_pert=0.03                # perturbation std (in units of lattice constant a)

# --------------------------------
# 4. Run
# --------------------------------
exe=./CCO.out

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Seed          : ${seed}"
echo "Threads       : ${threads}"
echo "Save prefix   : ${fname}"
echo "Nx, Ny        : ${Nx}, ${Ny}  (N = $((Nx*Ny)))"
echo "K             : ${K}"
echo "chi (target)  : ${chi}"
echo "sigma_pert    : ${sigma_pert}"

time_start=$(date +%s)

${exe} crystal ${timelimit} ${seed} ${Verbosity} <<< "${d} ${Nx} ${Ny} ${K} ${val} ${sigma} ${threads} \
${fname} ${Nc} ${max_steps} ${sigma_pert}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
