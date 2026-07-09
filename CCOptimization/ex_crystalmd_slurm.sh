#!/bin/bash
#SBATCH --job-name=CCO_crystalmd
#SBATCH --array=1-10               # 10 independent MD trajectories (fresh velocities per task)
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=20         # OpenMP threads per task; also sets $SLURM_CPUS_PER_TASK
#SBATCH --time=20:00:00            # wall-clock limit (hh:mm:ss); should match timelimit below
#SBATCH --mem=4G                   # memory per task (MD of 625 particles is light)
#SBATCH --output=logs/%x_%A_%a.out # stdout: jobname_arrayjobid_taskid.out
#SBATCH --error=logs/%x_%A_%a.err  # stderr

mkdir -p logs

# --------------------------------
# 1. Derived from SLURM environment
# --------------------------------
# Each array task gets a UNIQUE seed => different Maxwell-Boltzmann velocities
# => an independent MD trajectory. 10 tasks x 1000 samples = 10000 per chi.
seed_base=1000
seed=$(( seed_base * SLURM_ARRAY_TASK_ID ))
threads=${SLURM_CPUS_PER_TASK}
Verbosity=3
fname=crystalmd_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Parameters read from params.dat (written by setup_crystalmd_runs.py)
# --------------------------------
# line 1: N   = points per side (total points = N*N = 625 for N=25)
# line 2: chi = target stealthiness; the code derives K from the discrete k-count
N=`awk 'NR==1' params.dat`
chi=`awk 'NR==2' params.dat`

# --------------------------------
# 3. Fixed MD parameters
# --------------------------------
T_E=1e-6                       # MD temperature (dimensionless)
timestep=0.05                  # FIXED timestep (auto-tuning is off; ~6x below the
                               # stability edge ~0.3, stable and reproducible)
steps_per_sample=2000          # MD steps between saved snapshots (spacing)
num_samples=1000               # snapshots saved per task  (x10 tasks = 10000 per chi)
equil_samples=200              # warm-up rounds before sampling (each = steps_per_sample steps)
timelimit=20                   # simulation time limit in hours (matches --time above)

# --------------------------------
# 4. Run
# --------------------------------
exe=./CCO.out

echo "Array task ID     : ${SLURM_ARRAY_TASK_ID}"
echo "Seed              : ${seed}"
echo "Threads           : ${threads}"
echo "N (per side)      : ${N}   (total = $((N*N)))"
echo "chi (target)      : ${chi}"
echo "T_E               : ${T_E}"
echo "steps_per_sample  : ${steps_per_sample}"
echo "num_samples       : ${num_samples}"
echo "equil_samples     : ${equil_samples}"
echo "Save prefix       : ${fname}"

time_start=$(date +%s)

# stdin order matches GetCrystalMD:
# N chi T_E timestep steps_per_sample num_samples equil_samples threads savename
${exe} crystalmd ${timelimit} ${seed} ${Verbosity} <<< "${N} ${chi} ${T_E} ${timestep} \
${steps_per_sample} ${num_samples} ${equil_samples} ${threads} ${fname}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
