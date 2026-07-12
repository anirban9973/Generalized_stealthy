#!/bin/bash
#SBATCH --job-name=CCO_quench
#SBATCH --array=1-80               # 1:1 with the 80 thermal packs from ex_crystalmd_slurm.sh
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=20         # OpenMP threads per task; also sets $SLURM_CPUS_PER_TASK
#SBATCH --time=6:00:00             # wall-clock limit (hh:mm:ss); MUST match timelimit below
#SBATCH --mem=16G                  # memory per task
#SBATCH --output=log_quench/%x_%A_%a.out   # separate log dir from the MD run
#SBATCH --error=log_quench/%x_%A_%a.err

mkdir -p log_quench

# --------------------------------
# 1. Derived from SLURM environment
# --------------------------------
threads=${SLURM_CPUS_PER_TASK}
Verbosity=3
# task i quenches ALL the thermal snapshots of MD trajectory i.
#   loads  crystalmd_run<i>.ConfigPack                 (READ ONLY - the thermalized
#                                                       configs are never modified)
#   writes ground_crystalmd_run<i>_Success.ConfigPack  (distinct prefix, so the thermal
#                                                       packs can never be overwritten;
#                                                       configpack_to_h5 default mode
#                                                       already gathers *_Success)
loadname=crystalmd_run${SLURM_ARRAY_TASK_ID}
savename=ground_crystalmd_run${SLURM_ARRAY_TASK_ID}

# --------------------------------
# 2. Parameters read from params.dat (written by setup_crystalmd_runs.py)
# --------------------------------
# line 1: N (unused here - N comes from the loaded configs); line 2: chi
chi=`awk 'NR==2' params.dat`

# --------------------------------
# 3. Fixed quench parameters
# --------------------------------
max_steps=100000               # L-BFGS eval ceiling per config (they start near the minimum)
tolerance=1e-16                # accept a quenched config if Phi < tolerance
timelimit=6                    # hours (MUST match --time above)

# --------------------------------
# 4. Run
# --------------------------------
exe=./CCO.out

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Threads       : ${threads}"
echo "chi (target)  : ${chi}"
echo "Load          : ${loadname}.ConfigPack   (read only)"
echo "Save          : ${savename}_Success.ConfigPack"
echo "tolerance     : ${tolerance}"

time_start=$(date +%s)

# stdin order matches GetCrystalQuench: chi threads loadname savename max_steps tolerance
${exe} quench ${timelimit} 0 ${Verbosity} <<< "${chi} ${threads} ${loadname} ${savename} \
${max_steps} ${tolerance}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
