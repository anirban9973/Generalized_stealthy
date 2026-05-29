#!/bin/bash
#SBATCH --job-name=CCO_hole
#SBATCH --array=1               # array indices; each index = one independent seed
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=40          # OpenMP threads per task
#SBATCH --time=10:00:00            # wall-clock limit; should match timelimit below
#SBATCH --mem=32G                   # memory per task
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

mkdir -p logs

# --------------------------------
# 1. Derived from SLURM environment
# --------------------------------
seed_base=1000
seed=$(( seed_base * SLURM_ARRAY_TASK_ID ))
threads=${SLURM_CPUS_PER_TASK}
Verbosity=3
fname=hole_run${SLURM_ARRAY_TASK_ID}

# --------------------------------
# 2. Shell potential parameters
# --------------------------------
d=2                            # space dimension

# Single stealthy shell benchmark: S(k) = 0 for |k| in [0, K]
# K = 4*sqrt(pi*chi) at unit number density  (Torquato convention: chi = K^2 / (16*pi))
chi=0.3                        # stealthiness parameter

M=1                            # number of shells
S0s=(0.0)                      # stealthy (S0 = 0)
vareps0=0.0                    # soft-core repulsion strength (0 = off)
phi_fake=0.15                  # fictitious packing fraction (do not touch)
sigma=0.20                     # soft-core exclusion radius (unit number density)

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=10                   # simulation time limit in hours
N=400                          # number of particles

# Hole sweep parameters
Rf_start=1.0                   # initial exclusion radius (physical units)
delta_Rf=0.01                  # R_f step size (physical units)
N_trial=100                    # L-BFGS trials per R_f value
max_steps=100000               # L-BFGS steps per trial

# --------------------------------
# 4. Derived quantities (do not touch)
# --------------------------------
pi=`echo "4*a(1)" | bc -l`
if [ "$d" == 1 ]; then
    v=2.0
elif [ "$d" == 2 ]; then
    v=${pi}
elif [ "$d" == 3 ]; then
    v=`echo "4.*${pi}/3." | bc -l`
fi
a=`echo "e(l(${phi_fake}/${v})/${d})" | bc -l`

# K = 4*sqrt(pi*chi) at unit number density; K1=0 so delta=K
K=`echo "4*sqrt(${pi}*${chi})" | bc -l`
K1s=(0.0)
deltas=(${K})

echo "chi           : ${chi}"
echo "K (unit rho)  : ${K}"

# Build shell string: "M  K1a_0 deltaa_0 S0_0 ..."
shell_str="${M}"
for (( n=0; n<M; n++ )); do
    K1a=`echo "${K1s[$n]}*${a}"       | bc -l`
    deltaa=`echo "${deltas[$n]}*${a}" | bc -l`
    shell_str="${shell_str} ${K1a} ${deltaa} ${S0s[$n]}"
done

exe=./CCO.out

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Seed          : ${seed}"
echo "Threads       : ${threads}"
echo "Save prefix   : ${fname}"
echo "Rf_start      : ${Rf_start}"
echo "N_trial       : ${N_trial}"
echo "max_steps     : ${max_steps}"

time_start=$(date +%s)

${exe} hole ${timelimit} ${seed} ${Verbosity} <<< "${d} ${shell_str} \
${vareps0} ${sigma} ${phi_fake} ${threads} ${N} \
${fname} ${Rf_start} ${delta_Rf} ${N_trial} ${max_steps}"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
