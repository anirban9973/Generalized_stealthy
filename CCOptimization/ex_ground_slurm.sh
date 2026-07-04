#!/bin/bash
#SBATCH --job-name=CCO_ground
#SBATCH --array=1-10               # array indices; each index = one independent run
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
# seed = seed_base * array task ID  (each task gets a unique, reproducible seed)
seed_base=1000
seed=$(( seed_base * SLURM_ARRAY_TASK_ID ))

# threads = number of CPUs allocated to this task by SLURM
threads=${SLURM_CPUS_PER_TASK}

# each array task writes to its own save file to avoid conflicts
beg_idx=0
Verbosity=10
fname=GS_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Shell potential parameters
# --------------------------------
d=2                            # space dimension

# Single annular stealthy shell: S(k) = 0 for |k| in [k1, k1+delta]
# k1    = 2*sqrt(pi) * (chi_gen - chi0) / sqrt(chi0)   (unit number density)
# delta = 4*sqrt(pi * chi0)                             (unit number density)
# k2    = k1 + delta
# chi0 and chi_gen are read from params.dat (line 1: chi0, line 2: chi_gen)
chi0=`awk 'NR==1' params.dat`
chi_gen=`awk 'NR==2' params.dat`

M=1                            # number of shells
S0s=(0.0)                      # stealthy (S0 = 0)
vareps0=0.0                    # soft-core repulsion strength (0 = off)
phi_fake=0.15                  # fictitious packing fraction (do not touch)
sigma=0.20                     # soft-core exclusion radius (unit number density)

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=20                   # simulation time limit in hours
N=1600                         # number of particles
Nc=100000                      # large number; time limit is the actual stopper

# Minimization parameters
eps0="tolerance 1e-16"
max_eval="maxsteps 10000000"
algorithm="algorithm LBFGS"

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

# k1 = 2*sqrt(pi) * (chi_gen - chi0) / sqrt(chi0)
# delta = 4*sqrt(pi*chi0)
k1=`echo "2*sqrt(${pi})*(${chi_gen}-${chi0})/sqrt(${chi0})" | bc -l`
delta=`echo "4*sqrt(${pi}*${chi0})" | bc -l`
K1s=(${k1})
deltas=(${delta})

echo "chi_gen       : ${chi_gen}"
echo "chi0          : ${chi0}"
echo "k1 (unit rho) : ${k1}"
echo "delta (unit rho): ${delta}"

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
echo "Save prefix   : ${fname}_GS"
echo "N             : ${N}"

time_start=$(date +%s)

${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_GS \
ground $eps0 $max_eval $algorithm run"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
