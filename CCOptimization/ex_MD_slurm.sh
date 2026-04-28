#!/bin/bash
#SBATCH --job-name=CCO_MD
#SBATCH --array=1-20               # array indices; each index = one independent run
#SBATCH --ntasks=1                 # one process per array task
#SBATCH --cpus-per-task=4          # OpenMP threads per task; also sets $SLURM_CPUS_PER_TASK
#SBATCH --time=04:00:00            # wall-clock limit (hh:mm:ss); should match timelimit below
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

beg_idx=0
Verbosity=3
fname=MD_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Parameters of pair potential.
# --------------------------------
d=2                            # space dimension

# Generalized stealthy: M annular shells, each defined by (K1, delta) in unit number density.
# S(k) = 0 for k in union of [K1^(n), K1^(n) + delta^(n)], n = 1..M.
# To use a single standard stealthy shell: set M=1, K1s=(0.0), deltas=(your_delta).
M=2                            # number of stealthy shells
K1s=(0.0  0.5)                 # lower bounds k1^(n) for each shell (unit number density)
deltas=(0.3  0.2)              # widths delta^(n) for each shell (unit number density)

S0=0.                          # S0 > 0: equiluminous. S0 = 0.0: stealthy
vareps0=0.0                    # relative strength of soft-core repulsion. (0 means no soft-core repulsion.)
phi_fake=0.15                  # (Do not touch) fictitious packing fraction; used only for unit conversion
sigma=0.20                     # exclusion radius of soft-core repulsion (unit number density); must be > particle diameter

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=4                    # simulation time limit in hours; should match #SBATCH --time
N=1000                          # number of particles
Nc=1                          # number of configurations to collect

# --------------------------------
# 4. MD simulation parameters
# (inserted between "MD" and "run"; order- and case-insensitive; all optional)
# --------------------------------
T="temperature 1e-3"           # dimensionless temperature in units of v0. Default: dimension-dependent small value.
dt="timestep 0.05"             # MD time step. Default: 0.01.
int_samp="samplesteps 1000"    # MD steps between saved snapshots. Default: 1e5.
num_equi="numequilibration 10" # number of equilibration blocks (each = samplesteps steps). Default: 500.
restore="restore"              # save checkpoint to .MDDump; load it on restart. Remove if not needed.

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

# Build shell string: "M  K1a_0 deltaa_0  K1a_1 deltaa_1 ..."
shell_str="${M}"
for (( n=0; n<M; n++ )); do
    K1a=`echo "${K1s[$n]}*${a}"       | bc -l`
    deltaa=`echo "${deltas[$n]}*${a}" | bc -l`
    shell_str="${shell_str} ${K1a} ${deltaa}"
done
# ------------------------------------

exe=./CCO.out

echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Seed          : ${seed}"
echo "Threads       : ${threads}"
echo "Save prefix   : ${fname}_MD"

time_start=$(date +%s)

${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $S0 $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_MD \
MD $T $dt $int_samp $num_equi $restore run"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
