#!/bin/bash
#SBATCH --job-name=CCO_hole
#SBATCH --array=1-100              # array indices; each index = one independent R_c measurement
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
fname=hole_run${SLURM_ARRAY_TASK_ID}   # unique output prefix per array task

# --------------------------------
# 2. Parameters of pair potential.
# --------------------------------
d=2                            # space dimension

# Generalized stealthy: M annular shells, each defined by (K1, delta) in unit number density.
# S(k) = 0 for k in union of [K1^(n), K1^(n) + delta^(n)], n = 1..M.
# To use a single standard stealthy shell: set M=1, K1s=(0.0), deltas=(your_delta).
M=2                            # number of stealthy shells
K1s=(0.0  3.2354553131)        # lower bounds k1^(n) for each shell (unit number density)
deltas=(2.2419964866  2.2419964866)  # widths delta^(n) for each shell (unit number density)

S0=0.                          # S0 > 0: equiluminous. S0 = 0.0: stealthy
vareps0=0.0                    # relative strength of soft-core repulsion (0 = off)
phi_fake=0.15                  # (Do not touch) fictitious packing fraction; used only for unit conversion
sigma=0.20                     # exclusion radius of soft-core repulsion (unit number density)

# --------------------------------
# 3. Computational parameters
# --------------------------------
timelimit=4                    # simulation time limit in hours; should match #SBATCH --time
N=100                          # number of particles
Nc=1                           # placeholder; not used in hole mode

# --------------------------------
# 4. Hole scan parameters
# (inserted between "hole" and "run"; order- and case-insensitive; all optional)
# --------------------------------
R_min="rmin 0.0"               # starting probe radius (absolute units)
R_max="rmax 1.0"               # upper bound for scan (absolute units)
dR="dr 0.01"                   # step size
N_trial="ntrial 100"           # LBFGS attempts per R_f before declaring failure
tol="tolerance 1e-10"          # energy threshold for a feasible hole
max_steps="maxsteps 100000"    # max LBFGS evaluations per attempt

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
echo "Save prefix   : ${fname}_hole"

time_start=$(date +%s)

${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $S0 $vareps0 ${sigma} ${phi_fake} \
${threads} ${N} ${Nc} random ${fname}_hole \
hole $R_min $R_max $dR $N_trial $tol $max_steps run"

time_end=$(date +%s)
elapsed=$(( time_end - time_start ))
hrs=$(( elapsed / 3600 ))
mins=$(( (elapsed % 3600) / 60 ))
secs=$(( elapsed % 60 ))
echo "Total wall-clock time: ${hrs}h ${mins}m ${secs}s"
