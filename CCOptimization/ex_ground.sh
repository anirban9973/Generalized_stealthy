#!/bin/bash

# --------------------------------
# 1. commandline arguments
# --------------------------------
timelimit=1				    	# Time limit of simulation in hours
seed=1234						# random seed to generate initial conditions 
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=3						# verbosity of output messages  

# Interactive arguments.
fname=test		# save file

# --------------------------------
# 2. Parameters of pair potential.
# --------------------------------
d=1                            # space dimension

# Generalized stealthy: M annular shells, each defined by (K1, delta) in unit number density.
# S(k) = 0 for k in union of [K1^(n), K1^(n) + delta^(n)], n = 1..M.
# To use a single standard stealthy shell: set M=1, K1s=(0.0), deltas=(your_delta).
M=2                            # number of shells
K1s=(0.0  0.5)                 # lower bounds k1^(n) per shell (unit number density)
deltas=(0.3  0.2)              # widths delta^(n) per shell (unit number density)
S0s=(0.0  0.0)                 # S0 per shell: 0.0 = stealthy, >0 = equiluminous
vareps0=0.0			# relative strength of the soft-core repulsion. (0 means no soft-core repulsions.)
phi_fake=0.15                      # (Do not touch) fictitious packing fraction; the particle radius is computed from this value in the unit number density
sigma=0.20                     # Exclusion radius of the soft-core repulsion in unit number density; this value must be larger than the particle diameter
angle=90.0                     # Box angle in degrees (2D only): 90=square, 60=hexagonal; ignored for d!=2

# --------------------------------
# 3. Computational parameters
# --------------------------------
threads=2                      # Number of threads in openmp
N=100	                       # Number of particles
Nc=10                          # Number of configurations
initial=random                 # Choice of initial condition (random/input)

# Minimization parameters
# - These parameters can be skipped. They are also insensitive in order and case-insensitive.
eps0="tolerance 1e-14"		# energy tolerance for ground states in the unit of v0. Default value is 1e-14.
max_eval="maxsteps 100000" 	# the maximum number of evaluations before quitting calculations. Default value is 10000.
algorithm="algorithm LBFGS"	# Minimization algorithm. Default option is LBFGS. Other available options: LocalGradientDescent, ConjugateGradient, SteepestDescent, MINOP


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
    K1a=`echo "${K1s[$n]}*${a}"    | bc -l`
    deltaa=`echo "${deltas[$n]}*${a}" | bc -l`
    shell_str="${shell_str} ${K1a} ${deltaa} ${S0s[$n]}"
done
# ------------------------------------


# Filenames
exe=./CCO.out

echo "# Run the following command:"

# ground states from Random Initial Conditions
${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $vareps0 ${sigma} ${phi_fake} ${angle} \
${threads} ${N} ${Nc} random ${fname}_GS \
ground $eps0 $max_eval $algorithm run"


# ground states from Input Initial Conditions
#loadconfig=${fname}_GS # the name of loaded ConfigPack file.
#${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} <<< "${d} ${shell_str} $S0 $vareps0 ${sigma} ${phi_fake} \
#${threads} ${N} ${Nc} input ${loadconfig} ${fname}_GS2 \
#ground $eps0 $max_eval $algorithm run" > log3

