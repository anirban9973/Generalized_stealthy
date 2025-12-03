#!/bin/bash

# commandline arguments
timelimit=4				    	# Time limit of simulation in hours
seed=1234						# random seed to generate initial conditions 
beg_idx=0 						# the index at which the configurations start to be loaded.
Verbosity=2						# verbosity of output messages  
pot_shape=flat					# (optional) the shape of collective-coordinate potential. flat/power-law/overlap are possible.

# Interactive arguments.
# Parameters of a packing (random configrations)
dimension=2						# space dimension (1,2,3)
init=random
num_radii=2						# the number of different radii (>= 1)
N1=30; r1=1						# the number of particles of radius 1; the relative size of radius 1.
N2=70; r2=0.5					# the number of particles of radius 2; the relative size of radius 2. [You need these variables as many as ${num_radii}].
phi=0.304							# target packing fraction  (particle radii are scaled to meet this target packing fraction.)
val=100.0							# relative strength of soft-core contact potential.

# Parameters of stealthy potential
num_species=2					# the number of species (This value can be different from ${num_radii})
#N1=30							# particle number of species 1.
S0_1=0.0						# equiluminousity of species 1
type1=chi						# unit of constraints for species 1 (chi/wavenumber)
lower1=0.0						# lower bound
upper1=0.15						# upper bound

#N2=70							
S0_2=0.0
type2=chi
lower2=0.0
upper2=0.15

threads=4
savename=${dimension}D_bidisperse_${upper1}_${N1}_${N2}

Nc=1000 		# target number of configurations
mode=ground
tol="tolerance 1e-18"
step="maxsteps 10000"

EXC=./CCO.out

# random initial condition
$EXC multi $timelimit $seed $beg_idx $Verbosity $pot_shape <<< "$init 
$dimension 
$num_radii
$N1	$r1
$N2 $r2
$phi
$val
$num_species
$N1 $S0_1 $type1 $lower1 $upper1 
$N2 $S0_2 $type2 $lower2 $upper2 
$threads
$savename
$Nc
$mode 
$tol
$step
run
"

