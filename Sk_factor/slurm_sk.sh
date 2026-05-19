#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1                    # Request 1 GPU
#SBATCH --cpus-per-task=1               # CPUs for data loading
#SBATCH --mem=32G                       # Memory
#SBATCH --output=output_structure_factor.txt
#SBATCH --error=error_structure_factor.txt
#SBATCH --time=00:30:00
#SBATCH --job-name=structure_factor
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=am4197@princeton.edu

set -e

format_time() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    printf "%02d:%02d:%02d" $hours $minutes $seconds
}

module load anaconda3/2025.12
conda activate data_analysis

echo "Job ID: $SLURM_JOB_ID"
echo "Node: $HOSTNAME"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "Start time: $(date)"

which python
python --version

# Check GPU availability
nvidia-smi

echo "Starting structure factor calculation..."
compute_start=$(date +%s)

python Sk.py 

compute_end=$(date +%s)
compute_time=$((compute_end - compute_start))

echo "Computation time: $(format_time $compute_time)"
echo "Done!"
