"""
Set up and submit crystalline (chi > 1/2) stealthy ground-state runs.

Edit ONLY this file: set N (points per side) and the list of chi values.
For each chi it creates a folder, drops in CCO.out + ex_crystal_slurm.sh, writes
params.dat, and submits the sbatch array.

Per-folder file written:
    params.dat   <- line 1: N (points per side; total points = N*N)
                    line 2: chi (target stealthiness)

The slurm script reads N and chi from params.dat; everything else (Nc, sigma_pert,
max_steps, threads, timelimit, verbosity) is fixed inside ex_crystal_slurm.sh.
"""

import os
import shutil
import subprocess

# ============================================================
#  EDIT HERE
# ============================================================
N = 50                                              # points per side (total = N*N)
chi_values = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]   # target stealthiness values
# ============================================================

src_exe   = "CCO.out"
src_slurm = "ex_crystal_slurm.sh"


def fmt(v):
    """Format float without trailing zeros: 0.55 -> '0.55'."""
    return f"{v:g}"


for chi in chi_values:
    run_dir = f"chi_{fmt(chi)}"
    os.makedirs(run_dir, exist_ok=True)

    shutil.copy(src_exe,   run_dir)
    shutil.copy(src_slurm, run_dir)

    # params.dat: line 1 = N (per side), line 2 = chi
    with open(os.path.join(run_dir, "params.dat"), "w") as f:
        f.write(f"{N}\n{chi}\n")

    print(f"Submitting: {run_dir}  (N/side={N}, total={N*N}, chi={chi})")
    result = subprocess.run(
        ["sbatch", src_slurm],
        cwd=run_dir,
        capture_output=True, text=True
    )
    print(f"  {result.stdout.strip()}")
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
