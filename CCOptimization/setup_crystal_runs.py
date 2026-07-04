"""
Create folder structure for crystalline (chi > 1/2) ground-state runs,
copy CCO.out and ex_ground_slurm.sh into each chi folder, write params.dat,
and submit sbatch jobs.

Folder structure:
    chi_<chi>/
        CCO.out
        ex_ground_slurm.sh
        params.dat      <- line 1: chi (= chi0), line 2: chi (= chi_gen)
                           chi0 = chi_gen => k1 = 0 (standard stealthy disk [0, K])

chi values: 0.55, 0.60, 0.65, 0.70, 0.75, 0.80
"""

import os
import shutil
import subprocess

# ---- parameters ----
chi_values = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]

src_exe   = "CCO.out"
src_slurm = "ex_ground_slurm.sh"

# --------------------

def fmt(v):
    """Format float without trailing zeros: 0.55 -> '0.55'."""
    return f"{v:g}"

for chi in chi_values:
    run_dir = f"chi_{fmt(chi)}"
    os.makedirs(run_dir, exist_ok=True)

    # Copy executables
    shutil.copy(src_exe,   run_dir)
    shutil.copy(src_slurm, run_dir)

    # chi0 = chi_gen = chi => k1 = 0 (standard stealthy disk)
    with open(os.path.join(run_dir, "params.dat"), "w") as f:
        f.write(f"{chi}\n{chi}\n")

    print(f"Submitting: {run_dir}  (chi={chi})")
    result = subprocess.run(
        ["sbatch", src_slurm],
        cwd=run_dir,
        capture_output=True, text=True
    )
    print(f"  {result.stdout.strip()}")
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
