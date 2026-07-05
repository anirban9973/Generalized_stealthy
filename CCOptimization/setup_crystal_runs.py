"""
Create folder structure for crystalline (chi > 1/2) stealthy ground-state runs,
copy CCO.out and ex_crystal_slurm.sh into each chi folder, write params.dat,
and submit sbatch jobs.

Reads K_values.dat (produced by compute_K_crystal.py) with columns:
    chi_target  chi_actual  K

Folder structure:
    chi_<chi_target>/
        CCO.out
        ex_crystal_slurm.sh
        params.dat      <- line 1: K, line 2: chi_actual

Each chi submits an array of 10 independent tasks (see ex_crystal_slurm.sh).
"""

import os
import shutil
import subprocess

# ---- inputs ----
kfile     = "K_values.dat"
src_exe   = "CCO.out"
src_slurm = "ex_crystal_slurm.sh"

# --------------------

def fmt(v):
    """Format float without trailing zeros: 0.55 -> '0.55'."""
    return f"{v:g}"

# Read K values
rows = []
with open(kfile) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        chi_target, chi_actual, K = line.split()
        rows.append((float(chi_target), float(chi_actual), float(K)))

if not rows:
    raise SystemExit(f"No data found in {kfile}. Run compute_K_crystal.py first.")

for chi_target, chi_actual, K in rows:
    run_dir = f"chi_{fmt(chi_target)}"
    os.makedirs(run_dir, exist_ok=True)

    # Copy executables
    shutil.copy(src_exe,   run_dir)
    shutil.copy(src_slurm, run_dir)

    # Write params.dat: line1 = K, line2 = chi_actual
    with open(os.path.join(run_dir, "params.dat"), "w") as f:
        f.write(f"{K:.10f}\n{chi_actual:.10f}\n")

    print(f"Submitting: {run_dir}  (chi_actual={chi_actual:.6f}, K={K:.6f})")
    result = subprocess.run(
        ["sbatch", src_slurm],
        cwd=run_dir,
        capture_output=True, text=True
    )
    print(f"  {result.stdout.strip()}")
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
