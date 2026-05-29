"""
Create folder structure for hole statistics sweep over chi0 and chi_gen,
copy CCO.out and ex_hole_slurm.sh into each chi_gen folder, write params.dat,
and submit sbatch jobs.

Folder structure:
    chi0_<chi0>/chi_gen_<chi_gen>/
        CCO.out
        ex_hole_slurm.sh
        params.dat      <- line 1: chi0, line 2: chi_gen

chi0 values  : 0.1, 0.2, 0.3, 0.4
chi_gen range: from chi0 up to 0.49 in steps of 0.05 (last value always 0.49)
"""

import os
import shutil
import subprocess

# ---- parameters ----
chi0_values   = [0.1, 0.2, 0.3, 0.4]
chi_gen_pool  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.49]

src_exe    = "CCO.out"
src_slurm  = "ex_hole_slurm.sh"

# --------------------

def fmt(v):
    """Format float without trailing zeros: 0.10 -> '0.1', 0.15 -> '0.15'."""
    return f"{v:g}"

for chi0 in chi0_values:
    chi0_dir = f"chi0_{fmt(chi0)}"
    os.makedirs(chi0_dir, exist_ok=True)

    chi_gen_values = [cg for cg in chi_gen_pool if cg >= chi0 - 1e-9]

    for chi_gen in chi_gen_values:
        run_dir = os.path.join(chi0_dir, f"chi_gen_{fmt(chi_gen)}")
        os.makedirs(run_dir, exist_ok=True)

        # Copy executables
        shutil.copy(src_exe,   run_dir)
        shutil.copy(src_slurm, run_dir)

        # Write params.dat
        with open(os.path.join(run_dir, "params.dat"), "w") as f:
            f.write(f"{chi0}\n{chi_gen}\n")

        print(f"Submitting: {run_dir}")
        result = subprocess.run(
            ["sbatch", src_slurm],
            cwd=run_dir,
            capture_output=True, text=True
        )
        print(f"  {result.stdout.strip()}")
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr.strip()}")
