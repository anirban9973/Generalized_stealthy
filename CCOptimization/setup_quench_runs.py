"""
Submit the QUENCH pass over the thermalized MD snapshots.

Run this AFTER setup_crystalmd_runs.py has finished, i.e. once each chi folder
already contains the thermal packs crystalmd_run<i>.ConfigPack.

For each chi folder this script:
  - checks the folder and its params.dat exist (params.dat is NOT overwritten --
    the quench must use exactly the chi the MD used),
  - counts the thermal packs and skips the folder if there are none,
  - copies the freshly built CCO.out + ex_crystalquench_slurm.sh in,
  - submits the sbatch array (task i quenches crystalmd_run<i>.ConfigPack).

Each task READS crystalmd_run<i>.ConfigPack and WRITES
ground_crystalmd_run<i>_Success.ConfigPack -- the thermalized packs are never
modified, so they remain available for later analysis.

Edit ONLY the chi_values list below (it must match what the MD was run with).
"""

import os
import glob
import shutil
import subprocess

# ============================================================
#  EDIT HERE  (must match the chi values the MD was run with)
# ============================================================
chi_values = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]
# ============================================================

src_exe   = "CCO.out"
src_slurm = "ex_crystalquench_slurm.sh"


def fmt(v):
    """Format float without trailing zeros: 0.55 -> '0.55'."""
    return f"{v:g}"


for chi in chi_values:
    run_dir = f"chi_{fmt(chi)}"

    if not os.path.isdir(run_dir):
        print(f"SKIP {run_dir}: folder does not exist (run the MD first).")
        continue

    params = os.path.join(run_dir, "params.dat")
    if not os.path.isfile(params):
        print(f"SKIP {run_dir}: params.dat missing (the quench reads chi from it).")
        continue

    # thermal packs written by the MD: crystalmd_run<i>.ConfigPack
    packs = [p for p in glob.glob(os.path.join(run_dir, "crystalmd_run*.ConfigPack"))
             if "_Success" not in os.path.basename(p)]
    if not packs:
        print(f"SKIP {run_dir}: no crystalmd_run*.ConfigPack found (MD not done?).")
        continue

    shutil.copy(src_exe,   run_dir)      # fresh binary (must include the 'quench' mode)
    shutil.copy(src_slurm, run_dir)

    print(f"Submitting: {run_dir}  (chi={chi}, {len(packs)} thermal packs to quench)")
    result = subprocess.run(
        ["sbatch", src_slurm],
        cwd=run_dir,
        capture_output=True, text=True
    )
    print(f"  {result.stdout.strip()}")
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
