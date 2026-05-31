"""
Create folder structure for hole statistics sweep over chi0 and chi_gen,
copy CCO.out and ex_hole_slurm.sh into each chi_gen folder, write params.dat,
and submit sbatch jobs.

Folder structure:
    chi0_<chi0>/chi_gen_<chi_gen>/
        CCO.out
        ex_hole_slurm.sh
        params.dat      <- line 1: chi0, line 2: chi_gen, line 3: Rf_start

Rf_start is taken from R_c in ../N_400/chi0_X/chi_gen_Y/Rc_biased.dat if it
exists; otherwise defaults to 0.1.

chi0 values  : 0.1, 0.2, 0.3, 0.4
chi_gen range: from chi0 up to 0.49 in steps of 0.05 (last value always 0.49)
"""

import os
import shutil
import subprocess

# ---- parameters ----
chi0_values   = [0.1, 0.2, 0.3, 0.4]
chi_gen_pool  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.46, 0.47, 0.48, 0.49]

src_exe    = "CCO.out"
src_slurm  = "ex_hole_slurm.sh"

# Directory containing previous Rc_biased.dat results
ref_dir    = "../N_400"
rf_default = 0.1

# --------------------

def fmt(v):
    """Format float without trailing zeros: 0.10 -> '0.1', 0.15 -> '0.15'."""
    return f"{v:g}"

def read_Rc(path):
    """Return R_c (col 4) from first data line of Rc_biased.dat, or None."""
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    return float(line.split()[3])
    except (FileNotFoundError, IndexError, ValueError):
        pass
    return None

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

        # Look for Rc from a previous run
        ref_rc_path = os.path.join(ref_dir, chi0_dir,
                                   f"chi_gen_{fmt(chi_gen)}", "Rc_biased.dat")
        Rc = read_Rc(ref_rc_path)
        Rf_start = Rc if Rc is not None else rf_default
        if Rc is not None:
            print(f"  Found Rc={Rc} from {ref_rc_path} -> Rf_start={Rf_start}")
        else:
            print(f"  No Rc_biased.dat found -> Rf_start={rf_default} (default)")

        # Write params.dat: line1=chi0, line2=chi_gen, line3=Rf_start
        with open(os.path.join(run_dir, "params.dat"), "w") as f:
            f.write(f"{chi0}\n{chi_gen}\n{Rf_start}\n")

        print(f"Submitting: {run_dir}")
        result = subprocess.run(
            ["sbatch", src_slurm],
            cwd=run_dir,
            capture_output=True, text=True
        )
        print(f"  {result.stdout.strip()}")
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr.strip()}")
