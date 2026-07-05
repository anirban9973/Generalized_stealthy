"""
Set up and submit crystalline (chi > 1/2) stealthy ground-state runs.

Single knob: N (particles per side).  The simulation box is always a rhombic
(triangular) commensurate cell with Nx = Ny = N, so the total particle number is
N*N.  For each target chi this script computes the constraint radius K on the
discrete plateau nearest the target (the logic formerly in compute_K_crystal.py),
then creates one folder per chi, drops in CCO.out + ex_crystal_slurm.sh, and
submits an sbatch array.

Per-folder files written:
    N_value_parameter.dat   <- line 1: N        (read by the slurm script -> Nx, Ny)
    params.dat              <- line 1: K, line 2: chi_actual

Cell (unit number density rho = 1):
    a1 = a*(1, 0),  a2 = a*(1/2, sqrt(3)/2),  a = sqrt(2/sqrt(3))
    box  A1 = N*a1,  A2 = N*a2
Allowed constraint wavevectors are the box reciprocal lattice; matching the C++
code (GetKs keeps one of each +/-k pair),
    M_indep(K) = (# nonzero reciprocal points with |k| <= K) / 2,
    chi(K)     = M_indep(K) / (d * (Npart - 1)),   d = 2.
The triangular lattice is an exact Phi = 0 ground state only while K stays below
the first Bragg peak |b1| = 4*pi/(a*sqrt(3)); any violation is flagged.
"""

import os
import shutil
import subprocess
import numpy as np

# ---- user parameters ----
N           = 50            # particles per side; Nx = Ny = N, total = N*N
chi_targets = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]

src_exe   = "CCO.out"
src_slurm = "ex_crystal_slurm.sh"

# --------------------

def fmt(v):
    """Format float without trailing zeros: 0.55 -> '0.55'."""
    return f"{v:g}"

def compute_K_table(N, chi_targets):
    """Return (rows, k_bragg, Npart) where rows = [(chi_target, chi_actual, K), ...]."""
    Nx = Ny = N
    Npart = Nx * Ny
    d = 2
    a = np.sqrt(2.0 / np.sqrt(3.0))

    A1 = Nx * a * np.array([1.0, 0.0])
    A2 = Ny * a * np.array([0.5, np.sqrt(3.0) / 2.0])
    A  = np.column_stack([A1, A2])
    B  = 2.0 * np.pi * np.linalg.inv(A).T
    B1, B2 = B[:, 0], B[:, 1]
    k_bragg = 4.0 * np.pi / (a * np.sqrt(3.0))

    # enumerate reciprocal-lattice points out to a radius safely above the
    # largest target K (K ~ 4*sqrt(pi*chi)); |B1| = 2*pi/(N*a) sets the spacing
    K_max_est = 4.0 * np.sqrt(np.pi * max(chi_targets)) * 1.25
    spacing   = 2.0 * np.pi / (N * a)
    Mrange    = int(np.ceil(K_max_est / spacing)) * 2 + 5

    ms = np.arange(-Mrange, Mrange + 1)
    mm, nn = np.meshgrid(ms, ms)
    mm, nn = mm.ravel(), nn.ravel()
    kx = mm * B1[0] + nn * B2[0]
    ky = mm * B1[1] + nn * B2[1]
    kr = np.sqrt(kx * kx + ky * ky)
    kr = np.sort(kr[kr > 1e-9])

    # group into shells: cumulative independent count and chi per shell
    tol = 1e-6
    shell_r, shell_cum = [], []
    i, n = 0, len(kr)
    while i < n:
        r = kr[i]
        j = i
        while j < n and kr[j] - r < tol:
            j += 1
        shell_r.append(r)
        shell_cum.append(j)
        i = j
    shell_r = np.array(shell_r)
    M_indep = np.array(shell_cum) / 2.0
    chi_shell = M_indep / (d * (Npart - 1))

    rows = []
    print(f"N = {Npart}  (Nx=Ny={N}),  a = {a:.6f},  k_Bragg = {k_bragg:.6f}\n")
    print(f"{'chi_target':>10}  {'chi_actual':>10}  {'K':>10}  {'K/k_Bragg':>10}  {'M_indep':>8}")
    for chi_t in chi_targets:
        s = int(np.argmin(np.abs(chi_shell - chi_t)))
        if s + 1 < len(shell_r):
            K = 0.5 * (shell_r[s] + shell_r[s + 1])
        else:
            K = shell_r[s] * 1.0001
        chi_a = chi_shell[s]
        flag = "  <-- K >= k_Bragg!" if K >= k_bragg else ""
        print(f"{chi_t:>10.4f}  {chi_a:>10.6f}  {K:>10.6f}  "
              f"{K / k_bragg:>10.4f}  {int(M_indep[s]):>8d}{flag}")
        rows.append((chi_t, chi_a, K))
    print()
    return rows, k_bragg, Npart


if __name__ == "__main__":
    rows, k_bragg, Npart = compute_K_table(N, chi_targets)

    for chi_target, chi_actual, K in rows:
        run_dir = f"chi_{fmt(chi_target)}"
        os.makedirs(run_dir, exist_ok=True)

        shutil.copy(src_exe,   run_dir)
        shutil.copy(src_slurm, run_dir)

        with open(os.path.join(run_dir, "N_value_parameter.dat"), "w") as f:
            f.write(f"{N}\n")

        with open(os.path.join(run_dir, "params.dat"), "w") as f:
            f.write(f"{K:.17g}\n{chi_actual:.17g}\n")

        print(f"Submitting: {run_dir}  (chi_actual={chi_actual:.6f}, K={K:.6f})")
        result = subprocess.run(
            ["sbatch", src_slurm],
            cwd=run_dir,
            capture_output=True, text=True
        )
        print(f"  {result.stdout.strip()}")
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr.strip()}")
