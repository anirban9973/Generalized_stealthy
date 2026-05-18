import numpy as np
import cupy as cp
from pathlib import Path
import sys
import glob

# =============================================================================
# PARAMETERS
# =============================================================================
dim  = 2
nmax = 1000

p_values   = [1]
file_lists = {}

for p in p_values:
    file_lists[p] = glob.glob(f"../collected_data/config_no*.dat")


def read_basis(filepath, dim):
    """Read the d basis vectors from a .dat file (rows 4..4+dim-1)."""
    basis = []
    for i in range(dim):
        row = np.loadtxt(filepath, skiprows=4 + i, max_rows=1)
        basis.append(row)
    return np.array(basis)   # shape (dim, dim), rows = basis vectors


def reciprocal_basis(A):
    """Compute reciprocal basis B = 2pi * inv(A).T; rows of B are b_i."""
    return 2 * np.pi * np.linalg.inv(A).T


for p in p_values:
    files = file_lists[p]

    if len(files) == 0:
        print(f"No files found for p = {p}, skipping.")
        sys.stdout.flush()
        continue

    print(f"\n{'='*60}")
    print(f"Processing p = {p}  ({len(files)} configs)")
    print(f"{'='*60}")
    sys.stdout.flush()

    # --- read box geometry from first file ---
    first_file = files[0]
    N = int(np.loadtxt(first_file, skiprows=2, max_rows=1))

    A = read_basis(first_file, dim)          # real-space basis (rows = a_i)
    B = reciprocal_basis(A)                  # reciprocal basis (rows = b_i)
    k_unit = np.linalg.norm(B[0])            # |b1|: fundamental k-spacing

    print(f"N = {N}")
    print(f"Real-space basis vectors:")
    for i, ai in enumerate(A):
        print(f"  a{i+1} = {ai}")
    print(f"Reciprocal basis vectors:")
    for i, bi in enumerate(B):
        print(f"  b{i+1} = {bi}  (|b{i+1}| = {np.linalg.norm(bi):.6f})")
    print(f"k_unit = |b1| = {k_unit:.6f}")

    shell_width = k_unit / 5
    print(f"shell_width = k_unit/5 = {shell_width:.6f}")
    print(f"nmax = {nmax}\n")
    sys.stdout.flush()

    shell_centers = k_unit * np.arange(1, nmax + 1)

    # =========================================================================
    # GENERATE ALL K-VECTORS ON GPU — integer combinations of reciprocal basis
    # =========================================================================
    print("Generating k-vectors on GPU...")
    sys.stdout.flush()

    b1_gpu = cp.array(B[0])
    b2_gpu = cp.array(B[1])

    n_range          = cp.arange(-nmax, nmax + 1)
    nx_grid, ny_grid = cp.meshgrid(n_range, n_range, indexing='ij')
    nx_flat          = nx_grid.flatten()
    ny_flat          = ny_grid.flatten()

    mask_nonzero = (nx_flat != 0) | (ny_flat != 0)
    nx_flat      = nx_flat[mask_nonzero]
    ny_flat      = ny_flat[mask_nonzero]

    # k = n1*b1 + n2*b2  (works for any box geometry)
    k_vecs_gpu = nx_flat[:, None] * b1_gpu + ny_flat[:, None] * b2_gpu
    k_mags_gpu = cp.linalg.norm(k_vecs_gpu, axis=1)
    n_k        = len(k_vecs_gpu)

    print(f"Generated {n_k} k-vectors")
    sys.stdout.flush()

    # =========================================================================
    # PRE-FILTER K-VECTORS INTO SHELLS
    # =========================================================================
    print("Pre-filtering k-vectors into shells...")
    sys.stdout.flush()

    shell_k_indices = []
    shell_k_counts  = []

    for n in range(1, nmax + 1):
        k_center = k_unit * n
        k_lower  = k_center - shell_width / 2
        k_upper  = k_center + shell_width / 2

        mask    = (k_mags_gpu >= k_lower) & (k_mags_gpu < k_upper)
        indices = cp.where(mask)[0]

        shell_k_indices.append(indices)
        shell_k_counts.append(len(indices))

    # =========================================================================
    # LOOP OVER CONFIGS, ACCUMULATE S(k)
    # =========================================================================
    S_k_accum   = np.zeros(nmax)
    n_processed = 0

    for file in files:
        if not Path(file).exists():
            print(f"  File not found: {file}, skipping.")
            sys.stdout.flush()
            continue

        print(f"  Reading: {Path(file).name}")
        sys.stdout.flush()

        positions     = np.loadtxt(file, skiprows=7)
        positions_gpu = cp.array(positions)

        S_k = np.zeros(nmax)

        for shell_idx in range(nmax):
            indices = shell_k_indices[shell_idx]
            n_shell = len(indices)

            if n_shell == 0:
                S_k[shell_idx] = np.nan
                continue

            k_vecs_shell = k_vecs_gpu[indices]
            phases       = cp.dot(positions_gpu, k_vecs_shell.T)
            sum_exp      = cp.sum(cp.exp(1j * phases), axis=0)
            S_k_shell    = cp.abs(sum_exp)**2 / N
            S_k[shell_idx] = float(cp.mean(S_k_shell))

        S_k_accum += S_k
        n_processed += 1

    if n_processed == 0:
        print(f"No valid configs processed for p = {p}, skipping output.")
        sys.stdout.flush()
        continue

    S_k_avg = S_k_accum / n_processed

    print(f"\n  Averaged over {n_processed} configs.")
    sys.stdout.flush()

    # =========================================================================
    # SAVE
    # =========================================================================
    output_data = np.column_stack([
        shell_centers,
        S_k_avg,
        shell_k_counts
    ])

    outfile = f"structure_factor.dat"
    np.savetxt(outfile, output_data,
               header="k  S(k)  k_count",
               fmt='%.17g %.17g %d')

    print(f"Written: {outfile}")
    sys.stdout.flush()
