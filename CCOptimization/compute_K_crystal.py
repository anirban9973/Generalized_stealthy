"""
Tabulate the discrete stealthiness parameter chi(K) for a rhombic (triangular)
simulation cell and pick, for each target chi, the constraint radius K that
lands on the plateau nearest the target.

Cell (unit number density rho = 1):
    Triangular primitive vectors  a1 = a*(1, 0),  a2 = a*(1/2, sqrt(3)/2)
    with a = sqrt(2/sqrt(3))  so that the area per particle is 1.
    Simulation box                A1 = Nx*a1,  A2 = Ny*a2   (Nx*Ny particles).

The allowed constraint wavevectors are the reciprocal lattice of the box,
k = m*B1 + n*B2 with (m, n) != (0, 0).  Matching the C++ code (GetKs keeps one
of each +/-k pair), the number of independent constrained modes inside a disk of
radius K is  M_indep(K) = (# nonzero reciprocal points with |k| <= K) / 2,
and   chi(K) = M_indep(K) / (d * (N - 1)),  d = 2.

The triangular lattice is an exact Phi = 0 stealthy ground state as long as K is
below the first Bragg peak |b1| = 4*pi/(a*sqrt(3)); this script flags any K that
violates that bound.

Output: K_values.dat  with columns  chi_target  chi_actual  K
"""

import numpy as np

# ---- cell parameters ----
Nx, Ny = 50, 50
N = Nx * Ny            # 2500
d = 2
a = np.sqrt(2.0 / np.sqrt(3.0))      # triangular lattice constant at rho = 1

chi_targets = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]

# ---- box and reciprocal basis ----
A1 = Nx * a * np.array([1.0, 0.0])
A2 = Ny * a * np.array([0.5, np.sqrt(3.0) / 2.0])
A  = np.column_stack([A1, A2])       # columns are box basis vectors
B  = 2.0 * np.pi * np.linalg.inv(A).T
B1, B2 = B[:, 0], B[:, 1]            # reciprocal basis vectors

# first Bragg peak of the triangular particle lattice
k_bragg = 4.0 * np.pi / (a * np.sqrt(3.0))

# ---- enumerate reciprocal-lattice points and their radii ----
Mrange = 120                         # |m|,|n| range; |B1|~0.117 so covers K up to ~14
ms = np.arange(-Mrange, Mrange + 1)
mm, nn = np.meshgrid(ms, ms)
mm, nn = mm.ravel(), nn.ravel()
kx = mm * B1[0] + nn * B2[0]
ky = mm * B1[1] + nn * B2[1]
kr = np.sqrt(kx * kx + ky * ky)
kr = np.sort(kr[kr > 1e-9])          # drop the origin

# ---- group into shells; cumulative independent count and chi per shell ----
tol = 1e-6
shell_r, shell_cum = [], []
i, n = 0, len(kr)
while i < n:
    r = kr[i]
    j = i
    while j < n and kr[j] - r < tol:
        j += 1
    shell_r.append(r)
    shell_cum.append(j)              # total points with |k| <= r
    i = j
shell_r = np.array(shell_r)
M_indep = np.array(shell_cum) / 2.0  # +/- pairing -> independent modes
chi_shell = M_indep / (d * (N - 1))

# ---- pick K per target ----
print(f"N = {N}  (Nx={Nx}, Ny={Ny}),  a = {a:.6f},  k_Bragg = {k_bragg:.6f}\n")
print(f"{'chi_target':>10}  {'chi_actual':>10}  {'K':>10}  {'K/k_Bragg':>10}  {'M_indep':>8}")

rows = []
for chi_t in chi_targets:
    s = int(np.argmin(np.abs(chi_shell - chi_t)))
    # K on the plateau: midpoint between this shell and the next
    if s + 1 < len(shell_r):
        K = 0.5 * (shell_r[s] + shell_r[s + 1])
    else:
        K = shell_r[s] * 1.0001
    chi_a = chi_shell[s]
    flag = "  <-- K >= k_Bragg!" if K >= k_bragg else ""
    print(f"{chi_t:>10.4f}  {chi_a:>10.6f}  {K:>10.6f}  "
          f"{K / k_bragg:>10.4f}  {int(M_indep[s]):>8d}{flag}")
    rows.append((chi_t, chi_a, K))

# ---- write output ----
with open("K_values.dat", "w") as f:
    f.write("# chi_target\tchi_actual\tK\n")
    for chi_t, chi_a, K in rows:
        f.write(f"{chi_t}\t{chi_a:.10f}\t{K:.10f}\n")

print("\nWrote K_values.dat")
